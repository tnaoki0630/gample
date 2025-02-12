
#import "PIC.h"

// Metal Shaderコード
static NSString *const kMetalShaderSource = @R"(
#include <metal_stdlib>
using namespace metal;

struct Particle {
    float2 position;
    float2 velocity;
    float charge;
    float mass;
};

struct Field {
    float2 value;
};

// 粒子更新カーネル
kernel void updateParticles(device Particle* particles [[buffer(0)]],
                          device const Field* fields [[buffer(1)]],
                          const device uint& particleCount [[buffer(2)]],
                          const device float& deltaTime [[buffer(3)]],
                          const device float2& gridSize [[buffer(4)]],
                          uint id [[thread_position_in_grid]]) {
    if (id >= particleCount) return;
    
    Particle& p = particles[id];
    
    // 最近接格子点からの電場を計算
    uint2 gridPos = uint2(p.position);
    Field localField = fields[gridPos.y * uint(gridSize.x) + gridPos.x];
    
    // 運動方程式の解析
    float2 acceleration = (p.charge / p.mass) * localField.value;
    p.velocity += acceleration * deltaTime;
    p.position += p.velocity * deltaTime;
    
    // 周期的境界条件
    p.position = fmod(p.position + gridSize, gridSize);
}

// 電場計算カーネル
kernel void computeFields(device Field* fields [[buffer(0)]],
                        device const Particle* particles [[buffer(1)]],
                        const device uint& particleCount [[buffer(2)]],
                        const device float2& gridSize [[buffer(3)]],
                        uint2 pos [[thread_position_in_grid]]) {
    if (pos.x >= uint(gridSize.x) || pos.y >= uint(gridSize.y)) return;
    
    uint index = pos.y * uint(gridSize.x) + pos.x;
    float2 fieldValue = float2(0, 0);
    
    // 各粒子からのクーロン力を計算
    for (uint i = 0; i < particleCount; i++) {
        const Particle& p = particles[i];
        float2 r = p.position - float2(pos);
        float r2 = dot(r, r);
        if (r2 > 0.0001f) {  // 特異点を避ける
            fieldValue += p.charge * normalize(r) / r2;
        }
    }
    
    fields[index].value = fieldValue;
}
)";

@implementation PIC

- (instancetype)initWithDevice:(id<MTLDevice>)device
                particleCount:(NSUInteger)particleCount
                    gridSize:(NSUInteger)gridSize {
    self = [super init];
    if (self) {
        _device = device;
        _commandQueue = [device newCommandQueue];
        
        // コンピュートパイプラインの設定
        NSError *error = nil;
        id<MTLLibrary> library = [device newLibraryWithSource:kMetalShaderSource
                                                     options:nil
                                                       error:&error];
        if (!library) {
            NSLog(@"Failed to create Metal library: %@", error);
            return nil;
        }
        
        id<MTLFunction> updateParticlesFunction = [library newFunctionWithName:@"updateParticles"];
        id<MTLFunction> computeFieldsFunction = [library newFunctionWithName:@"computeFields"];
        
        _updateParticlesPipeline = [device newComputePipelineStateWithFunction:updateParticlesFunction
                                                                        error:&error];
        _computeFieldsPipeline = [device newComputePipelineStateWithFunction:computeFieldsFunction
                                                                      error:&error];
        
        // バッファの初期化
        _particleBuffer = [device newBufferWithLength:sizeof(Particle) * particleCount
                                            options:MTLResourceStorageModeShared];
        _fieldBuffer = [device newBufferWithLength:sizeof(Field) * gridSize * gridSize
                                         options:MTLResourceStorageModeShared];
        
        // 初期粒子分布の設定
        [self initializeParticles:particleCount gridSize:gridSize];
    }
    return self;
}

- (void)initializeParticles:(NSUInteger)particleCount gridSize:(NSUInteger)gridSize {
    Particle *particles = (Particle *)self.particleBuffer.contents;
    
    for (NSUInteger i = 0; i < particleCount; i++) {
        particles[i].position = (vector_float2){
            (float)(arc4random_uniform(gridSize)),
            (float)(arc4random_uniform(gridSize))
        };
        particles[i].velocity = (vector_float2){0.0f, 0.0f};
        particles[i].charge = (i % 2) * 2.0f - 1.0f;  // 正負交互に配置
        particles[i].mass = 1.0f;
    }
}

- (void)update:(float)deltaTime {
    id<MTLCommandBuffer> commandBuffer = [self.commandQueue commandBuffer];
    
    // 電場の計算
    {
        id<MTLComputeCommandEncoder> computeEncoder = [commandBuffer computeCommandEncoder];
        [computeEncoder setComputePipelineState:self.computeFieldsPipeline];
        [computeEncoder setBuffer:self.fieldBuffer offset:0 atIndex:0];
        [computeEncoder setBuffer:self.particleBuffer offset:0 atIndex:1];
        // ... その他のパラメータの設定 ...
        
        MTLSize gridSize = MTLSizeMake(32, 32, 1);
        MTLSize threadGroupSize = MTLSizeMake(8, 8, 1);
        [computeEncoder dispatchThreadgroups:gridSize
                    threadsPerThreadgroup:threadGroupSize];
        [computeEncoder endEncoding];
    }
    
    // 粒子の更新
    {
        id<MTLComputeCommandEncoder> computeEncoder = [commandBuffer computeCommandEncoder];
        [computeEncoder setComputePipelineState:self.updateParticlesPipeline];
        [computeEncoder setBuffer:self.particleBuffer offset:0 atIndex:0];
        [computeEncoder setBuffer:self.fieldBuffer offset:0 atIndex:1];
        // ... その他のパラメータの設定 ...
        
        MTLSize gridSize = MTLSizeMake(1024, 1, 1);
        MTLSize threadGroupSize = MTLSizeMake(256, 1, 1);
        [computeEncoder dispatchThreadgroups:gridSize
                    threadsPerThreadgroup:threadGroupSize];
        [computeEncoder endEncoding];
    }
    
    [commandBuffer commit];
}

@end
