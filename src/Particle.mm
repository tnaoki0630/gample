#import "Particle.h"

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
    
    // 参照ではなくポインタとして扱う
    device Particle* p = &particles[id];
    
    // 最近接格子点からの電場を計算
    uint2 gridPos = uint2(p->position);
    Field localField = fields[gridPos.y * uint(gridSize.x) + gridPos.x];
    
    // 運動方程式の解析
    float2 acceleration = (p->charge / p->mass) * localField.value;
    p->velocity += acceleration * deltaTime;
    p->position += p->velocity * deltaTime;
    
    // 周期的境界条件
    p->position = fmod(p->position + gridSize, gridSize);
}
)";

@implementation Particle
// 初期設定
- (instancetype)initWithDevice:(id<MTLDevice>)device
                initWithParam:(Init*) init{
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
        
        _updateParticlesPipeline = [device newComputePipelineStateWithFunction:updateParticlesFunction
                                                                        error:&error];
        
        // バッファの初期化
        // _particleBuffer = [device newBufferWithLength:sizeof(Particle) * particleCount
        //                                     options:MTLResourceStorageModeShared];
        
        // 初期粒子分布の設定
        // [self initializeParticles:particleCount gridSize:gridSize];
    }
    return self;
}
// 時間更新
- (void)update:(double)dt
        withEMField:(EMField*)fld{

}
@end