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
        // 参照ではなくポインタとして直接アクセス
        const device Particle* p = &particles[i];
        float2 r = p->position - float2(pos);
        float r2 = dot(r, r);
        if (r2 > 0.0001f) {  // 特異点を避ける
            fieldValue += p->charge * normalize(r) / r2;
        }
    }
    
    fields[index].value = fieldValue;
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
        
        // 初期粒子分布の設定

    }
    return self;
}
// 時間更新
- (void)update:(double)dt
        withEMField:(EMField*)fld{

}
@end