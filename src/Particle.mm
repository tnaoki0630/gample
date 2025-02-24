#import "Particle.h"
#import "Init.h"
#import "Constant.h"
#include <random>

static NSString *const kMetalShaderSource = @R"(
#include <metal_stdlib>
using namespace metal;

struct ParticleState {
    float x;
    float y;
    float vx;
    float vy;
    float vz;
};

// 粒子更新カーネル
kernel void updateParticles(device ParticleState* particles[[buffer(0)]],
                          uint id [[thread_position_in_grid]]) {
    if (id >= 10000) return;
    
    // 参照ではなくポインタとして扱う
    device ParticleState* p = &particles[id];
}
)";

@implementation Particle
// 初期設定
- (instancetype)initWithDevice:(id<MTLDevice>)device
                withParticleParam:(ParamForParticle) ParticleParam
                withFieldParam:(ParamForField) FieldParam{
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
        _particleBuffer = [device newBufferWithLength:sizeof(ParticleState) * ParticleParam.pNum
                                            options:MTLResourceStorageModeShared];
        
        // 初期粒子分布の設定
        [self generateParticles:ParticleParam withFieldParam:FieldParam];
    }
    return self;
}

- (void)generateParticles:(ParamForParticle)ParticleParam
            withFieldParam:(ParamForField)FieldParam{
    // 乱数の初期化
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    // 粒子の初期化
    ParticleState *particles = (ParticleState *)self.particleBuffer.contents;
    double Lx = FieldParam.dx * FieldParam.ngx;
    double Ly = FieldParam.dy * FieldParam.ngy;
    double vth_epi  = sqrt(2*kb*ParticleParam.initT/ParticleParam.m);
    std::normal_distribution norm_dist(0.0, vth_epi);
    for (int i = 0; i < ParticleParam.pNum; i++) {
        particles[i].x = (double)arc4random_uniform(Lx);
        particles[i].y = (double)arc4random_uniform(Ly);
        particles[i].vx = (double)norm_dist(engine);
        particles[i].vy = (double)norm_dist(engine);
        particles[i].vz = (double)norm_dist(engine);
    }
}
// 時間更新
- (void)update:(double)dt
        withEMField:(EMField*)fld{

}
@end