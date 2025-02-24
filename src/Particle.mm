#import "Particle.h"
#import "Init.h"
#import "Constant.h"
#include <random>

// Metal シェーダーのソースコード
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

// 電場データ構造体
struct EMFieldData {
    device float* Ex;
    device float* Ey;
    device float* Ez;
    device float* Bx;
    device float* By;
    device float* Bz;
    int ngx;
    int ngy;
    float dx;
    float dy;
};

// シミュレーションパラメータ構造体
struct SimulationParams {
    int particleCount;
    float charge;
    float mass;
    float weight;
    float dt;
};

// 粒子更新カーネル
kernel void updateParticles(device ParticleState* particles [[buffer(0)]],
                          device const float* Ex [[buffer(1)]],
                          device const float* Ey [[buffer(2)]],
                          device const float* Ez [[buffer(3)]],
                          constant SimulationParams& params [[buffer(4)]],
                          uint id [[thread_position_in_grid]]) {
    if (id >= params.particleCount) return;
    
    device ParticleState& p = particles[id];
    
    // 電場の補間や力の計算をここに実装
    // (実際の物理計算はプロジェクトの要件に応じて実装)
    
    // 位置と速度の更新（単純な例）
    p.x += p.vx * params.dt;
    p.y += p.vy * params.dt;
    
    // 加速度の計算と速度の更新
    float qm = params.charge / params.mass;
    p.vx += qm * Ex[0] * params.dt;  // 簡略化のため、最も近いグリッドポイントの値を使用
    p.vy += qm * Ey[0] * params.dt;
    p.vz += qm * Ez[0] * params.dt;
}
)";

@implementation Particle {
    // プライベートインスタンス変数
    id<MTLBuffer> _ExBuffer;
    id<MTLBuffer> _EyBuffer;
    id<MTLBuffer> _EzBuffer;
    id<MTLBuffer> _BxBuffer;
    id<MTLBuffer> _ByBuffer;
    id<MTLBuffer> _BzBuffer;
    id<MTLBuffer> _paramsBuffer;
}

// 初期設定
- (instancetype)initWithDevice:(id<MTLDevice>)device
                withParticleParam:(ParamForParticle)ParticleParam
                withFieldParam:(ParamForField)FieldParam {
    self = [super init];
    if (self) {
        _device = device;
        _commandQueue = [device newCommandQueue];
        
        // 物理パラメータの設定
        _charge = ParticleParam.q;
        _mass = ParticleParam.m;
        _weight = ParticleParam.w;
        
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
        
        // パーティクルバッファの初期化
        _particleBuffer = [device newBufferWithLength:sizeof(ParticleState) * ParticleParam.pNum
                                            options:MTLResourceStorageModeShared];
        
        // 電場バッファの初期化
        size_t fieldSize = sizeof(float) * FieldParam.ngx * FieldParam.ngy;
        _ExBuffer = [device newBufferWithLength:fieldSize options:MTLResourceStorageModeShared];
        _EyBuffer = [device newBufferWithLength:fieldSize options:MTLResourceStorageModeShared];
        _EzBuffer = [device newBufferWithLength:fieldSize options:MTLResourceStorageModeShared];
        // 磁場バッファの初期化
        _BxBuffer = [device newBufferWithLength:fieldSize options:MTLResourceStorageModeShared];
        _ByBuffer = [device newBufferWithLength:fieldSize options:MTLResourceStorageModeShared];
        _BzBuffer = [device newBufferWithLength:fieldSize options:MTLResourceStorageModeShared];
        
        // シミュレーションパラメータバッファの初期化
        _paramsBuffer = [device newBufferWithLength:sizeof(SimulationParams)
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
- (void)update:(double)dt withEMField:(EMField*)fld {
    // シミュレーションパラメータの更新
    SimulationParams* params = (SimulationParams*)_paramsBuffer.contents;
    params->particleCount = (int)(_particleBuffer.length / sizeof(ParticleState));
    params->charge = _charge;
    params->mass = _mass;
    params->weight = _weight;
    params->dt = dt;
    
    // EMFieldのバッファを直接使用
    id<MTLBuffer> ExBuffer = [fld ExBuffer];
    id<MTLBuffer> EyBuffer = [fld EyBuffer];
    id<MTLBuffer> EzBuffer = [fld EzBuffer];
    id<MTLBuffer> BxBuffer = [fld BxBuffer];
    id<MTLBuffer> ByBuffer = [fld ByBuffer];
    id<MTLBuffer> BzBuffer = [fld BzBuffer];
    
    // コマンドバッファとエンコーダの作成
    id<MTLCommandBuffer> commandBuffer = [_commandQueue commandBuffer];
    id<MTLComputeCommandEncoder> computeEncoder = [commandBuffer computeCommandEncoder];
    
    // パイプラインとバッファの設定
    [computeEncoder setComputePipelineState:_updateParticlesPipeline];
    [computeEncoder setBuffer:_particleBuffer offset:0 atIndex:0];
    [computeEncoder setBuffer:ExBuffer offset:0 atIndex:1];
    [computeEncoder setBuffer:EyBuffer offset:0 atIndex:2];
    [computeEncoder setBuffer:EzBuffer offset:0 atIndex:3];
    [computeEncoder setBuffer:BxBuffer offset:0 atIndex:1];
    [computeEncoder setBuffer:ByBuffer offset:0 atIndex:2];
    [computeEncoder setBuffer:BzBuffer offset:0 atIndex:3];
    [computeEncoder setBuffer:_paramsBuffer offset:0 atIndex:4];
    
    // グリッドとスレッドグループのサイズ設定
    NSUInteger gridSize = params->particleCount;
    MTLSize threadGroupSize = MTLSizeMake(256, 1, 1);
    MTLSize gridSizeMetalStyle = MTLSizeMake((gridSize + threadGroupSize.width - 1) / threadGroupSize.width, 1, 1);
    
    // ディスパッチ
    [computeEncoder dispatchThreadgroups:gridSizeMetalStyle
                  threadsPerThreadgroup:threadGroupSize];
    
    // エンコーダと実行
    [computeEncoder endEncoding];
    [commandBuffer commit];
}

- (void)integrateChargeDensity:(EMField*)fld{
    
};

@end