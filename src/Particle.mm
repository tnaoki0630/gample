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
};

// シミュレーションパラメータ構造体
struct SimulationParams {
    int particleCount;
    float charge;
    float mass;
    float weight;
    float dt;
    float constE;
    float constB;
    float constX;
    float constY;
    int ngx;
    int ngy;
    float dx;
    float dy;
};

// 粒子更新カーネル
kernel void updateParticles(device ParticleState* particles [[buffer(0)]],
                          constant SimulationParams& params [[buffer(1)]],
                          device const float* Ex [[buffer(2)]],
                          device const float* Ey [[buffer(3)]],
                          device const float* Ez [[buffer(4)]],
                          device const float* Bx [[buffer(5)]],
                          device const float* By [[buffer(6)]],
                          device const float* Bz [[buffer(7)]],
                          uint id [[thread_position_in_grid]]) {
    if (id >= params.particleCount) return;
    
    device ParticleState& p = particles[id];
    
    // electro-magnetic field on each particles
    float xh = p.x - 0.5;
    float yh = p.y - 0.5;
    int i1 = int(xh);
    int j1 = int(p.y);
    int i2 = int(p.x);
    int j2 = int(yh);
    float hv[2][2] ;
    hv[0][0] = xh  - float(i1);
    hv[0][1] = p.y - float(j1);
    hv[1][0] = p.x - float(i2);
    hv[1][1] = yh  - float(i2);
    
    // 5th-order weighting
    float sc;
    float sf[6][2][2];
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            sc = 2 + hv[i][j];
            sf[0][i][j] = 1/120 *pow(3-sc, 5);
            sc = 1 + hv[i][j];
            sf[1][i][j] = 1/120 *(51 +75*sc -210*pow(sc,2) +150*pow(sc,3) -45*pow(sc,4) +5*pow(sc,5));
            sc = hv[i][j];
            sf[2][i][j] = 1/60 *(33 -30*pow(sc,2) +15*pow(sc,4) -5*pow(sc,5));
            sc = 1 - hv[i][j];
            sf[3][i][j] = 1/60 *(33 -30*pow(sc,2) +15*pow(sc,4) -5*pow(sc,5));
            sc = 2 - hv[i][j];
            sf[4][i][j] = 1/120 *(51 +75*sc -210*pow(sc,2) +150*pow(sc,3) -45*pow(sc,4) +5*pow(sc,5));
            sc = 3 - hv[i][j];
            sf[5][i][j] = 1/120 *pow(3-sc, 5);
        }
    }

    // EMField for particle
    float Epx = 0.0;
    float Epy = 0.0;
    float Epz = 0.0;
    float Bpx = 0.0;
    float Bpy = 0.0;
    float Bpz = 0.0;
    int ii, jj;
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            ii = i1+(i-2);
            jj = j1+(j-2);
            Epx = Epx + sf[i][0][0]*sf[j][0][1]*Ex[ii+jj*params.ngy];
            ii = i2+(i-2);
            jj = j2+(j-2);
            Epy = Epy + sf[i][1][0]*sf[j][1][1]*Ey[ii+jj*params.ngy];
            ii = i1+(i-2);
            jj = j2+(j-2);
            Bpz = Bpz + sf[i][0][0]*sf[j][1][1]*Bz[ii+jj*params.ngy];
        }
    }
    
    // acceleration by electric field
    float umx = p.vx + params.constE*Epx;
    float umy = p.vy + params.constE*Epy;
    float umz = p.vz + params.constE*Epz;

    // preparing for rotation
    float btx = params.constB*Bpx;
    float bty = params.constB*Bpy;
    float btz = params.constB*Bpz;

    // vector product
    float v0x = umx + umy*btz - umz*bty;
    float v0y = umy + umz*btx - umx*btz;
    float v0z = umz + umx*bty - umy*btx;

    // rotation by magnetic field
    float btbt = btx*btx + bty*bty + btz*btz;
    float ssx = 2*btx/(1+btbt);
    float ssy = 2*bty/(1+btbt);
    float ssz = 2*btz/(1+btbt);
    float upx = umx + v0y*ssz - v0z*ssy;
    float upy = umy + v0z*ssx - v0x*ssz;
    float upz = umz + v0x*ssy - v0y*ssx;

    // acceleration by electric field
    p.vx = upx + params.constE*Epx;
    p.vy = upy + params.constE*Epy;
    p.vz = upz + params.constE*Epz;

    // updating position
    p.x = p.x + p.vx * params.constX;
    p.y = p.y + p.vy * params.constY;

    // periodic boundary
    p.x = fmod(p.x, float(params.ngx));
    p.y = fmod(p.y, float(params.ngy));
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
        _ngx = FieldParam.ngx;
        _ngy = FieldParam.ngy;
        _dx = FieldParam.dx;
        _dy = FieldParam.dy;
        
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
        // uniform distribution for position
        particles[i].x = (double)arc4random_uniform(Lx);
        particles[i].y = (double)arc4random_uniform(Ly);
        // shift from real coordinate to integer coodinate
        particles[i].x = particles[i].x/FieldParam.dx;
        particles[i].y = particles[i].y/FieldParam.dy;
        // Maxwellian for velocity 
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
    params->dt = dt;
    params->constE = 0.5*_charge*dt/_mass;
    params->constB = 0.5*_charge*dt/_mass/c;
    params->constX = dt/_dx;
    params->constY = dt/_dy;
    
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
    [computeEncoder setBuffer:_paramsBuffer offset:0 atIndex:1];
    [computeEncoder setBuffer:ExBuffer offset:0 atIndex:2];
    [computeEncoder setBuffer:EyBuffer offset:0 atIndex:3];
    [computeEncoder setBuffer:EzBuffer offset:0 atIndex:4];
    [computeEncoder setBuffer:BxBuffer offset:0 atIndex:5];
    [computeEncoder setBuffer:ByBuffer offset:0 atIndex:6];
    [computeEncoder setBuffer:BzBuffer offset:0 atIndex:7];
    
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

- (void)outputPhaseSpace:(int)i withParticleParam:(ParamForParticle)ParticleParam{
    // 粒子状態の内容にアクセス
    ParticleState *p = (ParticleState*)[_particleBuffer contents];
    // dx, dy の取得用にパラメータ取得
    SimulationParams *params = (SimulationParams*)_paramsBuffer.contents;

    // 位置を格納する変数
    float x,y;
    // 各粒子についてバイナリファイルに出力する 
    for (int idx = 0; idx < 20; idx++) {
        NSString *filePath = [NSString stringWithFormat:@"bin/PhaseSpace_%@_%d.bin", ParticleParam.pName, idx];
        // std::ofstream を利用してバイナリ書き出し（Objective-C++としてコンパイル）
        std::ofstream ofs([filePath UTF8String], std::ios::binary | std::ios::app);
        if (!ofs) {
            NSLog(@"Failed to open file: %@", filePath);
            continue;
        }

        // 位置を物理次元に戻す
        x = p[idx].x*_dx;
        y = p[idx].y*_dy;
        if (idx == 0) NSLog(@"x: %f", x);
        
        // phasespace を出力
        ofs.write(reinterpret_cast<const char*>(&i), sizeof(int));
        ofs.write(reinterpret_cast<const char*>(&x), sizeof(float));
        ofs.write(reinterpret_cast<const char*>(&y), sizeof(float));
        ofs.write(reinterpret_cast<const char*>(&p[idx].vx), sizeof(float));
        ofs.write(reinterpret_cast<const char*>(&p[idx].vy), sizeof(float));
        ofs.write(reinterpret_cast<const char*>(&p[idx].vz), sizeof(float));
        
        ofs.close();
    }
}

@end