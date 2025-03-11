#import "Particle.h"
#import "Init.h"
#import "Constant.h"
#include <random>
#import <string>

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
    uint particleCount;
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
};

// 粒子更新カーネル
kernel void updateParticles(
                        device ParticleState* particles     [[ buffer(0) ]],
                        constant SimulationParams& params   [[ buffer(1) ]],
                        device const float* Ex              [[ buffer(2) ]],
                        device const float* Ey              [[ buffer(3) ]],
                        device const float* Ez              [[ buffer(4) ]],
                        device const float* Bx              [[ buffer(5) ]],
                        device const float* By              [[ buffer(6) ]],
                        device const float* Bz              [[ buffer(7) ]],
                        device float* print                 [[ buffer(8) ]],
                        uint id                             [[ thread_position_in_grid ]]
                        ) {
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
    hv[1][1] = yh  - float(j2);
    
    // 5th-order weighting
    float sc;
    float sf[6][2][2];
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            sc = 2.0 + hv[i][j];
            sf[0][i][j] = 1.0/120.0 *pow(3.0-sc, 5);
            sc = 1.0 + hv[i][j];
            sf[1][i][j] = 1.0/120.0 *(51.0 +75.0*sc -210.0*pow(sc,2) +150.0*pow(sc,3) -45.0*pow(sc,4) +5.0*pow(sc,5));
            sc = hv[i][j];
            sf[2][i][j] = 1.0/60.0 *(33.0 -30.0*pow(sc,2) +15.0*pow(sc,4) -5.0*pow(sc,5));
            sc = 1.0 - hv[i][j];
            sf[3][i][j] = 1.0/60.0 *(33.0 -30.0*pow(sc,2) +15.0*pow(sc,4) -5.0*pow(sc,5));
            sc = 2.0 - hv[i][j];
            sf[4][i][j] = 1.0/120.0 *(51.0 +75.0*sc -210.0*pow(sc,2) +150.0*pow(sc,3) -45.0*pow(sc,4) +5.0*pow(sc,5));
            sc = 3.0 - hv[i][j];
            sf[5][i][j] = 1.0/120.0 *pow(3.0-sc, 5);
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
            Epx += sf[i][0][0]*sf[j][0][1]*Ex[ii+jj*params.ngy];
            ii = i2+(i-2);
            jj = j2+(j-2);
            Epy += sf[i][1][0]*sf[j][1][1]*Ey[ii+jj*params.ngy];
            ii = i1+(i-2);
            jj = j2+(j-2);
            Bpz += sf[i][0][0]*sf[j][1][1]*Bz[ii+jj*params.ngy];
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
    p.x = fmod(p.x + float(params.ngx) - 2.0, float(params.ngx)) + 2.0;
    p.y = fmod(p.y + float(params.ngy) - 2.0, float(params.ngy)) + 2.0;

    // debug print
    print[id] = Bpz;
}

// 電荷密度更新カーネル
kernel void integrateChargeDensity(
                        device ParticleState* particles             [[ buffer(0) ]],
                        constant SimulationParams& params           [[ buffer(1) ]],
                        device float* temp                          [[ buffer(2) ]],
                        device float* partial                       [[ buffer(3) ]],
                        device float* print                         [[ buffer(4) ]],
                        device int &arrSize                         [[ buffer(5) ]],
                        device int &chunkSize                       [[ buffer(6) ]],
                        device int &chunkOffset                     [[ buffer(7) ]],
                        device int &threadGroupSize                 [[ buffer(8) ]],
                        device float &constRho                      [[ buffer(9) ]],
                        uint gid                                    [[ thread_position_in_grid ]],
                        uint tid                                    [[ thread_index_in_threadgroup ]],
                        uint groupID                                [[ threadgroup_position_in_grid ]]
                        ) {
    // 粒子番号
    uint const pid = gid + chunkOffset;
    // dataCount を超えたら積分処理をスキップ
    if (pid >= params.particleCount) return;

    // initialize(各スレッドがアクセスし得る範囲を初期化)
    for (int i = 0; i < arrSize; i++){
        temp[gid + i*chunkSize] = 0.0f;
    }
    
    // 粒子を取得
    device ParticleState& p = particles[pid];

    // electro-magnetic field on each particles
    int i1 = int(p.x);
    int j1 = int(p.y);
    float hv[2];
    hv[0] = p.x - float(i1);
    hv[1] = p.y - float(j1);
    
    // 5th-order weighting
    float sc;
    float sf[6][2];
    for (int i = 0; i < 2; i++) {
        sc = 2.0 + hv[i];
        sf[0][i] = 1.0/120.0 *pow(3.0-sc, 5);
        sc = 1.0 + hv[i];
        sf[1][i] = 1.0/120.0 *(51.0 +75.0*sc -210.0*pow(sc,2) +150.0*pow(sc,3) -45.0*pow(sc,4) +5.0*pow(sc,5));
        sc = hv[i];
        sf[2][i] = 1.0/60.0 *(33.0 -30.0*pow(sc,2) +15.0*pow(sc,4) -5.0*pow(sc,5));
        sc = 1.0 - hv[i];
        sf[3][i] = 1.0/60.0 *(33.0 -30.0*pow(sc,2) +15.0*pow(sc,4) -5.0*pow(sc,5));
        sc = 2.0 - hv[i];
        sf[4][i] = 1.0/120.0 *(51.0 +75.0*sc -210.0*pow(sc,2) +150.0*pow(sc,3) -45.0*pow(sc,4) +5.0*pow(sc,5));
        sc = 3.0 - hv[i];
        sf[5][i] = 1.0/120.0 *pow(3.0-sc, 5);
    }
    
    // index
    int idx_in, idx_out, offset;

    // accumulation
    int ii, jj;
    const int ngx = (params.ngx + 2*2);
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            ii = i1+(i-2);
            jj = j1+(j-2);
            idx_out = gid + (ii+jj*ngx)*chunkSize;
            temp[idx_out] = sf[i][0]*sf[j][1]*constRho;
        }
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);

    // reduction
    for (int i = 0; i < arrSize; i++){
        for (uint stride = threadGroupSize / 2; stride > 0; stride /= 2) {
            if (tid < stride) {
                offset = groupID*threadGroupSize + i*chunkSize;
                idx_in = tid + stride + offset;
                idx_out = tid + offset;
                temp[idx_out] += temp[idx_in];
            }
            threadgroup_barrier(mem_flags::mem_threadgroup);
        }
    }
    
    // output partialSum
    if (tid == 0) {
        for (int i = 0; i < arrSize; i++){
            idx_in = 0 + groupID*threadGroupSize + i*chunkSize;
            idx_out = groupID + i*chunkSize/threadGroupSize;
            partial[idx_out] = temp[idx_in];
        }
    }
    
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
        _pNum = ParticleParam.pNum;
        _charge = ParticleParam.q;
        _mass = ParticleParam.m;
        _weight = ParticleParam.w;
        _ngx = FieldParam.ngx;
        _ngy = FieldParam.ngy;
        _dx = FieldParam.dx;
        _dy = FieldParam.dy;

        // 並列計算パラメータ
        _integrationChunkSize = 4096;
        _threadGroupSize = 256;
        
        // コンピュートパイプラインの設定
        NSError *error = nil;
        id<MTLLibrary> library = [device newLibraryWithSource:kMetalShaderSource
                                                    options:nil
                                                    error:&error];
        if (!library) {
            NSLog(@"Failed to create Metal library: %@", error);
            return nil;
        }
        // カーネルの取得
        id<MTLFunction> updateParticlesFunction = [library newFunctionWithName:@"updateParticles"];
        _updateParticlesPipeline = [device newComputePipelineStateWithFunction:updateParticlesFunction error:&error];
        id<MTLFunction> integrateChargeDensityFunction = [library newFunctionWithName:@"integrateChargeDensity"];
        _integrateChargeDensityPipeline = [device newComputePipelineStateWithFunction:integrateChargeDensityFunction error:&error];
        
        // バッファサイズ
        size_t buffSize;
        // 粒子バッファ
        buffSize = sizeof(ParticleState)*ParticleParam.pNum;
        _particleBuffer = [device newBufferWithLength:buffSize options:MTLResourceStorageModeShared];
        // 電磁場バッファ
        size_t EMbuffSize = sizeof(float)*(FieldParam.ngx+4)*(FieldParam.ngy+4);
        _ExBuffer = [device newBufferWithLength:EMbuffSize options:MTLResourceStorageModeShared];
        _EyBuffer = [device newBufferWithLength:EMbuffSize options:MTLResourceStorageModeShared];
        _EzBuffer = [device newBufferWithLength:EMbuffSize options:MTLResourceStorageModeShared];
        _BxBuffer = [device newBufferWithLength:EMbuffSize options:MTLResourceStorageModeShared];
        _ByBuffer = [device newBufferWithLength:EMbuffSize options:MTLResourceStorageModeShared];
        _BzBuffer = [device newBufferWithLength:EMbuffSize options:MTLResourceStorageModeShared];
        // 電荷密度用バッファ
        buffSize = EMbuffSize*_integrationChunkSize;
        _integrateTemporaryBuffer = [device newBufferWithLength:buffSize options:MTLResourceStorageModeShared];
        buffSize = EMbuffSize*(_integrationChunkSize/_threadGroupSize);
        _integratePartialBuffer = [device newBufferWithLength:buffSize options:MTLResourceStorageModeShared];
        // シミュレーションパラメータバッファ
        buffSize = sizeof(SimulationParams);
        _paramsBuffer = [device newBufferWithLength:buffSize options:MTLResourceStorageModeShared];
        // デバッグ出力用バッファ
        buffSize = sizeof(float)*ParticleParam.pNum;
        _printBuffer = [device newBufferWithLength:buffSize options:MTLResourceStorageModeShared];

        // 初期粒子分布の設定
        [self generateParticles:ParticleParam withFieldParam:FieldParam];

    }
    return self;
}

- (void)generateParticles:(ParamForParticle)ParticleParam
            withFieldParam:(ParamForField)FieldParam{
    // 乱数生成器
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    float Lx = FieldParam.dx * FieldParam.ngx;
    float Ly = FieldParam.dy * FieldParam.ngy;
    std::uniform_real_distribution<> unif_dist(0.0, 1.0);
    float vth_epi  = sqrt(2*kb*ParticleParam.initT/ParticleParam.m);
    std::normal_distribution<> norm_dist(0.0, vth_epi);
    // 粒子の初期化
    ParticleState *particles = (ParticleState *)self.particleBuffer.contents;
    for (int i = 0; i < ParticleParam.pNum; i++) {
        if ([ParticleParam.GenerateType isEqualToString:@"UniformGaussian"]){
            // uniform distribution for position
            particles[i].x = (float)unif_dist(engine)*Lx;
            particles[i].y = (float)unif_dist(engine)*Ly;
            // Maxwellian for velocity 
            particles[i].vx = (float)norm_dist(engine);
            particles[i].vy = (float)norm_dist(engine);
            particles[i].vz = (float)norm_dist(engine);
        } else if ([ParticleParam.GenerateType isEqualToString:@"UniformConstant"]){
            // uniform distribution for position
            particles[i].x = (float)unif_dist(engine)*Lx;
            particles[i].y = (float)unif_dist(engine)*Ly;
            // constant for velosity
            particles[i].vx = (float)ParticleParam.initU[0];
            particles[i].vy = (float)ParticleParam.initU[1];
            particles[i].vz = (float)ParticleParam.initU[2];
        }
        // shift from real coordinate to integer coodinate
        particles[i].x = particles[i].x/FieldParam.dx;
        particles[i].y = particles[i].y/FieldParam.dy;
        // shift origin for high-order weighting
        if (FieldParam.weightOrder == 5){
            particles[i].x = particles[i].x + 2.0;
            particles[i].y = particles[i].y + 2.0;
        }
        // check
        // NSLog(@"initial x[%d]: %f", i, particles[i].x);
    }
}


// 時間更新
- (void)update:(double)dt withEMField:(EMField*)fld {
    // シミュレーションパラメータの更新
    SimulationParams* params = (SimulationParams*)_paramsBuffer.contents;
    params->particleCount = _pNum;
    params->dt = dt;
    params->constE = 0.5*_charge*dt/_mass;
    params->constB = 0.5*_charge*dt/_mass/c;
    params->constX = dt/_dx;
    params->constY = dt/_dy;
    params->ngx = _ngx;
    params->ngy = _ngy;
    
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
    [computeEncoder setBuffer:_particleBuffer   offset:0 atIndex:0];
    [computeEncoder setBuffer:_paramsBuffer     offset:0 atIndex:1];
    [computeEncoder setBuffer:ExBuffer          offset:0 atIndex:2];
    [computeEncoder setBuffer:EyBuffer          offset:0 atIndex:3];
    [computeEncoder setBuffer:EzBuffer          offset:0 atIndex:4];
    [computeEncoder setBuffer:BxBuffer          offset:0 atIndex:5];
    [computeEncoder setBuffer:ByBuffer          offset:0 atIndex:6];
    [computeEncoder setBuffer:BzBuffer          offset:0 atIndex:7];
    [computeEncoder setBuffer:_printBuffer      offset:0 atIndex:8];
    
    // グリッドとスレッドグループのサイズ設定
    uint gridSize = params->particleCount;
    uint threadGroupNum = (gridSize + _threadGroupSize - 1) / _threadGroupSize;
    MTLSize threadGroupSize = MTLSizeMake(_threadGroupSize, 1, 1);
    MTLSize gridSizeMetalStyle = MTLSizeMake(threadGroupNum, 1, 1);
    
    // ディスパッチ
    [computeEncoder dispatchThreadgroups:gridSizeMetalStyle
                            threadsPerThreadgroup:threadGroupSize];
    
    // 粒子状態取得
    ParticleState* p = (ParticleState*)[_particleBuffer contents];
    for (int idx = 0; idx < 1; idx++){
        NSLog(@"before update: p.x[%d]: %f", idx, p[idx].x);
    }

    // エンコーディングと実行
    [computeEncoder endEncoding];
    [commandBuffer commit];
    [commandBuffer waitUntilCompleted];

    // 粒子状態取得
    p = (ParticleState*)[_particleBuffer contents];
    for (int idx = 0; idx < 1; idx++){
        NSLog(@"after update: p.x[%d]: %f", idx, p[idx].x);
    }
    // デバッグ出力
    float* prt = (float*)_printBuffer.contents;
    for (int idx = 0; idx < 1; idx++){
        NSLog(@"debug print: val[%d]: %f", idx, prt[idx]);
    }
}

- (void)integrateChargeDensity:(EMField*)fld{
    // 粒子情報取得
    SimulationParams* params = (SimulationParams*)_paramsBuffer.contents;
    // 格子情報取得
    int ng = (fld.ngx+2*fld.ngb)*(fld.ngy+2*fld.ngb);
    // constant 引数バッファ
    id<MTLBuffer> chunkSizeBuffer = [_device newBufferWithBytes:&_integrationChunkSize
                                                length:sizeof(int)
                                                options:MTLResourceStorageModeShared];
    id<MTLBuffer> arrSizeBuffer = [_device newBufferWithBytes:&ng
                                                length:sizeof(int)
                                                options:MTLResourceStorageModeShared];
    id<MTLBuffer> threadGroupSizeBuffer = [_device newBufferWithBytes:&_threadGroupSize
                                                length:sizeof(int)
                                                options:MTLResourceStorageModeShared];
    // constant for rho
    float constRho = _charge * _weight / (_dx * _dy);
    id<MTLBuffer> constRhoBuffer = [_device newBufferWithBytes:&constRho
                                                length:sizeof(int)
                                                options:MTLResourceStorageModeShared];
    // 分割して積分
    uint chunkNum = (params->particleCount + _integrationChunkSize - 1) / _integrationChunkSize;
    for (int chunk = 0; chunk < chunkNum; chunk++){
        // インデックスのオフセット計算
        uint chunkOffset = chunk * _integrationChunkSize;
        id<MTLBuffer> chunkOffsetBuffer = [_device newBufferWithBytes:&chunkOffset
                                                length:sizeof(uint)
                                                options:MTLResourceStorageModeShared];

        // コマンドバッファとエンコーダの作成
        id<MTLCommandBuffer> commandBuffer = [_commandQueue commandBuffer];
        id<MTLComputeCommandEncoder> computeEncoder = [commandBuffer computeCommandEncoder];
        
        // パイプラインとバッファの設定
        [computeEncoder setComputePipelineState:_integrateChargeDensityPipeline];
        [computeEncoder setBuffer:_particleBuffer           offset:0 atIndex:0];
        [computeEncoder setBuffer:_paramsBuffer             offset:0 atIndex:1];
        [computeEncoder setBuffer:_integrateTemporaryBuffer offset:0 atIndex:2];
        [computeEncoder setBuffer:_integratePartialBuffer   offset:0 atIndex:3];
        [computeEncoder setBuffer:_printBuffer              offset:0 atIndex:4];
        [computeEncoder setBuffer:arrSizeBuffer             offset:0 atIndex:5];
        [computeEncoder setBuffer:chunkSizeBuffer           offset:0 atIndex:6];
        [computeEncoder setBuffer:chunkOffsetBuffer         offset:0 atIndex:7];
        [computeEncoder setBuffer:threadGroupSizeBuffer     offset:0 atIndex:8];
        [computeEncoder setBuffer:constRhoBuffer            offset:0 atIndex:9];

        // グリッドとスレッドグループのサイズ設定
        uint threadGroupNum = _integrationChunkSize/_threadGroupSize;
        MTLSize gridSizeMetalStyle = MTLSizeMake(threadGroupNum, 1, 1);
        MTLSize threadGroupSize = MTLSizeMake(_threadGroupSize, 1, 1);

        // ディスパッチ
        [computeEncoder dispatchThreadgroups:gridSizeMetalStyle
                                threadsPerThreadgroup:threadGroupSize];

        // エンコーディングと実行
        [computeEncoder endEncoding];
        [commandBuffer commit];
        [commandBuffer waitUntilCompleted];

        // スレッドグループごとの部分和を加算
        float* partialSums = (float*)_integratePartialBuffer.contents;
        for (int i = 0; i < ng; i++){
            for (int j = 0; j < threadGroupNum; j++){
                fld.rho[i] += partialSums[j + i*threadGroupNum];
            }
        }
        // check
        int i,j,idx;
        i = 10;
        j = 10;
        idx = i + j *(fld.ngx+2*fld.ngb);
        // NSLog(@"integration: chunk: %d, fld.rho[%d,%d]: %f", chunk, i, j, fld.rho[idx]);
    }
    // 加算後
    int i,j,idx;
    i = 10;
    j = 10;
    idx = i + j *(fld.ngy+fld.ngb);
    NSLog(@"after integration: fld.rho[%d,%d]: %f", i, j, fld.rho[idx]);
};

- (void)outputPhaseSpace:(int)cycle withParticleParam:(ParamForParticle)ParticleParam{
    // 粒子状態の内容にアクセス
    ParticleState *p = (ParticleState*)[_particleBuffer contents];
    // dx, dy の取得用にパラメータ取得
    SimulationParams *params = (SimulationParams*)_paramsBuffer.contents;

    // 位置を格納する変数
    float x,y;
    // 各粒子についてバイナリファイルに出力する 
    for (int idx = 0; idx < 20; idx++) {
        NSString *filePath = [NSString stringWithFormat:@"bin/PhaseSpace_%@_%d.bin", ParticleParam.pName, idx];
        // バイナリ書き出し
        std::ofstream ofs([filePath UTF8String], std::ios::binary | std::ios::app);
        if (!ofs) {
            NSLog(@"Failed to open file: %@", filePath);
            continue;
        }

        // 位置を物理次元に戻す
        x = (p[idx].x - 2.0)*_dx;
        y = (p[idx].y - 2.0)*_dy;
        
        // phasespace を出力
        ofs.write(reinterpret_cast<const char*>(&cycle), sizeof(int));
        ofs.write(reinterpret_cast<const char*>(&x), sizeof(float));
        ofs.write(reinterpret_cast<const char*>(&y), sizeof(float));
        ofs.write(reinterpret_cast<const char*>(&p[idx].vx), sizeof(float));
        ofs.write(reinterpret_cast<const char*>(&p[idx].vy), sizeof(float));
        ofs.write(reinterpret_cast<const char*>(&p[idx].vz), sizeof(float));
        
        ofs.close();
    }
}

@end
