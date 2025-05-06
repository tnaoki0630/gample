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
    int piflag;
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
    uint pNum;
    float constE;
    float constB;
    float constX;
    float constY;
    int ngx;
    int ngy;
    int ngb;
    int BC_Xmin;
    int BC_Xmax;
    int BC_Ymin;
    int BC_Ymax;
};

// 粒子更新カーネル
kernel void updateParticles(
                        device ParticleState* ptcl          [[ buffer(0) ]],
                        constant SimulationParams& prm      [[ buffer(1) ]],
                        device const float* Ex              [[ buffer(2) ]],
                        device const float* Ey              [[ buffer(3) ]],
                        device const float* Ez              [[ buffer(4) ]],
                        device const float* Bx              [[ buffer(5) ]],
                        device const float* By              [[ buffer(6) ]],
                        device const float* Bz              [[ buffer(7) ]],
                        device float* print                 [[ buffer(8) ]],
                        uint id                             [[ thread_position_in_grid ]]
                        ) {
    if (id >= prm.pNum) return;
    
    device ParticleState& p = ptcl[id];

    // electro-magnetic field on each ptcl
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
            Epx += sf[i][0][0]*sf[j][0][1]*Ex[ii+jj*prm.ngy];
            ii = i2+(i-2);
            jj = j2+(j-2);
            Epy += sf[i][1][0]*sf[j][1][1]*Ey[ii+jj*prm.ngy];
            ii = i1+(i-2);
            jj = j2+(j-2);
            Bpz += sf[i][0][0]*sf[j][1][1]*Bz[ii+jj*prm.ngy];
        }
    }
    
    // acceleration by electric field
    float umx = p.vx + prm.constE*Epx;
    float umy = p.vy + prm.constE*Epy;
    float umz = p.vz + prm.constE*Epz;

    // preparing for rotation
    float btx = prm.constB*Bpx;
    float bty = prm.constB*Bpy;
    float btz = prm.constB*Bpz;

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
    p.vx = upx + prm.constE*Epx;
    p.vy = upy + prm.constE*Epy;
    p.vz = upz + prm.constE*Epz;

    // updating position
    p.x = p.x + p.vx * prm.constX;
    p.y = p.y + p.vy * prm.constY;

    // periodic boundary
    if (prm.BC_Xmin == 0) {
        p.x = fmod(p.x + float(prm.ngx - prm.ngb), float(prm.ngx)) + float(prm.ngb);
    }
    if (prm.BC_Ymin == 0) {
        p.y = fmod(p.y + float(prm.ngy - prm.ngb), float(prm.ngy)) + float(prm.ngb);
    }

    // delete boundary
    if (prm.BC_Xmin == 1 && p.x - prm.ngb < 0.0f){
        p.piflag = 1;
        p.x += 1.0f; // for check
    }
    if (prm.BC_Xmax == 1 && p.x - prm.ngb > float(prm.ngx)){
        p.piflag = 1;
        p.x -= 1.0f; // for check
    }
    if (prm.BC_Ymin == 1 && p.y - prm.ngb < 0.0f){
        p.piflag = 1;
        p.y += 1.0f; // for check
    }
    if (prm.BC_Ymax == 1 && p.y - prm.ngb > float(prm.ngy)){
        p.piflag = 1;
        p.y -= 1.0f; // for check
    }

    // debug print
    print[id] = p.piflag;
}

// 電荷密度更新カーネル
kernel void integrateChargeDensity(
                        device ParticleState* ptcl          [[ buffer(0) ]],
                        constant SimulationParams& prm      [[ buffer(1) ]],
                        device float* temp                  [[ buffer(2) ]],
                        device float* partial               [[ buffer(3) ]],
                        device float* print                 [[ buffer(4) ]],
                        device int &arrSize                 [[ buffer(5) ]],
                        device int &chunkSize               [[ buffer(6) ]],
                        device int &chunkOffset             [[ buffer(7) ]],
                        device int &threadGroupSize         [[ buffer(8) ]],
                        device float &constRho              [[ buffer(9) ]],
                        uint gid                            [[ thread_position_in_grid ]],
                        uint tid                            [[ thread_index_in_threadgroup ]],
                        uint groupID                        [[ threadgroup_position_in_grid ]]
                        ) {
    // 粒子番号
    uint const pid = gid + chunkOffset;
    // dataCount を超えたら積分処理をスキップ
    if (pid >= prm.pNum) return;

    // initialize(各スレッドがアクセスし得る範囲を初期化)
    for (int i = 0; i < arrSize; i++){
        temp[gid + i*chunkSize] = 0.0f;
    }
    
    // 粒子を取得
    device ParticleState& p = ptcl[pid];

    // electro-magnetic field on each ptcl
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
    const int nx = (prm.ngx + 2*prm.ngb);
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            ii = i1+(i-prm.ngb);
            jj = j1+(j-prm.ngb);
            idx_out = gid + (ii+jj*(nx+1))*chunkSize;
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

// プライベートインスタンス変数
@implementation Particle {
    NSString* _pName;
    float _m;
    float _q;
    float _w;
    id<MTLBuffer> _ExBuffer;
    id<MTLBuffer> _EyBuffer;
    id<MTLBuffer> _EzBuffer;
    id<MTLBuffer> _BxBuffer;
    id<MTLBuffer> _ByBuffer;
    id<MTLBuffer> _BzBuffer;
}

// 初期設定
- (instancetype)initWithDevice:(id<MTLDevice>)device withParam:(Init*)initParam specimen:(int)s{
    self = [super init];
    if (self) {

        // パラメータ取得
        NSArray *ParticleParams = [initParam getParamForParticle];
        NSValue *value = ParticleParams[s];
        struct ParamForParticle particleParam;
        [value getValue:&particleParam];
        struct ParamForField fieldParam = [initParam getParamForField];
        
        // 粒子パラメータ
        _pName = particleParam.pName;
        _m = particleParam.m;
        _q = particleParam.q;
        _w = particleParam.w;

        // 並列計算パラメータ
        _integrationChunkSize = 4096;
        _threadGroupSize = 256;
        
        // コマンドキューの作成
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
        // カーネルの取得
        id<MTLFunction> updateParticlesFunction = [library newFunctionWithName:@"updateParticles"];
        _updateParticlesPipeline = [device newComputePipelineStateWithFunction:updateParticlesFunction error:&error];
        id<MTLFunction> integrateChargeDensityFunction = [library newFunctionWithName:@"integrateChargeDensity"];
        _integrateChargeDensityPipeline = [device newComputePipelineStateWithFunction:integrateChargeDensityFunction error:&error];
        
        // バッファサイズ
        size_t buffSize;
        // シミュレーションパラメータバッファ
        buffSize = sizeof(SimulationParams);
        _paramsBuffer = [device newBufferWithLength:buffSize options:MTLResourceStorageModeShared];
        SimulationParams* prm = (SimulationParams*)_paramsBuffer.contents;        
        // パラメータ格納
        prm->pNum = particleParam.pNum;
        prm->ngx = fieldParam.ngx;
        prm->ngy = fieldParam.ngy;
        prm->ngb = fieldParam.ngb;
        // 境界条件格納
        NSArray* particleBCs = [initParam getParticleBoundaries];
        for (int pos = 0; pos < particleBCs.count; pos++){
            NSValue* value = particleBCs[pos];
            struct BoundaryConditionForParticle BC;
            [value getValue:&BC];
            if ([BC.position isEqualToString:@"Xmin"]){
                if ([BC.type isEqualToString:@"periodic"]){
                    prm->BC_Xmin = 0;
                }else if ([BC.type isEqualToString:@"Delete"]){
                    prm->BC_Xmin = 1;
                }else{
                    prm->BC_Xmin = -1; // error
                }
            }else if ([BC.position isEqualToString:@"Xmax"]){
                if ([BC.type isEqualToString:@"periodic"]){
                    prm->BC_Xmax = 0;
                }else if ([BC.type isEqualToString:@"Delete"]){
                    prm->BC_Xmax = 1;
                }else{
                    prm->BC_Xmax = -1; // error
                }
            }else if ([BC.position isEqualToString:@"Ymin"]){
                if ([BC.type isEqualToString:@"periodic"]){
                    prm->BC_Ymin = 0;
                }else if ([BC.type isEqualToString:@"Delete"]){
                    prm->BC_Ymin = 1;
                }else{
                    prm->BC_Ymin = -1; // error
                }
            }else if ([BC.position isEqualToString:@"Ymax"]){
                if ([BC.type isEqualToString:@"periodic"]){
                    prm->BC_Ymax = 0;
                }else if ([BC.type isEqualToString:@"Delete"]){
                    prm->BC_Ymax = 1;
                }else{
                    prm->BC_Ymax = -1; // error
                }
            }
        }

        // 粒子バッファ
        buffSize = sizeof(ParticleState)*particleParam.pNum;
        _particleBuffer = [device newBufferWithLength:buffSize options:MTLResourceStorageModeShared];
        // 初期粒子分布の設定
        [self generateParticles:particleParam withFieldParam:fieldParam];

        // 電磁場バッファ
        int nx = fieldParam.ngx + 2*fieldParam.ngb;
        int ny = fieldParam.ngy + 2*fieldParam.ngb;
        size_t EMbuffSize = sizeof(float)*(nx+1)*(ny+1);
        
        // 電荷密度用バッファ
        buffSize = EMbuffSize*_integrationChunkSize;
        _integrateTemporaryBuffer = [device newBufferWithLength:buffSize options:MTLResourceStorageModeShared];
        buffSize = EMbuffSize*(_integrationChunkSize/_threadGroupSize);
        _integratePartialBuffer = [device newBufferWithLength:buffSize options:MTLResourceStorageModeShared];

        // デバッグ出力用バッファ
        buffSize = sizeof(float)*particleParam.pNum;
        _printBuffer = [device newBufferWithLength:buffSize options:MTLResourceStorageModeShared];

    }
    return self;
}

- (void)generateParticles:(ParamForParticle)particleParam
            withFieldParam:(ParamForField)fieldParam{
    // 乱数生成器
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    float Lx = fieldParam.dx * fieldParam.ngx;
    float Ly = fieldParam.dy * fieldParam.ngy;
    std::uniform_real_distribution<> unif_dist(0.0, 1.0);
    float vth_epi  = sqrt(2*kb*particleParam.initT/particleParam.m);
    std::normal_distribution<> norm_dist(0.0, vth_epi);
    // 粒子の初期化
    ParticleState *ptcl = (ParticleState *)self.particleBuffer.contents;
    for (int i = 0; i < particleParam.pNum; i++) {
        if ([particleParam.GenerateType isEqualToString:@"UniformGaussian"]){
            // uniform distribution for position
            ptcl[i].x = (float)unif_dist(engine)*Lx;
            ptcl[i].y = (float)unif_dist(engine)*Ly;
            // Maxwellian for velocity 
            ptcl[i].vx = (float)norm_dist(engine);
            ptcl[i].vy = (float)norm_dist(engine);
            ptcl[i].vz = (float)norm_dist(engine);
        } else if ([particleParam.GenerateType isEqualToString:@"UniformConstant"]){
            // uniform distribution for position
            ptcl[i].x = (float)unif_dist(engine)*Lx;
            ptcl[i].y = (float)unif_dist(engine)*Ly;
            // constant for velosity
            ptcl[i].vx = (float)particleParam.initU[0];
            ptcl[i].vy = (float)particleParam.initU[1];
            ptcl[i].vz = (float)particleParam.initU[2];
        }
        // shift from real coordinate to integer coodinate
        ptcl[i].x = ptcl[i].x/fieldParam.dx;
        ptcl[i].y = ptcl[i].y/fieldParam.dy;
        // shift origin for high-order weighting
        if (fieldParam.weightOrder == 5){
            ptcl[i].x = ptcl[i].x + 2.0;
            ptcl[i].y = ptcl[i].y + 2.0;
        }
        // check
        // NSLog(@"initial x[%d]: %f", i, ptcl[i].x);
    }
}


// 時間更新
- (void)update:(double)dt withEMField:(EMField*)fld {
    // シミュレーションパラメータの更新
    SimulationParams* prm = (SimulationParams*)_paramsBuffer.contents;
    prm->constE = 0.5*_q*dt/_m;
    prm->constB = 0.5*_q*dt/_m/c;
    prm->constX = dt/fld.dx;
    prm->constY = dt/fld.dy;
    
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
    uint gridSize = prm->pNum;
    uint threadGroupNum = (gridSize + _threadGroupSize - 1) / _threadGroupSize;
    MTLSize threadGroupSize = MTLSizeMake(_threadGroupSize, 1, 1);
    MTLSize gridSizeMetalStyle = MTLSizeMake(threadGroupNum, 1, 1);
    
    // ディスパッチ
    [computeEncoder dispatchThreadgroups:gridSizeMetalStyle
                            threadsPerThreadgroup:threadGroupSize];
    
    // 粒子状態取得
    // ParticleState* p = (ParticleState*)[_particleBuffer contents];
    // for (int idx = 0; idx < 1; idx++){
    //     NSLog(@"before update: p.x[%d]: %f", idx, p[idx].x);
    // }

    // エンコーディングと実行
    [computeEncoder endEncoding];
    [commandBuffer commit];
    [commandBuffer waitUntilCompleted];

    //// 粒子状態取得
    ParticleState* p = (ParticleState*)[_particleBuffer contents];
    for (int idx = 0; idx < prm->pNum; idx++){
        if (p[idx].piflag == 1){
            NSLog(@"flowout: idx = %d, p.x[idx] = %f, p.y[idx] = %f, p.vx[idx] = %f", idx, p[idx].x, p[idx].y, p[idx].vx);
            p[idx].piflag = 0;
        }
    }
    // // デバッグ出力
    // float* prt = (float*)_printBuffer.contents;
    // for (int idx = 0; idx < 10; idx++){
    //     NSLog(@"debug print: val[%d]: %f", idx, prt[idx]);
    // }
}

- (void)integrateChargeDensity:(EMField*)fld{
    // prepare
    SimulationParams* prm = (SimulationParams*)_paramsBuffer.contents;
    int ng = (fld.nx+1)*(fld.ny+1);
    id<MTLBuffer> rhoBuffer = [fld rhoBuffer];
    float* rho = (float*)rhoBuffer.contents;

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
    // constant
    float constRho = _q * _w * _m / (fld.dx * fld.dx);
    id<MTLBuffer> constRhoBuffer = [_device newBufferWithBytes:&constRho
                                                length:sizeof(int)
                                                options:MTLResourceStorageModeShared];
    // 分割して積分
    uint chunkNum = (prm->pNum + _integrationChunkSize - 1) / _integrationChunkSize;
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
                rho[i] += partialSums[j + i*threadGroupNum];
            }
        }
    }
    // 加算後
    // int i,j,idx;
    // i = 10;
    // j = 10;
    // idx = i + j *(fld.ngy+fld.ngb);
    // NSLog(@"after integration: rho[%d,%d]: %f", i, j, fld.rho[idx]);
};

- (void)outputPhaseSpace:(int)cycle withEMField:(EMField*)fld{
    // prepare
    ParticleState *p = (ParticleState*)[_particleBuffer contents];
    SimulationParams *prm = (SimulationParams*)_paramsBuffer.contents;

    // 位置を格納する変数
    float x,y;
    // 各粒子についてバイナリファイルに出力する 
    for (int idx = 0; idx < 20; idx++) {
        NSString *filePath = [NSString stringWithFormat:@"bin/PhaseSpace_%@_%d.bin", _pName, idx];
        // バイナリ書き出し
        std::ofstream ofs([filePath UTF8String], std::ios::binary | std::ios::app);
        if (!ofs) {
            NSLog(@"Failed to open file: %@", filePath);
            continue;
        }

        // 位置を物理次元に戻す
        x = (p[idx].x - float(prm->ngb))*fld.dx;
        y = (p[idx].y - float(prm->ngb))*fld.dy;
        
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
