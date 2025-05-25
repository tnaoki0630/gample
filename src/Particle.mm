#import "Particle.h"
#import "Init.h"
#import "Constant.h"
#include <random>
#import <string>

# define PI 3.141592653589793

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
    const int nx = (prm.ngx + 2*prm.ngb);
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            ii = i1+(i-prm.ngb);
            jj = j1+(j-prm.ngb);
            Epx += sf[i][0][0]*sf[j][0][1]*Ex[ii+jj*(nx+1)];
            ii = i2+(i-prm.ngb);
            jj = j2+(j-prm.ngb);
            Epy += sf[i][1][0]*sf[j][1][1]*Ey[ii+jj*(nx+1)];
            ii = i1+(i-prm.ngb);
            jj = j2+(j-prm.ngb);
            Bpz += sf[i][0][0]*sf[j][1][1]*Bz[ii+jj*(nx+1)];
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
    }
    if (prm.BC_Xmax == 1 && p.x - prm.ngb > float(prm.ngx)){
        p.piflag = 2;
    }
    if (prm.BC_Ymin == 1 && p.y - prm.ngb < 0.0f){
        p.piflag = 3;
    }
    if (prm.BC_Ymax == 1 && p.y - prm.ngb > float(prm.ngy)){
        p.piflag = 4;
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
                        device int &pNumPerThread           [[ buffer(7) ]],
                        device int &threadGroupSize         [[ buffer(8) ]],
                        device float &constRho              [[ buffer(9) ]],
                        uint gid                            [[ thread_position_in_grid ]],
                        uint tid                            [[ thread_index_in_threadgroup ]],
                        uint groupID                        [[ threadgroup_position_in_grid ]]
                        ) {
    // initialize(各スレッドがアクセスし得る範囲を初期化)
    for (int i = 0; i < arrSize; i++){
        temp[gid + i*chunkSize] = 0.0f;
    }

    // 積分ループ
    for (int i = 0; i < pNumPerThread; i++){
        
        uint pid = i + pNumPerThread*gid;
        if (pid > prm.pNum) {
            // pNum を超えたら積分処理をスキップ
            break;
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

        // accumulation
        int ii, jj;
        const int nx = (prm.ngx + 2*prm.ngb);
        for (int i = 0; i < 6; i++) {
            for (int j = 0; j < 6; j++) {
                ii = i1+(i-prm.ngb);
                jj = j1+(j-prm.ngb);
                int idx_out = gid + (ii+jj*(nx+1))*chunkSize;
                temp[idx_out] += sf[i][0]*sf[j][1]*constRho;
            }
        }
    
    }

    // スレッドグループ内で同期
    threadgroup_barrier(mem_flags::mem_threadgroup);

    // reduction
    for (int i = 0; i < arrSize; i++){
        for (uint stride = threadGroupSize / 2; stride > 0; stride /= 2) {
            if (tid < stride) {
                uint offset = groupID*threadGroupSize + i*chunkSize;
                int idx_in = tid + stride + offset;
                int idx_out = tid + offset;
                temp[idx_out] += temp[idx_in];
            }
            threadgroup_barrier(mem_flags::mem_threadgroup);
        }
    }
    
    // output partialSum
    if (tid == 0) {
        for (int i = 0; i < arrSize; i++){
            int idx_in = 0 + groupID*threadGroupSize + i*chunkSize;
            int idx_out = groupID + i*chunkSize/threadGroupSize;
            partial[idx_out] = temp[idx_in];
        }
    }
    
}
)";

// プライベートインスタンス変数
@implementation Particle {
    NSString* _pName;
    int _pNumMax;
    int _pinum_Xmin;
    int _pinum_Xmax;
    int _pinum_Ymin;
    int _pinum_Ymax;
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
- (instancetype)initWithDevice:(id<MTLDevice>)device withParam:(Init*)initParam specimen:(int)s withLogger:(XmlLogger&)logger{
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
        _pNumMax = particleParam.pNumMax;
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
        SimulationParams* prm = (SimulationParams*)[_paramsBuffer contents];
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
        if(particleParam.pNum > particleParam.pNumMax){
            NSLog(@"pNum > pNumMax: pName = %@, pNum = %dd, pNumMax = %d", _pName, particleParam.pNum, particleParam.pNumMax);
            return nil;
        }
        buffSize = sizeof(ParticleState)*particleParam.pNumMax; // maxで初期化
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
        buffSize = sizeof(float)*particleParam.pNumMax;
        _printBuffer = [device newBufferWithLength:buffSize options:MTLResourceStorageModeShared];

    }
    return self;
}

- (void)generateParticles:(ParamForParticle)particleParam withFieldParam:(ParamForField)fieldParam{
    // 範囲指定
    float Xmin, Lx;
    float Ymin, Ly;
    if (particleParam.genX[0] < 0){
        Lx = fieldParam.dx * fieldParam.ngx;
        Xmin = 0.0;
    }else{
        Lx = particleParam.genX[1] - particleParam.genX[0];
        Xmin = particleParam.genX[0];
    }
    if (particleParam.genY[0] < 0){
        Ly = fieldParam.dy * fieldParam.ngy;
        Ymin = 0.0;
    }else{
        Ly = particleParam.genY[1] - particleParam.genY[0];
        Ymin = particleParam.genY[0];
    }
    // 乱数生成器
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_real_distribution<> unif_dist(0.0, 1.0);
    float vth_epi  = sqrt(2*kb*particleParam.genT/particleParam.m);
    std::normal_distribution<> norm_dist(0.0, vth_epi);
    // 粒子の初期化
    ParticleState *p = (ParticleState *)self.particleBuffer.contents;
    for (int i = 0; i < particleParam.pNum; i++) {
        if ([particleParam.genType isEqualToString:@"uniform-Gaussian"]){
            // uniform distribution for position
            p[i].x = Xmin + (float)unif_dist(engine)*Lx;
            p[i].y = Ymin + (float)unif_dist(engine)*Ly;
            // Maxwellian for velocity 
            p[i].vx = (float)particleParam.genU[0] + (float)norm_dist(engine);
            p[i].vy = (float)particleParam.genU[1] + (float)norm_dist(engine);
            p[i].vz = (float)particleParam.genU[2] + (float)norm_dist(engine);
        }
        // shift from real coordinate to integer coodinate
        p[i].x = p[i].x/fieldParam.dx;
        p[i].y = p[i].y/fieldParam.dy;
        // shift origin for high-order weighting
        p[i].x = p[i].x + (float)fieldParam.ngb;
        p[i].y = p[i].y + (float)fieldParam.ngb;
        // deletion flag
        p[i].piflag = 0;
    }
}

// 時間更新
- (void)update:(double)dt withEMField:(EMField*)fld  withLogger:(XmlLogger&)logger{
    // シミュレーションパラメータの更新
    SimulationParams* prm = (SimulationParams*)[_paramsBuffer contents];
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

    // エンコーディングと実行
    [computeEncoder endEncoding];
    [commandBuffer commit];
    [commandBuffer waitUntilCompleted];

    // 粒子状態取得
    // ParticleState* p = (ParticleState*)[_particleBuffer contents];
    // for (int idx = 0; idx < 1; idx++){
    //     NSLog(@"before update: p.x[%d]: %f", idx, p[idx].x);
    // }

    // // デバッグ出力
    // float* prt = (float*)_printBuffer.contents;
    // for (int idx = 0; idx < 10; idx++){
    //     NSLog(@"debug print: val[%d]: %f", idx, prt[idx]);
    // }

}

- (void)reduce:(XmlLogger&)logger{
    // obtain objects
    SimulationParams* prm = (SimulationParams*)[_paramsBuffer contents];
    ParticleState* p = (ParticleState*)[_particleBuffer contents];

    // initialize
    _pinum_Xmin = 0;
    _pinum_Xmax = 0;
    _pinum_Ymin = 0;
    _pinum_Ymax = 0;

    // reduce
    int kmax, kend, pulln;
    kmax = prm->pNum;
    pulln = 0;
    for (int k1 = 0; k1 < prm->pNum; k1++){
        // reached max
        if (k1 == kmax){
            if (p[k1].piflag == 1){
                pulln++;
            }
            break;
        }
        // delete
        if (p[k1].piflag > 0){
            // count particles for each boundary with charge
            if (p[k1].piflag == 1){
                _pinum_Xmin += (int)(_q/(float)ec);
            }else if(p[k1].piflag == 2){
                _pinum_Xmax += (int)(_q/(float)ec);
            }else if(p[k1].piflag == 3){
                _pinum_Ymin += (int)(_q/(float)ec);
            }else if(p[k1].piflag == 4){
                _pinum_Ymax += (int)(_q/(float)ec);
            }
            // keep kmax
            kend = kmax;
            // scan alive particle
            for (int k2 = kend; k2 > k1; k2--){
                kmax--;
                pulln++;
                // copy
                if (p[k2].piflag == 0){
                    p[k1] = p[k2];
                    break;
                }
            }
        }
    }
    // update
    prm->pNum -= pulln;
    // output log
    std::map<std::string, std::string>data ={
        {"particleNumber", std::to_string(prm->pNum)},
        {"pulledPtclNum", std::to_string(pulln)},
    };
    NSString* secName = [NSString stringWithFormat:@"flowout_%@", _pName];
    logger.logSection([secName UTF8String], data);
    // NSLog(@"flowout(%@): pulln = %d, pNum = %d", _pName, pulln, prm->pNum);
}

- (void)integrateChargeDensity:(EMField*)fld withLogger:(XmlLogger&)logger{
    // obtain objects
    SimulationParams* prm = (SimulationParams*)[_paramsBuffer contents];
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
    float constRho = _q * _w / (fld.dx * fld.dy);
    id<MTLBuffer> constRhoBuffer = [_device newBufferWithBytes:&constRho
                                                length:sizeof(int)
                                                options:MTLResourceStorageModeShared];
    // 分割して積分
    uint pNumPerThread = (prm->pNum + _integrationChunkSize - 1) / _integrationChunkSize;
    id<MTLBuffer> pNumPerThreadBuffer = [_device newBufferWithBytes:&pNumPerThread
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
    [computeEncoder setBuffer:pNumPerThreadBuffer       offset:0 atIndex:7];
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
};

- (void)outputPhaseSpace:(int)cycle withEMField:(EMField*)fld withLogger:(XmlLogger&)logger{
    // obtain objects
    ParticleState *p = (ParticleState*)[_particleBuffer contents];
    SimulationParams *prm = (SimulationParams*)[_paramsBuffer contents];

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

- (int)injection:(double)dt withParam:(Init*)initParam withCurrent:(int)current withLogger:(XmlLogger&)logger{
    // パラメータ取得
    struct ParamForField fieldParam = [initParam getParamForField];
    NSArray *Sources = [initParam getParticleSources];
    struct SourceForParticle source;
    // オブジェクト取得
    SimulationParams *prm = (SimulationParams*)[_paramsBuffer contents];
    // ログ出力用データセット
    std::map<std::string, std::string>data;
    
    for (int i = 0; i < Sources.count; i++){
        NSValue *value = Sources[i];
        [value getValue:&source];

        if([source.pName isEqualToString:_pName]){
            // 生成範囲
            float Xmin, Xmax, Lx;
            float Ymin, Ymax, Ly;
            if (source.genX[0] < 0){
                // auto
                Lx = fieldParam.dx * fieldParam.ngx;
                Xmin = 0.0;
                Xmax = Lx;
            }else{
                Lx = source.genX[1] - source.genX[0];
                Xmin = source.genX[0];
                Xmax = source.genX[1];
            }
            if (source.genY[0] < 0){
                // auto
                Ly = fieldParam.dy * fieldParam.ngy;
                Ymin = 0.0;
                Ymax = Ly;
            }else{
                Ly = source.genY[1] - source.genY[0];
                Ymin = source.genY[0];
                Ymax = source.genY[1];
            }
            // 乱数生成器
            std::random_device seed_gen;
            std::default_random_engine engine(seed_gen());
            std::uniform_real_distribution<> unif_dist(0.0, 1.0);
            float vth_epi  = sqrt(2*kb*source.genT/_m);
            std::normal_distribution<> norm_dist(0.0, vth_epi);
            // 粒子生成数
            int addn;
            double addn_d;
            if ([source.genType isEqualToString:@"hollow-cathode"]){
                addn = -current; // アノード電流の符号を反転した分だけ電子が流入
            }else{
                addn_d = source.src*dt;
                if(unif_dist(engine) < addn_d - (int)addn_d){
                    addn = (int)addn_d + 1;
                }else{
                    addn = (int)addn_d;
                }
            }
            // check
            if(prm->pNum+addn > _pNumMax){
                NSLog(@"pNum+addn > pNumMax: pName = %@, pNum+addn = %d, pNumMax = %d", _pName, prm->pNum+addn, _pNumMax);
                return 1;
            }
            // 粒子の追加
            ParticleState *p = (ParticleState*)[_particleBuffer contents];
            for (int i = prm->pNum; i < prm->pNum + addn; i++) {
                if ([source.genType isEqualToString:@"uniform-Gaussian"] || [source.genType isEqualToString:@"hollow-cathode"]){
                    // uniform distribution for position
                    p[i].x = Xmin + (float)unif_dist(engine)*Lx;
                    p[i].y = Ymin + (float)unif_dist(engine)*Ly;
                    // Maxwellian for velocity 
                    p[i].vx = (float)source.genU[0] + (float)norm_dist(engine);
                    p[i].vy = (float)source.genU[1] + (float)norm_dist(engine);
                    p[i].vz = (float)source.genU[2] + (float)norm_dist(engine);
                } else if ([source.genType isEqualToString:@"Xsinusoidal-Gaussian"]){
                    // uniform distribution for position
                    p[i].x = (Xmin+Xmax)/2 + Lx/PI*sin(2*unif_dist(engine)-1.0);
                    p[i].y = Ymin + (float)unif_dist(engine)*Ly;
                    // Maxwellian for velocity 
                    p[i].vx = (float)source.genU[0] + (float)norm_dist(engine);
                    p[i].vy = (float)source.genU[1] + (float)norm_dist(engine);
                    p[i].vz = (float)source.genU[2] + (float)norm_dist(engine);
                } else if ([source.genType isEqualToString:@"Ysinusoidal-Gaussian"]){
                    // uniform distribution for position
                    p[i].x = Xmin + (float)unif_dist(engine)*Lx;
                    p[i].y = (Ymin+Ymax)/2 + Ly/PI*sin(2*unif_dist(engine)-1.0);
                    // Maxwellian for velocity 
                    p[i].vx = (float)source.genU[0] + (float)norm_dist(engine);
                    p[i].vy = (float)source.genU[1] + (float)norm_dist(engine);
                    p[i].vz = (float)source.genU[2] + (float)norm_dist(engine);
                }
                // shift from real coordinate to integer coodinate
                p[i].x = p[i].x/fieldParam.dx;
                p[i].y = p[i].y/fieldParam.dy;
                // shift origin for high-order weighting
                p[i].x = p[i].x + (float)fieldParam.ngb;
                p[i].y = p[i].y + (float)fieldParam.ngb;
                // deletion flag
                p[i].piflag = 0;
            }
            // 粒子数更新
            prm->pNum += addn;
            data[[source.genType UTF8String]] = std::to_string(addn);
        }

    }
    // ログ出力
    NSString* secName = [NSString stringWithFormat:@"injection_%@", _pName];
    logger.logSection([secName UTF8String], data);
    return 0;
}

// アクセサ
- (NSString*)pName { return _pName; }
- (int)pinum_Xmin { return _pinum_Xmin; }
- (int)pinum_Xmax { return _pinum_Xmax; }
- (int)pinum_Ymin { return _pinum_Ymin; }
- (int)pinum_Ymax { return _pinum_Ymax; }

@end
