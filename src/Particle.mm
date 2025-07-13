#import "Particle.h"
#import "Init.h"
#import "Constant.h"
#import <random>
#import <string>
#import <iostream>

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
};

// 積分計算用構造体
struct integrationParams {
    uint pNum;
    int ngx;
    int ngy;
    int ngb;
    float scale;
};

// Function Constants
constant int BC_Xmin        [[function_constant(0)]];
constant int BC_Xmax        [[function_constant(1)]];
constant int BC_Ymin        [[function_constant(2)]];
constant int BC_Ymax        [[function_constant(3)]];
constant int weightOrder   [[function_constant(4)]];

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
    
    // get particle
    device ParticleState& p = ptcl[id];

    // electro-magnetic field on each ptcl
    float xh = p.x - 0.5;
    float yh = p.y - 0.5;
    int i1 = int(xh);
    int j1 = int(p.y);
    int i2 = int(p.x);
    int j2 = int(yh);
    // protect
    if (i1 < prm.ngb-1) i1 = prm.ngb-1;
    if (i1 > prm.ngb-1+prm.ngx) i1 = prm.ngb-1+prm.ngx;
    if (j1 < prm.ngb) j1 = prm.ngb;
    if (j1 >= prm.ngb+prm.ngy) j1 = prm.ngb+prm.ngy-1;
    if (i2 < prm.ngb) i2 = prm.ngb;
    if (i2 >= prm.ngb+prm.ngx) i2 = prm.ngb+prm.ngx-1;
    if (j2 < prm.ngb-1) j2 = prm.ngb-1;
    if (j2 > prm.ngb-1+prm.ngy) j2 = prm.ngb-1+prm.ngy;

    float hv[2][2] ;
    hv[0][0] = xh  - float(i1);
    hv[0][1] = p.y - float(j1);
    hv[1][0] = p.x - float(i2);
    hv[1][1] = yh  - float(j2);
    
    // weighting
    float sc;
    float sf[6][2][2];
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            // 1st-order weighting
            if (weightOrder == 1){
                sf[0][i][j] = hv[i][j];
                sf[1][i][j] = 1.0 - hv[i][j];
            // 5th-order weighting
            } else if (weightOrder == 5){
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
    for (int i = 0; i < weightOrder+1; i++) {
        for (int j = 0; j < weightOrder+1; j++) {
            ii = i1+(i-prm.ngb+1);
            jj = j1+(j-prm.ngb);
            Epx += sf[i][0][0]*sf[j][0][1]*Ex[ii+jj*(nx+2)];
            ii = i2+(i-prm.ngb);
            jj = j2+(j-prm.ngb+1);
            Epy += sf[i][1][0]*sf[j][1][1]*Ey[ii+jj*(nx+1)];
            ii = i1+(i-prm.ngb+1);
            jj = j2+(j-prm.ngb+1);
            Bpz += sf[i][0][0]*sf[j][1][1]*Bz[ii+jj*(nx+2)];
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
    if (BC_Xmin == 0) {
        p.x = fmod(p.x + float(prm.ngx - prm.ngb), float(prm.ngx)) + float(prm.ngb);
    }
    if (BC_Ymin == 0) {
        p.y = fmod(p.y + float(prm.ngy - prm.ngb), float(prm.ngy)) + float(prm.ngb);
    }

    // delete boundary
    if (BC_Xmin == 1){
        if(int(floor(p.x)) < prm.ngb){
            p.piflag = 1;
        }
    }
    if (BC_Xmax == 1){
        if(int(floor(p.x)) >= prm.ngb+prm.ngx){
            p.piflag = 2;
        }
    }
    if (BC_Ymin == 1){
        if(int(floor(p.y)) < prm.ngb){
            p.piflag = 3;
        }
    }
    if (BC_Ymax == 1){
        if(int(floor(p.y)) >= prm.ngb+prm.ngy){
            p.piflag = 4;
        }
    }

    // debug print
    if (!isfinite(p.x)) {
        print[id] = 1.0;
    }else{
        print[id] = 0.0;
    }
}

// 電荷密度更新カーネル
kernel void integrateChargeDensity(
                        device ParticleState* ptcl          [[ buffer(0) ]],
                        constant integrationParams& prm     [[ buffer(1) ]],
                        device atomic_float* rho            [[ buffer(2) ]],
                        device float* print                 [[ buffer(3) ]],
                        uint gid                            [[ thread_position_in_grid ]]
                        ) {
    
    // 粒子数を超えたらスキップ
    if(gid >= prm.pNum) return;

    // 変数定義
    const int nx = prm.ngx + 2*prm.ngb;
    const int ny = prm.ngy + 2*prm.ngb;
    const int ng = (nx+1)*(ny+1);
    int i1, j1;
    float hv[2];
    float sc;
    float sf[6][2];
    int ii, jj;

    // 粒子を取得
    device ParticleState& p = ptcl[gid];

    // electro-magnetic field on each ptcl
    i1 = int(p.x);
    j1 = int(p.y);

    // protect
    if (i1 < prm.ngb-1) i1 = prm.ngb-1;
    if (i1 > prm.ngb-1+prm.ngx) i1 = prm.ngb-1+prm.ngx;
    if (j1 < prm.ngb) j1 = prm.ngb;
    if (j1 >= prm.ngb+prm.ngy) j1 = prm.ngb+prm.ngy-1;
    
    hv[0] = p.x - float(i1);
    hv[1] = p.y - float(j1);
    
    // weighting
    for (int i = 0; i < 2; i++) {
        // 1st-order weighting
        if (weightOrder == 1){
            sf[0][i] = hv[i];
            sf[1][i] = 1.0 - hv[i];
        // 5th-order weighting
        } else if (weightOrder == 5){
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
    }

    // accumulation
    for (int i = 0; i < weightOrder+1; i++) {
        for (int j = 0; j < weightOrder+1; j++) {
            ii = i1+(i-prm.ngb);
            jj = j1+(j-prm.ngb);
            atomic_fetch_add_explicit(&(rho[ii+jj*(nx+1)]), sf[i][0]*sf[j][1]*prm.scale, memory_order_relaxed);
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
}

// 初期設定
- (instancetype)initWithDevice:(id<MTLDevice>)device withParam:(Init*)initParam specimen:(int)s withLogger:(XmlLogger&)logger{
    self = [super init];
    if (self) {

        // パラメータ取得
        struct ParamForField fieldParam = initParam.paramForField;
        struct ParamForComputing compParam = initParam.paramForComputing;
        std::vector<struct ParamForParticle> particleParam = initParam.paramForParticle;
        std::vector<struct BoundaryConditionForParticle> pBCs = initParam.particleBoundaries;
    
        // 粒子パラメータ
        _pName = particleParam[s].pName;
        _pNumMax = compParam.pNumMax;
        _m = particleParam[s].m;
        _q = particleParam[s].q;
        _w = particleParam[s].w;

        // 並列計算パラメータ
        _integrationChunkSize = compParam.integrationChunkSize;
        _threadGroupSize = compParam.threadGroupSize;
        
        // コマンドキューの作成
        _device = device;
        _commandQueue = [device newCommandQueue];

        // コンピュートパイプラインの設定
        NSError *error = nil;
        id<MTLLibrary> library = [device newLibraryWithSource:kMetalShaderSource options:nil error:&error];
        if (!library) {
            NSLog(@"Failed to create Metal library: %@", error);
            return nil;
        }
        
        // バッファサイズ
        size_t buffSize;
        // シミュレーションパラメータバッファ
        buffSize = sizeof(SimulationParams);
        _paramsBuffer = [device newBufferWithLength:buffSize options:MTLResourceStorageModeShared];
        SimulationParams* prm = (SimulationParams*)[_paramsBuffer contents];
        // パラメータ格納
        prm->pNum = particleParam[s].pNum;
        prm->ngx = fieldParam.ngx;
        prm->ngy = fieldParam.ngy;
        prm->ngb = fieldParam.ngb;

        // constant 引数付きでカーネルを生成
        MTLFunctionConstantValues *fc = [[MTLFunctionConstantValues alloc] init];
        // 境界条件格納
        int BC;
        for (int pos = 0; pos < pBCs.size(); pos++){
            if ([pBCs[pos].position isEqualToString:@"Xmin"]){
                if ([pBCs[pos].type isEqualToString:@"periodic"]){
                    BC = 0;
                }else if ([pBCs[pos].type isEqualToString:@"Delete"]){
                    BC = 1;
                }else{
                    BC = -1; // error
                }
                // function constant をセット
                [fc setConstantValue:&BC type:MTLDataTypeInt atIndex:0];
            }else if ([pBCs[pos].position isEqualToString:@"Xmax"]){
                if ([pBCs[pos].type isEqualToString:@"periodic"]){
                    BC = 0;
                }else if ([pBCs[pos].type isEqualToString:@"Delete"]){
                    BC = 1;
                }else{
                    BC = -1; // error
                }
                // function constant をセット
                [fc setConstantValue:&BC type:MTLDataTypeInt atIndex:1];
            }else if ([pBCs[pos].position isEqualToString:@"Ymin"]){
                if ([pBCs[pos].type isEqualToString:@"periodic"]){
                    BC = 0;
                }else if ([pBCs[pos].type isEqualToString:@"Delete"]){
                    BC = 1;
                }else{
                    BC = -1; // error
                }
                // function constant をセット
                [fc setConstantValue:&BC type:MTLDataTypeInt atIndex:2];
            }else if ([pBCs[pos].position isEqualToString:@"Ymax"]){
                if ([pBCs[pos].type isEqualToString:@"periodic"]){
                    BC = 0;
                }else if ([pBCs[pos].type isEqualToString:@"Delete"]){
                    BC = 1;
                }else{
                    BC = -1; // error
                }
                // function constant をセット
                [fc setConstantValue:&BC type:MTLDataTypeInt atIndex:3];
            }
            if (BC < 0){
                NSLog(@"[FatalError] invalid particle boundary");
                return nil;
            }
        }

        // weighting order
        int wo = fieldParam.weightOrder;
        [fc setConstantValue:&wo type:MTLDataTypeInt atIndex:4];

        // 引数付きでカーネルを作成
        id<MTLFunction> updateParticlesFunction = [library newFunctionWithName:@"updateParticles" constantValues:fc error:&error];
        _updateParticlesPipeline = [device newComputePipelineStateWithFunction:updateParticlesFunction error:&error];
        id<MTLFunction> integrateChargeDensityFunction = [library newFunctionWithName:@"integrateChargeDensity" constantValues:fc error:&error];
        _integrateChargeDensityPipeline = [device newComputePipelineStateWithFunction:integrateChargeDensityFunction error:&error];

        // 定数パラメータ格納バッファ
        buffSize = sizeof(integrationParams);
        _integrationParamsBuffer = [device newBufferWithLength:buffSize options:MTLResourceStorageModeShared];
        integrationParams* intgPrm = (integrationParams*)[_integrationParamsBuffer contents];
        intgPrm->ngx  = fieldParam.ngx;
        intgPrm->ngy  = fieldParam.ngy;
        intgPrm->ngb  = fieldParam.ngb;
        intgPrm->scale  = _q * _w / (fieldParam.dx * fieldParam.dy);
        intgPrm->pNum = 0; // 都度計算

        // 粒子バッファ
        if(particleParam[s].pNum > _pNumMax){
            NSLog(@"pNum > pNumMax: pName = %@, pNum = %d, pNumMax = %d", _pName, particleParam[s].pNum, _pNumMax);
            return nil;
        }
        buffSize = sizeof(ParticleState)*_pNumMax; // maxで初期化
        _particleBuffer = [device newBufferWithLength:buffSize options:MTLResourceStorageModeShared];
        // 初期粒子分布の設定
        float Xmin, Lx;
        float Ymin, Ly;
        if (particleParam[s].genX[0] < 0){
            Lx = fieldParam.dx * fieldParam.ngx;
            Xmin = 0.0;
        }else{
            Lx = particleParam[s].genX[1] - particleParam[s].genX[0];
            Xmin = particleParam[s].genX[0];
        }
        if (particleParam[s].genY[0] < 0){
            Ly = fieldParam.dy * fieldParam.ngy;
            Ymin = 0.0;
        }else{
            Ly = particleParam[s].genY[1] - particleParam[s].genY[0];
            Ymin = particleParam[s].genY[0];
        }

        // 乱数生成器
        std::random_device seed_gen;
        std::default_random_engine engine(seed_gen());
        std::uniform_real_distribution<> unif_dist(0.0, 1.0);
        float vth_epi  = sqrt(2*kb*particleParam[s].genT/particleParam[s].m);
        std::normal_distribution<> norm_dist(0.0, vth_epi);
        // 粒子の初期化
        ParticleState *p = (ParticleState *)self.particleBuffer.contents;
        for (int i = 0; i < particleParam[s].pNum; i++) {
            float rx, ry;
            if ([particleParam[s].genType isEqualToString:@"uniform-Gaussian"]){
                // uniform distribution for position
                rx = (float)unif_dist(engine)*Lx;
                ry = (float)unif_dist(engine)*Ly;
                p[i].x = Xmin + rx;
                p[i].y = Ymin + ry;
                // Maxwellian for velocity 
                p[i].vx = (float)particleParam[s].genU[0] + (float)norm_dist(engine);
                p[i].vy = (float)particleParam[s].genU[1] + (float)norm_dist(engine);
                p[i].vz = (float)particleParam[s].genU[2] + (float)norm_dist(engine);
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

        // 電磁場バッファ
        int nx = fieldParam.ngx + 2*fieldParam.ngb;
        int ny = fieldParam.ngy + 2*fieldParam.ngb;
        size_t EMbuffSize = sizeof(float)*(nx+1)*(ny+1);
    }
    return self;
}


// 時間更新
- (void)update:(double)dt withEMField:(EMField*)fld withMom:(Moment*)mom withLogger:(XmlLogger&)logger{
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
    
    // Moment のバッファを使用
    id<MTLBuffer> printBuffer = [mom printBuffer];

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
    [computeEncoder setBuffer:printBuffer       offset:0 atIndex:8];

    
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

    // デバッグ出力
    ParticleState* p = (ParticleState*)[_particleBuffer contents];
    float* prt = (float*)printBuffer.contents;
    float min_x = 1e20, max_x = -1e20;
    float min_y = 1e20, max_y = -1e20;
    float min_v = 1e20, max_v = -1e20;
    float min = 1e20, max = -1e20;
    for (int idx = 0; idx < prm->pNum; idx++){
        if(p[idx].piflag == 0){
            // debug print
            // if (min > prt[idx]){ min = prt[idx]; }
            // if (max < prt[idx]){ max = prt[idx]; }
            if (prt[idx] > 1e-20){ NSLog(@"[update] not finite: prt[%d] = %e",idx,prt[idx]); }
            // position
            if (min_x > p[idx].x){ min_x = p[idx].x; }
            if (max_x < p[idx].x){ max_x = p[idx].x; }
            if (min_y > p[idx].y){ min_y = p[idx].y; }
            if (max_y < p[idx].y){ max_y = p[idx].y; }
            // velocity
            float v = sqrt(p[idx].vx*p[idx].vx+p[idx].vy*p[idx].vy+p[idx].vz*p[idx].vz);
            if (v > c){ NSLog(@"[update] exceeded speed of light: p[%d].v = %e",idx,v); }
            if (min_v > v){ min_v = v; }
            if (max_v < v){ max_v = v; }
        }
    }
    // NSLog(@"debug print(%@): pNum = %d, min = %e, max = %e", _pName, prm->pNum, min, max);
    // NSLog(@"update(%@): pNum = %d, min_x = %e, max_x = %e, min_y = %e, max_y = %e, min_v = %e, max_v = %e", _pName, prm->pNum, min_x, max_x, min_y, max_y, min_v, max_v);

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
        if (k1 >= kmax){
            if (p[k1].piflag > 0){
                pulln++;
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
            for (int k2 = kend-1; k2 >= k1; k2--){
                kmax--;
                pulln++;
                // copy
                if (p[k2].piflag == 0){
                    p[k1] = p[k2];
                    break;
                }
                // count particles for each boundary with charge
                if (k2 > k1){
                    if (p[k2].piflag == 1){
                        _pinum_Xmin += (int)(_q/(float)ec);
                    }else if(p[k2].piflag == 2){
                        _pinum_Xmax += (int)(_q/(float)ec);
                    }else if(p[k2].piflag == 3){
                        _pinum_Ymin += (int)(_q/(float)ec);
                    }else if(p[k2].piflag == 4){
                        _pinum_Ymax += (int)(_q/(float)ec);
                    }
                }
            }
        }
    }
    // update
    prm->pNum -= pulln;
    // check reduction
    for (int k1 = 0; k1 < prm->pNum; k1++){
        if (p[k1].piflag > 0){
            NSLog(@"reduction failed: pName = %@, idx = %d", _pName, k1);
        }else if( int(p[k1].x) < prm->ngb
               || int(p[k1].x) > prm->ngb+prm->ngx
               || int(p[k1].y) < prm->ngb
               || int(p[k1].y) > prm->ngb+prm->ngy ){
            NSLog(@"spilled particle is detected: pName = %@, idx = %d, p.x = %e, p.y = %e", _pName, k1, p[k1].x, p[k1].y);
        }
    }
    // output log
    std::map<std::string, std::string>data ={
        {"particleNumber", std::to_string(prm->pNum)},
        {"pulledPtclNum", std::to_string(pulln)},
        {"Xmin", std::to_string(_pinum_Xmin)},
        {"Xmax", std::to_string(_pinum_Xmax)},
        {"Ymin", std::to_string(_pinum_Ymin)},
        {"Ymax", std::to_string(_pinum_Ymax)},
    };
    NSString* secName = [NSString stringWithFormat:@"flowout_%@", _pName];
    logger.logSection([secName UTF8String], data);
}

- (void)integrateChargeDensity:(EMField*)fld withMoment:(Moment*)mom withLogger:(XmlLogger&)logger{
    // obtain buffer
    id<MTLBuffer> rhoBuffer = [fld rhoBuffer];
    id<MTLBuffer> printBuffer = [mom printBuffer];
    SimulationParams* prm = (SimulationParams*)[_paramsBuffer contents];
    integrationParams* intgPrm = (integrationParams*)[_integrationParamsBuffer contents];
    intgPrm->pNum = prm->pNum;
    
    // コマンドバッファとエンコーダの作成
    id<MTLCommandBuffer> commandBuffer = [_commandQueue commandBuffer];
    id<MTLComputeCommandEncoder> computeEncoder = [commandBuffer computeCommandEncoder];
    
    // パイプラインとバッファの設定
    [computeEncoder setComputePipelineState:_integrateChargeDensityPipeline];
    [computeEncoder setBuffer:_particleBuffer           offset:0 atIndex:0];
    [computeEncoder setBuffer:_integrationParamsBuffer  offset:0 atIndex:1];
    [computeEncoder setBuffer:rhoBuffer                 offset:0 atIndex:2];
    [computeEncoder setBuffer:printBuffer               offset:0 atIndex:3];

    // グリッドとスレッドグループのサイズ設定
    uint threadGroupNum = (prm->pNum + _threadGroupSize - 1) / _threadGroupSize;
    MTLSize gridSizeMetalStyle = MTLSizeMake(threadGroupNum, 1, 1);
    MTLSize threadGroupSize = MTLSizeMake(_threadGroupSize, 1, 1);

    // ディスパッチ/エンコーディング/実行
    [computeEncoder dispatchThreadgroups:gridSizeMetalStyle threadsPerThreadgroup:threadGroupSize];
    [computeEncoder endEncoding];
    [commandBuffer commit];
    [commandBuffer waitUntilCompleted];

};

- (int)injection:(double)dt withParam:(Init*)initParam withCurrent:(int&)current withLogger:(XmlLogger&)logger{
    // パラメータ取得
    struct ParamForField fieldParam = initParam.paramForField;
    const std::vector<struct SourceForParticle> sources = initParam.particleSources;

    // オブジェクト取得
    SimulationParams *prm = (SimulationParams*)[_paramsBuffer contents];
    // ログ出力用データセット
    std::map<std::string, std::string>data;
    
    for (int i = 0; i < sources.size(); i++){
        if([sources[i].pName isEqualToString:_pName]){
            // 生成範囲
            float Xmin, Xmax, Lx;
            float Ymin, Ymax, Ly;
            if (sources[i].genX[0] < 0){
                // auto
                Lx = fieldParam.dx * fieldParam.ngx;
                Xmin = 0.0;
                Xmax = Lx;
            }else{
                Lx = sources[i].genX[1] - sources[i].genX[0];
                Xmin = sources[i].genX[0];
                Xmax = sources[i].genX[1];
            }
            if (sources[i].genY[0] < 0){
                // auto
                Ly = fieldParam.dy * fieldParam.ngy;
                Ymin = 0.0;
                Ymax = Ly;
            }else{
                Ly = sources[i].genY[1] - sources[i].genY[0];
                Ymin = sources[i].genY[0];
                Ymax = sources[i].genY[1];
            }
            // 乱数生成器
            std::random_device seed_gen;
            std::default_random_engine engine(seed_gen());
            std::uniform_real_distribution<> unif_dist(0.0, 1.0);
            float vth_epi  = sqrt(2*kb*sources[i].genT/_m);
            std::normal_distribution<> norm_dist(0.0, vth_epi);
            // 粒子生成数
            int addn;
            double addn_d;
            if ([sources[i].genType isEqualToString:@"hollow-cathode"]){
                if (current < 0){
                    addn = -current; // アノード電流の符号を反転した分だけ電子が流入
                    current = 0; // リセット
                }else{
                    addn = 0; // 追加しない
                    current = 0; // 自サイクルに引き継ぐならここはコメントアウトする。電子の方が無くなりやすいので、引き継がずにリセットした方が定常収束は早いはず。
                }
                data["keptCurrent"] = std::to_string(-current);
            }else{
                addn_d = sources[i].src/_w*dt;
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
            for (int idx = prm->pNum; idx < prm->pNum + addn; idx++) {
                if ([sources[i].genType isEqualToString:@"uniform-Gaussian"] || [sources[i].genType isEqualToString:@"hollow-cathode"]){
                    // uniform distribution for position
                    p[idx].x = Xmin + (float)unif_dist(engine)*Lx;
                    p[idx].y = Ymin + (float)unif_dist(engine)*Ly;
                    // Maxwellian for velocity 
                    p[idx].vx = (float)sources[i].genU[0] + (float)norm_dist(engine);
                    p[idx].vy = (float)sources[i].genU[1] + (float)norm_dist(engine);
                    p[idx].vz = (float)sources[i].genU[2] + (float)norm_dist(engine);
                } else if ([sources[i].genType isEqualToString:@"Xsinusoidal-Gaussian"]){
                    // uniform distribution for position
                    p[idx].x = (Xmin+Xmax)/2 + Lx/PI*asin(2*unif_dist(engine)-1.0);
                    p[idx].y = Ymin + (float)unif_dist(engine)*Ly;
                    // Maxwellian for velocity 
                    p[idx].vx = (float)sources[i].genU[0] + (float)norm_dist(engine);
                    p[idx].vy = (float)sources[i].genU[1] + (float)norm_dist(engine);
                    p[idx].vz = (float)sources[i].genU[2] + (float)norm_dist(engine);
                } else if ([sources[i].genType isEqualToString:@"Ysinusoidal-Gaussian"]){
                    // uniform distribution for position
                    p[idx].x = Xmin + (float)unif_dist(engine)*Lx;
                    p[idx].y = (Ymin+Ymax)/2 + Ly/PI*asin(2*unif_dist(engine)-1.0);
                    // Maxwellian for velocity 
                    p[idx].vx = (float)sources[i].genU[0] + (float)norm_dist(engine);
                    p[idx].vy = (float)sources[i].genU[1] + (float)norm_dist(engine);
                    p[idx].vz = (float)sources[i].genU[2] + (float)norm_dist(engine);
                }
                // shift from real coordinate to integer coodinate
                p[idx].x = p[idx].x/fieldParam.dx;
                p[idx].y = p[idx].y/fieldParam.dy;
                // shift origin for high-order weighting
                p[idx].x = p[idx].x + (float)fieldParam.ngb;
                p[idx].y = p[idx].y + (float)fieldParam.ngb;
                // initialize deletion flag
                p[idx].piflag = 0;
            }
            // 粒子数更新
            prm->pNum += addn;
            data[[sources[i].genType UTF8String]] = std::to_string(addn);
        }

    }
    // ログ出力
    NSString* secName = [NSString stringWithFormat:@"injection_%@", _pName];
    logger.logSection([secName UTF8String], data);
    return 0;
}

// アクセサ
- (NSString*)pName { return _pName; }
- (id<MTLBuffer>)paramsBuffer { return _paramsBuffer; }
- (id<MTLBuffer>)particleBuffer { return _particleBuffer;}
- (float)m {return _m; }
- (float)q {return _q; }
- (float)w {return _w; }
- (int)pinum_Xmin { return _pinum_Xmin; }
- (int)pinum_Xmax { return _pinum_Xmax; }
- (int)pinum_Ymin { return _pinum_Ymin; }
- (int)pinum_Ymax { return _pinum_Ymax; }

@end
