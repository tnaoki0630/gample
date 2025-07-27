#import "Moment.h"
#import "Init.h"
#import "Constant.h"
#import "Particle.h"

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

// 積分計算用構造体
struct integrationParams {
    uint pNum;
    int ngx;
    int ngy;
    int ngb;
    float scale;
};

// Function Constants
constant int weightOrder   [[function_constant(0)]];

// モーメント計算カーネル
kernel void integrateMoments(
                        device ParticleState* ptcl          [[ buffer(0) ]],
                        constant integrationParams& prm     [[ buffer(1) ]],
                        device atomic_float* n              [[ buffer(2) ]],
                        device atomic_float* ux             [[ buffer(3) ]],
                        device atomic_float* uy             [[ buffer(4) ]],
                        device atomic_float* uz             [[ buffer(5) ]],
                        device atomic_float* Pxx            [[ buffer(6) ]],
                        device atomic_float* Pxy            [[ buffer(7) ]],
                        device atomic_float* Pxz            [[ buffer(8) ]],
                        device atomic_float* Pyy            [[ buffer(9) ]],
                        device atomic_float* Pyz            [[ buffer(10) ]],
                        device atomic_float* Pzz            [[ buffer(11) ]],
                        device float* print                 [[ buffer(12) ]],
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
    int ii, jj, idx_out;

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
    
    // 5th-order weighting
    for (int i = 0; i < 2; i++) {
        // 1st-order weighting
        if (weightOrder == 1){
            sf[0][i] = 1.0 - hv[i];
            sf[1][i] = hv[i];
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
            idx_out = ii+jj*(nx+1);
            atomic_fetch_add_explicit(&( n[idx_out]   ), sf[i][0]*sf[j][1]*prm.scale , memory_order_relaxed);
            atomic_fetch_add_explicit(&( ux[idx_out]  ), sf[i][0]*sf[j][1]*p.vx      , memory_order_relaxed);
            atomic_fetch_add_explicit(&( uy[idx_out]  ), sf[i][0]*sf[j][1]*p.vy      , memory_order_relaxed);
            atomic_fetch_add_explicit(&( uz[idx_out]  ), sf[i][0]*sf[j][1]*p.vz      , memory_order_relaxed);
            atomic_fetch_add_explicit(&( Pxx[idx_out] ), sf[i][0]*sf[j][1]*p.vx*p.vx , memory_order_relaxed);
            atomic_fetch_add_explicit(&( Pxy[idx_out] ), sf[i][0]*sf[j][1]*p.vx*p.vy , memory_order_relaxed);
            atomic_fetch_add_explicit(&( Pxz[idx_out] ), sf[i][0]*sf[j][1]*p.vx*p.vz , memory_order_relaxed);
            atomic_fetch_add_explicit(&( Pyy[idx_out] ), sf[i][0]*sf[j][1]*p.vy*p.vy , memory_order_relaxed);
            atomic_fetch_add_explicit(&( Pyz[idx_out] ), sf[i][0]*sf[j][1]*p.vy*p.vz , memory_order_relaxed);
            atomic_fetch_add_explicit(&( Pzz[idx_out] ), sf[i][0]*sf[j][1]*p.vz*p.vz , memory_order_relaxed);
        }
    }
    if (!isfinite(p.x)) {
        print[gid] = 1.0;
    }else{
        print[gid] = 0.0;
    }
}
)";

@implementation Moment

- (instancetype)initWithDevice:(id<MTLDevice>)device withParam:(Init*)initParam withLogger:(XmlLogger&)logger{
    self = [super init];
    if (self) {

        // パラメータ取得
        struct ParamForField fieldParam = initParam.paramForField;
        struct ParamForComputing compParam = initParam.paramForComputing;
        std::vector<struct ParamForParticle> particleParam = initParam.paramForParticle;
        std::vector<struct BoundaryConditionForParticle> pBCs = initParam.particleBoundaries;

        // 並列計算パラメータ
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

        // constant 引数付きでカーネルを生成
        MTLFunctionConstantValues *fc = [[MTLFunctionConstantValues alloc] init];
        int wo = fieldParam.weightOrder;
        [fc setConstantValue:&wo type:MTLDataTypeInt atIndex:0];

        // 引数付きでカーネルを作成
        id<MTLFunction> integrateFunction = [library newFunctionWithName:@"integrateMoments" constantValues:fc error:&error];
        _integrateMomentsPipeline = [device newComputePipelineStateWithFunction:integrateFunction error:&error];

        // 定数パラメータ格納バッファ
        size_t buffSize = sizeof(integrationParams);
        _integrationParamsBuffer = [device newBufferWithLength:buffSize options:MTLResourceStorageModeShared];
        integrationParams* prm = (integrationParams*)[_integrationParamsBuffer contents];
        prm->ngx  = fieldParam.ngx;
        prm->ngy  = fieldParam.ngy;
        prm->ngb  = fieldParam.ngb;
        prm->pNum = 0; // ptclクラスごとに都度計算
        prm->scale = 1.0; // ptclクラスごとに都度計算

        // 積分計算用バッファ
        int nx = fieldParam.ngx + 2*fieldParam.ngb;
        int ny = fieldParam.ngy + 2*fieldParam.ngb;
        int ng = (nx+1)*(ny+1);
        _nBuffer   = [device newBufferWithLength:sizeof(float)*ng options:MTLResourceStorageModeShared];
        _uxBuffer  = [device newBufferWithLength:sizeof(float)*ng options:MTLResourceStorageModeShared];
        _uyBuffer  = [device newBufferWithLength:sizeof(float)*ng options:MTLResourceStorageModeShared];
        _uzBuffer  = [device newBufferWithLength:sizeof(float)*ng options:MTLResourceStorageModeShared];
        _PxxBuffer = [device newBufferWithLength:sizeof(float)*ng options:MTLResourceStorageModeShared];
        _PxyBuffer = [device newBufferWithLength:sizeof(float)*ng options:MTLResourceStorageModeShared];
        _PxzBuffer = [device newBufferWithLength:sizeof(float)*ng options:MTLResourceStorageModeShared];
        _PyyBuffer = [device newBufferWithLength:sizeof(float)*ng options:MTLResourceStorageModeShared];
        _PyzBuffer = [device newBufferWithLength:sizeof(float)*ng options:MTLResourceStorageModeShared];
        _PzzBuffer = [device newBufferWithLength:sizeof(float)*ng options:MTLResourceStorageModeShared];
        
        // デバッグ出力用バッファ
        buffSize = sizeof(float)*compParam.pNumMax;
        _printBuffer = [device newBufferWithLength:buffSize options:MTLResourceStorageModeShared];
    }
    return self;
}

- (void)integrateMoments:(Particle*)ptcl withEMField:(EMField*)fld withLogger:(XmlLogger&)logger{
    // initialize
    integrationParams* prm = (integrationParams*)[_integrationParamsBuffer contents];
    int ng = (prm->ngx + 2*prm->ngb + 1)*(prm->ngy + 2*prm->ngb + 1);
    float* n = (float*)[_nBuffer contents];
    float* ux = (float*)[_uxBuffer contents];
    float* uy = (float*)[_uyBuffer contents];
    float* uz = (float*)[_uzBuffer contents];
    float* Pxx = (float*)[_PxxBuffer contents];
    float* Pxy = (float*)[_PxyBuffer contents];
    float* Pxz = (float*)[_PxzBuffer contents];
    float* Pyy = (float*)[_PyyBuffer contents];
    float* Pyz = (float*)[_PyzBuffer contents];
    float* Pzz = (float*)[_PzzBuffer contents];
    for (int i = 0; i < ng; i++){
        n[i] = 0.0;
        ux[i] = 0.0;
        uy[i] = 0.0;
        uz[i] = 0.0;
        Pxx[i] = 0.0;
        Pxy[i] = 0.0;
        Pxz[i] = 0.0;
        Pyy[i] = 0.0;
        Pyz[i] = 0.0;
        Pzz[i] = 0.0;
    }
    // for (int i = 0; i < 10; i++){
    //     NSLog(@"Moment(before): n[%d] = %e",i,n[i]);
    // }
    
    // 粒子データは Particle クラスのバッファを使用
    id<MTLBuffer> particleBuffer = [ptcl particleBuffer];
    id<MTLBuffer> paramsBuffer = [ptcl paramsBuffer];
    SimulationParams* prmPtcl = (SimulationParams*)[paramsBuffer contents];
    // 定数パラメータ更新
    prm->pNum = prmPtcl->pNum;
    prm->scale = ptcl.w/(fld.dx*fld.dy);
    
    // グリッドとスレッドグループのサイズ設定
    uint threadGroupNum = (prm->pNum + _threadGroupSize - 1) / _threadGroupSize;
    MTLSize gridSizeMetalStyle = MTLSizeMake(threadGroupNum, 1, 1);
    MTLSize threadGroupSize = MTLSizeMake(_threadGroupSize, 1, 1);

    // コマンドバッファとエンコーダ
    id<MTLCommandBuffer> commandBuffer;
    id<MTLComputeCommandEncoder> computeEncoder;
    
    // 積分実行
    commandBuffer = [_commandQueue commandBuffer];
    computeEncoder = [commandBuffer computeCommandEncoder];
    [computeEncoder setComputePipelineState:_integrateMomentsPipeline];
    [computeEncoder setBuffer:particleBuffer            offset:0 atIndex:0];
    [computeEncoder setBuffer:_integrationParamsBuffer  offset:0 atIndex:1];
    [computeEncoder setBuffer:_nBuffer                  offset:0 atIndex:2];
    [computeEncoder setBuffer:_uxBuffer                 offset:0 atIndex:3];
    [computeEncoder setBuffer:_uyBuffer                 offset:0 atIndex:4];
    [computeEncoder setBuffer:_uzBuffer                 offset:0 atIndex:5];
    [computeEncoder setBuffer:_PxxBuffer                offset:0 atIndex:6];
    [computeEncoder setBuffer:_PxyBuffer                offset:0 atIndex:7];
    [computeEncoder setBuffer:_PxzBuffer                offset:0 atIndex:8];
    [computeEncoder setBuffer:_PyyBuffer                offset:0 atIndex:9];
    [computeEncoder setBuffer:_PyzBuffer                offset:0 atIndex:10];
    [computeEncoder setBuffer:_PzzBuffer                offset:0 atIndex:11];
    [computeEncoder setBuffer:_printBuffer              offset:0 atIndex:12];
    [computeEncoder dispatchThreadgroups:gridSizeMetalStyle threadsPerThreadgroup:threadGroupSize];
    [computeEncoder endEncoding];
    [commandBuffer commit];
    [commandBuffer waitUntilCompleted];

    // 積分量->モーメント量
    n = (float*)[_nBuffer contents];
    ux = (float*)[_uxBuffer contents];
    uy = (float*)[_uyBuffer contents];
    uz = (float*)[_uzBuffer contents];
    Pxx = (float*)[_PxxBuffer contents];
    Pxy = (float*)[_PxyBuffer contents];
    Pxz = (float*)[_PxzBuffer contents];
    Pyy = (float*)[_PyyBuffer contents];
    Pyz = (float*)[_PyzBuffer contents];
    Pzz = (float*)[_PzzBuffer contents];
    // cm/s*ptcl -> cm/s
    for (int i = 0; i < ng; i++){
        if(n[i] > 1e-20){
            ux[i] /= n[i]/prm->scale;
            uy[i] /= n[i]/prm->scale;
            uz[i] /= n[i]/prm->scale;
        }else{
            ux[i] = 0.0;
            uy[i] = 0.0;
            uz[i] = 0.0;
        }
    }
    // cm2/s2*ptcl -> g*cm2/s2*1/cm3 = 0.1 Pa
    // Cov[X,Y] = E[XY] - E[X]E[Y]
    for (int i = 0; i < ng; i++){
        if(n[i] > 1e-20){
            Pxx[i] = ptcl.m*(Pxx[i]/(n[i]/prm->scale) - ux[i]*ux[i])*n[i];
            Pxy[i] = ptcl.m*(Pxy[i]/(n[i]/prm->scale) - ux[i]*uy[i])*n[i];
            Pxz[i] = ptcl.m*(Pxz[i]/(n[i]/prm->scale) - ux[i]*uz[i])*n[i];
            Pyy[i] = ptcl.m*(Pyy[i]/(n[i]/prm->scale) - uy[i]*uy[i])*n[i];
            Pyz[i] = ptcl.m*(Pyz[i]/(n[i]/prm->scale) - uy[i]*uz[i])*n[i];
            Pzz[i] = ptcl.m*(Pzz[i]/(n[i]/prm->scale) - uz[i]*uz[i])*n[i];
        }else{
            Pxx[i] = 0.0;
            Pxy[i] = 0.0;
            Pxz[i] = 0.0;
            Pyy[i] = 0.0;
            Pyz[i] = 0.0;
            Pzz[i] = 0.0;
        }
    }

    // // debug print
    // float* prt = (float*)_printBuffer.contents;
    // for (int i = 0; i < prm->pNum; i++){
    //     if (prt[i] > 1e-20){
    //         NSLog(@"not finite: prt[%d] = %e",i,prt[i]);
    //     }
    //     // NSLog(@"Moment: n[%d] = %e",i,n[i]);
    // }

};

// アクセサ
- (id<MTLBuffer>)printBuffer { return _printBuffer; }
- (id<MTLBuffer>)nBuffer { return _nBuffer; }
- (id<MTLBuffer>)uxBuffer { return _uxBuffer; }
- (id<MTLBuffer>)uyBuffer { return _uyBuffer; }
- (id<MTLBuffer>)uzBuffer { return _uzBuffer; }
- (id<MTLBuffer>)PxxBuffer { return _PxxBuffer; }
- (id<MTLBuffer>)PxyBuffer { return _PxyBuffer; }
- (id<MTLBuffer>)PxzBuffer { return _PxzBuffer; }
- (id<MTLBuffer>)PyyBuffer { return _PyyBuffer; }
- (id<MTLBuffer>)PyzBuffer { return _PyzBuffer; }
- (id<MTLBuffer>)PzzBuffer { return _PzzBuffer; }

@end