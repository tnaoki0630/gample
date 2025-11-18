#import "Collision.h"
#import "Constant.h"
#import <random>
#import <string>
#import <iostream>

// Metal シェーダーのソースコード
static NSString *const kMetalShaderSource = @R"(
#include <metal_stdlib>
using namespace metal;

#define PI 3.141592653589793

// 乱数生成用構造体
struct RNGState { 
    ulong state;
    ulong inc; 
};

// PCG-XSH-RR 32bit
inline uint pcg32(thread RNGState &s) {
    ulong old = s.state;
    s.state = old * 6364136223846793005UL + s.inc;
    uint xorshifted = uint(((old >> 18u) ^ old) >> 27u);
    uint rot = uint(old >> 59u);
    return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
}
inline float rand01(thread RNGState &s) {
    return float(pcg32(s)) * (1.0f / 4294967296.0f);
}

// 粒子構造体
struct ParticleState {
    float x;
    float y;
    float vx;
    float vy;
    float vz;
    int piflag;
};

// 粒子生成用構造体
struct generationParams {
    uint addn;
    uint ele_pNum;
    uint ion_pNum;
    float Xmin;
    float Xmax;
    float Ymin;
    float Ymax;
    float ele_genU[3];
    float ele_vth;
    float ion_genU[3];
    float ion_vth;
    float dx;
    float dy;
    float ngb;
};

// Function Constants
constant int genType   [[function_constant(0)]];

kernel void artificialIonization(device ParticleState *ptcl_ele          [[ buffer(0) ]],
                                 device ParticleState *ptcl_ion          [[ buffer(1) ]],
                                 constant generationParams& prm        [[ buffer(2) ]],
                                 device RNGState *states            [[ buffer(3) ]],
                                 device float* print                 [[ buffer(4) ]],
                                 uint id [[thread_position_in_grid]])
{
    if(id >= prm.addn) return;
    
    // get rnd satate
    RNGState rng = states[id];
    
    // get electron
    device ParticleState &p_ele = ptcl_ele[prm.ele_pNum + id];
    
    // uniform distribution for position
    float rx = rand01(rng);
    float ry = rand01(rng);
    if (genType == 0 || genType == 3){
        // uniform
        p_ele.x = prm.Xmin + rx*(prm.Xmax-prm.Xmin);
        p_ele.y = prm.Ymin + ry*(prm.Ymax-prm.Ymin);
    }else if (genType == 1){
        // Xsinusoidal
        p_ele.x = (prm.Xmin+prm.Xmax)/2 + (prm.Xmax-prm.Xmin)/PI*asin(2*rx-1.0);
        p_ele.y = prm.Ymin + ry*(prm.Ymax-prm.Ymin);
    }else if (genType == 2){
        // Ysinusoidal
        p_ele.x = prm.Xmin + rx*(prm.Xmax-prm.Xmin);
        p_ele.y = (prm.Ymin+prm.Ymax)/2 + (prm.Ymax-prm.Ymin)/PI*asin(2*ry-1.0);
    }else if (genType == 4){
        // equal division
        int addn_x = (int)round(sqrt(prm.addn*(prm.Xmax-prm.Xmin)/(prm.Ymax-prm.Ymin))); // It is desirable that ppc be a perfect square.
        int addn_y = prm.addn/addn_x;
        float delta_x = (prm.Xmax-prm.Xmin)/addn_x;
        float delta_y = (prm.Ymax-prm.Ymin)/addn_y;
        p_ele.x = prm.Xmin + delta_x*(0.5+id%addn_x);
        p_ele.y = prm.Ymin + delta_y*(0.5+id/addn_x);
    }else if (genType == 5){
        // Xgaussian
        float rz = rand01(rng);
        p_ele.x = (prm.Xmin+prm.Xmax)/2 + (prm.Xmax-prm.Xmin)*sqrt(-log(rx))*cos(2.0*PI*rz);
        p_ele.y = prm.Ymin + ry*(prm.Ymax-prm.Ymin);
    }else if (genType == 6){
        // Ygaussian
        float rz = rand01(rng);
        p_ele.x = prm.Xmin + rx*(prm.Xmax-prm.Xmin);
        p_ele.y = (prm.Ymin+prm.Ymax)/2 + (prm.Ymax-prm.Ymin)*sqrt(-log(rz))*cos(2.0*PI*ry);
    }
    
    // shifted-Maxwellian for velocity 
    float rv1 = rand01(rng);
    float rv2 = rand01(rng);
    p_ele.vx = prm.ele_genU[0] + prm.ele_vth*sqrt(-log(rv1))*cos(2.0*PI*rv2);
    p_ele.vy = prm.ele_genU[1] + prm.ele_vth*sqrt(-log(rv1))*sin(2.0*PI*rv2);
    rv1 = rand01(rng);
    rv2 = rand01(rng);
    p_ele.vz = prm.ele_genU[2] + prm.ele_vth*sqrt(-log(rv1))*cos(2.0*PI*rv2);
    
    // shift from real coordinate to integer coodinate
    p_ele.x = p_ele.x/prm.dx;
    p_ele.y = p_ele.y/prm.dy;
    
    // shift origin for high-order weighting
    p_ele.x = p_ele.x + (float)prm.ngb;
    p_ele.y = p_ele.y + (float)prm.ngb;
    
    // deletion flag
    p_ele.piflag = 0;

    // debug print
    // if (sqrt(p_ele.vx*p_ele.vx+p_ele.vy*p_ele.vy+p_ele.vz*p_ele.vz) > 3e10) {
    if (!isfinite(p_ele.x)) {
        print[id] = 1.0;
    }else{
        print[id] = 0.0;
    }
    
    // escape for hollow-cathode
    if (genType == 3){
        states[id] = rng;
        return;
    }else{
    // add ion
        device ParticleState &p_ion = ptcl_ion[prm.ion_pNum + id];
        
        // position
        p_ion.x = p_ele.x;
        p_ion.y = p_ele.y;
        
        // shift particle position
        // p_ion.x -=  0.001/prm.dx;
        // if (p_ion.x > 90.0 && p_ion.x < 110.0){
        //     p_ion.x -= 0.001/prm.dx;
        // }
        
        // Maxwellian(Box-Muller)
        rv1 = rand01(rng);
        rv2 = rand01(rng);
        p_ion.vx = prm.ion_genU[0] + prm.ion_vth*sqrt(-log(rv1))*cos(2.0*PI*rv2);
        p_ion.vy = prm.ion_genU[1] + prm.ion_vth*sqrt(-log(rv1))*sin(2.0*PI*rv2);
        rv1 = rand01(rng);
        rv2 = rand01(rng);
        p_ion.vz = prm.ion_genU[2] + prm.ion_vth*sqrt(-log(rv1))*cos(2.0*PI*rv2);
        
        // deletion flag
        p_ion.piflag = 0;
        
        // update rnd state
        states[id] = rng;
    }
}
)";

@implementation Collision {
    size_t _initializedStates;
    uint64_t _masterSeed;
    float _srcValue;
    NSString* _srcName;
    int _srcValueForHC;
    std::default_random_engine _rndEngine;
    bool _existAI;
    bool _existHC;
}

// 初期設定
- (instancetype)initWithDevice:(id<MTLDevice>)device withParam:(Init*)initParam withParticles:(NSArray<Particle*>*)ptclArr withLogger:(XmlLogger&)logger{
    self = [super init];
    if (self) {

        // パラメータ取得
        struct ParamForField fieldParam = initParam.paramForField;
        struct ParamForComputing compParam = initParam.paramForComputing;
        std::vector<struct ParamForParticle> particleParam = initParam.paramForParticle;
        const std::vector<struct SourceForParticle> sources = initParam.particleSources;
    
        // 並列計算パラメータ
        _integrationChunkSize = compParam.integrationChunkSize;
        _threadGroupSize = compParam.threadGroupSize;
        // 粒子数制限
        _pNumMax = compParam.pNumMax;

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
        // 初期粒子生成
        int genType;
        if ([particleParam[0].genType isEqualToString:@"uniform-Gaussian"]){
            genType = 0;
        }else if([particleParam[0].genType isEqualToString:@"Xsinusoidal-Gaussian"]){
            genType = 1;
        }else if([particleParam[0].genType isEqualToString:@"Ysinusoidal-Gaussian"]){
            genType = 2;
        }else if([particleParam[0].genType isEqualToString:@"quietStart"]){
            genType = 4;
        }else if([particleParam[0].genType isEqualToString:@"Xgaussian-Gaussian"]){
            genType = 5;
        }else if([particleParam[0].genType isEqualToString:@"Ygaussian-Gaussian"]){
            genType = 6;
        }
        [fc setConstantValue:&genType type:MTLDataTypeInt atIndex:0];
        id<MTLFunction> function = [library newFunctionWithName:@"artificialIonization" constantValues:fc error:&error];
        _artificialIonizationPipeline = [device newComputePipelineStateWithFunction:function error:&error];

        // パラメータバッファ
        _paramsBuffer = [device newBufferWithLength:sizeof(generationParams) options:MTLResourceStorageModeShared];
        generationParams* prm = (generationParams*)[_paramsBuffer contents];
        // パラメータ格納
        prm->addn = particleParam[0].pNum;
        prm->ele_pNum = 0;
        prm->ion_pNum = 0;
        if (particleParam[0].genX[0] < 0){
            prm->Xmin = 0.0;
            prm->Xmax = fieldParam.dx*fieldParam.ngx;
        }else{
            prm->Xmin = particleParam[0].genX[0];
            prm->Xmax = particleParam[0].genX[1];
        }
        if (particleParam[0].genY[0] < 0){
            prm->Ymin = 0.0;
            prm->Ymax = fieldParam.dy*fieldParam.ngy;
        }else{
            prm->Ymin = particleParam[0].genY[0];
            prm->Ymax = particleParam[0].genY[1];
        }
        for (int s = 0; s < particleParam.size(); s++){
            if([particleParam[s].pName isEqualToString:@"electron"]){
                prm->ele_genU[0] = (float)particleParam[s].genU[0];
                prm->ele_genU[1] = (float)particleParam[s].genU[1];
                prm->ele_genU[2] = (float)particleParam[s].genU[2];
                prm->ele_vth = sqrt(2*kb*particleParam[s].genT/particleParam[s].m);
            }else{
                prm->ion_genU[0] = (float)particleParam[s].genU[0];
                prm->ion_genU[1] = (float)particleParam[s].genU[1];
                prm->ion_genU[2] = (float)particleParam[s].genU[2];
                prm->ion_vth = sqrt(2*kb*particleParam[s].genT/particleParam[s].m);
            }
        }
        prm->dx = fieldParam.dx;
        prm->dy = fieldParam.dy;
        prm->ngb = fieldParam.ngb;

        // デバッグ出力用バッファ
        _printBuffer = [device newBufferWithLength:sizeof(float)*compParam.pNumMax options:MTLResourceStorageModeShared];
        
        // 乱数状態
        _RNGStateBuffer = [device newBufferWithLength:sizeof(RNGState)*compParam.pNumMax options:MTLResourceStorageModeShared];
        RNGState *states = (RNGState*)_RNGStateBuffer.contents;

        // 固定シード or random_device()
        _masterSeed = 0ULL;

        // 未初期化状態を splitmix64 で初期化
        for (size_t i = 0; i < particleParam[0].pNum; ++i) {
            uint64_t s = splitmix64(_masterSeed);     // ① state 用
            uint64_t c = splitmix64(_masterSeed);     // ② inc 用
            states[i].state = s;
            states[i].inc   = (c | 1ULL);            // ビットOR演算で奇数化
        }
        // update count
        _initializedStates = particleParam[0].pNum;

        // CPU 側乱数
        std::random_device seed_gen;
        std::default_random_engine _rndEngine(seed_gen());

        // 粒子バッファ取得
        id<MTLBuffer> eleBuffer;
        id<MTLBuffer> ionBuffer;
        SimulationParams* eleParams;
        SimulationParams* ionParams;
        for (Particle* ptcl in ptclArr){
            if([ptcl.pName isEqualToString:@"electron"]){
                eleBuffer = [ptcl particleBuffer];
                eleParams = (SimulationParams*)[[ptcl paramsBuffer] contents];
            }else{
                ionBuffer = [ptcl particleBuffer];
                ionParams = (SimulationParams*)[[ptcl paramsBuffer] contents];
            }
        }
    
        // 粒子初期化
        id<MTLCommandBuffer> commandBuffer = [_commandQueue commandBuffer];
        id<MTLComputeCommandEncoder> computeEncoder = [commandBuffer computeCommandEncoder];
        [computeEncoder setComputePipelineState:_artificialIonizationPipeline];
        [computeEncoder setBuffer:eleBuffer             offset:0 atIndex:0];
        [computeEncoder setBuffer:ionBuffer             offset:0 atIndex:1];
        [computeEncoder setBuffer:_paramsBuffer         offset:0 atIndex:2];
        [computeEncoder setBuffer:_RNGStateBuffer       offset:0 atIndex:3];
        [computeEncoder setBuffer:_printBuffer          offset:0 atIndex:4];
        // グリッドとスレッドグループのサイズ設定
        uint threadGroupNum = (prm->addn + _threadGroupSize - 1) / _threadGroupSize;
        MTLSize gridSizeMetalStyle = MTLSizeMake(threadGroupNum, 1, 1);
        MTLSize threadGroupSize = MTLSizeMake(_threadGroupSize, 1, 1);
        // ディスパッチ/エンコーディング/実行
        [computeEncoder dispatchThreadgroups:gridSizeMetalStyle threadsPerThreadgroup:threadGroupSize];
        [computeEncoder endEncoding];
        [commandBuffer commit];
        [commandBuffer waitUntilCompleted];

        // 粒子数更新
        eleParams->pNum = prm->addn;
        ionParams->pNum = prm->addn;

        // debug print
        ParticleState* p;
        for (Particle* ptcl in ptclArr) {
            // electron の状態を取得
            if([ptcl.pName isEqualToString:@"electron"]){
                p = (ParticleState*)[[ptcl particleBuffer] contents];
                break;
            }
        }
        float* prt = (float*)_printBuffer.contents;
        for (int i = 0; i < prm->addn; i++){
            if (prt[i] > 1e-20){
                NSLog(@"[AI] detected not finite: prt[%d] = %e",i,prt[i]);
                // NSLog(@"[AI] detected over light speed: prt[%d] = %e",i,prt[i]);
            }
        }
        // // quiet start
        // int addn_x = (int)round(sqrt(prm->addn*(prm->Xmax-prm->Xmin)/(prm->Ymax-prm->Ymin))); // It is desirable that ppc be a perfect square.
        // int addn_y = prm->addn/addn_x;
        // float delta_x = (prm->Xmax-prm->Xmin)/addn_x;
        // float delta_y = (prm->Ymax-prm->Ymin)/addn_y;
        // for (int i = 0; i < prm->addn; i+=1){
        //     NSLog(@"[AI] p[pNum+%d].x = %e, p.x(@cpu) = %e, p[n].y = %e, p[n].vx = %e, p[n].vy = %e, p[n].vz = %e",i,p[i].x,(float)((prm->Xmin + delta_x*(0.5+i%addn_x))/prm->dx),p[i].y,p[i].vx,p[i].vy,p[i].vz);
        // }
        // NSLog(@"[AI] addn = %d, addn_x = %d, addn_y = %d, delta_x = %e, delta_y = %e",prm->addn,addn_x,addn_y,delta_x,delta_y);
        // NSLog(@"[AI] delX_gen = %e, delY_gen = %e",prm->Xmax-prm->Xmin, prm->Ymax-prm->Ymin);

        // 生成条件のカーネルをセット（元バージョンではそれぞれ1条件まで対応，2条件以上は pipeline を配列で管理する必要があるため後回し．）
        _existAI = false;
        _existHC = false;
        for (int i = 0; i < sources.size(); i++){
            if ([sources[i].genType isEqualToString:@"uniformIonization"]){
                genType = 0;
            }else if([sources[i].genType isEqualToString:@"XsinusoidalIonization"]){
                genType = 1;
            }else if([sources[i].genType isEqualToString:@"YsinusoidalIonization"]){
                genType = 2;
            }else if([sources[i].genType isEqualToString:@"hollow-cathode"]){
                genType = 3;
            }
            // artificial ionization
            if (genType == 0 || genType == 1 || genType == 2){
                // set exist flag
                _existAI = true;
                // set function constant
                [fc setConstantValue:&genType type:MTLDataTypeInt atIndex:0];
                function = [library newFunctionWithName:@"artificialIonization" constantValues:fc error:&error];
                _artificialIonizationPipeline = [device newComputePipelineStateWithFunction:function error:&error];
                
                // set param
                _srcName = sources[i].pName;
                _srcValue = sources[i].src;
                if (sources[i].genX[0] < 0){
                    prm->Xmin = 0.0;
                    prm->Xmax = fieldParam.dx*fieldParam.ngx;
                }else{
                    prm->Xmin = sources[i].genX[0];
                    prm->Xmax = sources[i].genX[1];
                }
                if (sources[i].genY[0] < 0){
                    prm->Ymin = 0.0;
                    prm->Ymax = fieldParam.dy*fieldParam.ngy;
                }else{
                    prm->Ymin = sources[i].genY[0];
                    prm->Ymax = sources[i].genY[1];
                }
                prm->ele_genU[0] = (float)sources[i].genU_ele[0];
                prm->ele_genU[1] = (float)sources[i].genU_ele[1];
                prm->ele_genU[2] = (float)sources[i].genU_ele[2];
                prm->ele_vth = sqrt(2*kb*sources[i].genT_ele/me);
                prm->ion_genU[0] = (float)sources[i].genU[0];
                prm->ion_genU[1] = (float)sources[i].genU[1];
                prm->ion_genU[2] = (float)sources[i].genU[2];
                for (int s = 0; s < particleParam.size(); s++){
                    if([particleParam[s].pName isEqualToString:sources[i].pName]){
                        prm->ion_vth = sqrt(2*kb*sources[i].genT/particleParam[s].m);
                    }
                }
            }
            // hollow-cathode
            if (genType == 3){
                // set exist flag
                _existHC = true;
                // set function constant
                [fc setConstantValue:&genType type:MTLDataTypeInt atIndex:0];
                function = [library newFunctionWithName:@"artificialIonization" constantValues:fc error:&error];
                _hollowCathodePipeline = [device newComputePipelineStateWithFunction:function error:&error];
                
                // params for hollowcathode
                _paramsForHCBuffer = [device newBufferWithLength:sizeof(generationParams) options:MTLResourceStorageModeShared];
                generationParams* prmHC = (generationParams*)[_paramsForHCBuffer contents];
        
                // set param
                _srcValueForHC = (int)sources[i].src; // 0: Xmin, 1: Xmax, 2: Ymin, 3: Ymax
                if (sources[i].genX[0] < 0){
                    prmHC->Xmin = 0.0;
                    prmHC->Xmax = fieldParam.dx*fieldParam.ngx;
                }else{
                    prmHC->Xmin = sources[i].genX[0];
                    prmHC->Xmax = sources[i].genX[1];
                }
                if (sources[i].genY[0] < 0){
                    prmHC->Ymin = 0.0;
                    prmHC->Ymax = fieldParam.dy*fieldParam.ngy;
                }else{
                    prmHC->Ymin = sources[i].genY[0];
                    prmHC->Ymax = sources[i].genY[1];
                }
                prmHC->ele_genU[0] = (float)sources[i].genU[0];
                prmHC->ele_genU[1] = (float)sources[i].genU[1];
                prmHC->ele_genU[2] = (float)sources[i].genU[2];
                prmHC->ele_vth = sqrt(2*kb*sources[i].genT/me);
                prmHC->ion_genU[0] = 0.0;
                prmHC->ion_genU[1] = 0.0;
                prmHC->ion_genU[2] = 0.0;
                prmHC->ion_vth = 0.0;
                prmHC->dx = fieldParam.dx;
                prmHC->dy = fieldParam.dy;
                prmHC->ngb = fieldParam.ngb;
            }
        }
        
    }
    return self;
}

- (bool)artificialIonization:(double)dt withParticles:(NSArray<Particle*>*)ptclArr withLogger:(XmlLogger&)logger{
        // logging
        std::map<std::string, std::string>loggingData;
        // check flag
        if(!_existAI){ return true; }

        // 粒子状態の取得とパラメータ更新
        generationParams* prm = (generationParams*)[_paramsBuffer contents];
        id<MTLBuffer> eleBuffer;
        id<MTLBuffer> ionBuffer;
        SimulationParams* eleParams;
        SimulationParams* ionParams;
        std::uniform_real_distribution<> unif_dist(0.0, 1.0);
        for (Particle* ptcl in ptclArr) {
            // 電子バッファを格納
            if([ptcl.pName isEqualToString:@"electron"]){
                eleBuffer = ptcl.particleBuffer;
                eleParams = (SimulationParams*)[[ptcl paramsBuffer] contents];
                double addn_d = _srcValue/ptcl.w*dt;
                if(unif_dist(_rndEngine) < addn_d - (int)addn_d){
                    prm->addn = (int)addn_d + 1;
                }else{
                    prm->addn = (int)addn_d;
                }
                prm->ele_pNum = eleParams->pNum;
            // イオンバッファを格納
            }else if([ptcl.pName isEqualToString:_srcName]){
                ionBuffer = ptcl.particleBuffer;
                ionParams = (SimulationParams*)[[ptcl paramsBuffer] contents];
                prm->ion_pNum = ionParams->pNum;
            }
        }

        // logging
        loggingData["addn_injected"] = std::to_string(prm->addn);
        logger.logSection("artificialIonization", loggingData);

        // check pNum
        if (eleParams->pNum+prm->addn > _pNumMax || ionParams->pNum+prm->addn > _pNumMax){
            NSLog(@"[FatalError] pNum+addn > pNumMax");
            return false;
        }

        // 乱数状態を取得
        RNGState* states = (RNGState*)[_RNGStateBuffer contents];
        // 未初期化状態があれば splitmix64 で初期化
        if(prm->addn > _initializedStates){
            for (size_t i = _initializedStates; i < prm->addn; ++i) {
                uint64_t s = splitmix64(_masterSeed);     // ① state 用
                uint64_t c = splitmix64(_masterSeed);     // ② inc 用
                states[i].state = s;
                states[i].inc   = (c | 1ULL);            // ビットOR演算で奇数化
                NSLog(@"[AI] states[%zu].state = %llu, states[%zu].state = %llu",i,s,i,c);
            }
            // update count
            _initializedStates = prm->addn;
        }

        // 粒子初期化
        id<MTLCommandBuffer> commandBuffer = [_commandQueue commandBuffer];
        id<MTLComputeCommandEncoder> computeEncoder = [commandBuffer computeCommandEncoder];
        [computeEncoder setComputePipelineState:_artificialIonizationPipeline];
        [computeEncoder setBuffer:eleBuffer             offset:0 atIndex:0];
        [computeEncoder setBuffer:ionBuffer             offset:0 atIndex:1];
        [computeEncoder setBuffer:_paramsBuffer         offset:0 atIndex:2];
        [computeEncoder setBuffer:_RNGStateBuffer       offset:0 atIndex:3];
        [computeEncoder setBuffer:_printBuffer          offset:0 atIndex:4];
        // グリッドとスレッドグループのサイズ設定
        uint threadGroupNum = (prm->addn + _threadGroupSize - 1) / _threadGroupSize;
        MTLSize gridSizeMetalStyle = MTLSizeMake(threadGroupNum, 1, 1);
        MTLSize threadGroupSize = MTLSizeMake(_threadGroupSize, 1, 1);
        // ディスパッチ/エンコーディング/実行
        [computeEncoder dispatchThreadgroups:gridSizeMetalStyle threadsPerThreadgroup:threadGroupSize];
        [computeEncoder endEncoding];
        [commandBuffer commit];
        [commandBuffer waitUntilCompleted];

        // 粒子数更新
        eleParams->pNum += prm->addn;
        ionParams->pNum += prm->addn;
        // NSLog(@"[AI] addn = %d, ele_pNum = %d, ion_pNum = %d",prm->addn,eleParams->pNum,ionParams->pNum);

        // debug print
        ParticleState* p;
        for (Particle* ptcl in ptclArr) {
            // electron の状態を取得
            if([ptcl.pName isEqualToString:@"electron"]){
                p = (ParticleState*)[[ptcl particleBuffer] contents];
                break;
            }
        }
        float* prt = (float*)_printBuffer.contents;
        for (int i = 0; i < prm->addn; i++){
            if (prt[i] > 1e-20){
                NSLog(@"[AI] detected not finite: prt[%d] = %e",i,prt[i]);
                // NSLog(@"[AI] detected over light speed: prt[%d] = %e",i,prt[i]);
            }
            int n = eleParams->pNum-prm->addn+i;
            // NSLog(@"[AI] p[pNum+%d].x = %e, p[n].y = %e, p[n].vx = %e, p[n].vy = %e, p[n].vz = %e",i,p[n].x,p[n].y,p[n].vx,p[n].vy,p[n].vz);
        }

        return true;
}

- (bool)hollowCathode:(double)dt withParticles:(NSArray<Particle*>*)ptclArr withLogger:(XmlLogger&)logger{
        // ログ出力用データセット
        std::map<std::string, std::string>loggingData;
        // check flag
        if(!_existHC){ return true; }
        
        // set params
        generationParams* prm = (generationParams*)[_paramsForHCBuffer contents];
        int addn = 0; // initialize
        id<MTLBuffer> eleBuffer;
        id<MTLBuffer> ionBuffer;
        SimulationParams* eleParams;
        SimulationParams* ionParams;
        for (Particle* ptcl in ptclArr) {
            if([ptcl.pName isEqualToString:@"electron"]){
                eleBuffer = ptcl.particleBuffer;
                eleParams = (SimulationParams*)[[ptcl paramsBuffer] contents];
                prm->ele_pNum = eleParams->pNum;
            }else{
                // ionization とカーネルを共有しているので取得しておく
                ionBuffer = ptcl.particleBuffer;
                ionParams = (SimulationParams*)[[ptcl paramsBuffer] contents];
                prm->ion_pNum = ionParams->pNum;
            }
            // accumlate boundary current
            if(_srcValueForHC == 0){
                addn -= ptcl.pinum_Xmin;
            }else if(_srcValueForHC == 1){
                addn -= ptcl.pinum_Xmax;
            }else if(_srcValueForHC == 2){
                addn -= ptcl.pinum_Ymin;
            }else if(_srcValueForHC == 3){
                addn -= ptcl.pinum_Ymax;
            }
        }

        // check addn
        if (eleParams->pNum+addn > _pNumMax){
            NSLog(@"[HC] exceeded maximum pNum");
            return false;
        }else if (addn <= 0){
            // neglect negative current
            loggingData["addn_neglected"] = std::to_string(-addn);
            loggingData["addn_injected"] = std::to_string(0);
            logger.logSection("hollowCathode", loggingData);
            return true;
        }else{
            // assign addn to prm
            loggingData["addn_neglected"] = std::to_string(0);
            loggingData["addn_injected"] = std::to_string(addn);
            logger.logSection("hollowCathode", loggingData);
            prm->addn = (uint)addn;
        }

        // 未初期化状態があれば splitmix64 で初期化
        if(prm->addn > _initializedStates){
            RNGState *states = (RNGState*)_RNGStateBuffer.contents;
            for (size_t i = _initializedStates; i < prm->addn; ++i) {
                uint64_t s = splitmix64(_masterSeed);     // ① state 用
                uint64_t c = splitmix64(_masterSeed);     // ② inc 用
                states[i].state = s;
                states[i].inc   = (c | 1ULL);            // ビットOR演算で奇数化
            }
            // update count
            _initializedStates = prm->addn;
        }

        // 粒子初期化
        id<MTLCommandBuffer> commandBuffer = [_commandQueue commandBuffer];
        id<MTLComputeCommandEncoder> computeEncoder = [commandBuffer computeCommandEncoder];
        [computeEncoder setComputePipelineState:_hollowCathodePipeline];
        [computeEncoder setBuffer:eleBuffer             offset:0 atIndex:0];
        [computeEncoder setBuffer:ionBuffer             offset:0 atIndex:1];
        [computeEncoder setBuffer:_paramsForHCBuffer    offset:0 atIndex:2];
        [computeEncoder setBuffer:_RNGStateBuffer       offset:0 atIndex:3];
        [computeEncoder setBuffer:_printBuffer          offset:0 atIndex:4];
        // グリッドとスレッドグループのサイズ設定
        uint threadGroupNum = (prm->addn + _threadGroupSize - 1) / _threadGroupSize;
        MTLSize gridSizeMetalStyle = MTLSizeMake(threadGroupNum, 1, 1);
        MTLSize threadGroupSize = MTLSizeMake(_threadGroupSize, 1, 1);
        // ディスパッチ/エンコーディング/実行
        [computeEncoder dispatchThreadgroups:gridSizeMetalStyle threadsPerThreadgroup:threadGroupSize];
        [computeEncoder endEncoding];
        [commandBuffer commit];
        [commandBuffer waitUntilCompleted];

        // 粒子数更新
        eleParams->pNum += prm->addn;
        // NSLog(@"[HC] addn = %d, pNum = %d",prm->addn,eleParams->pNum);

        // debug print
        ParticleState* p;
        for (Particle* ptcl in ptclArr) {
            // electron の状態を取得
            if([ptcl.pName isEqualToString:@"electron"]){
                p = (ParticleState*)[[ptcl particleBuffer] contents];
                break;
            }
        }
        float* prt = (float*)_printBuffer.contents;
        for (int i = 0; i < prm->addn; i++){
            if (prt[i] > 1e-20){
                NSLog(@"[HC] detected not finite: prt[%d] = %e",i,prt[i]);
                // NSLog(@"[AI] detected over light speed: prt[%d] = %e",i,prt[i]);
            }
            int n = eleParams->pNum-prm->addn+i;
            // NSLog(@"[HC] p[pNum+%d].x = %e, p[n].y = %e, p[n].vx = %e, p[n].vy = %e, p[n].vz = %e",i,p[n].x,p[n].y,p[n].vx,p[n].vy,p[n].vz);
        }

        return true;
}

@end