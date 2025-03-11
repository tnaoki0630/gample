#import <Foundation/Foundation.h>
#import "Init.h"
#import "Particle.h"
#import "Moment.h"
#import "EMField.h"
#import "DebugPrint.h"

int main(int argc, const char * argv[]) {
    @autoreleasepool {
        // デフォルトのMetalデバイスを取得
        id<MTLDevice> device = MTLCreateSystemDefaultDevice();
        if (!device) {
            NSLog(@"Metal is not supported on this device");
            return 1;
        }
        
        // 初期化用クラスの設定
        NSString *inputFilePath = @"data/condition.txt";
        Init *init = [[Init alloc]
                            parseInputFile:inputFilePath];

        // パース結果の出力(debugprint.mm)
        printInitContents(init);

        // 粒子の初期化
        NSArray *ParticleParams = [init getParamForParticle];
        NSMutableArray *ptclArr = [NSMutableArray arrayWithCapacity:ParticleParams.count];
        struct ParamForField FieldParam = [init getParamForField];
        for (int s = 0; s < ParticleParams.count; s++) {
            NSValue *value = ParticleParams[s];
            struct ParamForParticle ParticleParam;
            [value getValue:&ParticleParam];
            NSLog(@"initParticles: %@", ParticleParam.pName);
            Particle *ptcl = [[Particle alloc] initWithDevice:device
                                                withParticleParam:ParticleParam
                                                withFieldParam:FieldParam];
            [ptclArr addObject:ptcl];
        }
        // 場の初期化
        EMField *fld = [[EMField alloc] initWithDevice:device withParam:init];
        // モーメント量の初期化
        Moment *mom = [[Moment alloc] initialize];

        // シミュレーションループ
        struct ParamForTimeIntegration timeParams = [init getParamForTimeIntegration];
        int StartCycle = 0; // リスタート時は最終サイクルを引き継ぎたい
        double dt = timeParams.dt; // そのうち dt の更新が必要になるかも
        for (int i = StartCycle; i < timeParams.EndCycle; i++) {
            for (int s = 0; s < ParticleParams.count; s++) {
                Particle *ptcl = [ptclArr objectAtIndex:s];
                // 粒子の時間更新
                [ptcl update:dt withEMField:fld];
                // 電荷密度の更新
                [ptcl integrateChargeDensity:fld];
                // 粒子軌道の出力
                if (i%timeParams.OutputCycle == 0){
                    NSValue *value = ParticleParams[s];
                    struct ParamForParticle ParticleParam;
                    [value getValue:&ParticleParam];
                    [ptcl outputPhaseSpace:i withParticleParam:ParticleParam];
                }
            }
            // 電場の更新
            [fld solvePoisson];
            NSLog(@"Frame %d completed", i);
        }
    }
    return 0;
}