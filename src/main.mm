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
        Init *init = [[Init alloc] parseInputFile:inputFilePath];

        // パース結果の出力(debugprint.mm)
        printInitContents(init);

        struct FragForEquation EqFrags = [init getFragForEquation];
        
        // 粒子の初期化
        NSMutableArray *ptclArr = [NSMutableArray arrayWithCapacity:EqFrags.Particle];
        for (int s = 0; s < EqFrags.Particle; s++) {
            Particle *ptcl = [[Particle alloc] initWithDevice:device withParam:init specimen:s];
            [ptclArr addObject:ptcl];
        }
        // 場の初期化
        EMField *fld = [[EMField alloc] initWithDevice:device withParam:init];
        // モーメント量の初期化
        Moment *mom = [[Moment alloc] initialize];

        // 時間更新ループ
        struct ParamForTimeIntegration timeParams = [init getParamForTimeIntegration];
        int StartCycle = 1; // リスタート時は最終サイクルを引き継ぎたい
        double dt = timeParams.dt;
        for (int cycle = StartCycle; cycle <= timeParams.EndCycle; cycle++) {
            // 電荷密度の初期化
            [fld resetChargeDensity];
            // 粒子ループ
            for (int s = 0; s < EqFrags.Particle; s++) {
                Particle *ptcl = [ptclArr objectAtIndex:s];
                // 粒子の時間更新
                [ptcl update:dt withEMField:fld];
                // 流出粒子の処理
                [ptcl reduce];
                // 電荷密度の更新
                if (EqFrags.EMField == 1){
                    [ptcl integrateChargeDensity:fld];
                }
                // 粒子軌道の出力
                if (timeParams.ptclOutCycle != 0 && cycle%timeParams.ptclOutCycle == 0){
                    [ptcl outputPhaseSpace:cycle withEMField:fld];
                }
            }
            // 電場の更新
            if (EqFrags.EMField == 1){
                [fld solvePoisson];
                // 場の出力
                if (timeParams.fldOutCycle != 0 && cycle%timeParams.fldOutCycle == 0){
                    [fld outputField: cycle];
                }
            }
            // 粒子生成
            for (int s = 0; s < EqFrags.Particle; s++) {
                Particle *ptcl = [ptclArr objectAtIndex:s];
                [ptcl injection:dt withParam:init];
            }

            NSLog(@"Frame %d completed", cycle);
        }
    }
    return 0;
}