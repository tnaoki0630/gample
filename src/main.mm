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
        Particle *ptcl = [[Particle alloc] 
                            initWithDevice:device
                            initWithParam:init];
        // 場の初期化
        EMField *fld = [[EMField alloc]
                            initWithParam:init];
        // モーメント量の初期化
        Moment *mom = [[Moment alloc] initialize];

        // シミュレーションループ
        double dt = 1e-3;
        int maxCycle = 100;
        for (int i = 0; i < maxCycle; i++) {
            // 粒子の時間更新
            [ptcl update:dt withEMField:fld];
            // 電荷密度の更新
            [fld culcChargeDensity:ptcl];
            // 電場の更新
            [fld solvePoisson];
            NSLog(@"Frame %d completed", i);
        }
    }
    return 0;
}