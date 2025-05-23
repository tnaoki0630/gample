#include <chrono>
#include <iostream>
#include <numeric>
#import <map>
#import <string>
#import <Foundation/Foundation.h>
#import "Init.h"
#import "Particle.h"
#import "Moment.h"
#import "EMField.h"
#import "DebugPrint.h"

#define MEASURE(name, code)                                       \
    do {                                                           \
        auto _start = std::chrono::high_resolution_clock::now();  \
        code;                                                     \
        auto _end = std::chrono::high_resolution_clock::now();    \
        auto _time = std::chrono::duration_cast<std::chrono::microseconds>(_end - _start).count(); \
        NSLog(@"%s: %lld us", name, _time);                           \
    } while (0)

std::map<std::string, std::string> parseArgs(int argc, const char* argv[]) {
    std::map<std::string, std::string> options;
    for (int i = 1; i < argc - 1; ++i) {
        std::string key = argv[i];
        if (key[0] == '-') {
            std::string value = argv[i + 1];
            options[key] = value;
            ++i; // 値をスキップ
        }
    }
    return options;
}

int main(int argc, const char * argv[]) {
    @autoreleasepool {
        // デフォルトのMetalデバイスを取得
        id<MTLDevice> device = MTLCreateSystemDefaultDevice();
        if (!device) {
            NSLog(@"Metal is not supported on this device");
            return 1;
        }

        // I/O チェック
        auto args = parseArgs(argc, argv);
        NSString* inputPath;
        if (args.count("-i") == 1) {
            inputPath = [NSString stringWithUTF8String:args["-i"].c_str()];
            // ファイルの存在確認
            if (![[NSFileManager defaultManager] fileExistsAtPath:inputPath]) {
                NSLog(@"%@ is not found.", inputPath);
                return 1;
            }
        }else{
            NSLog(@"Usage: %s -i inputfile.txt -o outputfile.txt", argv[0]);
            return 1;
        }
        if (args.count("-o") == 1) {
            // NSLog 出力先を変更
            freopen(args["-o"].c_str(), "w", stderr);
        }

        // 初期化用クラスの設定
        // NSString *inputPath = @"data/condition.txt";
        Init *init = [[Init alloc] parseInputFile:inputPath];
        // 入力チェック、変数演算
        bool check = [init checkInput];
        if (!check){
            NSLog(@"invalid condition.");
            return 1;
        }

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

        // 初期場
        [fld resetChargeDensity];
        MEASURE("solvePoisson", [fld solvePoisson]);
        [fld outputField: 0];

        // 時間更新ループ
        struct ParamForTimeIntegration timeParams = [init getParamForTimeIntegration];
        int StartCycle = 1; // リスタート時は最終サイクルを引き継ぎたい
        double dt = timeParams.dt;
        int int_current;
        for (int cycle = StartCycle; cycle <= timeParams.EndCycle; cycle++) {
            // 電荷密度の初期化
            [fld resetChargeDensity];
            // 粒子ループ
            int_current = 0;
            for (int s = 0; s < EqFrags.Particle; s++) {
                Particle *ptcl = [ptclArr objectAtIndex:s];
                // 粒子の時間更新
                MEASURE("update", [ptcl update:dt withEMField:fld]);
                // 流出粒子の処理
                MEASURE("reduce", [ptcl reduce]);
                int_current += ptcl.pinum_Xmin;
                NSLog(@"flowout(#%d), Xmin = %d, Xmax = %d", s, ptcl.pinum_Xmin, ptcl.pinum_Xmax);
                // 電荷密度の更新
                if (EqFrags.EMField == 1){
                    MEASURE("integrateChargeDensity", [ptcl integrateChargeDensity:fld]);
                }
                // 粒子軌道の出力
                if (timeParams.ptclOutCycle != 0 && cycle%timeParams.ptclOutCycle == 0){
                    [ptcl outputPhaseSpace:cycle withEMField:fld];
                }
            }
            // 電場の更新
            if (EqFrags.EMField == 1){
                MEASURE("solvePoisson", [fld solvePoisson]);
                // 場の出力
                if (timeParams.fldOutCycle != 0 && cycle%timeParams.fldOutCycle == 0){
                    [fld outputField: cycle];
                }
            }
            // 粒子生成
            NSLog(@"int_current = %d", int_current);
            std::vector<int> ret(EqFrags.Particle);
            for (int s = 0; s < EqFrags.Particle; s++) {
                Particle *ptcl = [ptclArr objectAtIndex:s];
                MEASURE("injection", ret.push_back([ptcl injection:dt withParam:init withCurrent:int_current]));
            }
            // すべての injection が成功したら総和は0
            if (std::reduce(std::begin(ret), std::end(ret)) != 0){
                NSLog(@"injection failed.");
                return 1;
            }

            std::cout << "Frame " << cycle << "completed" << std::endl;
            NSLog(@"Frame %d completed", cycle);
        }
    }
    return 0;
}