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
#import "XmlLogger.h"

#define MEASURE(name, code, data)                                  \
    do {                                                           \
        auto _start = std::chrono::high_resolution_clock::now();  \
        code;                                                     \
        auto _end = std::chrono::high_resolution_clock::now();    \
        auto _time = std::chrono::duration_cast<std::chrono::microseconds>(_end - _start).count(); \
        data[name] = std::to_string(_time);\
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
        if (args.count("-i") == 1 && args.count("-o") == 1) {
            inputPath = [NSString stringWithUTF8String:args["-i"].c_str()];
            // ファイルの存在確認
            if (![[NSFileManager defaultManager] fileExistsAtPath:inputPath]) {
                NSLog(@"%@ is not found.", inputPath);
                return 1;
            }
        }else{
            NSLog(@"Usage: %s -i inputfile.txt -o outputfile.xml", argv[0]);
            return 1;
        }

        // logger の作成
        XmlLogger logger(args["-o"].c_str());
        // 所用時間格納
        std::map<std::string,std::string> dataElapsedTime;

        // 初期化パラメータクラス作成
        Init *init = [[Init alloc] parseInputFile:inputPath];
        // 入力チェック、変数演算
        bool check = [init checkInput];
        if (!check){
            NSLog(@"invalid condition.");
            return 1;
        }
        // パース結果の出力(debugprint.mm)
        printInitContents(init, logger);

        struct FragForEquation EqFrags = [init getFragForEquation];
        
        // 粒子の初期化
        NSMutableArray *ptclArr = [NSMutableArray arrayWithCapacity:EqFrags.Particle];
        for (int s = 0; s < EqFrags.Particle; s++) {
            Particle *ptcl = [[Particle alloc] initWithDevice:device withParam:init specimen:s withLogger:logger];
            [ptclArr addObject:ptcl];
        }
        // 場の初期化
        EMField *fld = [[EMField alloc] initWithDevice:device withParam:init withLogger:logger];
        // モーメント量の初期化
        Moment *mom = [[Moment alloc] initialize];

        // 時間更新ループ
        struct ParamForTimeIntegration timeParams = [init getParamForTimeIntegration];
        int StartCycle = 1; // リスタート時は最終サイクルを引き継ぎたい
        double dt = timeParams.dt;
        double time = 0.0;
        double comp = 0.0; // 保証項
        double y,t;
        int intCurrent = 0;
        for (int cycle = StartCycle; cycle <= timeParams.EndCycle; cycle++) {
            // 時間計算
            y = dt - comp;
            t = time + y;
            comp = (t - time) - y;
            time = t;
            logger.logCycleStart(cycle, time);

            // 電荷密度の初期化
            [fld resetChargeDensity];

            // 粒子ループ
            for (int s = 0; s < EqFrags.Particle; s++) {
                Particle *ptcl = [ptclArr objectAtIndex:s];
                std::string pName = [ptcl.pName UTF8String];
                // 粒子の時間更新
                MEASURE("update_"+pName, [ptcl update:dt withEMField:fld withLogger:logger], dataElapsedTime);
                // 流出粒子の処理
                MEASURE("reduce_"+pName, [ptcl reduce:logger], dataElapsedTime);
                intCurrent += ptcl.pinum_Xmin;
                std::map<std::string, std::string> data ={
                    {"Xmin", std::to_string(ptcl.pinum_Xmin)},
                    {"Xmax", std::to_string(ptcl.pinum_Xmax)},
                    {"Ymin", std::to_string(ptcl.pinum_Ymin)},
                    {"Ymax", std::to_string(ptcl.pinum_Ymax)},
                };
                logger.logSection("flowout_"+pName, data);
                // 電荷密度の更新
                if (EqFrags.EMField == 1){
                    MEASURE("integCDens_"+pName, [ptcl integrateChargeDensity:fld withLogger:logger], dataElapsedTime);
                }
                // 粒子軌道の出力
                if (timeParams.ptclOutCycle != 0 && cycle%timeParams.ptclOutCycle == 0){
                    [ptcl outputPhaseSpace:cycle withEMField:fld withLogger:logger];
                }
            }

            // 電場の更新
            if (EqFrags.EMField == 1){
                MEASURE("solvePoisson", [fld solvePoisson:logger], dataElapsedTime);
                // 場の出力
                if (timeParams.fldOutCycle != 0 && cycle%timeParams.fldOutCycle == 0){
                    [fld outputField:cycle withLogger:logger];
                }
            }

            // 粒子生成
            std::vector<int> ret(EqFrags.Particle);
            for (int s = 0; s < EqFrags.Particle; s++) {
                Particle *ptcl = [ptclArr objectAtIndex:s];
                std::string pName = [ptcl.pName UTF8String];
                MEASURE("injection_"+pName, ret.push_back([ptcl injection:dt withParam:init withCurrent:intCurrent withLogger:logger]), dataElapsedTime);
            }
            if (std::reduce(std::begin(ret), std::end(ret)) != 0){
                // すべての injection が成功したら総和は0
                NSLog(@"injection failed.");
                return 1;
            }

            NSLog(@"Frame %d completed", cycle);
            logger.logSection("elapsedTime", dataElapsedTime);
            logger.logCycleEnd();
        }
    }
    return 0;
}