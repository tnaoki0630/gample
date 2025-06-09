#import <iostream>
#import <numeric>
#import <Foundation/Foundation.h>
#import <mach/mach.h>
#import "Init.h"
#import "Particle.h"
#import "Moment.h"
#import "EMField.h"
#import "DebugPrint.h"
#import "XmlLogger.h"
#import "ResourceUtils.h"
#import "FileIO.h"

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
        // 所要時間計測
        auto start = std::chrono::high_resolution_clock::now();
        
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
            NSLog(@"Usage: %s -i inputfile.json", argv[0]);
            return 1;
        }

        // 初期化パラメータクラス作成
        Init *init = [[Init alloc] parseInputFile:inputPath];

        // logger の作成
        struct ParamForTimeIntegration timeParam = init.paramForTimeIntegration;
        XmlLogger logger([[NSString stringWithFormat:@"%@_log.xml", timeParam.ProjectName] UTF8String]);

        // パース結果の出力(debugprint.mm)
        printInitContents(init, logger);

        struct FlagForEquation EqFlags = init.flagForEquation;

        // 粒子クラスの初期化
        NSMutableArray *ptclArr = [NSMutableArray arrayWithCapacity:EqFlags.Particle];
        for (int s = 0; s < EqFlags.Particle; s++) {
            Particle *ptcl = [[Particle alloc] initWithDevice:device withParam:init specimen:s withLogger:logger];
            [ptclArr addObject:ptcl];
        }
        // 電磁場クラスの初期化
        EMField *fld = [[EMField alloc] initWithDevice:device withParam:init withLogger:logger];
        // 初期場の計算
        [fld resetChargeDensity];
        [fld solvePoisson:logger];
        // モーメント計算クラスの初期化
        Moment *mom = [[Moment alloc] initWithDevice:device withParam:init withLogger:logger];

        // 時間更新パラメータ
        int StartCycle = timeParam.Start;
        int intCurrent = 0;
        // リスタート
        if(StartCycle != 0){
            if(!loadProgress(StartCycle, ptclArr, fld, intCurrent, init)){
                NSLog(@"restart failed.");
                return 1;
            };
        }
        double dt = timeParam.TimeStep;
        double time = StartCycle*dt;
        double comp = 0.0; // 保証項
        double y,t;
        // ループ開始
        for (int cycle = StartCycle+1; cycle <= timeParam.End; cycle++) {
            @autoreleasepool {
                // ログ出力のオンオフ
                if(cycle%timeParam.LogOutput == 0){ logger.swichLog(true); 
                }else{ logger.swichLog(false); }
                // 時間計算
                y = dt - comp;
                t = time + y;
                comp = (t - time) - y;
                time = t;
                logger.logCycleStart(cycle, time);
                
                // 所要リソース格納辞書
                std::map<std::string,std::string> dataElapsedTime;
                std::map<std::string,std::string> dataMemUsage;

                // 電荷密度の初期化
                [fld resetChargeDensity];

                // 粒子ループ
                std::vector<struct ParamForParticle> particles = init.paramForParticle;
                for (int s = 0; s < EqFlags.Particle; s++) {
                    Particle *ptcl = [ptclArr objectAtIndex:s];
                    std::string pName = [particles[s].pName UTF8String];
                    // 粒子の時間更新
                    MEASURE("update_"+pName, [ptcl update:dt withEMField:fld withMom:mom withLogger:logger], dataElapsedTime);
                    // 流出粒子の処理
                    MEASURE("reduce_"+pName, [ptcl reduce:logger], dataElapsedTime);
                    // 流出電流の計算
                    intCurrent += ptcl.pinum_Xmin;
                    std::map<std::string, std::string> data ={
                        {"Xmin", std::to_string(ptcl.pinum_Xmin)},
                        {"Xmax", std::to_string(ptcl.pinum_Xmax)},
                        {"Ymin", std::to_string(ptcl.pinum_Ymin)},
                        {"Ymax", std::to_string(ptcl.pinum_Ymax)},
                    };
                    logger.logSection("flowout_"+pName, data);
                    // 電荷密度の更新
                    if (EqFlags.EMField == 1){
                        MEASURE("integCDens_"+pName, [ptcl integrateChargeDensity:fld  withMoment:mom withLogger:logger], dataElapsedTime);
                    }
                    // 粒子軌道の出力
                    if (timeParam.ParticleOutput != 0 && cycle%timeParam.ParticleOutput == 0){
                        outputPhaseSpace(cycle, ptcl, init, logger);
                    }
                    // モーメントの出力
                    if (timeParam.FieldOutput != 0 && cycle%timeParam.FieldOutput == 0){
                        [mom integrateMoments:ptcl withEMField:fld withLogger:logger];
                        outputMoments(cycle, mom, particles[s].pName, init, logger);
                    }
                }

                // 電場の更新
                if (EqFlags.EMField == 1){
                    MEASURE("solvePoisson", [fld solvePoisson:logger], dataElapsedTime);
                    // 場の出力
                    if (timeParam.FieldOutput != 0 && cycle%timeParam.FieldOutput == 0){
                        outputField(cycle, fld, init, logger);
                    }
                }

                // 粒子生成
                std::vector<int> ret(EqFlags.Particle);
                for (int s = 0; s < EqFlags.Particle; s++) {
                    Particle *ptcl = [ptclArr objectAtIndex:s];
                    std::string pName = [particles[s].pName UTF8String];
                    MEASURE("injection_"+pName, ret.push_back([ptcl injection:dt withParam:init withCurrent:intCurrent withLogger:logger]), dataElapsedTime);
                }
                if (std::reduce(std::begin(ret), std::end(ret)) != 0){
                    NSLog(@"injection failed.");
                    return 1;
                }

                // 途中経過出力
                if (timeParam.FieldOutput != 0 && cycle%timeParam.ProgressOutput == 0){
                    if(!saveProgress(cycle, ptclArr, fld, intCurrent, init)){
                        NSLog(@"output progress failed.");
                        return 1;
                    };
                }

                // リソース情報出力(unit: time[us], mem[kB])
                logger.logSection("elapsedTime", dataElapsedTime);
                size_t phys = MemoryUtils::currentPhysicalFootprint();
                size_t rss = MemoryUtils::currentRSS();
                dataMemUsage["physicalFootprint"] = std::to_string(phys/1024);
                dataMemUsage["residentSetSize"] = std::to_string(rss/1024);
                logger.logSection("memoryUsage", dataMemUsage);
                
                //サイクル終了
                NSLog(@"Frame %d completed", cycle);
                logger.logCycleEnd();
            }
        }

        // 計測終了
        auto end = std::chrono::high_resolution_clock::now();
        auto elapsedTime = std::chrono::duration_cast<std::chrono::seconds>(end - start).count(); \
        logger.logComment("culculation end.\ntotal elapsed time: "+std::to_string(elapsedTime)+" sec.\n");
        
    }
    return 0;
}