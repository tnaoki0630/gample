#import <iostream>
#import <numeric>
#import <Foundation/Foundation.h>
#import <mach/mach.h>
#import "Init.h"
#import "Particle.h"
#import "Moment.h"
#import "Collision.h"
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
        // flags
        const bool debug = false;

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

        // ログファイル名
        std::string timestamp;
        if (args.count("-o") == 1){
            timestamp = args["-o"].c_str();
        }else{
            // 無指定ならタイムスタンプ YYYYMMDDhhmmss を付与
            auto now = std::chrono::system_clock::now();
            auto nowt = std::chrono::system_clock::to_time_t(now);
            std::tm tm = *std::localtime(&nowt);
            std::ostringstream oss;
            oss << std::put_time(&tm, "%Y%m%d%H%M%S");
            timestamp = oss.str();
        }

        // 初期化パラメータクラス作成
        if (debug){ NSLog(@"parseInputFile"); }
        Init *init = [[Init alloc] parseInputFile:inputPath];
        if(!init){
            NSLog(@"[Fatal] parseInputFile Failed");
            return -1;
        }

        // logger の作成
        struct ParamForTimeIntegration timeParam = init.paramForTimeIntegration;
        NSLog(@"timeParam.TimeStep = %f",timeParam.TimeStep);
        XmlLogger logger([[NSString stringWithFormat:@"%@_%s_log.xml", timeParam.ProjectName, timestamp.c_str()] UTF8String]);
        
        // パース結果の出力(debugprint.mm)
        printInitContents(init, logger);

        // 粒子クラスの初期化
        struct FlagForEquation EqFlags = init.flagForEquation;
        NSLog(@"EqFlags.Particle = %d",EqFlags.Particle);
        NSLog(@"EqFlags.EMField = %d",EqFlags.EMField);
        NSLog(@"EqFlags.Collision = %d",EqFlags.MCCollision);

        NSMutableArray *ptclArr = [NSMutableArray arrayWithCapacity:EqFlags.Particle];
        for (int s = 0; s < EqFlags.Particle; s++) {
            if (debug){ NSLog(@"initParticle#%d",s); }
            Particle *ptcl = [[Particle alloc] initWithDevice:device withParam:init specimen:s withLogger:logger];
            if(!ptcl){
                NSLog(@"[Fatal] initParticle(#%d) Failed",s);
                return 1;
            }
            [ptclArr addObject:ptcl];
        }
        // 電磁場クラスの初期化
        if (debug){ NSLog(@"initField"); }
        EMField *fld = [[EMField alloc] initWithDevice:device withParam:init withLogger:logger];
        if(!fld){
            NSLog(@"[Fatal] initFld Failed");
            return 1;
        }

        // モーメント計算クラスの初期化
        Moment *mom = [[Moment alloc] initWithDevice:device withParam:init withLogger:logger];
        // 衝突計算クラスの初期化
        Collision *col = [[Collision alloc] initWithDevice:device withParam:init withParticles:ptclArr withLogger:logger];
        if(!col){
            NSLog(@"[Fatal] initCol Failed");
            return 1;
        }

        // 電荷密度の初期化
        [fld resetChargeDensity];

        // リスタート
        int StartCycle = timeParam.Start;
        if(StartCycle != 0){
            // 粒子状態読み込み
            if(!loadProgress(StartCycle, ptclArr, init)){
                NSLog(@"[Fatal] restart failed.");
                return 1;
            };
            // 電荷密度計算
            for (Particle* ptcl in ptclArr) {
                std::string pName = [ptcl.pName UTF8String];
                [ptcl integrateChargeDensity:fld withMoment:mom withLogger:logger];
            }
        }

        // 初期場計算
        [fld solvePoisson:logger];
        // 初期場出力
        outputField(StartCycle, fld, init, logger);
        for(Particle* ptcl in ptclArr){
            [mom integrateMoments:ptcl withEMField:fld withLogger:logger];
            outputMoments(StartCycle, mom, ptcl.pName, init, logger);
        }

        // 非定常ループ開始
        double dt = timeParam.TimeStep;
        double time = StartCycle*dt;
        double comp = 0.0; // 保証項
        double y,t;
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

                // 電場の更新
                if (EqFlags.EMField == 1){
                    // 電荷密度の初期化
                    [fld resetChargeDensity];
                    // 電荷密度の更新
                    if (debug){ NSLog(@"integCDens"); }
                    for (Particle* ptcl in ptclArr) {
                        std::string pName = [ptcl.pName UTF8String];
                        MEASURE("integCDens_"+pName, [ptcl integrateChargeDensity:fld withMoment:mom withLogger:logger], dataElapsedTime);
                    }
                    // 静電場計算
                    if (debug){ NSLog(@"poisson"); }
                    MEASURE("solvePoisson", [fld solvePoisson:logger], dataElapsedTime);
                    // 場の出力
                    if (timeParam.FieldOutput != 0 && cycle%timeParam.FieldOutput == 0){
                        outputField(cycle, fld, init, logger);
                    }
                }

                // 粒子更新
                std::vector<struct ParamForParticle> particles = init.paramForParticle;
                for (int s = 0; s < EqFlags.Particle; s++) {
                    Particle *ptcl = [ptclArr objectAtIndex:s];
                    std::string pName = [ptcl.pName UTF8String];
                    // 粒子の時間更新
                    if (debug){ NSLog(@"update"); }
                    MEASURE("update_"+pName, [ptcl update:dt withEMField:fld withMom:mom withLogger:logger], dataElapsedTime);
                    // 流出粒子の処理
                    if (debug){ NSLog(@"reduce"); }
                    MEASURE("reduce_"+pName, [ptcl reduce:logger], dataElapsedTime);
                    // 粒子軌道の出力
                    if (timeParam.ParticleOutput != 0 && cycle%timeParam.ParticleOutput == 0){
                        outputPhaseSpace(cycle, ptcl, init, logger);
                    }
                    // モーメントの出力
                    if (timeParam.FieldOutput != 0 && cycle%timeParam.FieldOutput == 0){
                        [mom integrateMoments:ptcl withEMField:fld withLogger:logger];
                        if (debug){ NSLog(@"moments"); }
                        outputMoments(cycle, mom, ptcl.pName, init, logger);
                    }
                }

                // 粒子生成
                bool ret;
                if (debug){ NSLog(@"AI"); }
                MEASURE("artificialIonization", ret = [col artificialIonization:dt withParticles:ptclArr withLogger:logger], dataElapsedTime);
                if (!ret){
                    NSLog(@"injection failed.");
                    return 1;
                }
                if (debug){ NSLog(@"HC"); }
                MEASURE("hollowCathode", ret = [col hollowCathode:dt withParticles:ptclArr withLogger:logger], dataElapsedTime);
                if (!ret){
                    NSLog(@"injection failed.");
                    return 1;
                }

                // 途中経過出力
                if (timeParam.FieldOutput != 0 && cycle%timeParam.ProgressOutput == 0){
                    if(!saveProgress(cycle, ptclArr, init)){
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
                if(cycle%timeParam.LogOutput == 0) NSLog(@"Frame %d completed", cycle);
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