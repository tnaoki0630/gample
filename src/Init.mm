#import "Init.h"
#import "Constant.h"
#import <string>
#import <fstream>
#import <iostream>
#import <sstream>
#import <vector>

@interface Init() {
    // Private properties to store parsed values
    struct FlagForEquation _flagEquation;
    struct ParamForTimeIntegration _timeIntegration;
    struct ParamForComputing _computing;
    std::vector<struct ParamForParticle> _particles;
    std::vector<struct BoundaryConditionForParticle> _particleBoundaries;
    std::vector<struct SourceForParticle> _particleSources;
    struct ParamForField _field;
    std::vector<struct BoundaryConditionForField> _fieldBoundaries;
}

@property (nonatomic, strong) NSString *filePath;

@end

@implementation Init

- (instancetype)parseInputFile:(NSString*)InputFilePath {
    self = [super init];
    if (self) {

        // ファイル読み込み
        NSError* __autoreleasing readError = nil;
        NSData *data = [NSData dataWithContentsOfFile:InputFilePath options:0 error:&readError];
        if (!data) {
            return nil;  // NSError に詳細が入っている
        }
        
        // JSON パース
        id jsonObj = [NSJSONSerialization JSONObjectWithData:data
                                                    options:NSJSONReadingMutableContainers
                                                    error:&readError];
        if (!jsonObj) {
            return nil;
        }
        if (![jsonObj isKindOfClass:[NSDictionary class]]) {
            if (readError) {
                readError = [NSError errorWithDomain:@"ConfigError" code:1001 userInfo:@{ NSLocalizedDescriptionKey : @"トップレベルが辞書ではありません" }];
            }
            return nil;
        }
        NSDictionary *root = (NSDictionary *)jsonObj;
        
        // FlagForEquation
        NSDictionary *flagD = root[@"FlagForEquation"];
        if (![flagD isKindOfClass:[NSDictionary class]]) {
            if (readError) {
                readError = [NSError errorWithDomain:@"ConfigError" code:1002 userInfo:@{ NSLocalizedDescriptionKey : @"FlagForEquation が辞書ではありません" }];
            }
            return nil;
        }
        NSMutableDictionary *mergedFlagD = [self FlagForEquationDefault];
        [mergedFlagD addEntriesFromDictionary:flagD];
        _flagEquation.Particle    = [mergedFlagD[@"Particle"]    integerValue];
        _flagEquation.EMField     = [mergedFlagD[@"EMField"]     integerValue];
        _flagEquation.MCCollision = [mergedFlagD[@"MCCollision"] integerValue];

        // ParamForTimeIntegration
        NSDictionary *timeD = root[@"ParamForTimeIntegration"];
        if (![timeD isKindOfClass:[NSDictionary class]]) {
            if (readError) {
                readError = [NSError errorWithDomain:@"ConfigError" code:1003 userInfo:@{ NSLocalizedDescriptionKey : @"ParamForTimeIntegration が辞書ではありません" }];
            }
            return nil;
        }
        NSMutableDictionary *mergedTimeD = [self ParamForTimeIntegrationDefault];
        [mergedTimeD addEntriesFromDictionary:timeD];
        _timeIntegration.Start          = [mergedTimeD[@"Start"] integerValue];
        _timeIntegration.End            = [mergedTimeD[@"End"] integerValue];
        _timeIntegration.TimeStep       = [mergedTimeD[@"TimeStep"] doubleValue];
        _timeIntegration.LogOutput      = [mergedTimeD[@"LogOutput"] integerValue];
        _timeIntegration.ParticleOutput = [mergedTimeD[@"ParticleOutput"] integerValue];
        _timeIntegration.FieldOutput    = [mergedTimeD[@"FieldOutput"] integerValue];
        _timeIntegration.ProgressOutput = [mergedTimeD[@"ProgressOutput"] integerValue];
        _timeIntegration.ProjectName    = mergedTimeD[@"ProjectName"];
        _timeIntegration.RestartName    = mergedTimeD[@"RestartName"];
        
        // ParamForComputing
        NSDictionary *compD = root[@"ParamForComputing"];
        if (![compD isKindOfClass:[NSDictionary class]]) {
            if (readError) {
                readError = [NSError errorWithDomain:@"ConfigError" code:1004 userInfo:@{ NSLocalizedDescriptionKey : @"ParamForComputing が辞書ではありません" }];
            }
            return nil;
        }
        NSMutableDictionary *mergedCompD = [self ParamForComputingDefault];
        [mergedCompD addEntriesFromDictionary:compD];
        _computing.pNumMax                  = [mergedCompD[@"MaximumParticleNumber"] integerValue];
        _computing.threadGroupSize          = [mergedCompD[@"ThreadGroupSize"] integerValue];
        _computing.integrationChunkSize     = [mergedCompD[@"IntegrationChunkSize"] integerValue];
        _computing.maxiter                  = [mergedCompD[@"MaxiterForPoisson"] integerValue];
        _computing.tolerance                = [mergedCompD[@"ToleranceForPoisson"] doubleValue];
        
        // ParamForField
        NSDictionary *fldD = root[@"ParamForField"];
        if (![fldD isKindOfClass:[NSDictionary class]]) {
            if (readError) {
                readError = [NSError errorWithDomain:@"ConfigError" code:1006 userInfo:@{ NSLocalizedDescriptionKey : @"ParamForField が辞書ではありません" }];
            }
            return nil;
        }
        NSMutableDictionary *mergedFldD = [self ParamForFieldDefault];
        [mergedFldD addEntriesFromDictionary:fldD];
        _field.ngx          = [mergedFldD[@"NumberOfGridX"] integerValue];
        _field.ngy          = [mergedFldD[@"NumberOfGridY"] integerValue];
        _field.dx           = [mergedFldD[@"GridSizeOfX"] doubleValue];
        _field.dy           = [mergedFldD[@"GridSizeOfY"] doubleValue];
        _field.InitTypeE    = mergedFldD[@"InitializeTypeOfE"];
        NSArray *ampArr; // vector<double> 読み込み
        ampArr = mergedFldD[@"AmplitudeOfE"];
        _field.ampE.clear();
        if ([ampArr isKindOfClass:[NSArray class]]) {
            for (id v in ampArr) _field.ampE.push_back([v doubleValue]);
        }
        _field.FilePathEx   = mergedFldD[@"FilePathOfEx"];
        _field.FilePathEy   = mergedFldD[@"FilePathOfEy"];
        _field.FilePathEz   = mergedFldD[@"FilePathOfEz"];
        _field.InitTypeB    = mergedFldD[@"InitializeTypeOfB"];
        ampArr = mergedFldD[@"AmplitudeOfB"];
        _field.ampB.clear();
        if ([ampArr isKindOfClass:[NSArray class]]) {
            for (id v in ampArr) _field.ampB.push_back([v doubleValue]);
        }
        _field.FilePathBx   = mergedFldD[@"FilePathOfBx"];
        _field.FilePathBy   = mergedFldD[@"FilePathOfBy"];
        _field.FilePathBz   = mergedFldD[@"FilePathOfBz"];
        _field.weightOrder  = [mergedFldD[@"WeightingOrder"] integerValue];
        _field.ngb  = _field.weightOrder/2;

        // BoundaryConditionForField
        NSArray* fbcArr = root[@"BoundaryConditionForField"];
        if (![fbcArr isKindOfClass:[NSArray class]]) {
            if (readError) {
                readError = [NSError errorWithDomain:@"ConfigError" code:1005 userInfo:@{ NSLocalizedDescriptionKey : @"BoundaryConditionForParticle が配列ではありません" }];
            }
            return nil;
        }
        _fieldBoundaries.clear();
        for (id elem in fbcArr) {
            if (![elem isKindOfClass:[NSDictionary class]]) continue;
            NSDictionary *d = (NSDictionary *)elem;
            NSMutableDictionary *mergedD = [self BoundaryCoditionDefault];
            [mergedD addEntriesFromDictionary:d];
            BoundaryConditionForField bc;
            bc.position = mergedD[@"RegionName"];
            bc.type     = mergedD[@"BCType"];
            bc.val      = [mergedD[@"Value"] doubleValue];
            _fieldBoundaries.push_back(bc);
        }

        // ParamForParticle (配列)
        NSArray *partArr = root[@"ParamForParticle"];
        if (![partArr isKindOfClass:[NSArray class]]) {
            if (readError) {
                readError = [NSError errorWithDomain:@"ConfigError" code:1004 userInfo:@{ NSLocalizedDescriptionKey : @"ParamForParticle が配列ではありません" }];
            }
            return nil;
        }
        _particles.clear();
        for (id elem in partArr) {
            if (![elem isKindOfClass:[NSDictionary class]]) continue;
            NSDictionary *d = (NSDictionary *)elem;
            NSMutableDictionary *mergedD = [self ParamForParticleDefault];
            [mergedD addEntriesFromDictionary:d];
            ParamForParticle p;
            p.pName             = mergedD[@"ParticleName"];
            if([mergedD[@"pNumSetMethod"] isEqualToString:@"InititalParticleNumber"]){
                p.pNum = [mergedD[@"pNumValue"] integerValue];
            }else if([mergedD[@"pNumSetMethod"] isEqualToString:@"PtclNumPerCell"]){
                p.pNum = [mergedD[@"pNumValue"] integerValue]*_field.ngx*_field.ngy;
            }else{
                NSLog(@"pNumSet failed.");
                return nil;
            }
            p.q = [mergedD[@"Charge"] doubleValue];
            p.m = [mergedD[@"Mass"] doubleValue];
            if([mergedD[@"WeightSetMethod"] isEqualToString:@"WeightValue"]){
                p.w = [mergedD[@"WeightValue"] doubleValue];
            }else if([mergedD[@"WeightSetMethod"] isEqualToString:@"InitialNumDens"]){
                double ppc = (double)(p.pNum)/(double)(_field.ngx*_field.ngy);
                p.w = [mergedD[@"WeightValue"] doubleValue]*_field.dx*_field.dy/ppc;
            }else{
                NSLog(@"WeightSet failed.");
                return nil;
            }
            p.genType = mergedD[@"GenerateType"];
            NSArray *velArr;    // vector<double> の読み込み
            velArr = mergedD[@"InitialPosX"];
            if ([velArr isKindOfClass:[NSArray class]]) {
                for (id v in velArr) p.genX.push_back([v doubleValue]);
            }
            velArr = mergedD[@"InitialPosY"];
            if ([velArr isKindOfClass:[NSArray class]]) {
                for (id v in velArr) p.genY.push_back([v doubleValue]);
            }
            velArr = mergedD[@"InitialVel"];
            if ([velArr isKindOfClass:[NSArray class]]) {
                for (id v in velArr) p.genU.push_back([v doubleValue]);
            }
            p.genT = [mergedD[@"InitialTemp"] doubleValue]*eVtoK;
            _particles.push_back(p);
        }
        
        // SourceForParticle
        NSArray *srcArr = root[@"SourceForParticle"];
        if (![srcArr isKindOfClass:[NSArray class]]) {
            if (readError) {
                readError = [NSError errorWithDomain:@"ConfigError" code:1004 userInfo:@{ NSLocalizedDescriptionKey : @"ParamForParticle が配列ではありません" }];
            }
            return nil;
        }
        _particleSources.clear();
        for (id elem in srcArr) {
            if (![elem isKindOfClass:[NSDictionary class]]) continue;
            NSDictionary *d = (NSDictionary *)elem;
            NSMutableDictionary *mergedD = [self SourceForParticleDefault];
            [mergedD addEntriesFromDictionary:d];
            SourceForParticle src;
            src.pName       = mergedD[@"ParticleName"];
            src.genType     = mergedD[@"GenerateType"];
            if([mergedD[@"SrcSetMethod"] isEqualToString:@"Source[1/cm/s]"]){
                src.src = [mergedD[@"SrcValue"] doubleValue];
            }else if([mergedD[@"SrcSetMethod"] isEqualToString:@"CurrentDensity[A/m2]"]){
                src.src = [mergedD[@"SrcValue"] doubleValue]*JtosJ*_field.ngy*_field.dy/ec; // 単一電荷を仮定しているので注意
            }else if([src.genType isEqualToString:@"hollow-cathode"]){
                src.src = 0;
            }else{
                NSLog(@"SrcSet failed.");
                return nil;
            }
            NSArray *velArr;
            velArr = mergedD[@"GeneratePosX"];
            if ([velArr isKindOfClass:[NSArray class]]) {
                for (id v in velArr) src.genX.push_back([v doubleValue]);
            }
            velArr = mergedD[@"GeneratePosY"];
            if ([velArr isKindOfClass:[NSArray class]]) {
                for (id v in velArr) src.genY.push_back([v doubleValue]);
            }
            velArr = mergedD[@"GenerateVel"];
            if ([velArr isKindOfClass:[NSArray class]]) {
                for (id v in velArr) src.genU.push_back([v doubleValue]);
            }
            src.genT              = [mergedD[@"GenerateTemp"] doubleValue]*eVtoK;
            _particleSources.push_back(src);
        }
        
        // BoundaryConditionForParticle
        NSArray *pbcArr = root[@"BoundaryConditionForParticle"];
        if (![pbcArr isKindOfClass:[NSArray class]]) {
            if (readError) {
                readError = [NSError errorWithDomain:@"ConfigError" code:1005 userInfo:@{ NSLocalizedDescriptionKey : @"BoundaryConditionForParticle が配列ではありません" }];
            }
            return nil;
        }
        _particleBoundaries.clear();
        for (id elem in pbcArr) {
            if (![elem isKindOfClass:[NSDictionary class]]) continue;
            NSDictionary *d = (NSDictionary *)elem;
            NSMutableDictionary *mergedD = [self BoundaryCoditionDefault];
            [mergedD addEntriesFromDictionary:d];
            BoundaryConditionForParticle bc;
            bc.position = d[@"RegionName"];
            bc.type     = d[@"BCType"];
            _particleBoundaries.push_back(bc);
        }
    }
    return self;
}

// default values
- (NSMutableDictionary*)FlagForEquationDefault{
    return [@{
        @"Particle":        @0,
        @"EMField":         @0,
        @"MCCollision":     @0,
    } mutableCopy];
}
- (NSMutableDictionary*)ParamForTimeIntegrationDefault{
    return [@{
        @"Start":           @0,
        @"End":             @0,
        @"TimeStep":        @0,
        @"LogOutput":       @0,
        @"ParticleOutput":  @0,
        @"FieldOutput":     @0,
        @"ProgressOutput":  @0,
        @"ProjectName":     @"undefined",
        @"RestartName":     @"undefined"
    } mutableCopy];
}
- (NSMutableDictionary*)ParamForComputingDefault{
    return [@{
        @"MaximumParticleNumber":   @0,
        @"ThreadGroupSize":         @256,
        @"IntegrationChunkSize":    @256,
        @"MaxiterForPoisson":       @200,
        @"ToleranceForPoisson":     @(1e-7)
    } mutableCopy];
}
- (NSMutableDictionary*)ParamForParticleDefault{
    return [@{
        @"ParticleName":            @"undefined",
        @"pNumSetMethod":           @"undefined",
        @"pNumValue":               @0,
        @"Charge":                  @0,
        @"Mass":                    @(-1),
        @"WeightSetMethod":         @"undefined",
        @"WeightValue":             @0,
        @"GenerateType":            @"undefined",
        @"InitialPosX":             @[@(-1), @(-1)],
        @"InitialPosY":             @[@(-1), @(-1)],
        @"InitialVel":              @[@0, @0, @0],
        @"InitialTemp":             @0
    } mutableCopy];
}
- (NSMutableDictionary*)BoundaryCoditionDefault{
    return [@{
        @"RegionName":  @"undefined",
        @"BCType":      @"undefined",
        @"Value":       @0
    } mutableCopy];
}
- (NSMutableDictionary*)SourceForParticleDefault{
    return [@{
        @"ParticleName":    @"undefined",
        @"GenerateType":    @"undefined",
        @"SrcSetMethod":    @"undefined",
        @"SrcVal":          @0,
        @"GeneratePosX":    @[@(-1), @(-1)],
        @"GeneratePosY":    @[@(-1), @(-1)],
        @"GenerateVel":     @[@0, @0, @0],
        @"GenerateTemp":    @0
    } mutableCopy];
}
- (NSMutableDictionary*)ParamForFieldDefault{
    return [@{
        @"NumberOfGridX":       @0,
        @"NumberOfGridY":       @0,
        @"GridSizeOfX":         @0,
        @"GridSizeOfY":         @0,
        @"InitializeTypeOfE":   @"undefined",
        @"AmplitudeOfE":        @[@0, @0, @0],
        @"FilePathOfEx":        @"undefined",
        @"FilePathOfEy":        @"undefined",
        @"FilePathOfEz":        @"undefined",
        @"InitializeTypeOfB":   @"undefined",
        @"AmplitudeOfB":        @[@0, @0, @0],
        @"FilePathOfBx":        @"undefined",
        @"FilePathOfBy":        @"undefined",
        @"FilePathOfBz":        @"undefined",
        @"WeightingOrder":      @5,
    } mutableCopy];
}

// Public accessor methods to get the parsed data
- (struct FlagForEquation)flagForEquation { return _flagEquation; }
- (struct ParamForTimeIntegration)paramForTimeIntegration { return _timeIntegration; }
- (struct ParamForComputing)paramForComputing { return _computing; }
- (std::vector<struct ParamForParticle>)paramForParticle { return _particles; }
- (std::vector<struct BoundaryConditionForParticle>)particleBoundaries { return _particleBoundaries; }
- (std::vector<struct SourceForParticle>)particleSources{ return _particleSources; }
- (struct ParamForField)paramForField { return _field; }
- (std::vector<struct BoundaryConditionForField>)fieldBoundaries { return _fieldBoundaries; }
@end