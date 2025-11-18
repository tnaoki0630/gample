// #import "ParticleParam.h"
#include <Foundation/Foundation.h>
#include <string>
#include <array>
#include <vector>

@interface Init : NSObject

struct FlagForEquation {
    int Particle;
    int EMField;
    int MCCollision;
};

struct ParamForTimeIntegration {
    int Start;
    int End;
    double TimeStep;
    int LogOutput;
    int ParticleOutput;
    int FieldOutput;
    int ProgressOutput;
    NSString* ProjectName;
    NSString* RestartName;
};

struct ParamForComputing {
    uint pNumMax;
    int threadGroupSize;
    int integrationChunkSize;
    float aggrThreshold;
    int chebDegree;
    int amgCycleType;
    int maxiter;
    float tolerance;
};

struct ParamForParticle {
    // general
    NSString* pName;
    uint pNum;
    double q;
    double m;
    double w;
    // generation
    NSString* genType;
    std::vector<double> genX;
    std::vector<double> genY;
    std::vector<double> genU;
    double genT;
    // output
    int outNum;
};

struct BoundaryConditionForParticle {
    NSString* position;
    NSString* type;
};

struct SourceForParticle {
    NSString* pName;
    NSString* genType;
    double src;
    std::vector<double> genX;
    std::vector<double> genY;
    std::vector<double> genU;
    std::vector<double> genU_ele;
    double genT;
    double genT_ele;
};

struct ParamForField {
    int ngx;
    int ngy;
    int ngb;
    double dx;
    double dy;
    NSString* InitTypeE;
    NSString* InitTypeB;
    std::vector<double> ampE;
    std::vector<double> ampB;
    NSString* FilePathEx;
    NSString* FilePathEy;
    NSString* FilePathEz;
    NSString* FilePathBx;
    NSString* FilePathBy;
    NSString* FilePathBz;
    int weightOrder;
    int BorisOrder;
    int checkCurDens;
    std::vector<int> outputPoint;
    int diff_ngb;
};

struct BoundaryConditionForField {
    NSString* position;
    NSString* type;
    double val;
};

- (instancetype)parseInputFile:(NSString*)inputFilePath;
// デフォルト辞書
- (NSMutableDictionary*)FlagForEquationDefault;
- (NSMutableDictionary*)ParamForTimeIntegrationDefault;
- (NSMutableDictionary*)ParamForParticleDefault;
- (NSMutableDictionary*)BoundaryCoditionDefault;
- (NSMutableDictionary*)SourceForParticleDefault;
- (NSMutableDictionary*)ParamForFieldDefault;
// アクセサ
- (struct FlagForEquation)flagForEquation;
- (struct ParamForTimeIntegration)paramForTimeIntegration;
- (struct ParamForComputing)paramForComputing;
- (std::vector<struct ParamForParticle>)paramForParticle;
- (std::vector<struct BoundaryConditionForParticle>)particleBoundaries;
- (std::vector<struct SourceForParticle>)particleSources;
- (struct ParamForField)paramForField;
- (std::vector<struct BoundaryConditionForField>)fieldBoundaries;

@end