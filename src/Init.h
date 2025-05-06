#include <Foundation/Foundation.h>
#include <string>
#include <array>

@interface Init : NSObject

struct FragForEquation {
    int Particle;
    int EMField;
    int MCCollision;
};

struct ParamForTimeIntegration {
    int EndCycle;
    int ptclOutCycle;
    int fldOutCycle;
    double dt;
};

struct ParamForParticle {
    NSString* pName;
    uint pNum;
    uint pNumMax;
    double q;
    double m;
    double w;
    NSString* GenerateType;
    double initX[2];
    double initY[2];
    double initU[3];
    double initT;
};

struct BoundaryConditionForParticle {
    NSString* position;
    NSString* type;
    double val;
};

struct ParamForField {
    int ngx;
    int ngy;
    int ngb;
    double dx;
    double dy;
    NSString* InitType;
    double ampE[3];
    double ampB[3];
    int weightOrder;
    int maxiter;
    float tolerance;
};

struct BoundaryConditionForField {
    NSString* position;
    NSString* type;
    double val;
};

- (instancetype)parseInputFile:(NSString*)inputFilePath;
- (BOOL)parseFile:(std::ifstream&) inputFile;
- (void)parseFlagForEquation:(const std::string&)line;
- (void)parseParamForTimeIntegration:(const std::string&)line;
- (void)parseParamForParticle:(const std::string&)line inputFile:(std::ifstream&)inputFile;
- (void)parseBoundaryConditionForParticle:(const std::string&)line inputFile:(std::ifstream&)inputFile;
- (void)parseParamForField:(const std::string&)line;
- (void)parseBoundaryConditionForField:(const std::string&)line inputFile:(std::ifstream&)inputFile;
- (struct FragForEquation)getFragForEquation;
- (struct ParamForTimeIntegration)getParamForTimeIntegration;
- (NSArray*)getParamForParticle;
- (NSArray*)getParticleBoundaries;
- (struct ParamForField)getParamForField;
- (NSArray*)getFieldBoundaries;

@end