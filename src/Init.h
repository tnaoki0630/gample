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
    NSString* genType;
    double genX[2];
    double genY[2];
    double genU[3];
    double genT;
};

struct BoundaryConditionForParticle {
    NSString* position;
    NSString* type;
    double val;
};

struct SourceForParticle {
    NSString* pName;
    NSString* genType;
    double src;
    double genX[2];
    double genY[2];
    double genU[3];
    double genT;
};

struct ParamForField {
    int ngx;
    int ngy;
    int ngb;
    double dx;
    double dy;
    NSString* InitTypeE;
    NSString* InitTypeB;
    double ampE[3];
    double ampB[3];
    NSString* FilePathEx;
    NSString* FilePathEy;
    NSString* FilePathEz;
    NSString* FilePathBx;
    NSString* FilePathBy;
    NSString* FilePathBz;
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
- (BOOL)checkInput;
- (void)parseFlagForEquation:(const std::string&)line;
- (void)parseParamForTimeIntegration:(const std::string&)line;
- (void)parseParamForParticle:(const std::string&)line inputFile:(std::ifstream&)inputFile;
- (void)parseBoundaryConditionForParticle:(const std::string&)line inputFile:(std::ifstream&)inputFile;
- (void)parseSourceForParticle:(const std::string&)line inputFile:(std::ifstream&)inputFile;
- (void)parseParamForField:(const std::string&)line;
- (void)parseBoundaryConditionForField:(const std::string&)line inputFile:(std::ifstream&)inputFile;
- (struct FragForEquation)getFragForEquation;
- (struct ParamForTimeIntegration)getParamForTimeIntegration;
- (NSArray*)getParamForParticle;
- (NSArray*)getParticleBoundaries;
- (NSArray*)getParticleSources;
- (struct ParamForField)getParamForField;
- (NSArray*)getFieldBoundaries;

@end