#import <Foundation/Foundation.h>
#import <string>

@interface Init : NSObject

struct FragForEquation {
    int Particle;
    int EMField;
    int MCCollision;
};

struct ParamForTimeIntegration {
    int StartCycle;
    int EndCycle;
    int OutputCycle;
    double dt;
};

struct ParamForParticle {
    char* pName;
    int pNum;
    int pNumMax;
    double q;
    double m;
    double w;
    double initalU;
    double initalV;
    double initalW;
    double initialT;
};

struct BoundaryConditionForParticle {
    char* position;
    char* type;
    double val;
};

struct ParamForField {
    int ngx;
    int ngy;
    double dx;
    double dy;
};

struct BoundaryConditionForField {
    char* position;
    char* type;
    double val;
};

- (instancetype)parseInputFile:(NSString*) InputFilePath;
- (BOOL)parseFile:(std::ifstream&) inputFile;
- (void)parseFlagForEquation:(const std::string&)line;
- (void)parseParamForTimeIntegration:(const std::string&)line;
- (void)parseParamForParticle:(const std::string&)line inputFile:(std::ifstream&)inputFile;
- (void)parseBoundaryConditionForParticle:(const std::string&)line inputFile:(std::ifstream&)inputFile;
- (void)parseParamForField:(const std::string&)line;
- (void)parseBoundaryConditionForField:(const std::string&)line inputFile:(std::ifstream&)inputFile;
- (struct FragForEquation)getFragForEquation;
- (struct ParamForTimeIntegration)getParamForTimeIntegration;
- (NSArray*)getParticles;
- (NSArray*)getParticleBoundaries;
- (struct ParamForField)getParamForField;
- (NSArray*)getFieldBoundaries;
- (void)dealloc;

@end