#import <Foundation/Foundation.h>

@interface Init : NSObject

struct ParamForTimeIntegration {
    int StartCycle;
    int EndCycle;
    int OutputCycle;
    double dt;
};

struct ParamForParticle {
    int NumberOfParticle;
    int NOPmax;
    double Charge;
    double Mass;
    double Weight;
};

struct ParamForField {
    int ngx;
    int ngy;
    double dx;
    double dy;
};

- (instancetype)initWithFile:(NSString*) InputFilePath;

@end