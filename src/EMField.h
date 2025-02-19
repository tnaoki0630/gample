#import <Foundation/Foundation.h>
#import "Init.h"
#import "Particle.h"
// コンパイルうまくいかないので前方宣言
@class Particle;

struct ElectroStatic {
    double rho;
    double phi;
};

struct ElectricField {
    double Ex;
    double Ey;
    double Ez;
};

struct MagneticField {
    double Bx;
    double By;
    double Bz;
};

@interface EMField : NSObject

- (instancetype)initWithParam:(Init*) init;

- (void)culcChargeDensity:(Particle*) ptcl;

- (void)solvePoisson;

@end