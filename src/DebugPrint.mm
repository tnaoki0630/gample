#import <Foundation/Foundation.h>
#import "Init.h"
#import "Particle.h"
#import "Moment.h"
#import "EMField.h"

void printInitContents(Init* init) {
    NSLog(@"==== パース結果出力開始 ====");
    
    // FlagForEquation の出力
    struct FragForEquation flags = [init getFragForEquation];
    NSLog(@"--- FlagForEquation ---");
    NSLog(@"Particle: %d", flags.Particle);
    NSLog(@"EMField: %d", flags.EMField);
    NSLog(@"MCCollision: %d", flags.MCCollision);
    
    // ParamForTimeIntegration の出力
    struct ParamForTimeIntegration timeParams = [init getParamForTimeIntegration];
    NSLog(@"--- ParamForTimeIntegration ---");
    NSLog(@"EndCycle: %d", timeParams.EndCycle);
    NSLog(@"OutputCycle: %d", timeParams.OutputCycle);
    NSLog(@"dt: %e", timeParams.dt);
    
    // ParamForParticle の出力
    NSArray *particles = [init getParamForParticle];
    NSLog(@"--- ParamForParticle (%lu entries) ---", (unsigned long)particles.count);
    for (NSUInteger i = 0; i < particles.count; i++) {
        NSValue *value = particles[i];
        struct ParamForParticle particle;
        [value getValue:&particle];
        
        NSLog(@"[Particle %lu]", (unsigned long)i + 1);
        NSLog(@"  Name: %@", particle.pName);
        NSLog(@"  InitialParticleNumber: %d", particle.pNum);
        NSLog(@"  MaxParticleNumber: %d", particle.pNumMax);
        NSLog(@"  Charge: %e", particle.q);
        NSLog(@"  Mass: %e", particle.m);
        NSLog(@"  Weight: %e", particle.w);
        NSLog(@"  Weight: %@", particle.GenerateType);
        NSLog(@"  initialU: %e", particle.initU[0]);
        NSLog(@"  initialV: %e", particle.initU[1]);
        NSLog(@"  initialW: %e", particle.initU[2]);
        NSLog(@"  initialT: %e", particle.initT);
    }
    
    // BoundaryConditionForParticle の出力
    NSArray *particleBoundaries = [init getParticleBoundaries];
    NSLog(@"--- BoundaryConditionForParticle (%lu entries) ---", (unsigned long)particleBoundaries.count);
    for (NSUInteger i = 0; i < particleBoundaries.count; i++) {
        NSValue *value = particleBoundaries[i];
        struct BoundaryConditionForParticle boundary;
        [value getValue:&boundary];
        
        NSLog(@"[Boundary %lu]", (unsigned long)i + 1);
        NSLog(@"  Position: %@", boundary.position);
        NSLog(@"  Type: %@", boundary.type);
        NSLog(@"  Value: %e", boundary.val);
    }
    
    // ParamForField の出力
    struct ParamForField fieldParams = [init getParamForField];
    NSLog(@"--- ParamForField ---");
    NSLog(@"NGX: %d", fieldParams.ngx);
    NSLog(@"NGY: %d", fieldParams.ngy);
    NSLog(@"DX: %e", fieldParams.dx);
    NSLog(@"DY: %e", fieldParams.dy);
    
    // BoundaryConditionForField の出力
    NSArray *fieldBoundaries = [init getFieldBoundaries];
    NSLog(@"--- BoundaryConditionForField (%lu entries) ---", (unsigned long)fieldBoundaries.count);
    for (NSUInteger i = 0; i < fieldBoundaries.count; i++) {
        NSValue *value = fieldBoundaries[i];
        struct BoundaryConditionForField boundary;
        [value getValue:&boundary];
        
        NSLog(@"[Boundary %lu]", (unsigned long)i + 1);
        NSLog(@"  Position: %@", boundary.position);
        NSLog(@"  Type: %@", boundary.type);
        NSLog(@"  Value: %e", boundary.val);
    }
    
    NSLog(@"==== パース結果出力終了 ====");
}