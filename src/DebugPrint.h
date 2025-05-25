#import <Foundation/Foundation.h>
#import "Init.h"
#import "XmlLogger.h"

void printInitContents(Init* init, XmlLogger& logger){
    NSMutableString *log = [NSMutableString string];

    // FlagForEquation の出力
    struct FragForEquation flags = [init getFragForEquation];
    [log appendFormat:@"--- FlagForEquation ---\n"];
    [log appendFormat:@"Particle: %d\n", flags.Particle];
    [log appendFormat:@"EMField: %d\n", flags.EMField];
    [log appendFormat:@"MCCollision: %d\n", flags.MCCollision];
    
    // ParamForTimeIntegration の出力
    struct ParamForTimeIntegration timeParams = [init getParamForTimeIntegration];
    [log appendFormat:@"--- ParamForTimeIntegration ---\n"];
    [log appendFormat:@"EndCycle: %d\n", timeParams.EndCycle];
    [log appendFormat:@"ptclOutCycle: %d\n", timeParams.ptclOutCycle];
    [log appendFormat:@"fldOutCycle: %d\n", timeParams.fldOutCycle];
    [log appendFormat:@"dt: %e\n", timeParams.dt];
    
    // ParamForParticle の出力
    NSArray *particles = [init getParamForParticle];
    [log appendFormat:@"--- ParamForParticle (%lu entries) ---\n", (unsigned long)particles.count];
    for (NSUInteger i = 0; i < particles.count; i++) {
        NSValue *value = particles[i];
        struct ParamForParticle particle;
        [value getValue:&particle];
        
        [log appendFormat:@"[Particle %lu]\n", (unsigned long)i + 1];
        [log appendFormat:@"  Name: %@\n", particle.pName];
        [log appendFormat:@"  InitialParticleNumber: %d\n", particle.pNum];
        [log appendFormat:@"  MaxParticleNumber: %d\n", particle.pNumMax];
        [log appendFormat:@"  Charge: %e\n", particle.q];
        [log appendFormat:@"  Mass: %e\n", particle.m];
        [log appendFormat:@"  Weight: %e\n", particle.w];
        [log appendFormat:@"  genType: %@\n", particle.genType];
        [log appendFormat:@"  initialU: %e\n", particle.genU[0]];
        [log appendFormat:@"  initialV: %e\n", particle.genU[1]];
        [log appendFormat:@"  initialW: %e\n", particle.genU[2]];
        [log appendFormat:@"  initialT: %e\n", particle.genT];
    }
    
    // BoundaryConditionForParticle の出力
    NSArray *particleBoundaries = [init getParticleBoundaries];
    [log appendFormat:@"--- BoundaryConditionForParticle (%lu entries) ---\n", (unsigned long)particleBoundaries.count];
    for (NSUInteger i = 0; i < particleBoundaries.count; i++) {
        NSValue *value = particleBoundaries[i];
        struct BoundaryConditionForParticle boundary;
        [value getValue:&boundary];
        
        [log appendFormat:@"[Boundary %lu]\n", (unsigned long)i + 1];
        [log appendFormat:@"  Position: %@\n", boundary.position];
        [log appendFormat:@"  Type: %@\n", boundary.type];
        [log appendFormat:@"  Value: %e\n", boundary.val];
    }
    
    // SourceForParticle の出力
    NSArray *particleSources = [init getParticleSources];
    [log appendFormat:@"--- SourceForParticle (%lu entries) ---\n", (unsigned long)particleSources.count];
    for (NSUInteger i = 0; i < particleSources.count; i++) {
        NSValue *value = particleSources[i];
        struct SourceForParticle source;
        [value getValue:&source];
        
        [log appendFormat:@"[Source %lu]\n", (unsigned long)i + 1];
        [log appendFormat:@"  pName: %@\n", source.pName];
        [log appendFormat:@"  genType: %@\n", source.genType];
        [log appendFormat:@"  src: %e\n", source.src];
        [log appendFormat:@"  genXmin: %e\n", source.genX[0]];
        [log appendFormat:@"  genXmax: %e\n", source.genX[1]];
        [log appendFormat:@"  genYmin: %e\n", source.genY[0]];
        [log appendFormat:@"  genYmax: %e\n", source.genY[1]];
        [log appendFormat:@"  genU: %e\n", source.genU[0]];
        [log appendFormat:@"  genV: %e\n", source.genU[1]];
        [log appendFormat:@"  genW: %e\n", source.genU[2]];
        [log appendFormat:@"  genT: %e\n", source.genT];
    }
    
    // ParamForField の出力
    struct ParamForField fieldParams = [init getParamForField];
    [log appendFormat:@"--- ParamForField ---\n"];
    [log appendFormat:@"NGX: %d\n", fieldParams.ngx];
    [log appendFormat:@"NGY: %d\n", fieldParams.ngy];
    [log appendFormat:@"DX: %e\n", fieldParams.dx];
    [log appendFormat:@"DY: %e\n", fieldParams.dy];
    
    // BoundaryConditionForField の出力
    NSArray *fieldBoundaries = [init getFieldBoundaries];
    [log appendFormat:@"--- BoundaryConditionForField (%lu entries) ---\n", (unsigned long)fieldBoundaries.count];
    for (NSUInteger i = 0; i < fieldBoundaries.count; i++) {
        NSValue *value = fieldBoundaries[i];
        struct BoundaryConditionForField boundary;
        [value getValue:&boundary];
        
        [log appendFormat:@"[Boundary %lu]\n", (unsigned long)i + 1];
        [log appendFormat:@"  Position: %@\n", boundary.position];
        [log appendFormat:@"  Type: %@\n", boundary.type];
        [log appendFormat:@"  Value: %e\n", boundary.val];
    }

    logger.logComment([log UTF8String]);
}