#import <Foundation/Foundation.h>
#import "Init.h"
#import "XmlLogger.h"

void printInitContents(Init* init, XmlLogger& logger){
    NSMutableString *log = [NSMutableString string];

    // FlagForEquation の出力
    struct FlagForEquation flags = init.flagForEquation;
    [log appendFormat:@"--- FlagForEquation ---\n"];
    [log appendFormat:@"Particle: %d\n", flags.Particle];
    [log appendFormat:@"EMField: %d\n", flags.EMField];
    [log appendFormat:@"MCCollision: %d\n", flags.MCCollision];
    
    // ParamForTimeIntegration の出力
    struct ParamForTimeIntegration timeParams = init.paramForTimeIntegration;
    [log appendFormat:@"--- ParamForTimeIntegration ---\n"];
    [log appendFormat:@"End: %d\n", timeParams.End];
    [log appendFormat:@"ParticleOutput: %d\n", timeParams.ParticleOutput];
    [log appendFormat:@"FieldOutput: %d\n", timeParams.FieldOutput];
    [log appendFormat:@"TimeStep: %e\n", timeParams.TimeStep];
    
    // ParamForComputing の出力
    struct ParamForComputing compParams = init.paramForComputing;
    [log appendFormat:@"--- ParamForComputing ---\n"];
    [log appendFormat:@"threadGroupSize: %d\n", compParams.threadGroupSize];
    [log appendFormat:@"integrationChunkSize: %d\n", compParams.integrationChunkSize];
    [log appendFormat:@"maxiter: %d\n", compParams.maxiter];
    [log appendFormat:@"tolerance: %e\n", compParams.tolerance];
    
    // ParamForParticle の出力
    std::vector<struct ParamForParticle> particles = init.paramForParticle;
    [log appendFormat:@"--- ParamForParticle (%zu entries) ---\n", (unsigned long)particles.size()];
    for (int i = 0; i < particles.size(); i++) {
        [log appendFormat:@"[Particle %lu]\n", (unsigned long)i + 1];
        [log appendFormat:@"  Name: %@\n", particles[i].pName];
        [log appendFormat:@"  InitialParticleNumber: %u\n", particles[i].pNum];
        [log appendFormat:@"  MaxParticleNumber: %u\n", particles[i].pNumMax];
        [log appendFormat:@"  Charge: %e\n", particles[i].q];
        [log appendFormat:@"  Mass: %e\n", particles[i].m];
        [log appendFormat:@"  Weight: %e\n", particles[i].w];
        [log appendFormat:@"  genType: %@\n", particles[i].genType];
        [log appendFormat:@"  initialT: %e\n", particles[i].genT];
    }
    
    // BoundaryConditionForParticle の出力
    std::vector<struct BoundaryConditionForParticle> pBCs = init.particleBoundaries;
    [log appendFormat:@"--- BoundaryConditionForParticle (%zu entries) ---\n", (unsigned long)pBCs.size()];
    for (int i = 0; i < pBCs.size(); i++) {
        [log appendFormat:@"[pBC %d]\n", i + 1];
        [log appendFormat:@"  Position: %@\n", pBCs[i].position];
        [log appendFormat:@"  Type: %@\n", pBCs[i].type];
    }

    // SourceForParticle の出力
    std::vector<struct SourceForParticle> sources = init.particleSources;
    [log appendFormat:@"--- SourceForParticle (%zu entries) ---\n", sources.size()];
    for (int i = 0; i < sources.size(); i++) {
        [log appendFormat:@"[Source %lu]\n", (unsigned long)i + 1];
        [log appendFormat:@"  pName: %@\n", sources[i].pName];
        [log appendFormat:@"  genType: %@\n", sources[i].genType];
        [log appendFormat:@"  src: %e\n", sources[i].src];
        [log appendFormat:@"  genXmin: %e\n", sources[i].genX[0]];
        [log appendFormat:@"  genXmax: %e\n", sources[i].genX[1]];
        [log appendFormat:@"  genYmin: %e\n", sources[i].genY[0]];
        [log appendFormat:@"  genYmax: %e\n", sources[i].genY[1]];
        [log appendFormat:@"  genU: %e\n", sources[i].genU[0]];
        [log appendFormat:@"  genV: %e\n", sources[i].genU[1]];
        [log appendFormat:@"  genW: %e\n", sources[i].genU[2]];
        [log appendFormat:@"  genT: %e\n", sources[i].genT];
    }
    
    // ParamForField の出力
    struct ParamForField fieldParams = init.paramForField;
    [log appendFormat:@"--- ParamForField ---\n"];
    [log appendFormat:@"ngx: %d\n", fieldParams.ngx];
    [log appendFormat:@"ngy: %d\n", fieldParams.ngy];
    [log appendFormat:@"ngb: %d\n", fieldParams.ngb];
    [log appendFormat:@"dx: %e\n", fieldParams.dx];
    [log appendFormat:@"dy: %e\n", fieldParams.dy];
    
    // BoundaryConditionForField の出力
    std::vector<struct BoundaryConditionForField> fBCs = init.fieldBoundaries;
    [log appendFormat:@"--- BoundaryConditionForField (%zu entries) ---\n", fBCs.size()];
    for (int i = 0; i < fBCs.size(); i++) {
        [log appendFormat:@"[Boundary %d]\n", i + 1];
        [log appendFormat:@"  Position: %@\n", fBCs[i].position];
        [log appendFormat:@"  Type: %@\n", fBCs[i].type];
        [log appendFormat:@"  Value: %e\n", fBCs[i].val];
    }

    logger.logComment([log UTF8String]);
}