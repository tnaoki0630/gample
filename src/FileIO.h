#import <Foundation/Foundation.h>
#import "Init.h"
#import "Particle.h"
#import "EMField.h"
#import "Moment.h"
#import "XmlLogger.h"
@class EMField;
@class Moment;

void outputPhaseSpace(int cycle, Particle* ptcl, Init* init, XmlLogger& logger);
void outputField(int cycle, EMField* fld, Init* init, XmlLogger& logger);
void outputMoments(int cycle, Moment* mom, NSString* pName, Init* init, XmlLogger& logger);
static void writeField(FILE* fp, const char* name, int type_id, float* array, int arrSize, float scale);
BOOL saveProgress(int cycle, NSMutableArray* ptclArr, EMField* fld, Init* init);
BOOL loadProgress(int cycle, NSMutableArray* ptclArr, EMField* fld, Init* init);
BOOL saveBuffer(id<MTLBuffer> buffer, uint  length, uint lengthMax, NSString* outPath);
BOOL loadBuffer(id<MTLBuffer> buffer, uint& length, uint lengthMax, NSString* inPath);