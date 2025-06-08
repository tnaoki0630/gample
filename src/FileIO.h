#import <Foundation/Foundation.h>
#import "Particle.h"
#import "EMField.h"
#import "XmlLogger.h"
@class EMField;
@class Moment;

BOOL saveProgress(int cycle, NSMutableArray* ptclArr, EMField* fld, int  current, Init* init);
BOOL loadProgress(int cycle, NSMutableArray* ptclArr, EMField* fld, int& current, Init* init);
BOOL saveBuffer(id<MTLBuffer> buffer, uint  length, uint lengthMax, NSString* outPath);
BOOL loadBuffer(id<MTLBuffer> buffer, uint& length, uint lengthMax, NSString* inPath);