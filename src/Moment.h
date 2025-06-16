#import <Foundation/Foundation.h>
#import <Metal/Metal.h>
#import "Init.h"
#import "EMField.h"
#import "Particle.h"
#import "XmlLogger.h"
@class Particle;
@class EMField;


// 積分計算用構造体
struct integrationParams {
    uint pNum;
    int ngx;
    int ngy;
    int ngb;
    float scale;
};

@interface Moment : NSObject

// metal関連
@property (nonatomic, strong) id<MTLDevice> device;
@property (nonatomic, strong) id<MTLCommandQueue> commandQueue;
@property (nonatomic, strong) id<MTLComputePipelineState> integrateMomentsPipeline;
@property (nonatomic, strong) id<MTLBuffer> integrationParamsBuffer;
@property (nonatomic, strong) id<MTLBuffer> nBuffer;
@property (nonatomic, strong) id<MTLBuffer> uxBuffer;
@property (nonatomic, strong) id<MTLBuffer> uyBuffer;
@property (nonatomic, strong) id<MTLBuffer> uzBuffer;
@property (nonatomic, strong) id<MTLBuffer> PxxBuffer;
@property (nonatomic, strong) id<MTLBuffer> PxyBuffer;
@property (nonatomic, strong) id<MTLBuffer> PxzBuffer;
@property (nonatomic, strong) id<MTLBuffer> PyyBuffer;
@property (nonatomic, strong) id<MTLBuffer> PyzBuffer;
@property (nonatomic, strong) id<MTLBuffer> PzzBuffer;
@property (nonatomic, strong) id<MTLBuffer> printBuffer;
@property (nonatomic) uint threadGroupSize;
@property (nonatomic) uint integrationChunkSize;

- (instancetype)initWithDevice:(id<MTLDevice>)device withParam:(Init*)initParam withLogger:(XmlLogger&)logger;
- (void)integrateMoments:(Particle*)ptcl withEMField:(EMField*)fld withLogger:(XmlLogger&)logger;
- (id<MTLBuffer>)nBuffer;
- (id<MTLBuffer>)uxBuffer;
- (id<MTLBuffer>)uyBuffer;
- (id<MTLBuffer>)uzBuffer;
- (id<MTLBuffer>)PxxBuffer;
- (id<MTLBuffer>)PxyBuffer;
- (id<MTLBuffer>)PxzBuffer;
- (id<MTLBuffer>)PyyBuffer;
- (id<MTLBuffer>)PyzBuffer;
- (id<MTLBuffer>)PzzBuffer;
- (id<MTLBuffer>)printBuffer;

@end