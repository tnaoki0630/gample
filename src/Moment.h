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
    uint ppt;
    int tgs;
    int ics;
};

@interface Moment : NSObject

// metal関連
@property (nonatomic, strong) id<MTLDevice> device;
@property (nonatomic, strong) id<MTLCommandQueue> commandQueue;
@property (nonatomic, strong) id<MTLComputePipelineState> integrateNumDensPipeline;
@property (nonatomic, strong) id<MTLComputePipelineState> integrateMeanVelXPipeline;
@property (nonatomic, strong) id<MTLComputePipelineState> integrateMeanVelYPipeline;
@property (nonatomic, strong) id<MTLComputePipelineState> integrateMeanVelZPipeline;
@property (nonatomic, strong) id<MTLComputePipelineState> integratePressureXXPipeline;
@property (nonatomic, strong) id<MTLComputePipelineState> integratePressureXYPipeline;
@property (nonatomic, strong) id<MTLComputePipelineState> integratePressureXZPipeline;
@property (nonatomic, strong) id<MTLComputePipelineState> integratePressureYYPipeline;
@property (nonatomic, strong) id<MTLComputePipelineState> integratePressureYZPipeline;
@property (nonatomic, strong) id<MTLComputePipelineState> integratePressureZZPipeline;
@property (nonatomic, strong) id<MTLBuffer> integrationParamsBuffer;
@property (nonatomic, strong) id<MTLBuffer> integrateTemporaryBuffer;
@property (nonatomic, strong) id<MTLBuffer> integratePartialBuffer;
@property (nonatomic, strong) id<MTLBuffer> printBuffer;
// スカラー変数
@property (nonatomic) uint threadGroupSize;
@property (nonatomic) uint integrationChunkSize;

- (instancetype)initWithDevice:(id<MTLDevice>)device withParam:(Init*)initParam withLogger:(XmlLogger&)logger;
- (void)integrateMoments:(Particle*)ptcl withEMField:(EMField*)fld withLogger:(XmlLogger&)logger;
- (void)outputMoments:(int)cycle withPtclName:(NSString*)pName withEMField:(EMField*)fld withLogger:(XmlLogger&)logger;
- (id<MTLBuffer>)integrateTemporaryBuffer;
- (id<MTLBuffer>)integratePartialBuffer;
- (id<MTLBuffer>)printBuffer;

@end