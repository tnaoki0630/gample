#import <Metal/Metal.h>
#import <simd/simd.h>
#import <fstream>
#import <string>
#import "Init.h"
#import "EMField.h"
#import "Moment.h"
#import "XmlLogger.h"
// コンパイルうまくいかないので前方宣言
@class EMField;
@class Moment;

// 粒子データ構造
struct ParticleState {
        float x;
        float y;
        float vx;
        float vy;
        float vz;
        int piflag;
};

// 電場データ構造体
struct EMFieldData {
        float* Ex;
        float* Ey;
        float* Ez;
        int ngx;
        int ngy;
        float dx;
        float dy;
};

// シミュレーションパラメータ構造体
struct SimulationParams {
        int pNum;
        float constE;
        float constB;
        float constX;
        float constY;
        int ngx;
        int ngy;
        int ngb;
};

/// \class Particle
/// \brief particle motion solver with Boris scheme
/// \ingroup solvers
@interface Particle : NSObject

// metal関連
@property (nonatomic, strong) id<MTLDevice> device;
@property (nonatomic, strong) id<MTLCommandQueue> commandQueue;
@property (nonatomic, strong) id<MTLComputePipelineState> updateParticlesPipeline;
@property (nonatomic, strong) id<MTLComputePipelineState> integrateChargeDensityPipeline;
@property (nonatomic, strong) id<MTLBuffer> particleBuffer;
@property (nonatomic, strong) id<MTLBuffer> paramsBuffer;
@property (nonatomic, strong) id<MTLBuffer> integrationParamsBuffer;
@property (nonatomic, strong) id<MTLBuffer> printBuffer;
@property (nonatomic) uint threadGroupSize;
@property (nonatomic) uint integrationChunkSize;
// 関数
- (instancetype)initWithDevice:(id<MTLDevice>)device withParam:(Init*)initParam specimen:(int)s withLogger:(XmlLogger&)logger;
- (void)update:(double)dt withEMField:(EMField*)fld withMom:(Moment*)mom withLogger:(XmlLogger&)logger;
- (void)reduce:(XmlLogger&)logger;
- (void)integrateChargeDensity:(EMField*)fld withMoment:(Moment*)mom withLogger:(XmlLogger&)logger;
- (int)injection:(double)dt withParam:(Init*)initParam withCurrent:(int&)current withLogger:(XmlLogger&)logger;
// アクセサ
- (NSString*)pName;
- (id<MTLBuffer>)paramsBuffer;
- (id<MTLBuffer>)particleBuffer;
- (float)m;
- (float)q;
- (float)w;
- (int)pinum_Xmin;
- (int)pinum_Xmax;
- (int)pinum_Ymin;
- (int)pinum_Ymax;
@end