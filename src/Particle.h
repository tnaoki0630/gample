#import <Metal/Metal.h>
#import <simd/simd.h>
#import <fstream>
#import <string>
#import "Init.h"
#import "EMField.h"
// コンパイルうまくいかないので前方宣言
@class EMField;

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
        uint pNum;
        float constE;
        float constB;
        float constX;
        float constY;
        int ngx;
        int ngy;
        int ngb;
        int BC_Xmin;
        int BC_Xmax;
        int BC_Ymin;
        int BC_Ymax;
};

@interface Particle : NSObject

// metal関連
@property (nonatomic, strong) id<MTLDevice> device;
@property (nonatomic, strong) id<MTLCommandQueue> commandQueue;
@property (nonatomic, strong) id<MTLComputePipelineState> updateParticlesPipeline;
@property (nonatomic, strong) id<MTLComputePipelineState> integrateChargeDensityPipeline;
@property (nonatomic, strong) id<MTLBuffer> particleBuffer;
@property (nonatomic, strong) id<MTLBuffer> integrateTemporaryBuffer;
@property (nonatomic, strong) id<MTLBuffer> integratePartialBuffer;
@property (nonatomic, strong) id<MTLBuffer> paramsBuffer;
@property (nonatomic, strong) id<MTLBuffer> printBuffer;
// スカラー変数
@property (nonatomic) uint threadGroupSize;
@property (nonatomic) uint integrationChunkSize;

// 初期設定
- (instancetype)initWithDevice:(id<MTLDevice>)device withParam:(Init*)initParam specimen:(int)s;
// 粒子生成
- (void)generateParticles:(ParamForParticle)ParticleParam withFieldParam:(ParamForField)FieldParam;
// 時間更新
- (void)update:(double)dt withEMField:(EMField*)fld;
// 粒子境界処理
- (void)reduce;
// 電荷密度更新
- (void)integrateChargeDensity:(EMField*)fld;
// 粒子軌道出力
- (void)outputPhaseSpace:(int)cycle withEMField:(EMField*)fld;
@end