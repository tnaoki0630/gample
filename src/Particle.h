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
        double x;
        double y;
        double vx;
        double vy;
        double vz;
};

@interface Particle : NSObject

// metal関連
@property (nonatomic, strong) id<MTLDevice> device;
@property (nonatomic, strong) id<MTLCommandQueue> commandQueue;
@property (nonatomic, strong) id<MTLBuffer> particleBuffer;
@property (nonatomic, strong) id<MTLComputePipelineState> updateParticlesPipeline;
// スカラー変数
@property (nonatomic) double charge;
@property (nonatomic) double mass;
@property (nonatomic) double weight;

// 初期設定
- (instancetype)initWithDevice:(id<MTLDevice>)device
                withParticleParam:(ParamForParticle)ParticleParam
                withFieldParam:(ParamForField)FieldParam;
// 粒子生成
- (void)generateParticles:(ParamForParticle)ParticleParam
                withFieldParam:(ParamForField)FieldParam;
// 時間更新
- (void)update:(double)dt
        withEMField:(EMField*)fld ;

@end