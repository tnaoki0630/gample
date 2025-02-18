#import <Metal/Metal.h>
#import <simd/simd.h>
#import <fstream>
#import <string>
#import "Init.h"
#import "EMField.h"

// 粒子データ構造
struct ParticleData {
    vector_double2 position;
    vector_double2 velocity;
    double charge;
    double mass;
    double weight;
};

@interface Particle : NSObject

@property (nonatomic, strong) id<MTLDevice> device;
@property (nonatomic, strong) id<MTLCommandQueue> commandQueue;
@property (nonatomic, strong) id<MTLBuffer> particleBuffer;
@property (nonatomic, strong) id<MTLComputePipelineState> updateParticlesPipeline;
// 初期設定
- (instancetype)initWithDevice:(id<MTLDevice>)device
                initWithParam:(Init*) init;
// 時間更新
- (void)update:(double)dt
        withEMField:(EMField*)fld ;

@end