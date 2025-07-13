#import <Metal/Metal.h>
#import <simd/simd.h>
#import <fstream>
#import <string>
#import "Init.h"
#import "Particle.h"
#import "XmlLogger.h"
@class Particle;

// --- 64bit ミキサ関数（PCG 公式が推奨）
static inline uint64_t splitmix64(uint64_t &seed) {
    uint64_t z = (seed += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}
struct RNGState { 
    unsigned long state;
    unsigned long inc;
};

// 粒子生成用構造体
struct generationParams {
    uint addn;
    uint ele_pNum;
    uint ion_pNum;
    float Xmin;
    float Xmax;
    float Ymin;
    float Ymax;
    float ele_genU[3];
    float ele_vth;
    float ion_genU[3];
    float ion_vth;
    float dx;
    float dy;
    float ngb;
};

@interface Collision : NSObject

// metal関連
@property (nonatomic, strong) id<MTLDevice> device;
@property (nonatomic, strong) id<MTLCommandQueue> commandQueue;
@property (nonatomic, strong) id<MTLComputePipelineState> artificialIonizationPipeline;
@property (nonatomic, strong) id<MTLComputePipelineState> hollowCathodePipeline;
@property (nonatomic, strong) id<MTLBuffer> RNGStateBuffer;
@property (nonatomic, strong) id<MTLBuffer> paramsBuffer;
@property (nonatomic, strong) id<MTLBuffer> paramsForHCBuffer;
@property (nonatomic, strong) id<MTLBuffer> printBuffer;
@property (nonatomic) int pNumMax;
@property (nonatomic) uint threadGroupSize;
@property (nonatomic) uint integrationChunkSize;

// 関数
- (instancetype)initWithDevice:(id<MTLDevice>)device withParam:(Init*)initParam withParticles:(NSArray<Particle*>*)ptclArr withLogger:(XmlLogger&)logger;
- (bool)artificialIonization:(double)dt withParticles:(NSArray<Particle*>*)ptclArr withLogger:(XmlLogger&)logger;
- (bool)hollowCathode:(double)dt withParticles:(NSArray<Particle*>*)ptclArr withLogger:(XmlLogger&)logger;
@end