#import <Metal/Metal.h>
#import <simd/simd.h>

// 粒子データ構造
struct Particle {
    vector_float2 position;
    vector_float2 velocity;
    float charge;
    float mass;
};

// 電場データ構造
struct Field {
    vector_float2 value;
};

@interface PIC : NSObject

@property (nonatomic, strong) id<MTLDevice> device;
@property (nonatomic, strong) id<MTLCommandQueue> commandQueue;
@property (nonatomic, strong) id<MTLBuffer> particleBuffer;
@property (nonatomic, strong) id<MTLBuffer> fieldBuffer;
@property (nonatomic, strong) id<MTLComputePipelineState> updateParticlesPipeline;
@property (nonatomic, strong) id<MTLComputePipelineState> computeFieldsPipeline;

- (instancetype)initWithDevice:(id<MTLDevice>)device
                particleCount:(NSUInteger)particleCount
                    gridSize:(NSUInteger)gridSize;
- (void)update:(float)deltaTime;

@end