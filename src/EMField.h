#import <Foundation/Foundation.h>
#import <Metal/Metal.h>
#import "Init.h"
#import "Particle.h"
@class Particle;

@interface EMField : NSObject

- (instancetype)initWithDevice:(id<MTLDevice>)device withParam:(Init*)initParam ;
- (void)solvePoisson;
- (void)resetChargeDensity;
- (void)outputField:(int)cycle;

// 電荷密度へのアクセサ
- (float*)rho;

// Metal バッファへのアクセサ
- (id<MTLBuffer>)ExBuffer;
- (id<MTLBuffer>)EyBuffer;
- (id<MTLBuffer>)EzBuffer;
- (id<MTLBuffer>)BxBuffer;
- (id<MTLBuffer>)ByBuffer;
- (id<MTLBuffer>)BzBuffer;

// グリッド情報へのアクセサ
- (int)ngx;
- (int)ngy;
- (int)ngb;
- (double)dx;
- (double)dy;

@end