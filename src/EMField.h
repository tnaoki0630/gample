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
- (bool)load1dField:(std::vector<float>&)field withFilePath:(NSString*)filePath;
- (bool)load2dField:(std::vector<float>&)field withFilePath:(NSString*)filePath;

// Metal バッファへのアクセサ
- (id<MTLBuffer>)rhoBuffer;
- (id<MTLBuffer>)ExBuffer;
- (id<MTLBuffer>)EyBuffer;
- (id<MTLBuffer>)EzBuffer;
- (id<MTLBuffer>)BxBuffer;
- (id<MTLBuffer>)ByBuffer;
- (id<MTLBuffer>)BzBuffer;

// グリッド情報へのアクセサ
- (int)ngx;
- (int)ngy;
- (int)nx;
- (int)ny;
- (int)ngy;
- (int)ngb;
- (double)dx;
- (double)dy;

@end