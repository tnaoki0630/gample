#import <Foundation/Foundation.h>
#import <Metal/Metal.h>
#import "Init.h"
#import "Particle.h"
#import "XmlLogger.h"
@class Particle;

@interface EMField : NSObject
// 関数
- (instancetype)initWithDevice:(id<MTLDevice>)device withParam:(Init*)initParam withLogger:(XmlLogger&)logger;
- (void)outputCSRmtx:(int)n row:(std::vector<int>)row_ptr collumn:(std::vector<int>)col_idx value:(std::vector<float>)val;
- (void)solvePoisson:(XmlLogger&)logger;
- (void)resetChargeDensity;
- (void)checkChargeDensity;
- (bool)load1dField:(std::vector<float>&)field withFilePath:(NSString*)filePath;
- (bool)load2dField:(std::vector<float>&)field withFilePath:(NSString*)filePath;
// アクセサ
- (float*)phi;
- (id<MTLBuffer>)rhoBuffer;
- (id<MTLBuffer>)ExBuffer;
- (id<MTLBuffer>)EyBuffer;
- (id<MTLBuffer>)EzBuffer;
- (id<MTLBuffer>)BxBuffer;
- (id<MTLBuffer>)ByBuffer;
- (id<MTLBuffer>)BzBuffer;
- (int)ngx;
- (int)ngy;
- (int)nx;
- (int)ny;
- (int)ngy;
- (int)ngb;
- (double)dx;
- (double)dy;

@end