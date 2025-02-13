#import <Foundation/Foundation.h>
#import "PIC.h"

int main(int argc, const char * argv[]) {
    @autoreleasepool {
        // デフォルトのMetalデバイスを取得
        id<MTLDevice> device = MTLCreateSystemDefaultDevice();
        if (!device) {
            NSLog(@"Metal is not supported on this device");
            return 1;
        }
        
        // シミュレーションの初期化
        NSUInteger particleCount = 1024;
        NSUInteger gridSize = 32;
        PIC *simulation = [[PIC alloc] 
                                      initWithDevice:device
                                      particleCount:particleCount
                                      gridSize:gridSize];
        
        // シミュレーションループ
        for (int i = 0; i < 100; i++) {  // 100フレーム分実行
            [simulation update:0.016f];
            if (i % 10 == 0) {  // 10ステップごとに出力
                [simulation writeVTKFile:@"plasma_particles" forTimestep:i];
                [simulation writeFieldVTKFile:@"plasma_field" forTimestep:i];
            }
            NSLog(@"Frame %d completed", i);
        }
    }
    return 0;
}