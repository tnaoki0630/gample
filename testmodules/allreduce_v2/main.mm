// main.mm
#import <Foundation/Foundation.h>
#import <Metal/Metal.h>
#include <cstdio>

int main(int argc, const char * argv[]) {
    // 出力先をテキストファイルに変更
    // freopen("output.log", "w", stdout);
    // freopen("output.log", "w", stderr);

    @autoreleasepool {
        // Metalデバイスの取得
        id<MTLDevice> device = MTLCreateSystemDefaultDevice();
        if (!device) {
            NSLog(@"Metal is not supported on this device");
            return -1;
        }
        
        // Metalシェーダーコードを文字列リテラルとして定義
        NSString *shaderSource = @R"(
        #include <metal_stdlib>
        using namespace metal;

        kernel void integration(    device const float *in      [[buffer(0)]],
                                    device atomic_float *out    [[buffer(1)]],
                                    constant uint      &inSize  [[buffer(2)]],
                                    constant uint      &outSize [[buffer(3)]],
                                    uint gid [[thread_position_in_grid]])
        {
            // break
            if(gid >= inSize) return;
            
            // atomic addition
            uint idx_out = gid%outSize;
            atomic_fetch_add_explicit(&(out[idx_out]), in[gid], memory_order_relaxed);
        }
        )";
        
        // シェーダーソースからライブラリをコンパイル
        NSError *error = nil;
        MTLCompileOptions *options = [[MTLCompileOptions alloc] init];
        id<MTLLibrary> library = [device newLibraryWithSource:shaderSource options:options error:&error];
        if (!library) {
            NSLog(@"Failed to compile Metal shader library: %@", error);
            return -1;
        }
        
        // カーネル関数を取得
        id<MTLFunction> kernelFunction = [library newFunctionWithName:@"integration"];
        if (!kernelFunction) {
            NSLog(@"Failed to get kernel function");
            return -1;
        }
        
        // コンピュートパイプラインの作成
        id<MTLComputePipelineState> pipelineState = [device newComputePipelineStateWithFunction:kernelFunction error:&error];
        if (!pipelineState) {
            NSLog(@"Failed to create pipeline state: %@", error);
            return -1;
        }

        // コマンドキューを起動
        id<MTLCommandQueue> commandQueue = [device newCommandQueue];

        // 最大共有メモリサイズ確認
        NSUInteger maxThreadgroupMemory = [device maxThreadgroupMemoryLength];
        NSLog(@"Max threadgroup memory: %lu bytes", (unsigned long)maxThreadgroupMemory);
        
        // 入力バッファ
        uint dataCount = 10000 -2;
        id<MTLBuffer> dataCountBuffer = [device newBufferWithBytes:&dataCount length:sizeof(int) options:MTLResourceStorageModeShared];
        id<MTLBuffer> inputBuffer = [device newBufferWithLength:sizeof(float)*dataCount options:MTLResourceStorageModeShared];
        // 出力バッファ
        uint arrSize = 99;
        id<MTLBuffer> arrSizeBuffer = [device newBufferWithBytes:&arrSize length:sizeof(int) options:MTLResourceStorageModeShared];
        id<MTLBuffer> outputBuffer = [device newBufferWithLength:sizeof(float)*arrSize options:MTLResourceStorageModeShared];
        
        // 初期化
        float* input = (float*)inputBuffer.contents;
        for (uint i = 0; i < dataCount; i++) {
            input[i] = 1.0f;
        }
        float* output = (float*)outputBuffer.contents;
        for (uint i = 0; i < arrSize; i++){
            output[i] = 0.0f;
        }
            
        // コマンドバッファの作成
        id<MTLCommandBuffer> commandBuffer = [commandQueue commandBuffer];
        id<MTLComputeCommandEncoder> encoder = [commandBuffer computeCommandEncoder];
            
        // カーネルのセットアップ
        [encoder setComputePipelineState:pipelineState];
        [encoder setBuffer:inputBuffer offset:0 atIndex:0];
        [encoder setBuffer:outputBuffer offset:0 atIndex:1];
        [encoder setBuffer:dataCountBuffer offset:0 atIndex:2];
        [encoder setBuffer:arrSizeBuffer offset:0 atIndex:3];
        
        // ディスパッチ設定
        const uint threadgroupSize = 128;
        uint numThreadgroups = (dataCount + threadgroupSize -1)/threadgroupSize;
        MTLSize MTLgridSize = MTLSizeMake(numThreadgroups, 1, 1);
        MTLSize MTLthreadgroupSize = MTLSizeMake(threadgroupSize, 1, 1);
        [encoder dispatchThreadgroups:MTLgridSize threadsPerThreadgroup:MTLthreadgroupSize];
        
        // 実行
        [encoder endEncoding];
        [commandBuffer commit];
        [commandBuffer waitUntilCompleted];
        
        NSLog(@"gridSize = %d",numThreadgroups);

        // 結果確認
        float sumsum = 0.0f;
        output = (float *)outputBuffer.contents;
        for (int i = 0; i < arrSize; i++){
            sumsum += output[i];
            NSLog(@"Final output[%d]: %f", i, output[i]);
        }
        NSLog(@"Final sumsum: %f", sumsum);
        
    }
    return 0;
}
