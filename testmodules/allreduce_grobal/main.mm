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

            kernel void reductionKernel(
                device float *input                 [[ buffer(0) ]],
                device float *temp                  [[ buffer(1) ]],
                device float *output                [[ buffer(2) ]],
                device int &arrSize                 [[ buffer(3) ]],
                device int &chunkSize               [[ buffer(4) ]],
                device int &chunkOffset             [[ buffer(5) ]],
                device int &threadGroupSize         [[ buffer(6) ]],
                device int &dataCount               [[ buffer(7) ]],
                uint gid                            [[ thread_position_in_grid ]],
                uint tid                            [[ thread_index_in_threadgroup ]],
                uint groupID                        [[ threadgroup_position_in_grid ]]
            ) { 
                // gid が dataCount を超えたら処理はスキップ
                if (gid + chunkOffset >= dataCount) return;

                // インデックス計算
                int idx_in, idx_out, offset;

                // 入力データ取得
                for (int i = 0; i < arrSize; i++){
                    idx_in = gid + chunkOffset;
                    idx_out = gid + i*chunkSize;
                    temp[idx_out] = input[idx_in]*i;
                }
                threadgroup_barrier(mem_flags::mem_threadgroup);
                
                // スレッドグループ内で段階的に和を計算（リダクション）
                for (int i = 0; i < arrSize; i++){
                    for (uint stride = threadGroupSize / 2; stride > 0; stride /= 2) {
                        if (tid < stride) {
                            offset = groupID*threadGroupSize + i*chunkSize;
                            idx_in = tid + stride + offset;
                            idx_out = tid + offset;
                            temp[idx_out] += temp[idx_in];
                        }
                        threadgroup_barrier(mem_flags::mem_threadgroup);
                    }
                }

                // スレッドグループ内の最初のスレッドが部分和を出力バッファに書き込む
                if (tid == 0) {
                    for (int i = 0; i < arrSize; i++){
                        idx_in = 0 + groupID*threadGroupSize + i*chunkSize;
                        idx_out = groupID + i*chunkSize/threadGroupSize;
                        output[idx_out] = temp[idx_in];
                    }
                }
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
        id<MTLFunction> kernelFunction = [library newFunctionWithName:@"reductionKernel"];
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
        
        // サンプル用入力データ
        // const int dataCount = pow(2, 15);   // dataCount%chunkSize == 0
        const int dataCount = pow(2, 15)-100;    // dataCount%chunkSize != 0
        float inputData[dataCount];
        for (int i = 0; i < dataCount; i++) {
            // inputData[i] = 1.0f + fmod(i,2);
            inputData[i] = 1.0f;
            // NSLog(@"inputdata: %f", inputData[i]);
        }

        // 総和計算用配列
        const int arrSize = 256;
        float Sums[arrSize];
        for (int i = 0; i < arrSize; i++){
            Sums[i] = 0.0f;
        }
        
        // 入力バッファ作成
        id<MTLBuffer> inputBuffer = [device newBufferWithBytes:inputData
                                                    length:sizeof(inputData)
                                                    options:MTLResourceStorageModeShared];

        // チャンクに分割して積分
        const int chunkSize = 1024;
        const int numChunks = (dataCount + chunkSize - 1)/chunkSize; // dataCount > 1 なら numChunk > 1, dataCount = 1 なら numChunk = 0
        // dataCount 用バッファ（constant 引数）
        id<MTLBuffer> dataCountBuffer = [device newBufferWithBytes:&dataCount
                                                    length:sizeof(int)
                                                    options:MTLResourceStorageModeShared];
        // chunkSize 用バッファ（constant 引数）
        id<MTLBuffer> chunkSizeBuffer = [device newBufferWithBytes:&chunkSize
                                                    length:sizeof(int)
                                                    options:MTLResourceStorageModeShared];
        // 和の配列サイズ
        id<MTLBuffer> arrSizeBuffer = [device newBufferWithBytes:&arrSize
                                                    length:sizeof(arrSize)
                                                    options:MTLResourceStorageModeShared];
                                                    
        // チャンクごとにカーネルを dispatch して部分和を計算
        for (int chunk = 0; chunk < numChunks; chunk++){
            uint chunkOffset = chunk * chunkSize;
            // chunkOffset 用バッファ（constant 引数）
            id<MTLBuffer> chunkOffsetBuffer = [device newBufferWithBytes:&chunkOffset
                                                    length:sizeof(uint)
                                                    options:MTLResourceStorageModeShared];

            // スレッドグループサイズとグループ数の設定
            const int threadgroupSize = 128;
            int numThreadgroups = chunkSize/threadgroupSize;
            id<MTLBuffer> threadGroupSizeBuffer = [device newBufferWithBytes:&threadgroupSize
                                                        length:sizeof(int)
                                                        options:MTLResourceStorageModeShared];

            // 結果を受け取るバッファを作成
            id<MTLBuffer> temporaryBuffer = [device newBufferWithLength:sizeof(float)*arrSize*chunkSize
                                                    options:MTLResourceStorageModeShared];
            id<MTLBuffer> outputBuffer = [device newBufferWithLength:sizeof(float)*arrSize*numThreadgroups
                                                    options:MTLResourceStorageModeShared];
            
            // コマンドバッファの作成
            id<MTLCommandBuffer> commandBuffer = [commandQueue commandBuffer];
            id<MTLComputeCommandEncoder> encoder = [commandBuffer computeCommandEncoder];
            
            // カーネルのセットアップ
            [encoder setComputePipelineState:pipelineState];
            [encoder setBuffer:inputBuffer offset:0 atIndex:0];
            [encoder setBuffer:temporaryBuffer offset:0 atIndex:1];
            [encoder setBuffer:outputBuffer offset:0 atIndex:2];
            [encoder setBuffer:arrSizeBuffer offset:0 atIndex:3];
            [encoder setBuffer:chunkSizeBuffer offset:0 atIndex:4];
            [encoder setBuffer:chunkOffsetBuffer offset:0 atIndex:5];
            [encoder setBuffer:threadGroupSizeBuffer offset:0 atIndex:6];
            [encoder setBuffer:dataCountBuffer offset:0 atIndex:7];
            
            // ディスパッチ設定
            MTLSize MTLgridSize = MTLSizeMake(chunkSize, 1, 1);
            MTLSize MTLthreadgroupSize = MTLSizeMake(threadgroupSize, 1, 1);
            [encoder dispatchThreads:MTLgridSize threadsPerThreadgroup:MTLthreadgroupSize];
            [encoder endEncoding];
            // 実行
            [commandBuffer commit];
            [commandBuffer waitUntilCompleted];
            
            // スレッドグループ毎に加算
            float* partialSums = (float *)outputBuffer.contents;
            for (int i = 0; i < arrSize; i++){
                for (int j = 0; j < numThreadgroups; j++) {
                    Sums[i] += partialSums[j+i*numThreadgroups];
                }
            }
        }
        
        // 結果確認
        float sumsum = 0.0f;
        for (int i = 0; i < arrSize; i++){
            sumsum += Sums[i];
            NSLog(@"Final Sums[%d]: %f", i, Sums[i]);
        }
        NSLog(@"Final sumsum: %f", sumsum);
        
    }
    return 0;
}
