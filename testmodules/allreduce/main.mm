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
                device float *partialSums           [[ buffer(1) ]],
                device int &arrSize                 [[ buffer(2) ]],
                device int &kGroupNum               [[ buffer(3) ]],
                device int &kGroupSize              [[ buffer(4) ]],
                device int &chunkSize               [[ buffer(5) ]],
                device int &chunkOffset             [[ buffer(6) ]],
                uint gid                            [[ thread_position_in_grid ]],
                uint tid                            [[ thread_index_in_threadgroup ]],
                uint groupID                        [[ threadgroup_position_in_grid ]],
                threadgroup float* localSums        [[ threadgroup(0) ]]
            ) { 
                // インデックス計算
                int idx_in, idx_out;
                
                // 部分和を初期化
                for (int i = tid; i < kGroupSize*arrSize; i++) {
                    localSums[i] = 0.0f;
                }
                threadgroup_barrier(mem_flags::mem_threadgroup);

                // 入力データ取得
                for (int i = 0; i < arrSize; i++){
                    idx_in = gid + chunkOffset;
                    idx_out = tid + i*kGroupSize;
                    localSums[idx_out] = input[idx_in];
                }
                threadgroup_barrier(mem_flags::mem_threadgroup);
                
                // スレッドグループ内で段階的に和を計算（リダクション）
                for (int i = 0; i < arrSize; i++){
                    for (uint stride = kGroupSize / 2; stride > 0; stride /= 2) {
                        if (tid < stride) {
                            idx_in = tid + stride + i*kGroupSize;
                            idx_out = tid + i*kGroupSize;
                            localSums[idx_out] += localSums[idx_in];
                        }
                        threadgroup_barrier(mem_flags::mem_threadgroup);
                    }
                }
                
                // スレッドグループ内の最初のスレッドが部分和を出力バッファに書き込む
                if (tid == 0) {
                    for (int i = 0; i < arrSize; i++){
                        idx_in = 0 + i*kGroupSize;
                        idx_out = groupID + i*kGroupNum;
                        partialSums[idx_out] = localSums[idx_in];
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
        
        // カーネル関数(reductionKernel)を取得
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
        const int dataCount = pow(2, 15);
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
        const int numChunks = 32;
        const int chunkSize = dataCount/numChunks;
        id<MTLBuffer> chunkSizeBuffer = [device newBufferWithBytes:&chunkSize
                                                    length:sizeof(chunkSize)
                                                    options:MTLResourceStorageModeShared];
        // 和の配列サイズ
        id<MTLBuffer> arrSizeBuffer = [device newBufferWithBytes:&arrSize
                                                    length:sizeof(arrSize)
                                                    options:MTLResourceStorageModeShared];
        // チャンクごとにカーネルを dispatch して部分和を計算
        for (int chunk = 0; chunk < numChunks; chunk++){
            uint chunkOffset = chunk * chunkSize;
            // chunkOffset と chunkSize 用のバッファ（constant 引数）
            id<MTLBuffer> chunkOffsetBuffer = [device newBufferWithBytes:&chunkOffset
                                                    length:sizeof(chunkOffset)
                                                    options:MTLResourceStorageModeShared];

            // スレッドグループサイズとグループ数の設定
            const int threadsPerThreadgroup = 32;
            int numThreadgroups = chunkSize/threadsPerThreadgroup;
            // バッファ作成
            id<MTLBuffer> kGroupNumBuffer = [device newBufferWithBytes:&numThreadgroups
                                                        length:sizeof(numThreadgroups)
                                                        options:MTLResourceStorageModeShared];
            id<MTLBuffer> kGroupSizeBuffer = [device newBufferWithBytes:&threadsPerThreadgroup
                                                        length:sizeof(threadsPerThreadgroup)
                                                        options:MTLResourceStorageModeShared];

            // 部分結果を受け取るバッファを作成
            id<MTLBuffer> partialSumsBuffer = [device newBufferWithLength:sizeof(float)*arrSize*numThreadgroups
                                                    options:MTLResourceStorageModeShared];
            
            // コマンドバッファの作成
            id<MTLCommandBuffer> commandBuffer = [commandQueue commandBuffer];
            id<MTLComputeCommandEncoder> encoder = [commandBuffer computeCommandEncoder];
            
            // カーネルのセットアップ
            [encoder setComputePipelineState:pipelineState];
            [encoder setBuffer:inputBuffer offset:0 atIndex:0];
            [encoder setBuffer:partialSumsBuffer offset:0 atIndex:1];
            [encoder setBuffer:arrSizeBuffer offset:0 atIndex:2];
            [encoder setBuffer:kGroupNumBuffer offset:0 atIndex:3];
            [encoder setBuffer:kGroupSizeBuffer offset:0 atIndex:4];
            [encoder setBuffer:chunkSizeBuffer offset:0 atIndex:5];
            [encoder setBuffer:chunkOffsetBuffer offset:0 atIndex:6];
            
            // threadgroup メモリのサイズ設定
            NSUInteger tgMemSize = sizeof(float)*threadsPerThreadgroup*arrSize;
            [encoder setThreadgroupMemoryLength:tgMemSize atIndex:0];
            
            // ディスパッチ設定
            MTLSize gridSize = MTLSizeMake(chunkSize, 1, 1);
            MTLSize threadgroupSize = MTLSizeMake(threadsPerThreadgroup, 1, 1);
            [encoder dispatchThreads:gridSize threadsPerThreadgroup:threadgroupSize];
            [encoder endEncoding];
            [commandBuffer commit];
            [commandBuffer waitUntilCompleted];
            
            // チャンク毎に加算
            float *partialSums = (float *)partialSumsBuffer.contents;
            for (int i = 0; i < arrSize; i++){
                for (int j = 0; j < numThreadgroups; j++) {
                    Sums[i] += partialSums[j+i*numThreadgroups];
                }
            }
            NSLog(@"chunk %d: partialSums[%d,%d]: %f", chunk, 0,0, partialSums[0+0*numThreadgroups]);
            NSLog(@"chunk %d: Sums[%d]: %f", chunk, 1, Sums[1]);
            NSLog(@"chunk = %d, chunkSize = %d, numThreadgroups = %d, threadsPerThreadgroup = %d", chunk, chunkSize, numThreadgroups, threadsPerThreadgroup);
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
