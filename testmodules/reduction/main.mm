// main.mm
#import <Foundation/Foundation.h>
#import <Metal/Metal.h>
#include <cstdio>

int main(int argc, const char * argv[]) {
    // リダクション演算(array)
    uint kGroupSize = 128;
    float localSums[256];
    for (int tid = 0; tid < 256; tid++){
        localSums[tid] = 1.0f;
    }
    NSLog(@"\nreduction(array)");
    int tid_shifted;
    for (int i = 0; i < 2; i++){
        NSLog(@"i = %d",i);
        for (uint stride = kGroupSize / 2; stride > 0; stride /= 2) {
            for (int tid = 0; tid < kGroupSize; tid++){
                if (tid < stride) {
                    tid_shifted = tid + i*kGroupSize;
                    NSLog(@"true:  tid = %d, tid + stride = %d", tid, tid + stride);
                    NSLog(@"true:  tid_shifted = %d, tid_shifted + stride = %d", tid_shifted, tid_shifted + stride);
                    NSLog(@"localSums[tid_shifted] = %f, localSums[tid_shifted + stride] = %f", localSums[tid_shifted], localSums[tid_shifted + stride]);
                    localSums[tid_shifted] += localSums[tid_shifted + stride];
                    NSLog(@"localSums[tid_shifted] = %f, localSums[tid_shifted + stride] = %f", localSums[tid_shifted], localSums[tid_shifted + stride]);
                } else {
                    // NSLog(@"false: tid = %d, tid + stride = %d", tid, tid + stride);
                }
            }
        }
    }
}