#import "Moment.h"
#import "Init.h"
#import "Constant.h"
#import "Particle.h"

// Metal シェーダーのソースコード
static NSString *const kMetalShaderSource = @R"(
#include <metal_stdlib>
using namespace metal;

struct ParticleState {
    float x;
    float y;
    float vx;
    float vy;
    float vz;
    int piflag;
};

// 積分計算用構造体
struct integrationParams {
    uint pNum;
    int ngx;
    int ngy;
    int ngb;
    uint ppt;
    int tgs;
    int ics;
};

// モーメント計算カーネル
kernel void integrateNumDens(
                        device ParticleState* ptcl          [[ buffer(0) ]],
                        constant integrationParams& prm     [[ buffer(1) ]],
                        device float* temp                  [[ buffer(2) ]],
                        device float* partial               [[ buffer(3) ]],
                        device float* print                 [[ buffer(4) ]],
                        uint gid                            [[ thread_position_in_grid ]],
                        uint tid                            [[ thread_index_in_threadgroup ]],
                        uint groupID                        [[ threadgroup_position_in_grid ]]
                        ) {
    // initialize(自スレッドがアクセスし得る範囲だけ)
    const int nx = prm.ngx + 2*prm.ngb;
    const int ny = prm.ngy + 2*prm.ngb;
    const int ng = (nx+1)*(ny+1);
    for (int i = 0; i < ng; i++){
        temp[i + gid*ng] = 0.0f;
    }

    // 変数定義
    int i1, j1;
    float hv[2];
    float sc;
    float sf[6][2];
    int ii, jj;

    // 積分ループ
    for (uint idx = 0; idx < prm.ppt; idx++){
        
        uint pid = idx + gid*prm.ppt;
        if (pid > prm.pNum) {
            // pNum を超えたら積分処理をスキップ
            break;
        }

        // 粒子を取得
        device ParticleState& p = ptcl[pid];

        // electro-magnetic field on each ptcl
        i1 = int(p.x);
        j1 = int(p.y);

        hv[0] = p.x - float(i1);
        hv[1] = p.y - float(j1);
        
        // 5th-order weighting
        for (int i = 0; i < 2; i++) {
            sc = 2.0 + hv[i];
            sf[0][i] = 1.0/120.0 *pow(3.0-sc, 5);
            sc = 1.0 + hv[i];
            sf[1][i] = 1.0/120.0 *(51.0 +75.0*sc -210.0*pow(sc,2) +150.0*pow(sc,3) -45.0*pow(sc,4) +5.0*pow(sc,5));
            sc = hv[i];
            sf[2][i] = 1.0/60.0 *(33.0 -30.0*pow(sc,2) +15.0*pow(sc,4) -5.0*pow(sc,5));
            sc = 1.0 - hv[i];
            sf[3][i] = 1.0/60.0 *(33.0 -30.0*pow(sc,2) +15.0*pow(sc,4) -5.0*pow(sc,5));
            sc = 2.0 - hv[i];
            sf[4][i] = 1.0/120.0 *(51.0 +75.0*sc -210.0*pow(sc,2) +150.0*pow(sc,3) -45.0*pow(sc,4) +5.0*pow(sc,5));
            sc = 3.0 - hv[i];
            sf[5][i] = 1.0/120.0 *pow(3.0-sc, 5);
        }

        // accumulation
        for (int i = 0; i < 6; i++) {
            for (int j = 0; j < 6; j++) {
                ii = i1+(i-prm.ngb);
                jj = j1+(j-prm.ngb);
                int idx_out = ii+jj*(nx+1) + gid*ng;
                temp[idx_out] += sf[i][0]*sf[j][1];
            }
        }
    
    }

    // スレッドグループ内で同期（未初期化値の参照を防ぐ）
    threadgroup_barrier(mem_flags::mem_threadgroup);

    // reduction
    for (int i = 0; i < ng; i++){
        for (uint stride = prm.tgs / 2; stride > 0; stride /= 2) {
            if (tid < stride) {
                uint offset = groupID*prm.tgs;
                int idx_in = i + (tid + stride + offset)*ng;
                int idx_out = i + (tid + offset)*ng;
                temp[idx_out] += temp[idx_in];
            }
            threadgroup_barrier(mem_flags::mem_threadgroup);
        }
    }
    
    // output partialSum
    if (tid == 0) {
        for (int i = 0; i < ng; i++){
            int idx_in = i + (0 + groupID*prm.tgs)*ng;
            int idx_out = i + groupID*ng;
            partial[idx_out] = temp[idx_in];
        }
    }
    
}

// モーメント計算カーネル
kernel void integrateMeanVelX(
                        device ParticleState* ptcl          [[ buffer(0) ]],
                        constant integrationParams& prm     [[ buffer(1) ]],
                        device float* temp                  [[ buffer(2) ]],
                        device float* partial               [[ buffer(3) ]],
                        device float* print                 [[ buffer(4) ]],
                        uint gid                            [[ thread_position_in_grid ]],
                        uint tid                            [[ thread_index_in_threadgroup ]],
                        uint groupID                        [[ threadgroup_position_in_grid ]]
                        ) {
    // initialize(自スレッドがアクセスし得る範囲だけ)
    const int nx = prm.ngx + 2*prm.ngb;
    const int ny = prm.ngy + 2*prm.ngb;
    const int ng = (nx+1)*(ny+1);
    for (int i = 0; i < ng; i++){
        temp[i + gid*ng] = 0.0f;
    }

    // 変数定義
    int i1, j1;
    float hv[2];
    float sc;
    float sf[6][2];
    int ii, jj;

    // 積分ループ
    for (uint idx = 0; idx < prm.ppt; idx++){
        
        uint pid = idx + gid*prm.ppt;
        if (pid > prm.pNum) {
            // pNum を超えたら積分処理をスキップ
            break;
        }

        // 粒子を取得
        device ParticleState& p = ptcl[pid];

        // electro-magnetic field on each ptcl
        i1 = int(p.x);
        j1 = int(p.y);

        hv[0] = p.x - float(i1);
        hv[1] = p.y - float(j1);
        
        // 5th-order weighting
        for (int i = 0; i < 2; i++) {
            sc = 2.0 + hv[i];
            sf[0][i] = 1.0/120.0 *pow(3.0-sc, 5);
            sc = 1.0 + hv[i];
            sf[1][i] = 1.0/120.0 *(51.0 +75.0*sc -210.0*pow(sc,2) +150.0*pow(sc,3) -45.0*pow(sc,4) +5.0*pow(sc,5));
            sc = hv[i];
            sf[2][i] = 1.0/60.0 *(33.0 -30.0*pow(sc,2) +15.0*pow(sc,4) -5.0*pow(sc,5));
            sc = 1.0 - hv[i];
            sf[3][i] = 1.0/60.0 *(33.0 -30.0*pow(sc,2) +15.0*pow(sc,4) -5.0*pow(sc,5));
            sc = 2.0 - hv[i];
            sf[4][i] = 1.0/120.0 *(51.0 +75.0*sc -210.0*pow(sc,2) +150.0*pow(sc,3) -45.0*pow(sc,4) +5.0*pow(sc,5));
            sc = 3.0 - hv[i];
            sf[5][i] = 1.0/120.0 *pow(3.0-sc, 5);
        }

        // accumulation
        for (int i = 0; i < 6; i++) {
            for (int j = 0; j < 6; j++) {
                ii = i1+(i-prm.ngb);
                jj = j1+(j-prm.ngb);
                int idx_out = ii+jj*(nx+1) + gid*ng;
                temp[idx_out] += sf[i][0]*sf[j][1]*p.vx;
            }
        }
    
    }

    // スレッドグループ内で同期（未初期化値の参照を防ぐ）
    threadgroup_barrier(mem_flags::mem_threadgroup);

    // reduction
    for (int i = 0; i < ng; i++){
        for (uint stride = prm.tgs / 2; stride > 0; stride /= 2) {
            if (tid < stride) {
                uint offset = groupID*prm.tgs;
                int idx_in = i + (tid + stride + offset)*ng;
                int idx_out = i + (tid + offset)*ng;
                temp[idx_out] += temp[idx_in];
            }
            threadgroup_barrier(mem_flags::mem_threadgroup);
        }
    }
    
    // output partialSum
    if (tid == 0) {
        for (int i = 0; i < ng; i++){
            int idx_in = i + (0 + groupID*prm.tgs)*ng;
            int idx_out = i + groupID*ng;
            partial[idx_out] = temp[idx_in];
        }
    }
    
}

// モーメント計算カーネル
kernel void integrateMeanVelY(
                        device ParticleState* ptcl          [[ buffer(0) ]],
                        constant integrationParams& prm     [[ buffer(1) ]],
                        device float* temp                  [[ buffer(2) ]],
                        device float* partial               [[ buffer(3) ]],
                        device float* print                 [[ buffer(4) ]],
                        uint gid                            [[ thread_position_in_grid ]],
                        uint tid                            [[ thread_index_in_threadgroup ]],
                        uint groupID                        [[ threadgroup_position_in_grid ]]
                        ) {
    // initialize(自スレッドがアクセスし得る範囲だけ)
    const int nx = prm.ngx + 2*prm.ngb;
    const int ny = prm.ngy + 2*prm.ngb;
    const int ng = (nx+1)*(ny+1);
    for (int i = 0; i < ng; i++){
        temp[i + gid*ng] = 0.0f;
    }

    // 変数定義
    int i1, j1;
    float hv[2];
    float sc;
    float sf[6][2];
    int ii, jj;

    // 積分ループ
    for (uint idx = 0; idx < prm.ppt; idx++){
        
        uint pid = idx + gid*prm.ppt;
        if (pid > prm.pNum) {
            // pNum を超えたら積分処理をスキップ
            break;
        }

        // 粒子を取得
        device ParticleState& p = ptcl[pid];

        // electro-magnetic field on each ptcl
        i1 = int(p.x);
        j1 = int(p.y);

        hv[0] = p.x - float(i1);
        hv[1] = p.y - float(j1);
        
        // 5th-order weighting
        for (int i = 0; i < 2; i++) {
            sc = 2.0 + hv[i];
            sf[0][i] = 1.0/120.0 *pow(3.0-sc, 5);
            sc = 1.0 + hv[i];
            sf[1][i] = 1.0/120.0 *(51.0 +75.0*sc -210.0*pow(sc,2) +150.0*pow(sc,3) -45.0*pow(sc,4) +5.0*pow(sc,5));
            sc = hv[i];
            sf[2][i] = 1.0/60.0 *(33.0 -30.0*pow(sc,2) +15.0*pow(sc,4) -5.0*pow(sc,5));
            sc = 1.0 - hv[i];
            sf[3][i] = 1.0/60.0 *(33.0 -30.0*pow(sc,2) +15.0*pow(sc,4) -5.0*pow(sc,5));
            sc = 2.0 - hv[i];
            sf[4][i] = 1.0/120.0 *(51.0 +75.0*sc -210.0*pow(sc,2) +150.0*pow(sc,3) -45.0*pow(sc,4) +5.0*pow(sc,5));
            sc = 3.0 - hv[i];
            sf[5][i] = 1.0/120.0 *pow(3.0-sc, 5);
        }

        // accumulation
        for (int i = 0; i < 6; i++) {
            for (int j = 0; j < 6; j++) {
                ii = i1+(i-prm.ngb);
                jj = j1+(j-prm.ngb);
                int idx_out = ii+jj*(nx+1) + gid*ng;
                temp[idx_out] += sf[i][0]*sf[j][1]*p.vy;
            }
        }
    
    }

    // スレッドグループ内で同期（未初期化値の参照を防ぐ）
    threadgroup_barrier(mem_flags::mem_threadgroup);

    // reduction
    for (int i = 0; i < ng; i++){
        for (uint stride = prm.tgs / 2; stride > 0; stride /= 2) {
            if (tid < stride) {
                uint offset = groupID*prm.tgs;
                int idx_in = i + (tid + stride + offset)*ng;
                int idx_out = i + (tid + offset)*ng;
                temp[idx_out] += temp[idx_in];
            }
            threadgroup_barrier(mem_flags::mem_threadgroup);
        }
    }
    
    // output partialSum
    if (tid == 0) {
        for (int i = 0; i < ng; i++){
            int idx_in = i + (0 + groupID*prm.tgs)*ng;
            int idx_out = i + groupID*ng;
            partial[idx_out] = temp[idx_in];
        }
    }
    
}

// モーメント計算カーネル
kernel void integrateMeanVelZ(
                        device ParticleState* ptcl          [[ buffer(0) ]],
                        constant integrationParams& prm     [[ buffer(1) ]],
                        device float* temp                  [[ buffer(2) ]],
                        device float* partial               [[ buffer(3) ]],
                        device float* print                 [[ buffer(4) ]],
                        uint gid                            [[ thread_position_in_grid ]],
                        uint tid                            [[ thread_index_in_threadgroup ]],
                        uint groupID                        [[ threadgroup_position_in_grid ]]
                        ) {
    // initialize(自スレッドがアクセスし得る範囲だけ)
    const int nx = prm.ngx + 2*prm.ngb;
    const int ny = prm.ngy + 2*prm.ngb;
    const int ng = (nx+1)*(ny+1);
    for (int i = 0; i < ng; i++){
        temp[i + gid*ng] = 0.0f;
    }

    // 変数定義
    int i1, j1;
    float hv[2];
    float sc;
    float sf[6][2];
    int ii, jj;

    // 積分ループ
    for (uint idx = 0; idx < prm.ppt; idx++){
        
        uint pid = idx + gid*prm.ppt;
        if (pid > prm.pNum) {
            // pNum を超えたら積分処理をスキップ
            break;
        }

        // 粒子を取得
        device ParticleState& p = ptcl[pid];

        // electro-magnetic field on each ptcl
        i1 = int(p.x);
        j1 = int(p.y);

        hv[0] = p.x - float(i1);
        hv[1] = p.y - float(j1);
        
        // 5th-order weighting
        for (int i = 0; i < 2; i++) {
            sc = 2.0 + hv[i];
            sf[0][i] = 1.0/120.0 *pow(3.0-sc, 5);
            sc = 1.0 + hv[i];
            sf[1][i] = 1.0/120.0 *(51.0 +75.0*sc -210.0*pow(sc,2) +150.0*pow(sc,3) -45.0*pow(sc,4) +5.0*pow(sc,5));
            sc = hv[i];
            sf[2][i] = 1.0/60.0 *(33.0 -30.0*pow(sc,2) +15.0*pow(sc,4) -5.0*pow(sc,5));
            sc = 1.0 - hv[i];
            sf[3][i] = 1.0/60.0 *(33.0 -30.0*pow(sc,2) +15.0*pow(sc,4) -5.0*pow(sc,5));
            sc = 2.0 - hv[i];
            sf[4][i] = 1.0/120.0 *(51.0 +75.0*sc -210.0*pow(sc,2) +150.0*pow(sc,3) -45.0*pow(sc,4) +5.0*pow(sc,5));
            sc = 3.0 - hv[i];
            sf[5][i] = 1.0/120.0 *pow(3.0-sc, 5);
        }

        // accumulation
        for (int i = 0; i < 6; i++) {
            for (int j = 0; j < 6; j++) {
                ii = i1+(i-prm.ngb);
                jj = j1+(j-prm.ngb);
                int idx_out = ii+jj*(nx+1) + gid*ng;
                temp[idx_out] += sf[i][0]*sf[j][1]*p.vz;
            }
        }
    
    }

    // スレッドグループ内で同期（未初期化値の参照を防ぐ）
    threadgroup_barrier(mem_flags::mem_threadgroup);

    // reduction
    for (int i = 0; i < ng; i++){
        for (uint stride = prm.tgs / 2; stride > 0; stride /= 2) {
            if (tid < stride) {
                uint offset = groupID*prm.tgs;
                int idx_in = i + (tid + stride + offset)*ng;
                int idx_out = i + (tid + offset)*ng;
                temp[idx_out] += temp[idx_in];
            }
            threadgroup_barrier(mem_flags::mem_threadgroup);
        }
    }
    
    // output partialSum
    if (tid == 0) {
        for (int i = 0; i < ng; i++){
            int idx_in = i + (0 + groupID*prm.tgs)*ng;
            int idx_out = i + groupID*ng;
            partial[idx_out] = temp[idx_in];
        }
    }
    
}

// モーメント計算カーネル
kernel void integratePressureXX(
                        device ParticleState* ptcl          [[ buffer(0) ]],
                        constant integrationParams& prm     [[ buffer(1) ]],
                        device float* temp                  [[ buffer(2) ]],
                        device float* partial               [[ buffer(3) ]],
                        device float* print                 [[ buffer(4) ]],
                        uint gid                            [[ thread_position_in_grid ]],
                        uint tid                            [[ thread_index_in_threadgroup ]],
                        uint groupID                        [[ threadgroup_position_in_grid ]]
                        ) {
    // initialize(自スレッドがアクセスし得る範囲だけ)
    const int nx = prm.ngx + 2*prm.ngb;
    const int ny = prm.ngy + 2*prm.ngb;
    const int ng = (nx+1)*(ny+1);
    for (int i = 0; i < ng; i++){
        temp[i + gid*ng] = 0.0f;
    }

    // 変数定義
    int i1, j1;
    float hv[2];
    float sc;
    float sf[6][2];
    int ii, jj;

    // 積分ループ
    for (uint idx = 0; idx < prm.ppt; idx++){
        
        uint pid = idx + gid*prm.ppt;
        if (pid > prm.pNum) {
            // pNum を超えたら積分処理をスキップ
            break;
        }

        // 粒子を取得
        device ParticleState& p = ptcl[pid];

        // electro-magnetic field on each ptcl
        i1 = int(p.x);
        j1 = int(p.y);

        hv[0] = p.x - float(i1);
        hv[1] = p.y - float(j1);
        
        // 5th-order weighting
        for (int i = 0; i < 2; i++) {
            sc = 2.0 + hv[i];
            sf[0][i] = 1.0/120.0 *pow(3.0-sc, 5);
            sc = 1.0 + hv[i];
            sf[1][i] = 1.0/120.0 *(51.0 +75.0*sc -210.0*pow(sc,2) +150.0*pow(sc,3) -45.0*pow(sc,4) +5.0*pow(sc,5));
            sc = hv[i];
            sf[2][i] = 1.0/60.0 *(33.0 -30.0*pow(sc,2) +15.0*pow(sc,4) -5.0*pow(sc,5));
            sc = 1.0 - hv[i];
            sf[3][i] = 1.0/60.0 *(33.0 -30.0*pow(sc,2) +15.0*pow(sc,4) -5.0*pow(sc,5));
            sc = 2.0 - hv[i];
            sf[4][i] = 1.0/120.0 *(51.0 +75.0*sc -210.0*pow(sc,2) +150.0*pow(sc,3) -45.0*pow(sc,4) +5.0*pow(sc,5));
            sc = 3.0 - hv[i];
            sf[5][i] = 1.0/120.0 *pow(3.0-sc, 5);
        }

        // accumulation
        for (int i = 0; i < 6; i++) {
            for (int j = 0; j < 6; j++) {
                ii = i1+(i-prm.ngb);
                jj = j1+(j-prm.ngb);
                int idx_out = ii+jj*(nx+1) + gid*ng;
                temp[idx_out] += sf[i][0]*sf[j][1]*p.vx*p.vx;
            }
        }
    
    }

    // スレッドグループ内で同期（未初期化値の参照を防ぐ）
    threadgroup_barrier(mem_flags::mem_threadgroup);

    // reduction
    for (int i = 0; i < ng; i++){
        for (uint stride = prm.tgs / 2; stride > 0; stride /= 2) {
            if (tid < stride) {
                uint offset = groupID*prm.tgs;
                int idx_in = i + (tid + stride + offset)*ng;
                int idx_out = i + (tid + offset)*ng;
                temp[idx_out] += temp[idx_in];
            }
            threadgroup_barrier(mem_flags::mem_threadgroup);
        }
    }
    
    // output partialSum
    if (tid == 0) {
        for (int i = 0; i < ng; i++){
            int idx_in = i + (0 + groupID*prm.tgs)*ng;
            int idx_out = i + groupID*ng;
            partial[idx_out] = temp[idx_in];
        }
    }
    
}

// モーメント計算カーネル
kernel void integratePressureXY(
                        device ParticleState* ptcl          [[ buffer(0) ]],
                        constant integrationParams& prm     [[ buffer(1) ]],
                        device float* temp                  [[ buffer(2) ]],
                        device float* partial               [[ buffer(3) ]],
                        device float* print                 [[ buffer(4) ]],
                        uint gid                            [[ thread_position_in_grid ]],
                        uint tid                            [[ thread_index_in_threadgroup ]],
                        uint groupID                        [[ threadgroup_position_in_grid ]]
                        ) {
    // initialize(自スレッドがアクセスし得る範囲だけ)
    const int nx = prm.ngx + 2*prm.ngb;
    const int ny = prm.ngy + 2*prm.ngb;
    const int ng = (nx+1)*(ny+1);
    for (int i = 0; i < ng; i++){
        temp[i + gid*ng] = 0.0f;
    }

    // 変数定義
    int i1, j1;
    float hv[2];
    float sc;
    float sf[6][2];
    int ii, jj;

    // 積分ループ
    for (uint idx = 0; idx < prm.ppt; idx++){
        
        uint pid = idx + gid*prm.ppt;
        if (pid > prm.pNum) {
            // pNum を超えたら積分処理をスキップ
            break;
        }

        // 粒子を取得
        device ParticleState& p = ptcl[pid];

        // electro-magnetic field on each ptcl
        i1 = int(p.x);
        j1 = int(p.y);

        hv[0] = p.x - float(i1);
        hv[1] = p.y - float(j1);
        
        // 5th-order weighting
        for (int i = 0; i < 2; i++) {
            sc = 2.0 + hv[i];
            sf[0][i] = 1.0/120.0 *pow(3.0-sc, 5);
            sc = 1.0 + hv[i];
            sf[1][i] = 1.0/120.0 *(51.0 +75.0*sc -210.0*pow(sc,2) +150.0*pow(sc,3) -45.0*pow(sc,4) +5.0*pow(sc,5));
            sc = hv[i];
            sf[2][i] = 1.0/60.0 *(33.0 -30.0*pow(sc,2) +15.0*pow(sc,4) -5.0*pow(sc,5));
            sc = 1.0 - hv[i];
            sf[3][i] = 1.0/60.0 *(33.0 -30.0*pow(sc,2) +15.0*pow(sc,4) -5.0*pow(sc,5));
            sc = 2.0 - hv[i];
            sf[4][i] = 1.0/120.0 *(51.0 +75.0*sc -210.0*pow(sc,2) +150.0*pow(sc,3) -45.0*pow(sc,4) +5.0*pow(sc,5));
            sc = 3.0 - hv[i];
            sf[5][i] = 1.0/120.0 *pow(3.0-sc, 5);
        }

        // accumulation
        for (int i = 0; i < 6; i++) {
            for (int j = 0; j < 6; j++) {
                ii = i1+(i-prm.ngb);
                jj = j1+(j-prm.ngb);
                int idx_out = ii+jj*(nx+1) + gid*ng;
                temp[idx_out] += sf[i][0]*sf[j][1]*p.vx*p.vy;
            }
        }
    
    }

    // スレッドグループ内で同期（未初期化値の参照を防ぐ）
    threadgroup_barrier(mem_flags::mem_threadgroup);

    // reduction
    for (int i = 0; i < ng; i++){
        for (uint stride = prm.tgs / 2; stride > 0; stride /= 2) {
            if (tid < stride) {
                uint offset = groupID*prm.tgs;
                int idx_in = i + (tid + stride + offset)*ng;
                int idx_out = i + (tid + offset)*ng;
                temp[idx_out] += temp[idx_in];
            }
            threadgroup_barrier(mem_flags::mem_threadgroup);
        }
    }
    
    // output partialSum
    if (tid == 0) {
        for (int i = 0; i < ng; i++){
            int idx_in = i + (0 + groupID*prm.tgs)*ng;
            int idx_out = i + groupID*ng;
            partial[idx_out] = temp[idx_in];
        }
    }
    
}

// モーメント計算カーネル
kernel void integratePressureXZ(
                        device ParticleState* ptcl          [[ buffer(0) ]],
                        constant integrationParams& prm     [[ buffer(1) ]],
                        device float* temp                  [[ buffer(2) ]],
                        device float* partial               [[ buffer(3) ]],
                        device float* print                 [[ buffer(4) ]],
                        uint gid                            [[ thread_position_in_grid ]],
                        uint tid                            [[ thread_index_in_threadgroup ]],
                        uint groupID                        [[ threadgroup_position_in_grid ]]
                        ) {
    // initialize(自スレッドがアクセスし得る範囲だけ)
    const int nx = prm.ngx + 2*prm.ngb;
    const int ny = prm.ngy + 2*prm.ngb;
    const int ng = (nx+1)*(ny+1);
    for (int i = 0; i < ng; i++){
        temp[i + gid*ng] = 0.0f;
    }

    // 変数定義
    int i1, j1;
    float hv[2];
    float sc;
    float sf[6][2];
    int ii, jj;

    // 積分ループ
    for (uint idx = 0; idx < prm.ppt; idx++){
        
        uint pid = idx + gid*prm.ppt;
        if (pid > prm.pNum) {
            // pNum を超えたら積分処理をスキップ
            break;
        }

        // 粒子を取得
        device ParticleState& p = ptcl[pid];

        // electro-magnetic field on each ptcl
        i1 = int(p.x);
        j1 = int(p.y);

        hv[0] = p.x - float(i1);
        hv[1] = p.y - float(j1);
        
        // 5th-order weighting
        for (int i = 0; i < 2; i++) {
            sc = 2.0 + hv[i];
            sf[0][i] = 1.0/120.0 *pow(3.0-sc, 5);
            sc = 1.0 + hv[i];
            sf[1][i] = 1.0/120.0 *(51.0 +75.0*sc -210.0*pow(sc,2) +150.0*pow(sc,3) -45.0*pow(sc,4) +5.0*pow(sc,5));
            sc = hv[i];
            sf[2][i] = 1.0/60.0 *(33.0 -30.0*pow(sc,2) +15.0*pow(sc,4) -5.0*pow(sc,5));
            sc = 1.0 - hv[i];
            sf[3][i] = 1.0/60.0 *(33.0 -30.0*pow(sc,2) +15.0*pow(sc,4) -5.0*pow(sc,5));
            sc = 2.0 - hv[i];
            sf[4][i] = 1.0/120.0 *(51.0 +75.0*sc -210.0*pow(sc,2) +150.0*pow(sc,3) -45.0*pow(sc,4) +5.0*pow(sc,5));
            sc = 3.0 - hv[i];
            sf[5][i] = 1.0/120.0 *pow(3.0-sc, 5);
        }

        // accumulation
        for (int i = 0; i < 6; i++) {
            for (int j = 0; j < 6; j++) {
                ii = i1+(i-prm.ngb);
                jj = j1+(j-prm.ngb);
                int idx_out = ii+jj*(nx+1) + gid*ng;
                temp[idx_out] += sf[i][0]*sf[j][1]*p.vx*p.vz;
            }
        }
    
    }

    // スレッドグループ内で同期（未初期化値の参照を防ぐ）
    threadgroup_barrier(mem_flags::mem_threadgroup);

    // reduction
    for (int i = 0; i < ng; i++){
        for (uint stride = prm.tgs / 2; stride > 0; stride /= 2) {
            if (tid < stride) {
                uint offset = groupID*prm.tgs;
                int idx_in = i + (tid + stride + offset)*ng;
                int idx_out = i + (tid + offset)*ng;
                temp[idx_out] += temp[idx_in];
            }
            threadgroup_barrier(mem_flags::mem_threadgroup);
        }
    }
    
    // output partialSum
    if (tid == 0) {
        for (int i = 0; i < ng; i++){
            int idx_in = i + (0 + groupID*prm.tgs)*ng;
            int idx_out = i + groupID*ng;
            partial[idx_out] = temp[idx_in];
        }
    }
    
}

// モーメント計算カーネル
kernel void integratePressureYY(
                        device ParticleState* ptcl          [[ buffer(0) ]],
                        constant integrationParams& prm     [[ buffer(1) ]],
                        device float* temp                  [[ buffer(2) ]],
                        device float* partial               [[ buffer(3) ]],
                        device float* print                 [[ buffer(4) ]],
                        uint gid                            [[ thread_position_in_grid ]],
                        uint tid                            [[ thread_index_in_threadgroup ]],
                        uint groupID                        [[ threadgroup_position_in_grid ]]
                        ) {
    // initialize(自スレッドがアクセスし得る範囲だけ)
    const int nx = prm.ngx + 2*prm.ngb;
    const int ny = prm.ngy + 2*prm.ngb;
    const int ng = (nx+1)*(ny+1);
    for (int i = 0; i < ng; i++){
        temp[i + gid*ng] = 0.0f;
    }

    // 変数定義
    int i1, j1;
    float hv[2];
    float sc;
    float sf[6][2];
    int ii, jj;

    // 積分ループ
    for (uint idx = 0; idx < prm.ppt; idx++){
        
        uint pid = idx + gid*prm.ppt;
        if (pid > prm.pNum) {
            // pNum を超えたら積分処理をスキップ
            break;
        }

        // 粒子を取得
        device ParticleState& p = ptcl[pid];

        // electro-magnetic field on each ptcl
        i1 = int(p.x);
        j1 = int(p.y);

        hv[0] = p.x - float(i1);
        hv[1] = p.y - float(j1);
        
        // 5th-order weighting
        for (int i = 0; i < 2; i++) {
            sc = 2.0 + hv[i];
            sf[0][i] = 1.0/120.0 *pow(3.0-sc, 5);
            sc = 1.0 + hv[i];
            sf[1][i] = 1.0/120.0 *(51.0 +75.0*sc -210.0*pow(sc,2) +150.0*pow(sc,3) -45.0*pow(sc,4) +5.0*pow(sc,5));
            sc = hv[i];
            sf[2][i] = 1.0/60.0 *(33.0 -30.0*pow(sc,2) +15.0*pow(sc,4) -5.0*pow(sc,5));
            sc = 1.0 - hv[i];
            sf[3][i] = 1.0/60.0 *(33.0 -30.0*pow(sc,2) +15.0*pow(sc,4) -5.0*pow(sc,5));
            sc = 2.0 - hv[i];
            sf[4][i] = 1.0/120.0 *(51.0 +75.0*sc -210.0*pow(sc,2) +150.0*pow(sc,3) -45.0*pow(sc,4) +5.0*pow(sc,5));
            sc = 3.0 - hv[i];
            sf[5][i] = 1.0/120.0 *pow(3.0-sc, 5);
        }

        // accumulation
        for (int i = 0; i < 6; i++) {
            for (int j = 0; j < 6; j++) {
                ii = i1+(i-prm.ngb);
                jj = j1+(j-prm.ngb);
                int idx_out = ii+jj*(nx+1) + gid*ng;
                temp[idx_out] += sf[i][0]*sf[j][1]*p.vy*p.vy;
            }
        }
    
    }

    // スレッドグループ内で同期（未初期化値の参照を防ぐ）
    threadgroup_barrier(mem_flags::mem_threadgroup);

    // reduction
    for (int i = 0; i < ng; i++){
        for (uint stride = prm.tgs / 2; stride > 0; stride /= 2) {
            if (tid < stride) {
                uint offset = groupID*prm.tgs;
                int idx_in = i + (tid + stride + offset)*ng;
                int idx_out = i + (tid + offset)*ng;
                temp[idx_out] += temp[idx_in];
            }
            threadgroup_barrier(mem_flags::mem_threadgroup);
        }
    }
    
    // output partialSum
    if (tid == 0) {
        for (int i = 0; i < ng; i++){
            int idx_in = i + (0 + groupID*prm.tgs)*ng;
            int idx_out = i + groupID*ng;
            partial[idx_out] = temp[idx_in];
        }
    }
    
}

// モーメント計算カーネル
kernel void integratePressureYZ(
                        device ParticleState* ptcl          [[ buffer(0) ]],
                        constant integrationParams& prm     [[ buffer(1) ]],
                        device float* temp                  [[ buffer(2) ]],
                        device float* partial               [[ buffer(3) ]],
                        device float* print                 [[ buffer(4) ]],
                        uint gid                            [[ thread_position_in_grid ]],
                        uint tid                            [[ thread_index_in_threadgroup ]],
                        uint groupID                        [[ threadgroup_position_in_grid ]]
                        ) {
    // initialize(自スレッドがアクセスし得る範囲だけ)
    const int nx = prm.ngx + 2*prm.ngb;
    const int ny = prm.ngy + 2*prm.ngb;
    const int ng = (nx+1)*(ny+1);
    for (int i = 0; i < ng; i++){
        temp[i + gid*ng] = 0.0f;
    }

    // 変数定義
    int i1, j1;
    float hv[2];
    float sc;
    float sf[6][2];
    int ii, jj;

    // 積分ループ
    for (uint idx = 0; idx < prm.ppt; idx++){
        
        uint pid = idx + gid*prm.ppt;
        if (pid > prm.pNum) {
            // pNum を超えたら積分処理をスキップ
            break;
        }

        // 粒子を取得
        device ParticleState& p = ptcl[pid];

        // electro-magnetic field on each ptcl
        i1 = int(p.x);
        j1 = int(p.y);

        hv[0] = p.x - float(i1);
        hv[1] = p.y - float(j1);
        
        // 5th-order weighting
        for (int i = 0; i < 2; i++) {
            sc = 2.0 + hv[i];
            sf[0][i] = 1.0/120.0 *pow(3.0-sc, 5);
            sc = 1.0 + hv[i];
            sf[1][i] = 1.0/120.0 *(51.0 +75.0*sc -210.0*pow(sc,2) +150.0*pow(sc,3) -45.0*pow(sc,4) +5.0*pow(sc,5));
            sc = hv[i];
            sf[2][i] = 1.0/60.0 *(33.0 -30.0*pow(sc,2) +15.0*pow(sc,4) -5.0*pow(sc,5));
            sc = 1.0 - hv[i];
            sf[3][i] = 1.0/60.0 *(33.0 -30.0*pow(sc,2) +15.0*pow(sc,4) -5.0*pow(sc,5));
            sc = 2.0 - hv[i];
            sf[4][i] = 1.0/120.0 *(51.0 +75.0*sc -210.0*pow(sc,2) +150.0*pow(sc,3) -45.0*pow(sc,4) +5.0*pow(sc,5));
            sc = 3.0 - hv[i];
            sf[5][i] = 1.0/120.0 *pow(3.0-sc, 5);
        }

        // accumulation
        for (int i = 0; i < 6; i++) {
            for (int j = 0; j < 6; j++) {
                ii = i1+(i-prm.ngb);
                jj = j1+(j-prm.ngb);
                int idx_out = ii+jj*(nx+1) + gid*ng;
                temp[idx_out] += sf[i][0]*sf[j][1]*p.vy*p.vz;
            }
        }
    
    }

    // スレッドグループ内で同期（未初期化値の参照を防ぐ）
    threadgroup_barrier(mem_flags::mem_threadgroup);

    // reduction
    for (int i = 0; i < ng; i++){
        for (uint stride = prm.tgs / 2; stride > 0; stride /= 2) {
            if (tid < stride) {
                uint offset = groupID*prm.tgs;
                int idx_in = i + (tid + stride + offset)*ng;
                int idx_out = i + (tid + offset)*ng;
                temp[idx_out] += temp[idx_in];
            }
            threadgroup_barrier(mem_flags::mem_threadgroup);
        }
    }
    
    // output partialSum
    if (tid == 0) {
        for (int i = 0; i < ng; i++){
            int idx_in = i + (0 + groupID*prm.tgs)*ng;
            int idx_out = i + groupID*ng;
            partial[idx_out] = temp[idx_in];
        }
    }
    
}

// モーメント計算カーネル
kernel void integratePressureZZ(
                        device ParticleState* ptcl          [[ buffer(0) ]],
                        constant integrationParams& prm     [[ buffer(1) ]],
                        device float* temp                  [[ buffer(2) ]],
                        device float* partial               [[ buffer(3) ]],
                        device float* print                 [[ buffer(4) ]],
                        uint gid                            [[ thread_position_in_grid ]],
                        uint tid                            [[ thread_index_in_threadgroup ]],
                        uint groupID                        [[ threadgroup_position_in_grid ]]
                        ) {
    // initialize(自スレッドがアクセスし得る範囲だけ)
    const int nx = prm.ngx + 2*prm.ngb;
    const int ny = prm.ngy + 2*prm.ngb;
    const int ng = (nx+1)*(ny+1);
    for (int i = 0; i < ng; i++){
        temp[i + gid*ng] = 0.0f;
    }

    // 変数定義
    int i1, j1;
    float hv[2];
    float sc;
    float sf[6][2];
    int ii, jj;

    // 積分ループ
    for (uint idx = 0; idx < prm.ppt; idx++){
        
        uint pid = idx + gid*prm.ppt;
        if (pid > prm.pNum) {
            // pNum を超えたら積分処理をスキップ
            break;
        }

        // 粒子を取得
        device ParticleState& p = ptcl[pid];

        // electro-magnetic field on each ptcl
        i1 = int(p.x);
        j1 = int(p.y);

        hv[0] = p.x - float(i1);
        hv[1] = p.y - float(j1);
        
        // 5th-order weighting
        for (int i = 0; i < 2; i++) {
            sc = 2.0 + hv[i];
            sf[0][i] = 1.0/120.0 *pow(3.0-sc, 5);
            sc = 1.0 + hv[i];
            sf[1][i] = 1.0/120.0 *(51.0 +75.0*sc -210.0*pow(sc,2) +150.0*pow(sc,3) -45.0*pow(sc,4) +5.0*pow(sc,5));
            sc = hv[i];
            sf[2][i] = 1.0/60.0 *(33.0 -30.0*pow(sc,2) +15.0*pow(sc,4) -5.0*pow(sc,5));
            sc = 1.0 - hv[i];
            sf[3][i] = 1.0/60.0 *(33.0 -30.0*pow(sc,2) +15.0*pow(sc,4) -5.0*pow(sc,5));
            sc = 2.0 - hv[i];
            sf[4][i] = 1.0/120.0 *(51.0 +75.0*sc -210.0*pow(sc,2) +150.0*pow(sc,3) -45.0*pow(sc,4) +5.0*pow(sc,5));
            sc = 3.0 - hv[i];
            sf[5][i] = 1.0/120.0 *pow(3.0-sc, 5);
        }

        // accumulation
        for (int i = 0; i < 6; i++) {
            for (int j = 0; j < 6; j++) {
                ii = i1+(i-prm.ngb);
                jj = j1+(j-prm.ngb);
                int idx_out = ii+jj*(nx+1) + gid*ng;
                temp[idx_out] += sf[i][0]*sf[j][1]*p.vz*p.vz;
            }
        }
    
    }

    // スレッドグループ内で同期（未初期化値の参照を防ぐ）
    threadgroup_barrier(mem_flags::mem_threadgroup);

    // reduction
    for (int i = 0; i < ng; i++){
        for (uint stride = prm.tgs / 2; stride > 0; stride /= 2) {
            if (tid < stride) {
                uint offset = groupID*prm.tgs;
                int idx_in = i + (tid + stride + offset)*ng;
                int idx_out = i + (tid + offset)*ng;
                temp[idx_out] += temp[idx_in];
            }
            threadgroup_barrier(mem_flags::mem_threadgroup);
        }
    }
    
    // output partialSum
    if (tid == 0) {
        for (int i = 0; i < ng; i++){
            int idx_in = i + (0 + groupID*prm.tgs)*ng;
            int idx_out = i + groupID*ng;
            partial[idx_out] = temp[idx_in];
        }
    }
    
}
)";

// プライベートインスタンス変数
@implementation Moment {
    float* _n;
    float* _ux;
    float* _uy;
    float* _uz;
    float* _Pxx;
    float* _Pxy;
    float* _Pxz;
    float* _Pyy;
    float* _Pyz;
    float* _Pzz;
}

- (instancetype)initWithDevice:(id<MTLDevice>)device withParam:(Init*)initParam withLogger:(XmlLogger&)logger{
    self = [super init];
    if (self) {

        // パラメータ取得
        struct ParamForField fieldParam = initParam.paramForField;
        struct ParamForComputing compParam = initParam.paramForComputing;
        std::vector<struct ParamForParticle> particleParam = initParam.paramForParticle;
        std::vector<struct BoundaryConditionForParticle> pBCs = initParam.particleBoundaries;

        // 並列計算パラメータ
        _integrationChunkSize = compParam.integrationChunkSize;
        _threadGroupSize = compParam.threadGroupSize;
        
        // コマンドキューの作成
        _device = device;
        _commandQueue = [device newCommandQueue];

        // コンピュートパイプラインの設定
        NSError *error = nil;
        id<MTLLibrary> library = [device newLibraryWithSource:kMetalShaderSource options:nil error:&error];
        if (!library) {
            NSLog(@"Failed to create Metal library: %@", error);
            return nil;
        }

        // カーネルの取得
        id<MTLFunction> integrateNumDensFunction = [library newFunctionWithName:@"integrateNumDens"];
        _integrateNumDensPipeline = [device newComputePipelineStateWithFunction:integrateNumDensFunction error:&error];
        id<MTLFunction> integrateMeanVelXFunction = [library newFunctionWithName:@"integrateMeanVelX"];
        _integrateMeanVelXPipeline = [device newComputePipelineStateWithFunction:integrateMeanVelXFunction error:&error];
        id<MTLFunction> integrateMeanVelYFunction = [library newFunctionWithName:@"integrateMeanVelY"];
        _integrateMeanVelYPipeline = [device newComputePipelineStateWithFunction:integrateMeanVelYFunction error:&error];
        id<MTLFunction> integrateMeanVelZFunction = [library newFunctionWithName:@"integrateMeanVelZ"];
        _integrateMeanVelZPipeline = [device newComputePipelineStateWithFunction:integrateMeanVelZFunction error:&error];
        id<MTLFunction> integratePressureXXFunction = [library newFunctionWithName:@"integratePressureXX"];
        _integratePressureXXPipeline = [device newComputePipelineStateWithFunction:integratePressureXXFunction error:&error];
        id<MTLFunction> integratePressureXYFunction = [library newFunctionWithName:@"integratePressureXY"];
        _integratePressureXYPipeline = [device newComputePipelineStateWithFunction:integratePressureXYFunction error:&error];
        id<MTLFunction> integratePressureXZFunction = [library newFunctionWithName:@"integratePressureXZ"];
        _integratePressureXZPipeline = [device newComputePipelineStateWithFunction:integratePressureXZFunction error:&error];
        id<MTLFunction> integratePressureYYFunction = [library newFunctionWithName:@"integratePressureYY"];
        _integratePressureYYPipeline = [device newComputePipelineStateWithFunction:integratePressureYYFunction error:&error];
        id<MTLFunction> integratePressureYZFunction = [library newFunctionWithName:@"integratePressureYZ"];
        _integratePressureYZPipeline = [device newComputePipelineStateWithFunction:integratePressureYZFunction error:&error];
        id<MTLFunction> integratePressureZZFunction = [library newFunctionWithName:@"integratePressureZZ"];
        _integratePressureZZPipeline = [device newComputePipelineStateWithFunction:integratePressureZZFunction error:&error];

        // 定数パラメータ格納バッファ
        size_t buffSize = sizeof(integrationParams);
        _integrationParamsBuffer = [device newBufferWithLength:buffSize options:MTLResourceStorageModeShared];
        integrationParams* prm = (integrationParams*)[_integrationParamsBuffer contents];
        prm->ngx  = fieldParam.ngx;
        prm->ngy  = fieldParam.ngy;
        prm->ngb  = fieldParam.ngb;
        prm->tgs  = _threadGroupSize;
        prm->ics  = _integrationChunkSize;
        prm->pNum = 0; // ptclクラスごとに都度計算
        prm->ppt  = 0; // ptclクラスごとに都度計算

        // 積分計算用バッファ
        int nx = fieldParam.ngx + 2*fieldParam.ngb;
        int ny = fieldParam.ngy + 2*fieldParam.ngb;
        int ng = (nx+1)*(ny+1);
        buffSize = sizeof(float)*ng*_integrationChunkSize;
        _integrateTemporaryBuffer = [device newBufferWithLength:buffSize options:MTLResourceStorageModePrivate];
        buffSize = sizeof(float)*ng*(_integrationChunkSize/_threadGroupSize);
        _integratePartialBuffer = [device newBufferWithLength:buffSize options:MTLResourceStorageModeShared];

        // デバッグ出力用バッファ
        buffSize = sizeof(float)*compParam.pNumMax;
        _printBuffer = [device newBufferWithLength:buffSize options:MTLResourceStorageModeShared];

        // malloc
        _n = (float *)malloc(sizeof(float)*ng);
        _ux = (float *)malloc(sizeof(float)*ng);
        _uy = (float *)malloc(sizeof(float)*ng);
        _uz = (float *)malloc(sizeof(float)*ng);
        _Pxx = (float *)malloc(sizeof(float)*ng);
        _Pxy = (float *)malloc(sizeof(float)*ng);
        _Pxz = (float *)malloc(sizeof(float)*ng);
        _Pyy = (float *)malloc(sizeof(float)*ng);
        _Pyz = (float *)malloc(sizeof(float)*ng);
        _Pzz = (float *)malloc(sizeof(float)*ng);
    }
    return self;
}

- (void)integrateMoments:(Particle*)ptcl withEMField:(EMField*)fld withLogger:(XmlLogger&)logger{
    // initialize
    integrationParams* prm = (integrationParams*)[_integrationParamsBuffer contents];
    int ng = (prm->ngx + 2*prm->ngb + 1)*(prm->ngy + 2*prm->ngb + 1);
    for (int i = 0; i < ng; i++){
        _n[i] = 0.0;
        _ux[i] = 0.0;
        _uy[i] = 0.0;
        _uz[i] = 0.0;
        _Pxx[i] = 0.0;
        _Pxy[i] = 0.0;
        _Pxz[i] = 0.0;
        _Pyy[i] = 0.0;
        _Pyz[i] = 0.0;
        _Pzz[i] = 0.0;
    }
    
    // 粒子データは Particle クラスのバッファを使用
    id<MTLBuffer> particleBuffer = [ptcl particleBuffer];
    id<MTLBuffer> paramsBuffer = [ptcl paramsBuffer];
    SimulationParams* prmPtcl = (SimulationParams*)[paramsBuffer contents];
    
    // 定数パラメータ更新
    prm->pNum = prmPtcl->pNum;
    uint pNumPerThread = (prmPtcl->pNum + _integrationChunkSize - 1) / _integrationChunkSize;
    prm->ppt  = pNumPerThread;

    // グリッドとスレッドグループのサイズ設定
    uint threadGroupNum = _integrationChunkSize/_threadGroupSize;
    MTLSize gridSizeMetalStyle = MTLSizeMake(threadGroupNum, 1, 1);
    MTLSize threadGroupSize = MTLSizeMake(_threadGroupSize, 1, 1);
    
    // コマンドバッファとエンコーダ
    id<MTLCommandBuffer> commandBuffer;
    id<MTLComputeCommandEncoder> computeEncoder;
    
    // 部分和取得
    float* partialSums;

    // 積分実行: n
    commandBuffer = [_commandQueue commandBuffer];
    computeEncoder = [commandBuffer computeCommandEncoder];
    [computeEncoder setComputePipelineState:_integrateNumDensPipeline];
    [computeEncoder setBuffer:particleBuffer            offset:0 atIndex:0];
    [computeEncoder setBuffer:_integrationParamsBuffer  offset:0 atIndex:1];
    [computeEncoder setBuffer:_integrateTemporaryBuffer offset:0 atIndex:2];
    [computeEncoder setBuffer:_integratePartialBuffer   offset:0 atIndex:3];
    [computeEncoder setBuffer:_printBuffer              offset:0 atIndex:4];
    [computeEncoder dispatchThreadgroups:gridSizeMetalStyle threadsPerThreadgroup:threadGroupSize];
    [computeEncoder endEncoding];
    [commandBuffer commit];
    [commandBuffer waitUntilCompleted];
    partialSums = (float*)_integratePartialBuffer.contents;
    for (int j = 0; j < threadGroupNum; j++){
        for (int i = 0; i < ng; i++){
            _n[i] += partialSums[i + j*ng];
        }
    }
    
    // 積分実行: ux
    commandBuffer = [_commandQueue commandBuffer];
    computeEncoder = [commandBuffer computeCommandEncoder];
    [computeEncoder setComputePipelineState:_integrateMeanVelXPipeline];
    [computeEncoder setBuffer:particleBuffer            offset:0 atIndex:0];
    [computeEncoder setBuffer:_integrationParamsBuffer  offset:0 atIndex:1];
    [computeEncoder setBuffer:_integrateTemporaryBuffer offset:0 atIndex:2];
    [computeEncoder setBuffer:_integratePartialBuffer   offset:0 atIndex:3];
    [computeEncoder setBuffer:_printBuffer              offset:0 atIndex:4];
    [computeEncoder dispatchThreadgroups:gridSizeMetalStyle threadsPerThreadgroup:threadGroupSize];
    [computeEncoder endEncoding];
    [commandBuffer commit];
    [commandBuffer waitUntilCompleted];
    partialSums = (float*)_integratePartialBuffer.contents;
    for (int j = 0; j < threadGroupNum; j++){
        for (int i = 0; i < ng; i++){
            _ux[i] += partialSums[i + j*ng];
        }
    }

    // 積分実行: uy
    commandBuffer = [_commandQueue commandBuffer];
    computeEncoder = [commandBuffer computeCommandEncoder];
    [computeEncoder setComputePipelineState:_integrateMeanVelYPipeline];
    [computeEncoder setBuffer:particleBuffer            offset:0 atIndex:0];
    [computeEncoder setBuffer:_integrationParamsBuffer  offset:0 atIndex:1];
    [computeEncoder setBuffer:_integrateTemporaryBuffer offset:0 atIndex:2];
    [computeEncoder setBuffer:_integratePartialBuffer   offset:0 atIndex:3];
    [computeEncoder setBuffer:_printBuffer              offset:0 atIndex:4];
    [computeEncoder dispatchThreadgroups:gridSizeMetalStyle threadsPerThreadgroup:threadGroupSize];
    [computeEncoder endEncoding];
    [commandBuffer commit];
    [commandBuffer waitUntilCompleted];
    partialSums = (float*)_integratePartialBuffer.contents;
    for (int j = 0; j < threadGroupNum; j++){
        for (int i = 0; i < ng; i++){
            _uy[i] += partialSums[i + j*ng];
        }
    }

    // 積分実行: uz
    commandBuffer = [_commandQueue commandBuffer];
    computeEncoder = [commandBuffer computeCommandEncoder];
    [computeEncoder setComputePipelineState:_integrateMeanVelZPipeline];
    [computeEncoder setBuffer:particleBuffer            offset:0 atIndex:0];
    [computeEncoder setBuffer:_integrationParamsBuffer  offset:0 atIndex:1];
    [computeEncoder setBuffer:_integrateTemporaryBuffer offset:0 atIndex:2];
    [computeEncoder setBuffer:_integratePartialBuffer   offset:0 atIndex:3];
    [computeEncoder setBuffer:_printBuffer              offset:0 atIndex:4];
    [computeEncoder dispatchThreadgroups:gridSizeMetalStyle threadsPerThreadgroup:threadGroupSize];
    [computeEncoder endEncoding];
    [commandBuffer commit];
    [commandBuffer waitUntilCompleted];
    partialSums = (float*)_integratePartialBuffer.contents;
    for (int j = 0; j < threadGroupNum; j++){
        for (int i = 0; i < ng; i++){
            _uz[i] += partialSums[i + j*ng];
        }
    }

    // 積分実行: Pxx
    commandBuffer = [_commandQueue commandBuffer];
    computeEncoder = [commandBuffer computeCommandEncoder];
    [computeEncoder setComputePipelineState:_integratePressureXXPipeline];
    [computeEncoder setBuffer:particleBuffer            offset:0 atIndex:0];
    [computeEncoder setBuffer:_integrationParamsBuffer  offset:0 atIndex:1];
    [computeEncoder setBuffer:_integrateTemporaryBuffer offset:0 atIndex:2];
    [computeEncoder setBuffer:_integratePartialBuffer   offset:0 atIndex:3];
    [computeEncoder setBuffer:_printBuffer              offset:0 atIndex:4];
    [computeEncoder dispatchThreadgroups:gridSizeMetalStyle threadsPerThreadgroup:threadGroupSize];
    [computeEncoder endEncoding];
    [commandBuffer commit];
    [commandBuffer waitUntilCompleted];
    partialSums = (float*)_integratePartialBuffer.contents;
    for (int j = 0; j < threadGroupNum; j++){
        for (int i = 0; i < ng; i++){
            _Pxx[i] += partialSums[i + j*ng];
        }
    }

    // 積分実行: Pxy
    commandBuffer = [_commandQueue commandBuffer];
    computeEncoder = [commandBuffer computeCommandEncoder];
    [computeEncoder setComputePipelineState:_integratePressureXYPipeline];
    [computeEncoder setBuffer:particleBuffer            offset:0 atIndex:0];
    [computeEncoder setBuffer:_integrationParamsBuffer  offset:0 atIndex:1];
    [computeEncoder setBuffer:_integrateTemporaryBuffer offset:0 atIndex:2];
    [computeEncoder setBuffer:_integratePartialBuffer   offset:0 atIndex:3];
    [computeEncoder setBuffer:_printBuffer              offset:0 atIndex:4];
    [computeEncoder dispatchThreadgroups:gridSizeMetalStyle threadsPerThreadgroup:threadGroupSize];
    [computeEncoder endEncoding];
    [commandBuffer commit];
    [commandBuffer waitUntilCompleted];
    partialSums = (float*)_integratePartialBuffer.contents;
    for (int j = 0; j < threadGroupNum; j++){
        for (int i = 0; i < ng; i++){
            _Pxy[i] += partialSums[i + j*ng];
        }
    }

    // 積分実行: Pxz
    commandBuffer = [_commandQueue commandBuffer];
    computeEncoder = [commandBuffer computeCommandEncoder];
    [computeEncoder setComputePipelineState:_integratePressureXZPipeline];
    [computeEncoder setBuffer:particleBuffer            offset:0 atIndex:0];
    [computeEncoder setBuffer:_integrationParamsBuffer  offset:0 atIndex:1];
    [computeEncoder setBuffer:_integrateTemporaryBuffer offset:0 atIndex:2];
    [computeEncoder setBuffer:_integratePartialBuffer   offset:0 atIndex:3];
    [computeEncoder setBuffer:_printBuffer              offset:0 atIndex:4];
    [computeEncoder dispatchThreadgroups:gridSizeMetalStyle threadsPerThreadgroup:threadGroupSize];
    [computeEncoder endEncoding];
    [commandBuffer commit];
    [commandBuffer waitUntilCompleted];
    partialSums = (float*)_integratePartialBuffer.contents;
    for (int j = 0; j < threadGroupNum; j++){
        for (int i = 0; i < ng; i++){
            _Pxz[i] += partialSums[i + j*ng];
        }
    }

    // 積分実行: Pyy
    commandBuffer = [_commandQueue commandBuffer];
    computeEncoder = [commandBuffer computeCommandEncoder];
    [computeEncoder setComputePipelineState:_integratePressureYYPipeline];
    [computeEncoder setBuffer:particleBuffer            offset:0 atIndex:0];
    [computeEncoder setBuffer:_integrationParamsBuffer  offset:0 atIndex:1];
    [computeEncoder setBuffer:_integrateTemporaryBuffer offset:0 atIndex:2];
    [computeEncoder setBuffer:_integratePartialBuffer   offset:0 atIndex:3];
    [computeEncoder setBuffer:_printBuffer              offset:0 atIndex:4];
    [computeEncoder dispatchThreadgroups:gridSizeMetalStyle threadsPerThreadgroup:threadGroupSize];
    [computeEncoder endEncoding];
    [commandBuffer commit];
    [commandBuffer waitUntilCompleted];
    partialSums = (float*)_integratePartialBuffer.contents;
    for (int j = 0; j < threadGroupNum; j++){
        for (int i = 0; i < ng; i++){
            _Pyy[i] += partialSums[i + j*ng];
        }
    }

    // 積分実行: Pyz
    commandBuffer = [_commandQueue commandBuffer];
    computeEncoder = [commandBuffer computeCommandEncoder];
    [computeEncoder setComputePipelineState:_integratePressureYZPipeline];
    [computeEncoder setBuffer:particleBuffer            offset:0 atIndex:0];
    [computeEncoder setBuffer:_integrationParamsBuffer  offset:0 atIndex:1];
    [computeEncoder setBuffer:_integrateTemporaryBuffer offset:0 atIndex:2];
    [computeEncoder setBuffer:_integratePartialBuffer   offset:0 atIndex:3];
    [computeEncoder setBuffer:_printBuffer              offset:0 atIndex:4];
    [computeEncoder dispatchThreadgroups:gridSizeMetalStyle threadsPerThreadgroup:threadGroupSize];
    [computeEncoder endEncoding];
    [commandBuffer commit];
    [commandBuffer waitUntilCompleted];
    partialSums = (float*)_integratePartialBuffer.contents;
    for (int j = 0; j < threadGroupNum; j++){
        for (int i = 0; i < ng; i++){
            _Pyz[i] += partialSums[i + j*ng];
        }
    }

    // 積分実行: Pzz
    commandBuffer = [_commandQueue commandBuffer];
    computeEncoder = [commandBuffer computeCommandEncoder];
    [computeEncoder setComputePipelineState:_integratePressureZZPipeline];
    [computeEncoder setBuffer:particleBuffer            offset:0 atIndex:0];
    [computeEncoder setBuffer:_integrationParamsBuffer  offset:0 atIndex:1];
    [computeEncoder setBuffer:_integrateTemporaryBuffer offset:0 atIndex:2];
    [computeEncoder setBuffer:_integratePartialBuffer   offset:0 atIndex:3];
    [computeEncoder setBuffer:_printBuffer              offset:0 atIndex:4];
    [computeEncoder dispatchThreadgroups:gridSizeMetalStyle threadsPerThreadgroup:threadGroupSize];
    [computeEncoder endEncoding];
    [commandBuffer commit];
    [commandBuffer waitUntilCompleted];
    partialSums = (float*)_integratePartialBuffer.contents;
    for (int j = 0; j < threadGroupNum; j++){
        for (int i = 0; i < ng; i++){
            _Pzz[i] += partialSums[i + j*ng];
        }
    }

    // 積分量->モーメント量
    // ptcl -> 1/cm3
    float constN = ptcl.w/(fld.dx*fld.dy);
    for (int i = 0; i < ng; i++){
        _n[i] *= constN;
    }
    // cm/s*ptcl -> cm/s
    for (int i = 0; i < ng; i++){
        if(_n[i] > 1e-20){
            _ux[i] /= _n[i]/constN;
            _uy[i] /= _n[i]/constN;
            _uz[i] /= _n[i]/constN;
        }else{
            _ux[i] = 0.0;
            _uy[i] = 0.0;
            _uz[i] = 0.0;
        }
    }
    // cm2/s2*ptcl -> g*cm2/s2*1/cm3 = 0.1 Pa
    // Cov[X,Y] = E[XY] - E[X]E[Y]
    for (int i = 0; i < ng; i++){
        if(_n[i] > 1e-20){
            _Pxx[i] = ptcl.m*(_Pxx[i]/(_n[i]/constN) - _ux[i]*_ux[i])*_n[i];
            _Pxy[i] = ptcl.m*(_Pxy[i]/(_n[i]/constN) - _ux[i]*_uy[i])*_n[i];
            _Pxz[i] = ptcl.m*(_Pxz[i]/(_n[i]/constN) - _ux[i]*_uz[i])*_n[i];
            _Pyy[i] = ptcl.m*(_Pyy[i]/(_n[i]/constN) - _uy[i]*_uy[i])*_n[i];
            _Pyz[i] = ptcl.m*(_Pyz[i]/(_n[i]/constN) - _uy[i]*_uz[i])*_n[i];
            _Pzz[i] = ptcl.m*(_Pzz[i]/(_n[i]/constN) - _uz[i]*_uz[i])*_n[i];
        }else{
            _Pxx[i] = 0.0;
            _Pxy[i] = 0.0;
            _Pxz[i] = 0.0;
            _Pyy[i] = 0.0;
            _Pyz[i] = 0.0;
            _Pzz[i] = 0.0;
        }
    }

};

static void writeField(FILE* fp, const char* name, int type_id, float* array, int arrSize, float scale) {
    // 変数名（固定長32文字）
    char name_buf[32] = {};
    strncpy(name_buf, name, sizeof(name_buf) - 1); // 末尾NULL保護
    fwrite(name_buf, sizeof(char), 32, fp);
    
    // タイプID（int）
    // 0: node-center value,    1: left-shifted value
    // 2: bottom-shifted value, 3: left&bottom-shifted value
    // 4: potential
    fwrite(&type_id, sizeof(int), 1, fp);
    
    // 配列本体をスケーリングして出力
    float* output = (float *)malloc(sizeof(float)*arrSize);
    for (int i = 0; i < arrSize; ++i) {
        output[i] = array[i] * scale;
    }
    fwrite(output, sizeof(float), arrSize, fp);
    // ここで必ず free する
    free(output);
}

// 積分計算用バッファへのアクセサ
- (id<MTLBuffer>)integrateTemporaryBuffer { return _integrateTemporaryBuffer; }
- (id<MTLBuffer>)integratePartialBuffer { return _integratePartialBuffer; }
- (id<MTLBuffer>)printBuffer { return _printBuffer; }
// モーメント配列へのアクセサ
- (float*)n { return _n; }
- (float*)ux { return _ux; }
- (float*)uy { return _uy; }
- (float*)uz { return _uz; }
- (float*)Pxx { return _Pxx; }
- (float*)Pxy { return _Pxy; }
- (float*)Pxz { return _Pxz; }
- (float*)Pyy { return _Pyy; }
- (float*)Pyz { return _Pyz; }
- (float*)Pzz { return _Pzz; }

@end