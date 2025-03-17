#import "EMField.h"
#import <Metal/Metal.h>
#include <amgcl/make_solver.hpp>
#include <amgcl/solver/bicgstab.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/adapter/crs_tuple.hpp>

// amgcl 
typedef amgcl::backend::builtin<float> Backend;
typedef amgcl::make_solver<
    amgcl::amg<
        Backend,
        amgcl::coarsening::smoothed_aggregation,
        amgcl::relaxation::spai0
    >,
    amgcl::solver::bicgstab<Backend>
> Solver;

@interface EMField () {
    id<MTLDevice> _device;
    float* _rho;                 // 電荷密度
    float* _phi;                 // ポテンシャル
    id<MTLBuffer> _ExBuffer;     // 電場 x成分
    id<MTLBuffer> _EyBuffer;     // 電場 y成分
    id<MTLBuffer> _EzBuffer;     // 電場 z成分
    id<MTLBuffer> _BxBuffer;     // 磁場 x成分
    id<MTLBuffer> _ByBuffer;     // 磁場 y成分
    id<MTLBuffer> _BzBuffer;     // 磁場 z成分
    
    int _ngx;                    // x方向グリッド数
    int _ngy;                    // y方向グリッド数
    int _ngb;                    // バッファ領域グリッド数
    float _dx;                  // x方向グリッド幅
    float _dy;                  // y方向グリッド幅
    
    int _weightOrder;            // weighting のオーダー
}
@end

// プライベートインスタンス変数
@implementation EMField {
    // number of grid for buffer
    int _ngbrho[2];
    int _ngbphi[2];
    int _ngbEx[2];
    int _ngbEy[2];
    int _ngbEz[2];
    int _ngbBx[2];
    int _ngbBy[2];
    int _ngbBz[2];
    // BoundaryCoditions
    NSString* _BC_Xmin;
    NSString* _BC_Xmax;
    NSString* _BC_Ymin;
    NSString* _BC_Ymax;
    std::vector<float> _phi_Xmin;
    std::vector<float> _phi_Xmax;
    std::vector<float> _phi_Ymin;
    std::vector<float> _phi_Ymax;
    std::vector<float> _dphidx_Xmin;
    std::vector<float> _dphidx_Xmax;
    std::vector<float> _dphidy_Ymin;
    std::vector<float> _dphidy_Ymax;
    int _nkx;
    int _nky;
    int _nk;
    float _Cd;
    float _Cx;
    float _Cy;
    // smart pointer
    std::unique_ptr<Solver> _solver;
}

- (instancetype)initWithDevice:(id<MTLDevice>)device withParam:(Init*)initParam {
    self = [super init];
    if (self) {
        _device = device;
        
        // フィールドパラメータの取得
        struct ParamForField fieldParam = [initParam getParamForField];
        _ngx = fieldParam.ngx;
        _ngy = fieldParam.ngy;
        _dx = fieldParam.dx;
        _dy = fieldParam.dy;
        _weightOrder = fieldParam.weightOrder;

        // 小行列の配列サイズ（A[1][1] が解析領域で、それ以外は領域境界）
        if (_weightOrder == 5){
            // general
            _ngb = 2;
            // rho,phi,Ez
            _ngbrho[0] = 2;
            _ngbphi[0] = 2;
            _ngbEz[0] = 2;
            _ngbrho[1] = 2;
            _ngbphi[1] = 2;
            _ngbEz[1] = 2;
            // Ex, By
            _ngbEx[0] = 3;
            _ngbBy[0] = 3;
            _ngbEx[1] = 2;
            _ngbBy[1] = 2;
            // Ey, Bx
            _ngbEy[0] = 2;
            _ngbBx[0] = 2;
            _ngbEy[1] = 3;
            _ngbBx[1] = 3;
            // Bz
            _ngbBz[0] = 3;
            _ngbBz[0] = 3;
        } else {
            _ngb = 0;
        }

        // バッファサイズの計算
        NSUInteger gridSize = (_ngx + 2*_ngb) * (_ngy + 2*_ngb) ;
        NSUInteger bufferSize = sizeof(float) * gridSize;
        
        // Metal バッファの初期化
        _ExBuffer = [device newBufferWithLength:bufferSize options:MTLResourceStorageModeShared];
        _EyBuffer = [device newBufferWithLength:bufferSize options:MTLResourceStorageModeShared];
        _EzBuffer = [device newBufferWithLength:bufferSize options:MTLResourceStorageModeShared];
        _BxBuffer = [device newBufferWithLength:bufferSize options:MTLResourceStorageModeShared];
        _ByBuffer = [device newBufferWithLength:bufferSize options:MTLResourceStorageModeShared];
        _BzBuffer = [device newBufferWithLength:bufferSize options:MTLResourceStorageModeShared];
        
        // バッファの初期化
        float* Ex = (float*)_ExBuffer.contents;
        float* Ey = (float*)_EyBuffer.contents;
        float* Ez = (float*)_EzBuffer.contents;
        float* Bx = (float*)_BxBuffer.contents;
        float* By = (float*)_ByBuffer.contents;
        float* Bz = (float*)_BzBuffer.contents;
        
        // malloc
        _rho = (float *)malloc(bufferSize);
        _phi = (float *)malloc(bufferSize);
        
        // initialize(Uniform)
        for (int i = 0; i < gridSize; i++) {
            _rho[i] = 0.0f;
            _phi[i] = 0.0f;
            Ex[i] = fieldParam.ampE[0];
            Ey[i] = fieldParam.ampE[1];
            Ez[i] = fieldParam.ampE[2];
            Bx[i] = fieldParam.ampB[0];
            By[i] = fieldParam.ampB[1];
            Bz[i] = fieldParam.ampB[2];
        }

        // construct Poisson solver
        float val_Xmin, val_Xmax, val_Ymin, val_Ymax;
        NSArray* fieldBCs = [initParam getFieldBoundaries];
        for (int pos = 0; pos < fieldBCs.count; pos++){
            NSValue* value = fieldBCs[pos];
            struct BoundaryConditionForField BC;
            [value getValue:&BC];
            if ([BC.position isEqualToString:@"Xmin"]){
                _BC_Xmin = BC.type;
                if ([BC.type isEqualToString:@"Dirichlet"]){
                    _phi_Xmin.assign(_ngy, BC.val);
                }else if ([BC.type isEqualToString:@"Neumann"]){
                    _dphidx_Xmin.assign(_ngy, BC.val);
                }
            }else if ([BC.position isEqualToString:@"Xmax"]){
                _BC_Xmax = BC.type;
                if ([BC.type isEqualToString:@"Dirichlet"]){
                    _phi_Xmax.assign(_ngy, BC.val);
                }else if ([BC.type isEqualToString:@"Neumann"]){
                    _dphidx_Xmax.assign(_ngy, BC.val);
                }
            }else if ([BC.position isEqualToString:@"Ymin"]){
                _BC_Ymin = BC.type;
                if ([BC.type isEqualToString:@"Dirichlet"]){
                    _phi_Ymin.assign(_ngx, BC.val);
                }else if ([BC.type isEqualToString:@"Neumann"]){
                    _dphidy_Ymin.assign(_ngx, BC.val);
                }
            }else if ([BC.position isEqualToString:@"Ymax"]){
                _BC_Ymax = BC.type;
                if ([BC.type isEqualToString:@"Dirichlet"]){
                    _phi_Ymax.assign(_ngx, BC.val);
                }else if ([BC.type isEqualToString:@"Neumann"]){
                    _dphidy_Ymax.assign(_ngx, BC.val);
                }
            }
        }
        // 定数係数
        _Cd = -2*(1.0/_dx/_dx + 1.0/_dy/_dy);
        _Cx = 1.0/_dx/_dx;
        _Cy = 1.0/_dy/_dy;
        
        // CSR 行列
        std::vector<int> ptr;
        std::vector<int> col;
        std::vector<float> val;
        std::vector<float> rhs;

        // 解くべき行列サイズ
        _nkx = _ngx-2;
        if ([_BC_Xmin isEqualToString:@"Neumann"]){
            _nkx++;
        }
        if ([_BC_Xmax isEqualToString:@"Neumann"]){
            _nkx++;
        }
        if ([_BC_Xmin isEqualToString:@"periodic"]){
            _nkx++;
        }
        _nky = _ngy-2;
        if ([_BC_Ymin isEqualToString:@"Neumann"]){
            _nky++;
        }
        if ([_BC_Ymax isEqualToString:@"Neumann"]){
            _nky++;
        }
        if ([_BC_Ymin isEqualToString:@"periodic"]){
            _nky++;
        }
        _nk = _nkx*_nky;

        // construct CSR matrix
        ptr.push_back(0);
        for (int j = 0; j < _nky; ++j) {
            for (int i = 0; i < _nkx; ++i) {
                int k = i + j*_nkx;
                int row_entries = 0;

                // Neumann(前進/後進 差分)
                if (i == 0 && [_BC_Xmin isEqualToString:@"Neumann"]){
                    col.push_back(k);
                    val.push_back(-1.0/_dx);
                    row_entries++;
                    col.push_back(k+1);
                    val.push_back(1.0/_dx);
                    row_entries++;
                    rhs.push_back(_dphidx_Xmin[j]);
                    // ptr を更新してループを抜ける
                    ptr.push_back(ptr.back() + row_entries);
                    continue;
                }else if(i == _nkx-1 && [_BC_Xmax isEqualToString:@"Neumann"]){
                    col.push_back(k);
                    val.push_back(1.0/_dx);
                    row_entries++;
                    col.push_back(k-1);
                    val.push_back(-1.0/_dx);
                    row_entries++;
                    rhs.push_back(_dphidx_Xmax[j]);
                    // ptr を更新してループを抜ける
                    ptr.push_back(ptr.back() + row_entries);
                    continue;
                }else if (j == 0 && [_BC_Ymin isEqualToString:@"Neumann"]){
                    col.push_back(k);
                    val.push_back(-1.0/_dy);
                    row_entries++;
                    col.push_back(k+_nkx);
                    val.push_back(1.0/_dy);
                    row_entries++;
                    rhs.push_back(_dphidy_Ymin[i]);
                    // ptr を更新してループを抜ける
                    ptr.push_back(ptr.back() + row_entries);
                    continue;
                }else if(j == _nky-1 && [_BC_Ymax isEqualToString:@"Neumann"]){
                    col.push_back(k);
                    val.push_back(1.0/_dx);
                    row_entries++;
                    col.push_back(k-_nkx);
                    val.push_back(-1.0/_dx);
                    row_entries++;
                    rhs.push_back(_dphidy_Ymax[i]);
                    // ptr を更新してループを抜ける
                    ptr.push_back(ptr.back() + row_entries);
                    continue;
                }

                // diagonal
                col.push_back(k);
                val.push_back(_Cd);
                row_entries++;
                rhs.push_back(_rho[k]);

                // imin側
                if (i == 0){
                    if ([_BC_Xmin isEqualToString:@"Dirichlet"]){
                        // ディリクレ境界に接してるなら b に足す
                        rhs[k] += -_Cx*_phi_Xmin[j];
                    }else if ([_BC_Xmin isEqualToString:@"periodic"]){
                        // 要素の追加
                        col.push_back(k+(_nkx-1));
                        val.push_back(_Cx);
                        row_entries++;
                    }
                }else{
                    // 領域内なら val を追加
                    col.push_back(k-1);
                    val.push_back(_Cx);
                    row_entries++;
                }

                // imax側
                if (i == _nkx-1){
                    if ([_BC_Xmax isEqualToString:@"Dirichlet"]){
                        rhs[k] += -_Cx*_phi_Xmax[j];
                    }else if ([_BC_Xmax isEqualToString:@"periodic"]){
                        // 要素の追加
                        col.push_back(k-(_nkx-1));
                        val.push_back(_Cx);
                        row_entries++;
                    }
                }else{
                    col.push_back(k+1);
                    val.push_back(_Cx);
                    row_entries++;
                }

                // jmin側
                if (j == 0){
                    if ([_BC_Ymin isEqualToString:@"Dirichlet"]){
                        rhs[k] += -_Cy*_phi_Ymin[j];
                    }else if ([_BC_Ymin isEqualToString:@"periodic"]){
                        // 要素の追加
                        col.push_back(k+(_nky-1));
                        val.push_back(_Cy);
                        row_entries++;
                    }
                }else{
                    col.push_back(k-_nkx);
                    val.push_back(_Cy);
                    row_entries++;
                }

                // jmax側
                if (j == _nky-1){
                    if ([_BC_Ymax isEqualToString:@"Dirichlet"]){
                        rhs[k] += -_Cy*_phi_Ymax[i];
                    }else if ([_BC_Ymax isEqualToString:@"periodic"]){
                        // 要素の追加
                        col.push_back(k-(_nky-1)*_nkx);
                        val.push_back(_Cy);
                        row_entries++;
                    }
                }else{
                    // 要素の追加
                    col.push_back(k+_nkx);
                    val.push_back(_Cy);
                    row_entries++;
                }

                // CSR ポインタ更新
                ptr.push_back(ptr.back() + row_entries);
            }
        }
        // --- amgcl の設定 ---
        Solver::params prm;
        prm.solver.maxiter = fieldParam.maxiter;
        prm.solver.tol = fieldParam.tolerance;
        // _solver を生成してインスタンス変数として確保
        _solver = std::unique_ptr<Solver>(
            new Solver( std::tie(_nk, ptr, col, val), prm )
        );
    }
    return self;
}

- (void)solvePoisson{

}

- (void)resetChargeDensity{
    // 電荷密度の初期化
    int ng = (_ngx + 2*_ngb) * (_ngy + 2*_ngb);
    for (int i = 0; i < ng; i++){
        _rho[i] = 0.0f;
    }
}

- (void)outputField:(int)cycle{
    // サイクル番号を用いてファイル名を作成 (例: field_0001.bin)
    NSString *fileName = [NSString stringWithFormat:@"bin/field_%08d.bin", cycle];
    const char *filePath = [fileName UTF8String];
    
    FILE *fp = fopen(filePath, "wb");
    if (!fp) {
        NSLog(@"Error: Unable to open file %s for writing", filePath);
        return;
    }
    
    // ヘッダ情報としてメッシュ情報を出力
    fwrite(&_ngx, sizeof(int), 1, fp);
    fwrite(&_ngy, sizeof(int), 1, fp);
    fwrite(&_ngb, sizeof(int), 1, fp);
    fwrite(&_dx, sizeof(float), 1, fp);
    fwrite(&_dy, sizeof(float), 1, fp);
    
    // 全グリッドサイズ (境界含む) の計算
    int ng = (_ngx + 2 * _ngb) * (_ngy + 2 * _ngb);
    
    // 順にフィールドデータを書き出す: _rho, _phi, Ex, Ey, Bz
    fwrite(_rho, sizeof(float), ng, fp);
    fwrite(_phi, sizeof(float), ng, fp);
    
    float *Ex = (float *)_ExBuffer.contents;
    float *Ey = (float *)_EyBuffer.contents;
    float *Bz = (float *)_BzBuffer.contents;
    
    fwrite(Ex, sizeof(float), ng, fp);
    fwrite(Ey, sizeof(float), ng, fp);
    fwrite(Bz, sizeof(float), ng, fp);
    
    fclose(fp);
    NSLog(@"Field data successfully written to %s", filePath);
}

// 電荷密度へのアクセサ
- (float*)rho { return _rho; }

// 電場バッファへのアクセサ
- (id<MTLBuffer>)ExBuffer { return _ExBuffer; }
- (id<MTLBuffer>)EyBuffer { return _EyBuffer; }
- (id<MTLBuffer>)EzBuffer { return _EzBuffer; }

// 磁場バッファへのアクセサ
- (id<MTLBuffer>)BxBuffer { return _BxBuffer; }
- (id<MTLBuffer>)ByBuffer { return _ByBuffer; }
- (id<MTLBuffer>)BzBuffer { return _BzBuffer; }

// グリッド情報へのアクセサ
- (int)ngx { return _ngx; }
- (int)ngy { return _ngy; }
- (int)ngb { return _ngb; }
- (double)dx { return _dx; }
- (double)dy { return _dy; }

@end

#import "EMField.h"