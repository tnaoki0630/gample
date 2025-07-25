#import <Metal/Metal.h>
#import <Accelerate/Accelerate.h>
#import "EMField.h"
#import "Constant.h"
#import "XmlLogger.h"
#include <amgcl/make_solver.hpp>
#include <amgcl/backend/builtin.hpp>
#include <amgcl/solver/bicgstab.hpp>
#include <amgcl/solver/cg.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/relaxation/chebyshev.hpp>
#include <amgcl/adapter/crs_tuple.hpp>

# define PI 3.141592653589793

// amgcl 
typedef amgcl::backend::builtin<double> Backend;
typedef amgcl::make_solver<
    amgcl::amg<
        Backend,
        amgcl::coarsening::smoothed_aggregation,
        amgcl::relaxation::spai0
    >,
    amgcl::solver::cg<Backend>
>  Solver;

@interface EMField () {
    id<MTLDevice> _device;
    id<MTLBuffer> _rhoBuffer;    // 電荷密度
    id<MTLBuffer> _atomicRhoBuffer;    // 電荷密度
    double* _phi;                 // 電位
    id<MTLBuffer> _ExBuffer;     // 電場 x成分
    id<MTLBuffer> _EyBuffer;     // 電場 y成分
    id<MTLBuffer> _EzBuffer;     // 電場 z成分
    id<MTLBuffer> _BxBuffer;     // 磁場 x成分
    id<MTLBuffer> _ByBuffer;     // 磁場 y成分
    id<MTLBuffer> _BzBuffer;     // 磁場 z成分
    
    int _ngx;                    // x方向グリッド数
    int _ngy;                    // y方向グリッド数
    int _ngb;                    // バッファ領域グリッド数
    int _nx;                    // x方向グリッド数（バッファ込み）
    int _ny;                    // y方向グリッド数（バッファ込み）
    float _dx;                  // x方向グリッド幅
    float _dy;                  // y方向グリッド幅
    
    int _weightOrder;            // weighting のオーダー
}
@end

// プライベートインスタンス変数
@implementation EMField {
    // BoundaryCoditions
    NSString* _BC_Xmin;
    NSString* _BC_Xmax;
    NSString* _BC_Ymin;
    NSString* _BC_Ymax;
    std::vector<double> _phi_Xmin;
    std::vector<double> _phi_Xmax;
    std::vector<double> _phi_Ymin;
    std::vector<double> _phi_Ymax;
    std::vector<double> _dphidx_Xmin;
    std::vector<double> _dphidx_Xmax;
    std::vector<double> _dphidy_Ymin;
    std::vector<double> _dphidy_Ymax;
    int _nkx;
    int _ikmin;
    int _nky;
    int _jkmin;
    double _Cd;
    double _Cx;
    double _Cy;
    std::vector<double> _phi_sol;
    std::vector<double> _rhs_BC;
    std::unique_ptr<Solver> _solver;    // smart pointer
    double _cathodePos;
}

- (instancetype)initWithDevice:(id<MTLDevice>)device withParam:(Init*)initParam  withLogger:(XmlLogger&)logger{
    self = [super init];
    if (self) {
        _device = device;
        
        // フィールドパラメータの取得
        struct ParamForField fieldParam = initParam.paramForField;
        struct ParamForComputing compParam = initParam.paramForComputing;
        std::vector<struct BoundaryConditionForField> fieldBCs = initParam.fieldBoundaries;

        // プライベート変数に格納
        _ngx = fieldParam.ngx;
        _ngy = fieldParam.ngy;
        _ngb = fieldParam.ngb;
        _dx = fieldParam.dx;
        _dy = fieldParam.dy;
        _weightOrder = fieldParam.weightOrder;

        // 配列サイズの計算(rho,Ex,Ey)
        _nx = _ngx + 2*_ngb;
        _ny = _ngy + 2*_ngb;
        NSUInteger arrSize = (_nx+1)*(_ny+1);
        NSUInteger bufferSize = sizeof(float) * arrSize;
        
        // Metal バッファの作成
        _rhoBuffer = [device newBufferWithLength:sizeof(float)*(_nx+1)*(_ny+1) options:MTLResourceStorageModeShared];
        _ExBuffer  = [device newBufferWithLength:sizeof(float)*(_nx+2)*(_ny+1) options:MTLResourceStorageModeShared];
        _EyBuffer  = [device newBufferWithLength:sizeof(float)*(_nx+1)*(_ny+2) options:MTLResourceStorageModeShared];
        _EzBuffer  = [device newBufferWithLength:sizeof(float)*(_nx+1)*(_ny+1) options:MTLResourceStorageModeShared];
        _BxBuffer  = [device newBufferWithLength:sizeof(float)*(_nx+1)*(_ny+2) options:MTLResourceStorageModeShared];
        _ByBuffer  = [device newBufferWithLength:sizeof(float)*(_nx+2)*(_ny+1) options:MTLResourceStorageModeShared];
        _BzBuffer  = [device newBufferWithLength:sizeof(float)*(_nx+2)*(_ny+2) options:MTLResourceStorageModeShared];
        
        // malloc
        _phi = (double*)malloc(sizeof(double) * (_nx+3)*(_ny+3) );
        
        // initialize E
        float* Ex = (float*)_ExBuffer.contents;
        float* Ey = (float*)_EyBuffer.contents;
        float* Ez = (float*)_EzBuffer.contents;
        if ([fieldParam.InitTypeE isEqualToString:@"Uniform"]){
            for (int i = 0; i < (_nx+2)*(_ny+1); i++) {
                Ex[i] = fieldParam.ampE[0];
            }
            for (int i = 0; i < (_nx+1)*(_ny+2); i++) {
                Ey[i] = fieldParam.ampE[1];
            }
            for (int i = 0; i < (_nx+1)*(_ny+1); i++) {
                Ez[i] = fieldParam.ampE[2];
            }
        }
        // initialize B
        float* Bx = (float*)_BxBuffer.contents;
        float* By = (float*)_ByBuffer.contents;
        float* Bz = (float*)_BzBuffer.contents;
        if ([fieldParam.InitTypeB isEqualToString:@"Uniform"]){
            for (int i = 0; i < (_nx+1)*(_ny+2); i++) {
                Bx[i] = fieldParam.ampB[0];
            }
            for (int i = 0; i < (_nx+1)*(_ny+2); i++) {
                By[i] = fieldParam.ampB[1];
            }
            for (int i = 0; i < (_nx+2)*(_ny+2); i++) {
                Bz[i] = fieldParam.ampB[2];
            }
        }else if ([fieldParam.InitTypeB isEqualToString:@"From1dXFile"]){
            std::vector<float> Bx_input;
            std::vector<float> By_input;
            std::vector<float> Bz_input;
            [self load1dField:Bx_input withFilePath:fieldParam.FilePathBx];
            [self load1dField:By_input withFilePath:fieldParam.FilePathBy];
            [self load1dField:Bz_input withFilePath:fieldParam.FilePathBz];
            for (int i = 0; i <= _nx; i++) {
                for (int j = 0; j <= _ny+1; j++) {
                    Bx[i+j*(_nx+1)] = Bx_input[i]*TtoG;
                }
            }
            for (int i = 0; i <= _nx+1; i++) {
                for (int j = 0; j <= _ny; j++) {
                    By[i+j*(_nx+2)] = By_input[i]*TtoG;
                }
            }
            for (int i = 0; i <= _nx+1; i++) {
                for (int j = 0; j <= _ny+1; j++) {
                    Bz[i+j*(_nx+2)] = Bz_input[i]*TtoG;
                }
            }
        }

        // construct Poisson solver
        float val_Xmin, val_Xmax, val_Ymin, val_Ymax;
        for (int pos = 0; pos < fieldBCs.size(); pos++){
            if ([fieldBCs[pos].position isEqualToString:@"Xmin"]){
                _BC_Xmin = fieldBCs[pos].type;
                if ([fieldBCs[pos].type isEqualToString:@"Dirichlet"]){
                    _phi_Xmin.assign(_ngy+1, fieldBCs[pos].val*VtosV);
                }else if ([fieldBCs[pos].type isEqualToString:@"Neumann"]){
                    _dphidx_Xmin.assign(_ngy+1, fieldBCs[pos].val);
                }
            }else if ([fieldBCs[pos].position isEqualToString:@"Xmax"]){
                _BC_Xmax = fieldBCs[pos].type;
                if ([fieldBCs[pos].type isEqualToString:@"Dirichlet"]){
                    _phi_Xmax.assign(_ngy+1, fieldBCs[pos].val*VtosV);
                }else if ([fieldBCs[pos].type isEqualToString:@"Neumann"]){
                    _dphidx_Xmax.assign(_ngy+1, fieldBCs[pos].val);
                }
            }else if ([fieldBCs[pos].position isEqualToString:@"Ymin"]){
                _BC_Ymin = fieldBCs[pos].type;
                if ([fieldBCs[pos].type isEqualToString:@"Dirichlet"]){
                    _phi_Ymin.assign(_ngx+1, fieldBCs[pos].val*VtosV);
                }else if ([fieldBCs[pos].type isEqualToString:@"Neumann"]){
                    _dphidy_Ymin.assign(_ngx+1, fieldBCs[pos].val);
                }
            }else if ([fieldBCs[pos].position isEqualToString:@"Ymax"]){
                _BC_Ymax = fieldBCs[pos].type;
                if ([fieldBCs[pos].type isEqualToString:@"Dirichlet"]){
                    _phi_Ymax.assign(_ngx+1, fieldBCs[pos].val*VtosV);
                }else if ([fieldBCs[pos].type isEqualToString:@"Neumann"]){
                    _dphidy_Ymax.assign(_ngx+1, fieldBCs[pos].val);
                }
            }else if ([fieldBCs[pos].position isEqualToString:@"hollow-cathode"]){
                _cathodePos = fieldBCs[pos].val;
            }
        }

        // 定数係数
        _Cd = -2*(1.0/_dx/_dx + 1.0/_dy/_dy);
        _Cx = 1.0/_dx/_dx;
        _Cy = 1.0/_dy/_dy;
        
        // CSR 行列
        std::vector<int> ptr;
        std::vector<int> col;
        std::vector<double> val;

        // 解くべき行列サイズ(境界条件に依存)
        _nkx = _ngx;
        _ikmin = _ngb;
        if ([_BC_Xmin isEqualToString:@"Dirichlet"]){
            _nkx--;
            _ikmin++;
        }
        if ([_BC_Xmax isEqualToString:@"Dirichlet"]){
            _nkx--;
        }
        if ([_BC_Xmin isEqualToString:@"periodic"]){
            _nkx--;
        }
        _nky = _ngy;
        _jkmin = _ngb;
        if ([_BC_Ymin isEqualToString:@"Dirichlet"]){
            _nky--;
            _jkmin++;
        }
        if ([_BC_Ymax isEqualToString:@"Dirichlet"]){
            _nky--;
        }
        if ([_BC_Ymin isEqualToString:@"periodic"]){
            _nky--;
        }
        arrSize = (_nkx+1)*(_nky+1);

        // construct CSR matrix
        ptr.push_back(0);
        for (int j = 0; j <= _nky; ++j) {
            for (int i = 0; i <= _nkx; ++i) {

                int k = i + j*(_nkx+1);
                int row_entries = 0;

                // Neumann(前進/後進 差分)
                if (i == 0 && [_BC_Xmin isEqualToString:@"Neumann"]){
                    col.push_back(k);
                    val.push_back(-1.0/_dx);
                    row_entries++;
                    col.push_back(k+1);
                    val.push_back(1.0/_dx);
                    row_entries++;
                    // 右辺ベクトルに対する境界条件の寄与
                    _rhs_BC.push_back(_dphidx_Xmin[j]);
                    // ptr を更新してループを抜ける
                    ptr.push_back(ptr.back() + row_entries);
                    continue;
                }else if(i == _nkx && [_BC_Xmax isEqualToString:@"Neumann"]){
                    col.push_back(k);
                    val.push_back(1.0/_dx);
                    row_entries++;
                    col.push_back(k-1);
                    val.push_back(-1.0/_dx);
                    row_entries++;
                    _rhs_BC.push_back(_dphidx_Xmax[j]);
                    // ptr を更新してループを抜ける
                    ptr.push_back(ptr.back() + row_entries);
                    continue;
                }else if (j == 0 && [_BC_Ymin isEqualToString:@"Neumann"]){
                    col.push_back(k);
                    val.push_back(-1.0/_dy);
                    row_entries++;
                    col.push_back(k+(_nkx+1));
                    val.push_back(1.0/_dy);
                    row_entries++;
                    _rhs_BC.push_back(_dphidy_Ymin[i]);
                    // ptr を更新してループを抜ける
                    ptr.push_back(ptr.back() + row_entries);
                    continue;
                }else if(j == _nky && [_BC_Ymax isEqualToString:@"Neumann"]){
                    col.push_back(k);
                    val.push_back(1.0/_dx);
                    row_entries++;
                    col.push_back(k-(_nkx+1));
                    val.push_back(-1.0/_dx);
                    row_entries++;
                    _rhs_BC.push_back(_dphidy_Ymax[i]);
                    // ptr を更新してループを抜ける
                    ptr.push_back(ptr.back() + row_entries);
                    continue;
                }

                // diagonal
                col.push_back(k);
                val.push_back(_Cd);
                row_entries++;
                _rhs_BC.push_back(0.0);

                // imin側
                if (i == 0){
                    if ([_BC_Xmin isEqualToString:@"Dirichlet"]){
                        // ディリクレ境界に接してるなら b に足す
                        _rhs_BC[k] += -_Cx*_phi_Xmin[j];
                    }else if ([_BC_Xmin isEqualToString:@"periodic"]){
                        // 要素の追加
                        col.push_back(k+_nkx);
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
                if (i == _nkx){
                    if ([_BC_Xmax isEqualToString:@"Dirichlet"]){
                        _rhs_BC[k] += -_Cx*_phi_Xmax[j];
                    }else if ([_BC_Xmax isEqualToString:@"periodic"]){
                        // 要素の追加
                        col.push_back(k-_nkx);
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
                        _rhs_BC[k] += -_Cy*_phi_Ymin[j];
                    }else if ([_BC_Ymin isEqualToString:@"periodic"]){
                        // 要素の追加
                        col.push_back(k+_nky*(_nkx+1));
                        val.push_back(_Cy);
                        row_entries++;
                    }
                }else{
                    col.push_back(k-(_nkx+1));
                    val.push_back(_Cy);
                    row_entries++;
                }

                // jmax側
                if (j == _nky){
                    if ([_BC_Ymax isEqualToString:@"Dirichlet"]){
                        _rhs_BC[k] += -_Cy*_phi_Ymax[i];
                    }else if ([_BC_Ymax isEqualToString:@"periodic"]){
                        // 要素の追加
                        col.push_back(k-_nky*(_nkx+1));
                        val.push_back(_Cy);
                        row_entries++;
                    }
                }else{
                    // 要素の追加
                    col.push_back(k+(_nkx+1));
                    val.push_back(_Cy);
                    row_entries++;
                }

                // CSR ポインタ更新
                ptr.push_back(ptr.back() + row_entries);
            }
        }
        // csr 行列の出力（対称性チェック用）
        // [self outputCSRmtx:arrSize row:ptr collumn:col value:val];

        // 初期解
        _phi_sol.assign(arrSize, 0.0);

        // --- amgcl の設定 ---
        if (![_BC_Xmin isEqualToString:@"Neumann"] && ![_BC_Xmin isEqualToString:@"Neumann"]
        &&  ![_BC_Ymin isEqualToString:@"Neumann"] && ![_BC_Ymin isEqualToString:@"Neumann"]){
            // 1) Prepare the single params object:
            Solver::params prm;

            // 2) AMG (preconditioner) parameters live in prm.precond:
            //   2a) Smoothed aggregation settings:
            prm.precond.coarsening.aggr.eps_strong  = compParam.aggrThreshold;  // 強結合判定しきい値（大きくすると粗グリッドが粗くなる）
            //   2b) AMG Cycle
            prm.precond.ncycle                      = compParam.amgCycleType; // サイクルタイプ（1: V-cycle, 2: W-cycle）

            // 3) Krylov‐solver (CG) parameters live in prm.solver:
            prm.solver.maxiter                = compParam.maxiter;
            prm.solver.tol                    = compParam.tolerance;

            // 4) Construct the solver once (setup), reuse later in your time‐loop:
            _solver = std::make_unique<Solver>(
                std::tie(arrSize, ptr, col, val),
                prm
            );
        }else{
            // Neumann 境界が一つでもあればエラー終了。_solver を使い回すので、メンバ変数の型定義を変更する必要がある。面倒なのでいったん非対応。
            NSLog(@"Neumann boundary condition is included. Switch solver typedef from CG to BiCGStab.");
            return nil;
            // Neumann 境界が一つでもあれば BiCGStab を使用（node-center での離散化なので）
            // Solver::params prm;
            // prm.solver.maxiter = compParam.maxiter;
            // prm.solver.tol = compParam.tolerance;
            // // _solver を生成してインスタンス変数として確保
            // _solver = std::unique_ptr<Solver>(
            //     new Solver( std::tie(arrSize, ptr, col, val), prm )
            // );
        }
    }
    return self;
}

- (bool)load1dField:(std::vector<float>&)field withFilePath:(NSString*)filePath{
    FILE* fp = std::fopen([filePath UTF8String], "rb");
    if (!fp) {
        NSLog(@"Error: Unable to open file %@ for reading", filePath);
        return false;
    }

    // ヘッダ読み込み
    int n;
    float d;
    if (std::fread(&n,  sizeof(int),   1, fp) != 1 ||
        std::fread(&d,  sizeof(float), 1, fp) != 1)
    {
        NSLog(@"Error: Unable to open file %@ for reading", filePath);
        std::fclose(fp);
        return false;
    }

    // データ本体を読み込み
    int arrSize = n+1;
    field.resize(arrSize);
    if (std::fread(field.data(), sizeof(float), arrSize, fp) != arrSize) {
        NSLog(@"Error: Unable to open file %@ for reading", filePath);
        std::fclose(fp);
        return false;
    }

    std::fclose(fp);
    return true;
}

- (bool)load2dField:(std::vector<float>&)field withFilePath:(NSString*)filePath{
    FILE* fp = std::fopen([filePath UTF8String], "rb");
    if (!fp) {
        NSLog(@"Error: Unable to open file %@ for reading", filePath);
        return false;
    }

    // ヘッダ読み込み
    int nx,ny;
    float dx,dy;
    if (std::fread(&nx,  sizeof(int),   1, fp) != 1 ||
        std::fread(&ny,  sizeof(int),   1, fp) != 1 ||
        std::fread(&dx,  sizeof(float), 1, fp) != 1 ||
        std::fread(&dy,  sizeof(float), 1, fp) != 1)
    {
        NSLog(@"Error: Unable to open file %@ for reading", filePath);
        std::fclose(fp);
        return false;
    }

    // データ本体を読み込み
    int arrSize = (nx+1)*(ny+1);
    field.resize(arrSize);
    if (std::fread(field.data(), sizeof(float), arrSize, fp) != arrSize) {
        std::cerr << "Error: Failed to read field data\n";
        std::fclose(fp);
        return false;
    }

    std::fclose(fp);
    return true;
}

- (void)outputCSRmtx:(int)n row:(std::vector<int>)row_ptr collumn:(std::vector<int>)col_idx value:(std::vector<float>)val{
    // MM 形式でファイルに書き出し (coordinate 形式)
    std::ofstream ofs("poisson.mtx");
    ofs << "%%MatrixMarket matrix coordinate real general\n";
    int nz = 0;
    for (int i = 0; i < n; ++i) {
        for (int idx = row_ptr[i]; idx < row_ptr[i+1]; ++idx) {
            ++nz;
        }
    }
    ofs << n << " " << n << " " << nz << "\n";

    for (int i = 0; i < n; ++i) {
        for (int idx = row_ptr[i]; idx < row_ptr[i+1]; ++idx) {
            int j = col_idx[idx];
            float v = val[idx];
            // MM 形式は 1-based index
            ofs << (i+1) << " " << (j+1) << " " << v << "\n";
        }
    }
    ofs.close();
}

- (void)solvePoisson:(XmlLogger&)logger{
    // 電荷密度の取得
    float* rho = (float*)_rhoBuffer.contents;

    // 周期境界の処理
    int idx_in, idx_out;
    if ([_BC_Xmin isEqualToString:@"periodic"]){
        for (int j = 0; j <= _ny; j++){
            // min側（境界の値を持つ）
            for (int i = 0; i < _ngb; i++){
                idx_in = i + j*(_nx+1);
                idx_out = _ngx+i + j*(_nx+1); // _ngb+_ngx-_ngb+i のように、オフセット分が打ち消しあう
                rho[idx_out] += rho[idx_in];
            }
            // max側（境界の値は持たない）
            for (int i = 0; i <= _ngb; i++){
                idx_in = _ngb+_ngx+i + j*(_nx+1);
                idx_out = _ngb+i + j*(_nx+1);
                rho[idx_out] += rho[idx_in];
            }
        }
    }
    if ([_BC_Ymin isEqualToString:@"periodic"]){
        for (int i = 0; i <= _nx; i++){
            // min側（境界の値を持つ）
            for (int j = 0; j < _ngb; j++){
                idx_in = i + j*(_nx+1);
                idx_out = i + (_ngy+j)*(_nx+1);
                rho[idx_out] += rho[idx_in];
            }
            // max側（境界の値は持たない）
            for (int j = 0; j <= _ngb; j++){
                idx_in = i + (_ngb+_ngy+j)*(_nx+1);
                idx_out = i + (_ngb+j)*(_nx+1);
                rho[idx_out] += rho[idx_in];
            }
        }
    }
    // Derichlet, Neumann の時は領域端の電荷密度を参照しないので無視。
    // 電荷密度保存を考えると領域内に足し込んだ方がいいのかもしれないが、いったん放置。

    // 右辺ベクトルの更新
    int arrSize = (_nkx+1)*(_nky+1);
    std::vector<double> rhs;
    for (int j = 0; j <= _nky; j++){
        for (int i = 0; i <= _nkx; i++){
            int k = i + j*(_nkx+1);
            rhs.push_back(_rhs_BC[k] -4*PI*static_cast<double>(rho[(i+_ikmin)+(j+_jkmin)*(_nx+1)]));
        }
    }

    // solve（solver, phi_sol は使い回す）
    int iters;
    double error;
    std::tie(iters, error) = (*_solver)(rhs, _phi_sol);

    // index
    bool isLeft, isRight, isBottom, isTop;
    int idx_in_i, idx_in_j;
    double dphidx, dphidy;

    // phi の全セルを走査して埋める
    for (int j = 0; j <= _ny+2; ++j) {
        // Ey11[0:1,0:2]
        if ([_BC_Ymin isEqualToString:@"Dirichlet"]){
            isBottom = (j <= _ngb+2);
            idx_in_j = j-(_ngb+2)-1;
        } else if ([_BC_Ymin isEqualToString:@"Neumann"]){
            isBottom = (j < _ngb+2);
            idx_in_j = j-(_ngb+2);
        } else if ([_BC_Ymin isEqualToString:@"periodic"]){
            isBottom = false;
            isTop = false;
            idx_in_j = (j-(_ngb+2)+_ngy)%_ngy;
        }
        // Ey13[0:1,0:2]
        if ([_BC_Ymax isEqualToString:@"Dirichlet"]){
            isTop = (j >= _ngb+2+_ngy);
        } else if ([_BC_Ymax isEqualToString:@"Neumann"]){
            isTop = (j > _ngb+2+_ngy);
        }

        for (int i = 0; i <= _nx+2; ++i) {
            // Ex11[0:2,0:1]
            if ([_BC_Xmin isEqualToString:@"Dirichlet"]){
                isLeft = (i <= _ngb+2);
                idx_in_i = i-(_ngb+2)-1;
            } else if ([_BC_Xmin isEqualToString:@"Neumann"]){
                isLeft = (i < _ngb+2);
                idx_in_i = i-(_ngb+2);
            } else if ([_BC_Xmin isEqualToString:@"periodic"]){
                // idx を調整するので無効化
                isLeft = false;
                isRight = false;
                idx_in_i = (i-(_ngb+2)+_ngx)%_ngx;
            }
            // Ex13[0:2,0:1]
            if ([_BC_Xmax isEqualToString:@"Dirichlet"]){
                isRight = (i >= _ngb+2+_ngx);
            } else if ([_BC_Xmax isEqualToString:@"Neumann"]){
                isRight = (i > _ngb+2+_ngx);
            }

            // 出力先アドレス
            idx_out = i + j*(_nx+3);
            // 領域内の値をコピー
            if (!isLeft && !isRight && !isBottom && !isTop){
                idx_in = idx_in_i + idx_in_j*(_nkx+1);
                _phi[idx_out] = _phi_sol[idx_in];
            }
            // 境界の値を線形近似(Left)
            if (isLeft && !isRight && !isBottom && !isTop){
                if ([_BC_Xmin isEqualToString:@"Dirichlet"]){
                    dphidx = _phi_sol[0 + idx_in_j*(_nkx+1)] - _phi_Xmin[idx_in_j];
                    _phi[idx_out] = _phi_Xmin[idx_in_j] + (i - (_ngb+2))*dphidx;
                } else if ([_BC_Xmin isEqualToString:@"Neumann"]){
                    dphidx = _dphidx_Xmin[idx_in_j];
                    _phi[idx_out] = _phi_sol[0 + idx_in_j*(_nkx+1)] + (i - (_ngb+2))*dphidx;
                }
            }// 境界の値を線形近似(Right)
            if (!isLeft && isRight && !isBottom && !isTop){
                if ([_BC_Xmax isEqualToString:@"Dirichlet"]){
                    dphidx = _phi_Xmax[idx_in_j] - _phi_sol[_nkx + idx_in_j*(_nkx+1)];
                    _phi[idx_out] = _phi_Xmax[idx_in_j] + (i - (_ngb+2+_ngx))*dphidx;
                } else if ([_BC_Xmax isEqualToString:@"Neumann"]){
                    dphidx = _dphidx_Xmax[idx_in_j];
                    _phi[idx_out] = _phi_sol[_nkx + idx_in_j*(_nkx+1)] + (i - (_ngb+2+_ngx))*dphidx;
                }
            }
            // 境界の値を線形近似(Bottom)
            if (!isLeft && !isRight && isBottom && !isTop){
                if ([_BC_Ymin isEqualToString:@"Dirichlet"]){
                    dphidy = _phi_sol[idx_in_i + 0*(_nkx+1)] - _phi_Ymin[idx_in_i];
                    _phi[idx_out] = _phi_Ymin[idx_in_i] + (j - (_ngb+2))*dphidy;
                } else if ([_BC_Ymin isEqualToString:@"Neumann"]){
                    dphidy = _dphidy_Ymin[idx_in_i];
                    _phi[idx_out] = _phi_sol[idx_in_i + 0*(_nkx+1)] + (j - (_ngb+2))*dphidy;
                }
            }
            // 境界の値を線形近似(Top)
            if (!isLeft && !isRight && !isBottom && isTop){
                if ([_BC_Ymax isEqualToString:@"Dirichlet"]){
                    dphidy = _phi_Ymax[idx_in_i] - _phi_sol[idx_in_i + _nky*(_nkx+1)];
                    _phi[idx_out] = _phi_Ymax[idx_in_i] + (j-(_ngb+2+_ngy))*dphidy;
                } else if ([_BC_Ymax isEqualToString:@"Neumann"]){
                    dphidy = _dphidy_Ymax[idx_in_i];
                    _phi[idx_out] = _phi_sol[idx_in_i + _nky*(_nkx+1)] + (j-(_ngb+2+_ngy))*dphidy;
                }
            }
            // 境界の値を線形近似(LeftBottom)
            if (isLeft && !isRight && isBottom && !isTop){
                if ([_BC_Xmin isEqualToString:@"Dirichlet"] && [_BC_Ymin isEqualToString:@"Dirichlet"]){
                    //// 境界条件との差を使用（微分値が大きくなりすぎることがあるので不採用）
                    // dphidx = _phi_sol[0 + 0*(_nkx+1)] - _phi_Xmin[0];
                    // dphidy = _phi_sol[0 + 0*(_nkx+1)] - _phi_Ymin[0];
                    // 解析領域内の差を使用
                    dphidx = _phi_sol[1 + 0*(_nkx+1)] - _phi_sol[0 + 0*(_nkx+1)];
                    dphidy = _phi_sol[0 + 1*(_nkx+1)] - _phi_sol[0 + 0*(_nkx+1)];
                    _phi[idx_out] = _phi_Xmin[0] + (i-(_ngb+2))*dphidx + (j-(_ngb+2))*dphidy;
                } else if ([_BC_Xmin isEqualToString:@"Neumann"] && [_BC_Ymin isEqualToString:@"Dirichlet"]){
                    dphidx = _dphidx_Xmin[0];
                    // dphidy = _phi_sol[0 + 0*(_nkx+1)] - _phi_Ymin[0];
                    dphidy = _phi_sol[0 + 1*(_nkx+1)] - _phi_sol[0 + 0*(_nkx+1)];
                    _phi[idx_out] = _phi_Ymin[0] + (i-(_ngb+2))*dphidx + (j-(_ngb+2))*dphidy;
                } else if ([_BC_Xmin isEqualToString:@"Dirichlet"] && [_BC_Ymin isEqualToString:@"Neumann"]){
                    // dphidx = _phi_sol[0 + 0*(_nkx+1)] - _phi_Xmin[0];
                    dphidx = _phi_sol[1 + 0*(_nkx+1)] - _phi_sol[0 + 0*(_nkx+1)];
                    dphidy = _dphidy_Ymin[0];
                    _phi[idx_out] = _phi_Xmin[0] + (i-(_ngb+2))*dphidx + (j-(_ngb+2))*dphidy;
                } else if ([_BC_Xmin isEqualToString:@"Neumann"] && [_BC_Ymin isEqualToString:@"Neumann"]){
                    dphidx = _dphidx_Xmin[0];
                    dphidy = _dphidy_Ymin[0];
                    _phi[idx_out] = _phi_sol[0 + 0*(_nkx+1)] + (i-(_ngb+2))*dphidx + (j-(_ngb+2))*dphidy;
                }
            }
            // 境界の値を線形近似(RightBottom)
            if (!isLeft && isRight && isBottom && !isTop){
                if ([_BC_Xmax isEqualToString:@"Dirichlet"] && [_BC_Ymin isEqualToString:@"Dirichlet"]){
                    // dphidx = _phi_Xmax[0] - _phi_sol[_nkx + 0*(_nkx+1)];
                    // dphidy = _phi_sol[_nkx + 0*(_nkx+1)] - _phi_Ymin[_nkx];
                    dphidx = _phi_sol[_nkx + 0*(_nkx+1)] - _phi_sol[_nkx-1 + 0*(_nkx+1)];
                    dphidy = _phi_sol[_nkx + 1*(_nkx+1)] - _phi_sol[_nkx + 0*(_nkx+1)];
                    _phi[idx_out] = _phi_Xmax[0] + (i-(_ngb+2+_ngx))*dphidx + (j-(_ngb+2))*dphidy;
                } else if ([_BC_Xmax isEqualToString:@"Neumann"] && [_BC_Ymin isEqualToString:@"Dirichlet"]){
                    dphidx = _dphidx_Xmax[0];
                    // dphidy = _phi_sol[_nkx + 0*(_nkx+1)] - _phi_Ymin[_nkx];
                    dphidy = _phi_sol[_nkx + 1*(_nkx+1)] - _phi_sol[_nkx + 0*(_nkx+1)];
                    _phi[idx_out] = _phi_Ymin[0] + (i-(_ngb+2+_ngx))*dphidx + (j-(_ngb+2))*dphidy;
                } else if ([_BC_Xmax isEqualToString:@"Dirichlet"] && [_BC_Ymin isEqualToString:@"Neumann"]){
                    // dphidx = _phi_Xmax[0] - _phi_sol[_nkx + 0*(_nkx+1)];
                    dphidx = _phi_sol[_nkx + 0*(_nkx+1)] - _phi_sol[_nkx-1 + 0*(_nkx+1)];
                    dphidy = _dphidy_Ymin[_nkx];
                    _phi[idx_out] = _phi_Xmax[0] + (i-(_ngb+2+_ngx))*dphidx + (j-(_ngb+2))*dphidy;
                } else if ([_BC_Xmax isEqualToString:@"Neumann"] && [_BC_Ymin isEqualToString:@"Neumann"]){
                    dphidx = _dphidx_Xmax[0];
                    dphidy = _dphidy_Ymin[_nkx];
                    _phi[idx_out] = _phi_sol[_nkx + 0*(_nkx+1)] + (i-(_ngb+2+_ngx))*dphidx + (j-(_ngb+2))*dphidy;
                }
            }
            // 境界の値を線形近似(LeftTop)
            if (isLeft && !isRight && !isBottom && isTop){
                if ([_BC_Xmin isEqualToString:@"Dirichlet"] && [_BC_Ymax isEqualToString:@"Dirichlet"]){
                    // dphidx = _phi_Xmin[_nky] - _phi_sol[0 + _nky*(_nkx+1)];
                    // dphidy = _phi_sol[0 + _nky*(_nkx+1)] - _phi_Ymax[0];
                    dphidx = _phi_sol[1 + _nky*(_nkx+1)] - _phi_sol[0 + _nky*(_nkx+1)];
                    dphidy = _phi_sol[0 + _nky*(_nkx+1)] - _phi_sol[0 + (_nky-1)*(_nkx+1)];
                    _phi[idx_out] = _phi_Xmin[_nky] + (i-(_ngb+2))*dphidx + (j-(_ngb+2+_ngy))*dphidy;
                } else if ([_BC_Xmin isEqualToString:@"Neumann"] && [_BC_Ymax isEqualToString:@"Dirichlet"]){
                    dphidx = _dphidx_Xmin[_nky];
                    // dphidy = _phi_sol[0 + _nky*(_nkx+1)] - _phi_Ymax[0];
                    dphidy = _phi_sol[0 + _nky*(_nkx+1)] - _phi_sol[0 + (_nky-1)*(_nkx+1)];
                    _phi[idx_out] = _phi_Ymax[_nky] + (i-(_ngb+2))*dphidx + (j-(_ngb+2+_ngy))*dphidy;
                } else if ([_BC_Xmin isEqualToString:@"Dirichlet"] && [_BC_Ymax isEqualToString:@"Neumann"]){
                    // dphidx = _phi_Xmin[_nky] - _phi_sol[0 + _nky*(_nkx+1)];
                    dphidx = _phi_sol[1 + _nky*(_nkx+1)] - _phi_sol[0 + _nky*(_nkx+1)];
                    dphidy = _dphidy_Ymax[0];
                    _phi[idx_out] = _phi_Xmin[_nky] + (i-(_ngb+2))*dphidx + (j-(_ngb+2+_ngy))*dphidy;
                } else if ([_BC_Xmin isEqualToString:@"Neumann"] && [_BC_Ymax isEqualToString:@"Neumann"]){
                    dphidx = _dphidx_Xmin[_nky];
                    dphidy = _dphidy_Ymax[0];
                    _phi[idx_out] = _phi_sol[0 + _nky*(_nkx+1)] + (i-(_ngb+2))*dphidx + (j-(_ngb+2+_ngy))*dphidy;
                }
            }
            // 境界の値を線形近似(RightTop)
            if (!isLeft && isRight && !isBottom && isTop){
                if ([_BC_Xmax isEqualToString:@"Dirichlet"] && [_BC_Ymax isEqualToString:@"Dirichlet"]){
                    // dphidx = _phi_Xmax[_nky] - _phi_sol[_nkx + _nky*(_nkx+1)];
                    // dphidy = _phi_sol[_nkx + _nky*(_nkx+1)] - _phi_Ymax[_nkx];
                    dphidx = _phi_sol[_nkx + _nky*(_nkx+1)] - _phi_sol[_nkx-1 + _nky*(_nkx+1)];
                    dphidy = _phi_sol[_nkx + _nky*(_nkx+1)] - _phi_sol[_nkx + (_nky-1)*(_nkx+1)];
                    _phi[idx_out] = _phi_Xmax[_nky] + (i-(_ngb+2+_ngx))*dphidx + (j-(_ngb+2+_ngy))*dphidy;
                } else if ([_BC_Xmax isEqualToString:@"Neumann"] && [_BC_Ymax isEqualToString:@"Dirichlet"]){
                    dphidx = _dphidx_Xmax[_nky];
                    // dphidy = _phi_sol[_nkx + _nky*(_nkx+1)] - _phi_Ymax[_nkx];
                    dphidy = _phi_sol[_nkx + _nky*(_nkx+1)] - _phi_sol[_nkx + (_nky-1)*(_nkx+1)];
                    _phi[idx_out] = _phi_Ymax[_nky] + (i-(_ngb+2+_ngx))*dphidx + (j-(_ngb+2+_ngy))*dphidy;
                } else if ([_BC_Xmax isEqualToString:@"Dirichlet"] && [_BC_Ymax isEqualToString:@"Neumann"]){
                    // dphidx = _phi_Xmax[_nky] - _phi_sol[_nkx + _nky*(_nkx+1)];
                    dphidx = _phi_sol[_nkx + _nky*(_nkx+1)] - _phi_sol[_nkx-1 + _nky*(_nkx+1)];
                    dphidy = _dphidy_Ymax[_nkx];
                    _phi[idx_out] = _phi_Xmax[_nky] + (i-(_ngb+2+_ngx))*dphidx + (j-(_ngb+2+_ngy))*dphidy;
                } else if ([_BC_Xmax isEqualToString:@"Neumann"] && [_BC_Ymax isEqualToString:@"Neumann"]){
                    dphidx = _dphidx_Xmax[_nky];
                    dphidy = _dphidy_Ymax[_nkx];
                    _phi[idx_out] = _phi_sol[_nkx + _nky*(_nkx+1)] + (i-(_ngb+2+_ngx))*dphidx + (j-(_ngb+2+_ngy))*dphidy;
                }
            }
        }
    }

    // adjust potential at hollow-cathode
    double mean = 0.0;
    for (int j = 0; j <= _nky; j++){
        idx_in = (int)(_cathodePos/_dx) + j*(_nkx+1);
        mean += _phi_sol[idx_in];
    }
    mean /= _nky+1;
    for (int j = 0; j <= _ny+2; j++){
        for (int i = 0; i <= _nx+2; i++){
            int k = i + j*(_nx+3);
            _phi[k] -= mean*((i-_ngb-2)*_dx)/_cathodePos;
        }
    }

    // 収束状況
    std::map<std::string, std::string> data ={
        {"iteration", std::to_string(iters)},
        {"error", fmtSci(error, 6)},
        {"meanCathode", fmtSci(mean*sVtoV, 6)},
    };
    logger.logSection("solvePoisson", data);
    // NSLog(@"mean = %e",mean*sVtoV);

    // metalバッファの取得
    float* Ex = (float*)_ExBuffer.contents;
    float* Ey = (float*)_EyBuffer.contents;
    // 勾配計算
    for (int i = 0; i <= _nx+1; ++i) {
        for (int j = 0; j <= _ny; ++j) {
            Ex[i + j*(_nx+2)] = -(_phi[i+1 + (j+1)*(_nx+3)] - _phi[i   + (j+1)*(_nx+3)])/_dx;
        }
    }
    for (int i = 0; i <= _nx; ++i) {
        for (int j = 0; j <= _ny+1; ++j) {
            Ey[i + j*(_nx+1)] = -(_phi[i+1 + (j+1)*(_nx+3)] - _phi[i+1 + (j  )*(_nx+3)])/_dy;
        }
    }

}

- (void)resetChargeDensity{
    // 電荷密度の取得
    float* rho = (float*)_rhoBuffer.contents;
    // 電荷密度の初期化
    int arrSize = (_nx+1) * (_ny+1);
    for (int i = 0; i < arrSize; i++){
        rho[i] = 0.0f;
    }
}

- (void)checkChargeDensity{
    // 電荷密度の取得
    float* rho = (float*)_rhoBuffer.contents;
    float* arho = (float*)_atomicRhoBuffer.contents;
    // 電荷密度の差分
    float diff, min = 1e20, max = -1e20;
    int arrSize = (_nx+1) * (_ny+1);
    for (int i = 0; i < arrSize; i++){
        diff = (rho[i] - arho[i])*(rho[i] - arho[i]);
        if (min > diff){ min = diff; }
        if (max < diff){ max = diff; }
    }
    NSLog(@"checkChargeDensity: min = %e, max = %e", min, max);
}

// 電位へのアクセサ（double->float)
- (float*)phi { 
    static std::vector<float> cache;
    size_t arrSize = (_nx+3) * (_ny+3);
    if (cache.size() != arrSize) cache.resize(arrSize);
    for (size_t i = 0; i < arrSize; ++i)
        cache[i] = static_cast<float>(_phi[i]);
    return cache.data();
}
// 電荷密度へのアクセサ
- (id<MTLBuffer>)rhoBuffer { return _rhoBuffer; }
- (id<MTLBuffer>)atomicRhoBuffer { return _atomicRhoBuffer; }

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
- (int)nx { return _nx; }
- (int)ny { return _ny; }
- (int)ngb { return _ngb; }
- (double)dx { return _dx; }
- (double)dy { return _dy; }

@end