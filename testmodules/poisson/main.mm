// main.mm
#include <amgcl/make_solver.hpp>
#include <amgcl/solver/bicgstab.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/adapter/crs_tuple.hpp>

int main(int argc, const char * argv[]) {
    // 格子点数
    int nx = 256;
    int ny = 256;
    int nk, nkx, nky;
    // 格子サイズ
    float dx = 7.8125e-3;
    float dy = 7.8125e-3;
    // 定数係数
    float Cd = -2*(1.0/dx/dx + 1.0/dy/dy);
    float Cx = 1.0/dx/dx;
    float Cy = 1.0/dy/dy;
    
    // CSR 行列
    std::vector<int> ptr;
    std::vector<int> col;
    std::vector<float> val;
    std::vector<float> rhs;

    // 境界条件
    std::string BC_imin = "Dirichlet";
    std::vector<float> phi_imin(ny, 200.0);
    std::vector<float> dphidx_imin(ny, 0.0);
    std::string BC_imax = "Dirichlet";
    std::vector<float> phi_imax(ny, 0.0);
    std::vector<float> dphidx_imax(ny, 0.0);
    std::string BC_jmin = "Neumann";
    std::vector<float> phi_jmin(nx, 0.0);
    std::vector<float> dphidy_jmin(nx, 0.0);
    std::string BC_jmax = "Neumann";
    std::vector<float> phi_jmax(nx, 0.0);
    std::vector<float> dphidy_jmax(nx, 0.0);
    
    // とくべき行列サイズ
    nkx = nx-2;
    if (BC_imin == "Neumann"){
        nkx++;
    }
    if (BC_imax == "Neumann"){
        nkx++;
    }
    if (BC_imin == "periodic"){
        nkx++;
    }
    nky = ny-2;
    if (BC_jmin == "Neumann"){
        nky++;
    }
    if (BC_jmax == "Neumann"){
        nky++;
    }
    if (BC_jmin == "periodic"){
        nky++;
    }
    nk = nkx*nky;

    // 初期条件
    std::vector<float> phi(nk, 0.0);
    // 電荷密度
    std::vector<float> rho(nk, 0.0);

    // construct CSR matrix
    ptr.push_back(0);
    for (int j = 0; j < nky; ++j) {
        for (int i = 0; i < nkx; ++i) {
            int k = i + j*nkx;
            int row_entries = 0;

            // Neumann(前進/後進 差分)
            if (i == 0 && BC_imin == "Neumann"){
                col.push_back(k);
                val.push_back(-1.0/dx);
                row_entries++;
                col.push_back(k+1);
                val.push_back(1.0/dx);
                row_entries++;
                rhs.push_back(dphidx_imin[j]);
                // ptr を更新してループを抜ける
                ptr.push_back(ptr.back() + row_entries);
                continue;
            }else if(i == nkx-1 && BC_imax == "Neumann"){
                col.push_back(k);
                val.push_back(1.0/dx);
                row_entries++;
                col.push_back(k-1);
                val.push_back(-1.0/dx);
                row_entries++;
                rhs.push_back(dphidx_imax[j]);
                // ptr を更新してループを抜ける
                ptr.push_back(ptr.back() + row_entries);
                continue;
            }else if (j == 0 && BC_jmin == "Neumann"){
                col.push_back(k);
                val.push_back(-1.0/dy);
                row_entries++;
                col.push_back(k+nkx);
                val.push_back(1.0/dy);
                row_entries++;
                rhs.push_back(dphidy_jmin[i]);
                // ptr を更新してループを抜ける
                ptr.push_back(ptr.back() + row_entries);
                continue;
            }else if(j == nky-1 && BC_jmax == "Neumann"){
                col.push_back(k);
                val.push_back(1.0/dx);
                row_entries++;
                col.push_back(k-nkx);
                val.push_back(-1.0/dx);
                row_entries++;
                rhs.push_back(dphidy_jmax[i]);
                // ptr を更新してループを抜ける
                ptr.push_back(ptr.back() + row_entries);
                continue;
            }

            // diagonal
            col.push_back(k);
            val.push_back(Cd);
            row_entries++;
            rhs.push_back(rho[k]);

            // imin側
            if (i == 0){
                if (BC_imin == "Dirichlet"){
                    // ディリクレ境界に接してるなら b に足す
                    rhs[k] += -Cx*phi_imin[j];
                }else if (BC_imin == "periodic"){
                    // 要素の追加
                    col.push_back(k+(nkx-1));
                    val.push_back(Cx);
                    row_entries++;
                }
            }else{
                // 領域内なら val を追加
                col.push_back(k-1);
                val.push_back(Cx);
                row_entries++;
            }

            // imax側
            if (i == nkx-1){
                if (BC_imax == "Dirichlet"){
                    rhs[k] += -Cx*phi_imax[j];
                }else if (BC_imax == "periodic"){
                    // 要素の追加
                    col.push_back(k-(nkx-1));
                    val.push_back(Cx);
                    row_entries++;
                }
            }else{
                col.push_back(k+1);
                val.push_back(Cx);
                row_entries++;
            }

            // jmin側
            if (j == 0){
                if (BC_jmin == "Dirichlet"){
                    rhs[k] += -Cy*phi_jmin[j];
                }else if (BC_jmin == "periodic"){
                    // 要素の追加
                    col.push_back(k+(nky-1));
                    val.push_back(Cx);
                    row_entries++;
                }
            }else{
                col.push_back(k-nkx);
                val.push_back(Cx);
                row_entries++;
            }

            // jmax側
            if (j == nky-1){
                if (BC_jmax == "Dirichlet"){
                    rhs[k] += -Cy*phi_jmax[i];
                }else if (BC_jmax == "periodic"){
                    // 要素の追加
                    col.push_back(k-(nky-1)*nkx);
                    val.push_back(Cx);
                    row_entries++;
                }
            }else{
                // 要素の追加
                col.push_back(k+nkx);
                val.push_back(Cx);
                row_entries++;
            }

            // CSR ポインタ更新
            ptr.push_back(ptr.back() + row_entries);

        }
    }

    // --- amgcl の設定 ---
    // バックエンドとソルバの型定義
    typedef amgcl::backend::builtin<float> Backend;
    typedef amgcl::make_solver<
    // Use AMG as preconditioner:
    amgcl::amg<
        Backend,
        amgcl::coarsening::smoothed_aggregation,
        amgcl::relaxation::spai0
        >,
    // And BiCGStab as iterative solver:
    amgcl::solver::bicgstab<Backend>
    > Solver;

    // ソルバの生成
    Solver::params prm;
    prm.solver.tol = 1e-5;
    Solver solve( std::tie(nk, ptr, col, val), prm );

    // 初期解（ゼロ初期化）
    std::vector<float> x(nk, 0.0);

    int iters;
    float error;

    // 連立方程式 A*x = rhs を解く
    std::tie(iters, error) = solve(rhs, x);

    // 解法結果の出力
    std::cout << "反復回数: " << iters << std::endl;
    std::cout << "最終誤差: " << error << std::endl;

    // fld 出力
    const char *filePath = "result.bin";
    FILE *fp = fopen(filePath, "wb");
    if (!fp) {
        std::cout << "Error" << std::endl;
        return 0;
    }
    
    // ヘッダ情報としてメッシュ情報を出力
    fwrite(&nx, sizeof(int), 1, fp);
    fwrite(&ny, sizeof(int), 1, fp);
    fwrite(&dx, sizeof(float), 1, fp);
    fwrite(&dy, sizeof(float), 1, fp);
    
    // 結果を書き出す
    fwrite(x.data(), sizeof(float), nk, fp);
    
    fclose(fp);
    std::cout << "done." << std::endl;

    return 0;
}
