// main.mm
#include <amgcl/make_solver.hpp>
#include <amgcl/solver/bicgstab.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/adapter/crs_tuple.hpp>

int main(int argc, const char * argv[]) {
    // 格子数
    int nx = 128;
    int ny = 64;
    int nb = 2; // 5th-order
    int arrSize, nkx, nky;
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
    std::vector<float> phi_imin(ny+1, 200.0);
    std::vector<float> dphidx_imin(ny+1, 0.0);
    std::string BC_imax = "Dirichlet";
    std::vector<float> phi_imax(ny+1, 0.0);
    std::vector<float> dphidx_imax(ny+1, 0.0);
    std::string BC_jmin = "Dirichlet";
    // std::string BC_jmin = "Neumann";
    std::vector<float> phi_jmin(nx+1, 0.0);
    std::vector<float> dphidy_jmin(nx+1, 0.0);
    std::string BC_jmax = "Dirichlet";
    // std::string BC_jmax = "Neumann";
    std::vector<float> phi_jmax(nx+1, 0.0);
    std::vector<float> dphidy_jmax(nx+1, 0.0);
    
    // とくべき行列サイズ
    nkx = nx-2; // 両方ディリクレなら自由度が2下がる
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
    arrSize = (nkx+1)*(nky+1);

    // 電荷密度
    std::vector<float> rho(arrSize, 0.0);

    // construct CSR matrix
    ptr.push_back(0);
    for (int j = 0; j <= nky; ++j) {
        for (int i = 0; i <= nkx; ++i) {
            int k = i + j*(nkx+1);
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
                col.push_back(k+(nkx+1));
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
                col.push_back(k-(nkx+1));
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
                    // col.push_back(k+(nkx-1)); // imax = nkx-1 で定義してたvar
                    col.push_back(k+nkx);
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
            // if (i == nkx-1){
            if (i == nkx){
                if (BC_imax == "Dirichlet"){
                    rhs[k] += -Cx*phi_imax[j];
                }else if (BC_imax == "periodic"){
                    // 要素の追加
                    // col.push_back(k-(nkx-1));
                    col.push_back(k-nkx);
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
                    // col.push_back(k+(nky-1));
                    col.push_back(k+nky);
                    val.push_back(Cx);
                    row_entries++;
                }
            }else{
                // col.push_back(k-nkx);
                col.push_back(k-(nkx+1));
                val.push_back(Cx);
                row_entries++;
            }

            // jmax側
            // if (j == nky-1){
            if (j == nky){
                if (BC_jmax == "Dirichlet"){
                    rhs[k] += -Cy*phi_jmax[i];
                }else if (BC_jmax == "periodic"){
                    // 要素の追加
                    // col.push_back(k-(nky-1)*nkx);
                    col.push_back(k-nky*(nkx+1));
                    val.push_back(Cx);
                    row_entries++;
                }
            }else{
                // 要素の追加
                // col.push_back(k+nkx);
                col.push_back(k+(nkx+1));
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
    Solver solve( std::tie(arrSize, ptr, col, val), prm );

    // 初期解（ゼロ初期化）
    std::vector<float> x(arrSize, 0.0);

    int iters;
    float error;

    // 連立方程式 A*x = rhs を解く
    std::tie(iters, error) = solve(rhs, x);

    // 収束状況の出力
    std::cout << "反復回数: " << iters << std::endl;
    std::cout << "最終誤差: " << error << std::endl;

    // 粒子用 fld 配列
    float *phi = (float *)malloc(sizeof(float)*(nx+2*nb+1+1)*(ny+2*nb+1+1));
    float *Ex = (float *)malloc(sizeof(float)*(nx+2*nb+1)*(ny+2*nb+1));
    float *Ey = (float *)malloc(sizeof(float)*(nx+2*nb+1)*(ny+2*nb+1));

    // index
    bool isLeft, isRight, isBottom, isTop;
    int idx_in, idx_out; // 入出力インデックス
    int idx_in_i, idx_in_j;
    float dphidx, dphidy;

    // phi の全セルを走査して埋める
    for (int i = 0; i <= nx+2*nb+1; ++i) {
        // Ex11[0:2,0:1]
        if (BC_imin == "Dirichlet"){
            isLeft = (i <= nb+1);
            idx_in_i = i-(nb+1)-1;
        } else if (BC_imin == "Neumann"){
            isLeft = (i < nb+1);
            idx_in_i = i-(nb+1);
        } else if (BC_imin == "periodic"){
            // idx を調整するので無効化
            isLeft = false;
            isRight = false;
            idx_in_i = (i-(nb+1)+(nx+1))%(nx+1);
        }
        // Ex13[0:1,0:1]
        if (BC_imax == "Dirichlet"){
            isRight = (i >= nb+1+nx);
        } else if (BC_imax == "Neumann"){
            isRight = (i > nb+1+nx);
        }
        for (int j = 0; j <= ny+2*nb+1; ++j) {
            // Ey11[0:1,0:2]
            if (BC_jmin == "Dirichlet"){
                isBottom = (j <= nb+1);
                idx_in_j = j-(nb+1)-1;
            } else if (BC_jmin == "Neumann"){
                isBottom = (j < nb+1);
                idx_in_j = j-(nb+1);
            } else if (BC_jmin == "periodic"){
                // idx を調整するので無効化
                isBottom = false;
                isTop = false;
                idx_in_j = (j-(nb+1)+(ny+1))%(ny+1);
            }
            // Ey13[0:1,0:1]
            if (BC_jmax == "Dirichlet"){
                isTop = (j >= nb+1+ny);
            } else if (BC_jmax == "Neumann"){
                isTop = (j > nb+1+ny);
            }

            // 出力先アドレス
            idx_out = i + j*(nx+2*nb+1+1);

            // 領域内の値をコピー
            if (!isLeft && !isRight && !isBottom && !isTop){
                idx_in = idx_in_i + idx_in_j*(nkx+1);
                phi[idx_out] = x[idx_in];
            }
            // 境界の値を扇形近似(Left)
            if (isLeft && !isRight && !isBottom && !isTop){
                if (BC_imin == "Dirichlet"){
                    dphidx = x[0 + idx_in_j*(nkx+1)] - phi_imin[idx_in_j];
                    phi[idx_out] = phi_imin[idx_in_j] + (i - (nb+1))*dphidx;
                } else if (BC_imin == "Neumann"){
                    dphidx = dphidx_imin[idx_in_j];
                    phi[idx_out] = x[0 + idx_in_j*(nkx+1)] + (i - (nb+1))*dphidx;
                }
            }
            // 境界の値を扇形近似(Right)
            if (!isLeft && isRight && !isBottom && !isTop){
                if (BC_imax == "Dirichlet"){
                    dphidx = phi_imax[idx_in_j] - x[nkx + idx_in_j*(nkx+1)];
                    phi[idx_out] = phi_imax[idx_in_j] + (i - (nb+1+nx))*dphidx;
                } else if (BC_imax == "Neumann"){
                    dphidx = dphidx_imax[idx_in_j];
                    phi[idx_out] = x[nkx + idx_in_j*(nkx+1)] + (i - (nb+1+nx))*dphidx;
                }
            }
            // 境界の値を扇形近似(Bottom)
            if (!isLeft && !isRight && isBottom && !isTop){
                if (BC_jmin == "Dirichlet"){
                    dphidy = x[idx_in_i + 0*(nkx+1)] - phi_jmin[idx_in_i];
                    phi[idx_out] = phi_jmin[idx_in_i] + (j - (nb+1))*dphidy;
                } else if (BC_jmin == "Neumann"){
                    dphidy = dphidy_jmin[idx_in_i];
                    phi[idx_out] = x[idx_in_i + 0*(nkx+1)] + (j - (nb+1))*dphidy;
                }
            }
            // 境界の値を扇形近似(Top)
            if (!isLeft && !isRight && !isBottom && isTop){
                if (BC_jmax == "Dirichlet"){
                    dphidy = phi_jmax[idx_in_i] - x[idx_in_i + nky*(nkx+1)];
                    phi[idx_out] = phi_jmax[idx_in_i] + (j - (nb+1+ny))*dphidy;
                } else if (BC_jmax == "Neumann"){
                    dphidy = dphidy_jmax[idx_in_i];
                    phi[idx_out] = x[idx_in_i + nky*(nkx+1)] + (j - (nb+1+ny))*dphidy;
                }
            }
            //// 境界チェック <- ok
            // std::cout
            // << "i = " << i
            // << ", j = " << j
            // << ", isLeft = " << isLeft
            // << ", isRight = " << isRight
            // << ", isBottom = " << isBottom
            // << ", isTop = " << isTop
            // << std::endl;
        }
    }
    // phi の全セルを走査して埋める
    for (int i = 0; i <= nx+2*nb+1; ++i) {
        // Ex11[0:2,0:1]
        if (BC_imin == "Dirichlet"){
            isLeft = (i <= nb+1);
            idx_in_i = i-(nb+1)-1;
        } else if (BC_imin == "Neumann"){
            isLeft = (i < nb+1);
            idx_in_i = i-(nb+1);
        } else if (BC_imin == "periodic"){
            // idx を調整するので無効化
            isLeft = false;
            isRight = false;
            idx_in_i = (i-(nb+1)+(nx+1))%(nx+1);
        }
        // Ex13[0:1,0:1]
        if (BC_imax == "Dirichlet"){
            isRight = (i >= nb+1+nx);
        } else if (BC_imax == "Neumann"){
            isRight = (i > nb+1+nx);
        }
        for (int j = 0; j <= ny+2*nb+1; ++j) {
            // Ey11[0:1,0:2]
            if (BC_jmin == "Dirichlet"){
                isBottom = (j <= nb+1);
                idx_in_j = j-(nb+1)-1;
            } else if (BC_jmin == "Neumann"){
                isBottom = (j < nb+1);
                idx_in_j = j-(nb+1);
            } else if (BC_jmin == "periodic"){
                // idx を調整するので無効化
                isBottom = false;
                isTop = false;
                idx_in_j = (j-(nb+1)+(ny+1))%(ny+1);
            }
            // Ey13[0:1,0:1]
            if (BC_jmax == "Dirichlet"){
                isTop = (j >= nb+1+ny);
            } else if (BC_jmax == "Neumann"){
                isTop = (j > nb+1+ny);
            }

            // 出力先アドレス
            idx_out = i + j*(nx+2*nb+1+1);

            // 境界の値を扇形近似(LeftBottom)
            if (isLeft && !isRight && isBottom && !isTop){
                // if (BC_imin == "Dirichlet"){
                //     delphi = x[0 + 0*(nkx+1)] - phi_imin[0];
                //     phi[idx_out] = phi_imin[0] + (i - (nb+1))*delphi;
                // } else if (BC_imin == "Neumann"){
                //     delphi = dphidx_imin[0];
                //     phi[idx_out] = x[0 + 0*(nkx+1)] + (i - (nb+1))*delphi;
                // }
                if (BC_imin == "Dirichlet"){
                    delphi = x[0 + idx_in_j*(nkx+1)] - phi_imin[idx_in_j];
                    phi[idx_out] = phi_imin[idx_in_j] + (i - (nb+1))*delphi;
                } else if (BC_imin == "Neumann"){
                    delphi = dphidx_imin[idx_in_j];
                    phi[idx_out] = x[0 + idx_in_j*(nkx+1)] + (i - (nb+1))*delphi;
                }
                delphi = phi[idx_in];
                phi[idx_out] = x[0 + idx_in_j*(nkx+1)] + (i - (nb+1))*delphi;
            }
            // 境界の値を扇形近似(RightBottom)
            if (!isLeft && isRight && isBottom && !isTop){
                // if (BC_jmin == "Dirichlet"){
                //     delphi = x[nkx + 0*(nkx+1)] - phi_jmin[nkx];
                //     phi[idx_out] = phi_jmin[nkx] + (j - (nb+1))*delphi;
                // } else if (BC_jmin == "Neumann"){
                //     delphi = dphidy_jmin[nkx];
                //     phi[idx_out] = x[nkx + 0*(nkx+1)] + (j - (nb+1))*delphi;
                // }
                if (BC_imax == "Dirichlet"){
                    delphi = phi_imax[idx_in_j] - x[nkx + idx_in_j*(nkx+1)];
                    phi[idx_out] = phi_imax[idx_in_j] + (i - (nb+1+nx))*delphi;
                } else if (BC_imax == "Neumann"){
                    delphi = dphidx_imax[idx_in_j];
                    phi[idx_out] = x[nkx + idx_in_j*(nkx+1)] + (i - (nb+1+nx))*delphi;
                }
            }
            // 境界の値を扇形近似(LeftTop)
            if (isLeft && !isRight && !isBottom && isTop){
                // if (BC_jmax == "Dirichlet"){
                //     delphi = phi_jmax[0] - x[0 + nky*(nkx+1)];
                //     phi[idx_out] = phi_jmax[0] + (j - (nb+1+ny))*delphi;
                // } else if (BC_jmax == "Neumann"){
                //     delphi = dphidy_jmax[0];
                //     phi[idx_out] = x[0 + nky*(nkx+1)] + (j - (nb+1+ny))*delphi;
                // }
                if (BC_imax == "Dirichlet"){
                    delphi = phi_imax[idx_in_j] - x[nkx + idx_in_j*(nkx+1)];
                    phi[idx_out] = phi_imax[idx_in_j] + (i - (nb+1+nx))*delphi;
                } else if (BC_imax == "Neumann"){
                    delphi = dphidx_imax[idx_in_j];
                    phi[idx_out] = x[nkx + idx_in_j*(nkx+1)] + (i - (nb+1+nx))*delphi;
                }
            }
            // 境界の値を扇形近似(RightTop)
            if (!isLeft && isRight && !isBottom && isTop){
                // if (BC_imax == "Dirichlet"){
                //     delphi = phi_imax[nky] - x[nkx + nky*(nkx+1)];
                //     phi[idx_out] = phi_imax[nky] + (i - (nb+1+nx))*delphi;
                // } else if (BC_imax == "Neumann"){
                //     delphi = dphidx_imax[nky];
                //     phi[idx_out] = x[nkx + nky*(nkx+1)] + (i - (nb+1+nx))*delphi;
                // }
                if (BC_jmax == "Dirichlet"){
                    delphi = phi_jmax[idx_in_i] - x[idx_in_i + nky*(nkx+1)];
                    phi[idx_out] = phi_jmax[idx_in_i] + (j - (nb+1+ny))*delphi;
                } else if (BC_jmax == "Neumann"){
                    delphi = dphidy_jmax[idx_in_i];
                    phi[idx_out] = x[idx_in_i + nky*(nkx+1)] + (j - (nb+1+ny))*delphi;
                }
            }
        }
    }

    // 勾配計算
    for (int i = 0; i <= nx+2*nb; ++i) {
        for (int j = 0; j <= ny+2*nb; ++j) {
            Ex[i + j*(nx+2*nb+1)] = (phi[i+1 + (j+1)*(nx+2*nb+1+1)] - phi[i   + (j+1)*(nx+2*nb+1+1)])/dx;
            Ey[i + j*(nx+2*nb+1)] = (phi[i+1 + (j+1)*(nx+2*nb+1+1)] - phi[i+1 + (j  )*(nx+2*nb+1+1)])/dy;
        }
    }
    
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
    fwrite(&nb, sizeof(int), 1, fp);
    fwrite(&dx, sizeof(float), 1, fp);
    fwrite(&dy, sizeof(float), 1, fp);
    
    // 結果を書き出す
    fwrite(phi, sizeof(float), (nx+2*nb+1+1)*(ny+2*nb+1+1), fp);
    fwrite(Ex, sizeof(float), (nx+2*nb+1)*(ny+2*nb+1), fp);
    fwrite(Ey, sizeof(float), (nx+2*nb+1)*(ny+2*nb+1), fp);
    
    fclose(fp);
    std::cout << "phi: " << phi[nb+1 + (nb+1+2)*(nx+2*nb+1+1)] << std::endl;
    std::cout << "done." << std::endl;

    return 0;
}
