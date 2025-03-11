import numpy as np
import matplotlib.pyplot as plt

cycle = 20
filename = f"bin/field_{cycle:08}.bin"

with open(filename, 'rb') as f:
    # ヘッダ情報の読み込み
    ngx = np.frombuffer(f.read(4), dtype=np.int32)[0]
    ngy = np.frombuffer(f.read(4), dtype=np.int32)[0]
    ngb = np.frombuffer(f.read(4), dtype=np.int32)[0]
    dx = np.frombuffer(f.read(4), dtype=np.float32)[0]
    dy = np.frombuffer(f.read(4), dtype=np.float32)[0]
    
    # 全グリッドサイズ（境界含む）
    total_nx = ngx + 2 * ngb
    total_ny = ngy + 2 * ngb
    total_cells = total_nx * total_ny
    
    # 各フィールドのデータ読み込み (float は 4 バイト)
    rho = np.frombuffer(f.read(total_cells * 4), dtype=np.float32).reshape((total_ny, total_nx))
    phi = np.frombuffer(f.read(total_cells * 4), dtype=np.float32).reshape((total_ny, total_nx))
    Ex = np.frombuffer(f.read(total_cells * 4), dtype=np.float32).reshape((total_ny, total_nx))
    Ey = np.frombuffer(f.read(total_cells * 4), dtype=np.float32).reshape((total_ny, total_nx))
    Bz = np.frombuffer(f.read(total_cells * 4), dtype=np.float32).reshape((total_ny, total_nx))

# 物理領域のみ抽出（境界部分を除去）
rho_phys = rho[ngb:ngb+ngy, ngb:ngb+ngx]
phi_phys = phi[ngb:ngb+ngy, ngb:ngb+ngx]
Ex_phys  = Ex[ngb:ngb+ngy, ngb:ngb+ngx]
Ey_phys  = Ey[ngb:ngb+ngy, ngb:ngb+ngx]
Bz_phys  = Bz[ngb:ngb+ngy, ngb:ngb+ngx]

# 物理座標の生成
x = np.linspace(0, ngx*dx, ngx)
y = np.linspace(0, ngy*dy, ngy)
X, Y = np.meshgrid(x, y)

print(dx)
print(ngx)
print(rho)

# プロット用の関数
def plot_field(field, title):
    plt.figure()
    # pcolormesh により実スケール (x, y 座標) で描画
    plt.pcolormesh(X, Y, field, shading='auto')
    plt.colorbar()
    plt.xlabel('x (units)')
    plt.ylabel('y (units)')
    plt.title(title)
    plt.tight_layout()
    plt.show()

# 各フィールドをプロット
plot_field(rho_phys, 'Charge Density (rho)')
# plot_field(phi_phys, 'Potential (phi)')
# plot_field(Ex_phys,  'Electric Field Ex')
# plot_field(Ey_phys,  'Electric Field Ey')
# plot_field(Bz_phys,  'Magnetic Field Bz')
