import numpy as np
import matplotlib.pyplot as plt

filename = f"result.bin"

with open(filename, 'rb') as f:
    # ヘッダ情報の読み込み
    ngx = np.frombuffer(f.read(4), dtype=np.int32)[0]
    ngy = np.frombuffer(f.read(4), dtype=np.int32)[0]
    dx = np.frombuffer(f.read(4), dtype=np.float32)[0]
    dy = np.frombuffer(f.read(4), dtype=np.float32)[0]
    
    # 全グリッドサイズ（境界含む）
    # ngx = ngx - 2
    # ngy = ngy - 2
    ngx = ngx - 2
    ngy = ngy - 1
    total_cells = ngx * ngy
    
    # 各フィールドのデータ読み込み (float は 4 バイト)
    phi = np.frombuffer(f.read(total_cells * 4), dtype=np.float32).reshape((ngy, ngx))
    
# 物理座標の生成
x = np.linspace(0, ngx*dx, ngx)
y = np.linspace(0, ngy*dy, ngy)
X, Y = np.meshgrid(x, y)

# プロット用の関数
def plot_field(field, title, figname):
    fig, ax = plt.subplots()
    im = ax.pcolormesh(X, Y, field, shading='auto')
    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.04)
    cbar.set_label('Field value')
    ax.set_xlabel('x (units)')
    ax.set_ylabel('y (units)')
    ax.set_title(title)
    ax.set_aspect('equal', adjustable='box')  # 軸のスケールを等しくする
    fig.tight_layout()
    fig.savefig(figname)
    plt.close(fig)

# 各フィールドをプロット
plot_field(phi, 'Potential (phi)', 'result.png')
