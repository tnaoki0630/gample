import numpy as np
import matplotlib.pyplot as plt

filename = f"result.bin"

with open(filename, 'rb') as f:
    # ヘッダ情報の読み込み
    nx = np.frombuffer(f.read(4), dtype=np.int32)[0]
    ny = np.frombuffer(f.read(4), dtype=np.int32)[0]
    nb = np.frombuffer(f.read(4), dtype=np.int32)[0]
    dx = np.frombuffer(f.read(4), dtype=np.float32)[0]
    dy = np.frombuffer(f.read(4), dtype=np.float32)[0]
    
    # phi の読み込み (float は 4 バイト) 
    nkx = nx+2*nb+1
    nky = ny+2*nb+1
    arrSize = (nkx+1)*(nky+1)
    phi_raw = np.frombuffer(f.read(arrSize * 4), dtype=np.float32)
    phi = phi_raw.reshape((nky+1, nkx+1))
    # Ex, Ey の読み込み (float は 4 バイト) 
    nkx = nx+2*nb
    nky = ny+2*nb
    arrSize = (nkx+1)*(nky+1)
    Ex = np.frombuffer(f.read(arrSize * 4), dtype=np.float32).reshape((nky+1, nkx+1))
    Ey = np.frombuffer(f.read(arrSize * 4), dtype=np.float32).reshape((nky+1, nkx+1))

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

# プロット
x = np.linspace(0, (nx+2*nb+1)*dx, nx+2*nb+1)
y = np.linspace(0, (ny+2*nb+1)*dy, ny+2*nb+1)
X, Y = np.meshgrid(x, y)
plot_field(phi[0:ny+2*nb+1,0:nx+2*nb+1], 'Potential (phi)', 'phi.png')
x = np.linspace(0, (nx+2*nb)*dx, nx+2*nb)
y = np.linspace(0, (ny+2*nb)*dy, ny+2*nb)
X, Y = np.meshgrid(x, y)
plot_field(Ex [0:ny+2*nb,0:nx+2*nb], 'Electric Field (Ex)', 'Ex.png')
plot_field(Ey [0:ny+2*nb,0:nx+2*nb], 'Electric Field (Ey)', 'Ey.png')

# 標準出力
np.set_printoptions(threshold=np.inf)
print(phi_raw[nb+1+(nb+2)*(nx+2*nb+2)])
print(phi[nb+2, nb+1])