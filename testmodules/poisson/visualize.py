import numpy as np
import matplotlib.pyplot as plt

filename = f"result.bin"

with open(filename, 'rb') as f:
    # グリッド情報の読み込み
    nx = np.frombuffer(f.read(4), dtype=np.int32)[0]
    ny = np.frombuffer(f.read(4), dtype=np.int32)[0]
    nb = np.frombuffer(f.read(4), dtype=np.int32)[0]
    dx = np.frombuffer(f.read(4), dtype=np.float32)[0]
    dy = np.frombuffer(f.read(4), dtype=np.float32)[0]
    # phi の読み込み 
    nkx = nx+2*nb+1
    nky = ny+2*nb+1
    arrSize = (nkx+1)*(nky+1)
    phi_raw = np.frombuffer(f.read(arrSize * 4), dtype=np.float32)
    phi = phi_raw.reshape((nky+1, nkx+1))
    # Ex, Ey の読み込み 
    nkx = nx+2*nb
    nky = ny+2*nb
    arrSize = (nkx+1)*(nky+1)
    Ex = np.frombuffer(f.read(arrSize * 4), dtype=np.float32).reshape((nky+1, nkx+1))
    Ey = np.frombuffer(f.read(arrSize * 4), dtype=np.float32).reshape((nky+1, nkx+1))

# プロット用の関数
def plot_field(field, title, figname, flag_shift, bool_buff):
    if flag_shift == 0:
        ## node-center
        buff = nb+1
        x = np.linspace(-buff*dx, (nx+nb)*dx, nx+nb+buff+1)
        y = np.linspace(-buff*dy, (ny+nb)*dy, ny+nb+buff+1)
    elif flag_shift == 1:
        ## shifted to left
        buff = nb
        x = np.linspace(-(buff+0.5)*dx, (nx+nb-0.5)*dx, nx+nb+buff+1)
        y = np.linspace(-buff*dy, (ny+nb)*dy, ny+nb+buff+1)
    elif flag_shift == 2:
        ## shifted to bottom
        buff = nb
        x = np.linspace(-buff*dx, (nx+nb)*dx, nx+nb+buff+1)
        y = np.linspace(-(buff+0.5)*dy, (ny+nb-0.5)*dy, ny+nb+buff+1)
    elif flag_shift == 3:
        ## shifted to left and bottom
        buff = nb
        x = np.linspace(-(buff+0.5)*dx, (nx+nb-0.5)*dx, nx+nb+buff+1)
        y = np.linspace(-(buff+0.5)*dy, (ny+nb-0.5)*dy, ny+nb+buff+1)
    X, Y = np.meshgrid(x, y)

    fig, ax = plt.subplots()
    if bool_buff:
        im = ax.scatter(X[0:nx+buff+nb+1, 0:ny+buff+nb+1], Y[0:nx+buff+nb+1, 0:ny+buff+nb+1], c=field[0:nx+buff+nb+1, 0:ny+buff+nb+1], s=20)
    else:
        im = ax.scatter(X[buff:nx+buff+1, buff:ny+buff+1], Y[buff:nx+buff+1, buff:ny+buff+1], c=field[buff:nx+buff+1, buff:ny+buff+1], s=20)
    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.04)
    cbar.set_label('Field value')
    ax.set_xlabel('x (units)')
    ax.set_ylabel('y (units)')
    ax.set_xlim(-(nb+2)*dx, (nx+nb+2)*dx)
    ax.set_ylim(-(nb+2)*dy, (ny+nb+2)*dy)
    ax.set_title(title)
    ax.set_aspect('equal', adjustable='box')  # 軸のスケールを等しくする
    # 格子線を明示的に描く
    for i in range(nx+1):
        ax.axvline(i*dx, color='gray', linestyle=':', linewidth=0.5)
    for j in range(ny+1):
        ax.axhline(j*dy, color='gray', linestyle=':', linewidth=0.5)
    fig.tight_layout()
    fig.savefig(figname)
    plt.close(fig)

# プロット
plot_field(phi  ,'Potential (phi)'      , 'phi_dirichlet_wo_buff.png' , 0 , False)
plot_field(Ex   ,'Electric Field (Ex)'  , 'Ex_dirichlet_wo_buff.png'  , 1 , False)
plot_field(Ey   ,'Electric Field (Ey)'  , 'Ey_dirichlet_wo_buff.png'  , 2 , False)

# 標準出力
# np.set_printoptions(threshold=np.inf)
# print(phi_raw[nb+1+(nb+2)*(nx+2*nb+2)])
# print(phi[nb+2, nb+1])