import numpy as np
import matplotlib.pyplot as plt

# プロット関数
def plotField2d(field, title, figname, type_id, bool_buff):
    # 解析領域の ikmin, jkmin を復元
    if type_id == 4:
        buff = ngb+1
    else:
        buff = ngb
    # 定義点を調整
    if type_id == 0 or type_id == 4:
        x = np.linspace(-buff*dx, (ngx+ngb)*dx, ngx+ngb+buff+1)
        y = np.linspace(-buff*dy, (ngy+ngb)*dy, ngy+ngb+buff+1)
    elif type_id == 1:
        ## shifted to left
        x = np.linspace(-(buff+0.5)*dx, (ngx+ngb-0.5)*dx, ngx+ngb+buff+2)
        y = np.linspace(-buff*dy, (ngy+ngb)*dy, ngy+ngb+buff+1)
    elif type_id == 2:
        ## shifted to bottom
        x = np.linspace(-buff*dx, (ngx+ngb)*dx, ngx+ngb+buff+1)
        y = np.linspace(-(buff+0.5)*dy, (ngy+ngb-0.5)*dy, ngy+ngb+buff+2)
    elif type_id == 3:
        ## shifted to left and bottom
        x = np.linspace(-(buff+0.5)*dx, (ngx+ngb-0.5)*dx, ngx+ngb+buff+2)
        y = np.linspace(-(buff+0.5)*dy, (ngy+ngb-0.5)*dy, ngy+ngb+buff+2)
    # メッシュ作成
    X, Y = np.meshgrid(x, y)

    # プロット
    fig, ax = plt.subplots()
    if (nx < 50 or ny < 50):
        if bool_buff:
            im = ax.scatter(X[0:ny+1, 0:nx+1], Y[0:ny+1, 0:nx+1], c=field[0:ny+1, 0:nx+1], s=20)
        else:
            im = ax.scatter(X[buff:ngy+buff+1, buff:ngx+buff+1], Y[buff:ngy+buff+1, buff:ngx+buff+1], c=field[buff:ngy+buff+1, buff:ngx+buff+1], s=20)
    else:
        if bool_buff:
            im = ax.pcolormesh(X[0:ny+1, 0:nx+1], Y[0:ny+1, 0:nx+1], field[0:ny+1, 0:nx+1], shading='auto')
        else:
            im = ax.pcolormesh(X[buff:ngy+buff+1, buff:ngx+buff+1], Y[buff:ngy+buff+1, buff:ngx+buff+1], field[buff:ngy+buff+1, buff:ngx+buff+1], shading='auto')
    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.04)
    cbar.set_label('Field value')
    ax.set_xlabel('x (units)')
    ax.set_ylabel('y (units)')
    ax.set_xlim(-(ngb+2)*dx, (ngx+ngb+2)*dx)
    ax.set_ylim(-(ngb+2)*dy, (ngy+ngb+2)*dy)
    ax.set_title(title)
    ax.set_aspect('equal', adjustable='box')
    # メッシュの描画
    if (nx < 50 or ny < 50):
        for i in range(ngx+1):
            ax.axvline(i*dx, color='gray', linestyle=':', linewidth=0.5)
        for j in range(ngy+1):
            ax.axhline(j*dy, color='gray', linestyle=':', linewidth=0.5)
    fig.tight_layout()
    # fig.savefig(figname)
    # plt.close(fig)
    plt.show()

# プロット関数
def plotField1dx(field, title, figname, type_id, j):
    # 解析領域の ikmin, jkmin を復元
    if type_id == 4:
        buff = ngb+2
    else:
        buff = ngb
    # 定義点を調整
    if type_id == 0 or type_id == 4:
        x = np.linspace(-buff*dx, (ngx+ngb)*dx, ngx+ngb+buff+1)
    elif type_id == 1:
        ## shifted to left
        x = np.linspace(-(buff+0.5)*dx, (ngx+ngb-0.5)*dx, ngx+ngb+buff+2)
    elif type_id == 2:
        ## shifted to bottom
        x = np.linspace(-buff*dx, (ngx+ngb)*dx, ngx+ngb+buff+1)
    elif type_id == 3:
        ## shifted to left and bottom
        x = np.linspace(-(buff+0.5)*dx, (ngx+ngb-0.5)*dx, ngx+ngb+buff+2)

    # プロット
    fig, ax = plt.subplots()
    ax.plot(x[0:nx+1], field[j, 0:nx+1], label = "j = "+str(j))
    ax.plot(x[0:nx+1], field[j+1, 0:nx+1], label = "j = "+str(j+1))
    ax.plot(x[0:nx+1], field[j+2, 0:nx+1], label = "j = "+str(j+2))
    ax.legend(loc=1)
    ax.set_xlabel('x (units)')
    ax.set_ylabel('y (units)')
    ax.set_xlim(-(ngb+2)*dx, (ngx+ngb+2)*dx)
    ax.set_title(title)
    # メッシュの描画
    ax.axvline(0, color='gray', linestyle=':', linewidth=0.5)
    ax.axvline(ngx*dx, color='gray', linestyle=':', linewidth=0.5)
    fig.tight_layout()
    # fig.savefig(figname)
    # plt.close(fig)
    plt.show()

if __name__ == '__main__':

    cycle = 100
    filename = f"bin/atomic_EMField_{cycle:08}.bin"
    #filename = f"bin/org_Moments_electron_{cycle:08}.bin"
    #filename = f"bin/org_Moments_ion_Xe1_{cycle:08}.bin"

    with open(filename, 'rb') as f:
        # ヘッダ情報の読み込み
        ngx = np.frombuffer(f.read(4), dtype=np.int32)[0]
        ngy = np.frombuffer(f.read(4), dtype=np.int32)[0]
        ngb = np.frombuffer(f.read(4), dtype=np.int32)[0]
        dx = np.frombuffer(f.read(4), dtype=np.float32)[0]
        dy = np.frombuffer(f.read(4), dtype=np.float32)[0]

        fields = []
        while True:
            # --- 1. 変数名（固定長32バイト）
            name_bytes = f.read(32)
            if not name_bytes:
                break  # EOF
            # 終端のNULL文字を削除して文字列として格
            name = name_bytes.decode("utf-8").rstrip('\x00')

            # --- 2. タイプID（int）
            type_id = int.from_bytes(f.read(4), byteorder="little")

            # 配列サイズの決定
            nx = ngx+2*ngb
            ny = ngy+2*ngb
            if type_id == 1 or type_id ==3:
                nx +=1
            if type_id == 2 or type_id ==3:
                ny +=1
            if type_id == 4:
                nx +=2
                ny +=2
            
            # --- 3. 配列
            arrSize = (nx+1)*(ny+1)
            arr = np.frombuffer(f.read(arrSize * 4), dtype=np.float32).reshape((ny+1, nx+1))

            # タプルに格納
            fields.append((name, type_id, arr))

    print(f"Mesh info: ngx={ngx}, ngy={ngy}, ngb={ngb}, dx={dx}, dy={dy}")

    # プロット
    for name, type_id, arr in fields:
        if name == "rho": rho = arr
        if name == "arho": arho = arr
        print(f"{name}: type_id={type_id}, shape={arr.shape}, min={arr.min()}, max={arr.max()}")
        plotField2d(arr, name, f"fig/{name}_wBuff.png" , type_id , True)
        plotField2d(arr, name, f"fig/{name}.png" , type_id , False)
        plotField1dx(arr, name, f"fig/{name}_1d_min.png" , type_id , 2)

    diff = rho-arho
    print(f"diff: min={diff.min()}, max={diff.max()}")
    plotField2d(diff, "diff", f"fig/diff.png", 0, False)