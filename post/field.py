import numpy as np
import matplotlib.pyplot as plt

# プロット関数
def plot_field(field, title, figname, type_id, bool_buff):
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
        x = np.linspace(-(buff+0.5)*dx, (ngx+ngb-0.5)*dx, ngx+ngb+buff+1)
        y = np.linspace(-buff*dy, (ngy+ngb)*dy, ngy+ngb+buff+1)
    elif type_id == 2:
        ## shifted to bottom
        x = np.linspace(-buff*dx, (ngx+ngb)*dx, ngx+ngb+buff+1)
        y = np.linspace(-(buff+0.5)*dy, (ngy+ngb-0.5)*dy, ngy+ngb+buff+1)
    elif type_id == 3:
        ## shifted to left and bottom
        x = np.linspace(-(buff+0.5)*dx, (ngx+ngb-0.5)*dx, ngx+ngb+buff+1)
        y = np.linspace(-(buff+0.5)*dy, (ngy+ngb-0.5)*dy, ngy+ngb+buff+1)
    # メッシュ作成
    X, Y = np.meshgrid(x, y)

    # プロット
    fig, ax = plt.subplots()
    if bool_buff:
        im = ax.scatter(X[0:nx+buff+ngb+1, 0:ny+buff+ngb+1], Y[0:nx+buff+ngb+1, 0:ny+buff+ngb+1], c=field[0:nx+buff+ngb+1, 0:ny+buff+ngb+1], s=20)
    else:
        im = ax.scatter(X[buff:nx+buff+1, buff:ny+buff+1], Y[buff:nx+buff+1, buff:ny+buff+1], c=field[buff:nx+buff+1, buff:ny+buff+1], s=20)
    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.04)
    cbar.set_label('Field value')
    ax.set_xlabel('x (units)')
    ax.set_ylabel('y (units)')
    ax.set_xlim(-(ngb+2)*dx, (ngx+ngb+2)*dx)
    ax.set_ylim(-(ngb+2)*dy, (ngy+ngb+2)*dy)
    ax.set_title(title)
    ax.set_aspect('equal', adjustable='box')  # 軸のスケールを等しくする
    # メッシュの描画
    for i in range(nx+1):
        ax.axvline(i*dx, color='gray', linestyle=':', linewidth=0.5)
    for j in range(ny+1):
        ax.axhline(j*dy, color='gray', linestyle=':', linewidth=0.5)
    fig.tight_layout()
    fig.savefig(figname)
    plt.close(fig)

if __name__ == '__main__':

    cycle = 20
    filename = f"bin/field_{cycle:08}.bin"

    with open(filename, 'rb') as f:
        # ヘッダ情報の読み込み
        ngx = np.frombuffer(f.read(4), dtype=np.int32)[0]
        ngy = np.frombuffer(f.read(4), dtype=np.int32)[0]
        ngb = np.frombuffer(f.read(4), dtype=np.int32)[0]
        dx = np.frombuffer(f.read(4), dtype=np.float32)[0]
        dy = np.frombuffer(f.read(4), dtype=np.float32)[0]

        print(f"Mesh info: ngx={ngx}, ngy={ngy}, ngb={ngb}, dx={dx}, dy={dy}")

        fields = []
        while True:
            # --- 1. 変数名（固定長32バイト）
            name_bytes = f.read(32)
            if not name_bytes:
                break  # EOF
            # 終端のNULL文字を削除して文字列として格納
            name = name_bytes.decode("utf-8").rstrip('\x00')

            # --- 2. タイプID（int）
            type_id = int.from_bytes(f.read(4), byteorder="little")

            # 配列サイズの決定
            if type_id == 4:
                nx = ngx+2*ngb+1
                ny = ngy+2*ngb+1
            else:
                nx = ngx+2*ngb
                ny = ngy+2*ngb

            # --- 3. 配列
            arrSize = (nx+1)*(ny+1)
            arr = np.frombuffer(f.read(arrSize * 4), dtype=np.float32).reshape((ny+1, nx+1))

            # タプルに格納
            fields.append((name, type_id, arr))

    # プロット
    for name, type_id, arr in fields:
        print(f"{name} (type_id={type_id}) shape={arr.shape}")
        plot_field(arr, name, f"fig/{name}.png" , type_id , False)