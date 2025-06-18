import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter

def getField(filename):
    with open(filename, 'rb') as f:
        # ヘッダ情報の読み込み
        ngx = np.frombuffer(f.read(4), dtype=np.int32)[0]
        ngy = np.frombuffer(f.read(4), dtype=np.int32)[0]
        ngb = np.frombuffer(f.read(4), dtype=np.int32)[0]
        dx = np.frombuffer(f.read(4), dtype=np.float32)[0]
        dy = np.frombuffer(f.read(4), dtype=np.float32)[0]
        # 辞書に格納
        meshinfo = {"ngx":ngx, "ngy":ngy, "ngb":ngb, "dx":dx, "dy":dy}

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
    return meshinfo, fields

# プロット関数
def plotField2d(field, title, figname, type_id, bool_buff):
    # 解析領域の ikmin, jkmin を復元
    if type_id == 4:
        buff = mesh["ngb"]+1
    else:
        buff = mesh["ngb"]
    # 定義点を調整
    if type_id == 0 or type_id == 4:
        x = np.linspace(-buff*mesh["dx"], (mesh["ngx"]+mesh["ngb"])*mesh["dx"], mesh["ngx"]+mesh["ngb"]+buff+1)
        y = np.linspace(-buff*mesh["dy"], (mesh["ngy"]+mesh["ngb"])*mesh["dy"], mesh["ngy"]+mesh["ngb"]+buff+1)
    elif type_id == 1:
        ## shifted to left
        x = np.linspace(-(buff+0.5)*mesh["dx"], (mesh["ngx"]+mesh["ngb"]-0.5)*mesh["dx"], mesh["ngx"]+mesh["ngb"]+buff+2)
        y = np.linspace(-buff*mesh["dy"], (mesh["ngy"]+mesh["ngb"])*mesh["dy"], mesh["ngy"]+mesh["ngb"]+buff+1)
    elif type_id == 2:
        ## shifted to bottom
        x = np.linspace(-buff*mesh["dx"], (mesh["ngx"]+mesh["ngb"])*mesh["dx"], mesh["ngx"]+mesh["ngb"]+buff+1)
        y = np.linspace(-(buff+0.5)*mesh["dy"], (mesh["ngy"]+mesh["ngb"]-0.5)*mesh["dy"], mesh["ngy"]+mesh["ngb"]+buff+2)
    elif type_id == 3:
        ## shifted to left and bottom
        x = np.linspace(-(buff+0.5)*mesh["dx"], (mesh["ngx"]+mesh["ngb"]-0.5)*mesh["dx"], mesh["ngx"]+mesh["ngb"]+buff+2)
        y = np.linspace(-(buff+0.5)*mesh["dy"], (mesh["ngy"]+mesh["ngb"]-0.5)*mesh["dy"], mesh["ngy"]+mesh["ngb"]+buff+2)
    # メッシュ作成
    X, Y = np.meshgrid(x, y)

    # プロット
    fig, ax = plt.subplots()
    if (nx < 50 or ny < 50):
        if bool_buff:
            im = ax.scatter(X[0:ny+1, 0:nx+1], Y[0:ny+1, 0:nx+1], c=field[0:ny+1, 0:nx+1], s=20)
        else:
            im = ax.scatter(X[buff:mesh["ngy"]+buff+1, buff:mesh["ngx"]+buff+1], Y[buff:mesh["ngy"]+buff+1, buff:mesh["ngx"]+buff+1], c=field[buff:mesh["ngy"]+buff+1, buff:mesh["ngx"]+buff+1], s=20)
    else:
        if bool_buff:
            im = ax.pcolormesh(X[0:ny+1, 0:nx+1], Y[0:ny+1, 0:nx+1], field[0:ny+1, 0:nx+1], shading='auto')
        else:
            im = ax.pcolormesh(X[buff:mesh["ngy"]+buff+1, buff:mesh["ngx"]+buff+1], Y[buff:mesh["ngy"]+buff+1, buff:mesh["ngx"]+buff+1], field[buff:mesh["ngy"]+buff+1, buff:mesh["ngx"]+buff+1], shading='auto')
    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.04)
    cbar.set_label('Field value')
    ax.set_xlabel('x (units)')
    ax.set_ylabel('y (units)')
    ax.set_xlim(-(mesh["ngb"]+2)*mesh["dx"], (mesh["ngx"]+mesh["ngb"]+2)*mesh["dx"])
    ax.set_ylim(-(mesh["ngb"]+2)*mesh["dy"], (mesh["ngy"]+mesh["ngb"]+2)*mesh["dy"])
    ax.set_title(title)
    ax.set_aspect('equal', adjustable='box')
    # メッシュの描画
    if (nx < 50 or ny < 50):
        for i in range(mesh["ngx"]+1):
            ax.axvline(i*mesh["dx"], color='gray', linestyle=':', linewidth=0.5)
        for j in range(mesh["ngy"]+1):
            ax.axhline(j*mesh["dy"], color='gray', linestyle=':', linewidth=0.5)
    fig.tight_layout()
    # fig.savefig(figname)
    # plt.close(fig)
    plt.show()

# プロット関数
def plotField1dx(field, title, figname, type_id, j):
    # 解析領域の ikmin, jkmin を復元
    if type_id == 4:
        buff = mesh["ngb"]+2
    else:
        buff = mesh["ngb"]
    # 定義点を調整
    if type_id == 0 or type_id == 4:
        x = np.linspace(-buff*mesh["dx"], (mesh["ngx"]+mesh["ngb"])*mesh["dx"], mesh["ngx"]+mesh["ngb"]+buff+1)
    elif type_id == 1:
        ## shifted to left
        x = np.linspace(-(buff+0.5)*mesh["dx"], (mesh["ngx"]+mesh["ngb"]-0.5)*mesh["dx"], mesh["ngx"]+mesh["ngb"]+buff+2)
    elif type_id == 2:
        ## shifted to bottom
        x = np.linspace(-buff*mesh["dx"], (mesh["ngx"]+mesh["ngb"])*mesh["dx"], mesh["ngx"]+mesh["ngb"]+buff+1)
    elif type_id == 3:
        ## shifted to left and bottom
        x = np.linspace(-(buff+0.5)*mesh["dx"], (mesh["ngx"]+mesh["ngb"]-0.5)*mesh["dx"], mesh["ngx"]+mesh["ngb"]+buff+2)

    # プロット
    fig, ax = plt.subplots()
    ax.plot(x[0:nx+1], field[j, 0:nx+1], label = "j = "+str(j))
    ax.plot(x[0:nx+1], field[j+1, 0:nx+1], label = "j = "+str(j+1))
    ax.plot(x[0:nx+1], field[j+2, 0:nx+1], label = "j = "+str(j+2))
    ax.legend(loc=1)
    ax.set_xlabel('x (units)')
    ax.set_ylabel('y (units)')
    ax.set_xlim(-(mesh["ngb"]+2)*mesh["dx"], (mesh["ngx"]+mesh["ngb"]+2)*mesh["dx"])
    ax.set_title(title)
    # メッシュの描画
    ax.axvline(0, color='gray', linestyle=':', linewidth=0.5)
    ax.axvline(mesh["ngx"]*mesh["dx"], color='gray', linestyle=':', linewidth=0.5)
    fig.tight_layout()
    # fig.savefig(figname)
    # plt.close(fig)
    plt.show()

def save_field_video(
    bin_pattern="bin/long_Moments_ion_Xe1_{cycle:08}.bin",
    name="ion_Xe1_n",
    type_id=0,
    dtoutput=10000,
    Sframe=1,
    Eframe=2,
    outfile="ion_Xe1_n.mp4",
    fps=5,
    dpi=150
):
    """
    bin_pattern : 読み込みファイル名のフォーマット文字列（cycle を受け取る）
    name        : fields 内で一致させるフィールド名
    type_id     : 同じく型 ID
    dtoutput    : ファイル名の cycle 間隔
    Sframe      : 開始フレーム
    Eframe      : 終了フレーム
    outfile     : 出力動画ファイル名
    fps         : 動画のフレームレート
    dpi         : 出力解像度
    """
    # 最初のフレームでメッシュとデータを取得し、範囲とカラー軸を決める
    mesh0, fields0 = getField(bin_pattern.format(cycle=Sframe*dtoutput))
    arr0 = next(arr for nm, tid, arr in fields0 if nm == name and tid == type_id)

    ngb, dx, dy = mesh0["ngb"], mesh0["dx"], mesh0["dy"]
    buff = ngb + 1 if type_id == 4 else ngb
    nrows, ncols = arr0.shape

    # ピクセル中心座標をそのまま格子間隔 dx,dy として仮定
    xmin, xmax = -buff * dx, -buff * dx + dx * (ncols - 1)
    ymin, ymax = -buff * dy, -buff * dy + dy * (nrows - 1)

    vmin, vmax = arr0.min(), arr0.max()

    fig, ax = plt.subplots()
    im = ax.imshow(
        arr0,
        origin="lower",
        extent=(xmin, xmax, ymin, ymax),
        vmin=vmin,
        vmax=vmax,
        interpolation="nearest",
    )
    ax.set_title(name)
    fig.colorbar(im, ax=ax, label=name)

    writer = FFMpegWriter(fps=fps, metadata={"title": name})
    with writer.saving(fig, outfile, dpi=dpi):
        for i in range(Sframe, Eframe):
            cycle = i * dtoutput
            mesh, fields = getField(bin_pattern.format(cycle=cycle))
            arr = next(arr for nm, tid, arr in fields if nm == name and tid == type_id)

            im.set_data(arr)
            # 必要なら軸やタイトルを更新
            ax.set_title(f"{name}  cycle={cycle}")
            writer.grab_frame()

    plt.close(fig)
    print(f"Saved video to {outfile}")

if __name__ == '__main__':

    cycle = 20
    mesh, fields = getField(f"bin/org_EMField_{cycle:08}.bin")
    # mesh, fields = getField(f"bin/long_Moments_ion_Xe1_{cycle:08}.bin")

    print(f"Mesh info: ngx={mesh["ngx"]}, ngy={mesh["ngy"]}, ngb={mesh["ngb"]}, dx={mesh["dx"]}, dy={mesh["dy"]}")
    nx = mesh["ngx"]+2*mesh["ngb"]
    ny = mesh["ngx"]+2*mesh["ngb"]

    # プロット
    for name, type_id, arr in fields:
        print(f"{name}: type_id={type_id}, shape={arr.shape}, min={arr.min()}, max={arr.max()}")
        #plotField2d(arr, name, f"fig/{name}_wBuff.png" , type_id , True)
        plotField2d(arr, name, f"fig/{name}.png" , type_id , False)
        plotField1dx(arr, name, f"fig/{name}_1d_min.png" , type_id , 2)

    # 動画保存
    # save_field_video("bin/long_Moments_ion_Xe1_{cycle:08}.bin", "ion_Xe1_n", 0, 10000, 1, 47, "long_ion_Xe1_n.mp4")