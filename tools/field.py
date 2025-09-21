import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter
import sys
import math

# constants
eps0 = 8.85418782e-12 # [F/m]
ec = 1.60217663e-19 # [C]

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
def plotField2d(field, title, figname, type_id, bool_buff, cmin=None, cmax=None):
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

    # カラーバー
    vmin = field.min() if cmin is None else cmin
    vmax = field.max() if cmax is None else cmax
    cmap = 'bwr' if abs(vmin+vmax) < 1e-20 else 'plasma'
    # プロット
    fig, ax = plt.subplots()
    # if (nx < 50 or ny < 50):
    #     if bool_buff:
    #         im = ax.scatter(X[0:ny+1, 0:nx+1], Y[0:ny+1, 0:nx+1], c=field[0:ny+1, 0:nx+1], s=20,marker='s',vmin=vmin,vmax=vmax,cmap=cmap)
    #     else:
    #         im = ax.scatter(X[buff:mesh["ngy"]+buff+1, buff:mesh["ngx"]+buff+1], Y[buff:mesh["ngy"]+buff+1, buff:mesh["ngx"]+buff+1], c=field[buff:mesh["ngy"]+buff+1, buff:mesh["ngx"]+buff+1], s=30,marker='s',vmin=vmin,vmax=vmax,cmap=cmap)
    # else:
    #     if bool_buff:
    #         im = ax.pcolormesh(X[0:ny+1, 0:nx+1], Y[0:ny+1, 0:nx+1], field[0:ny+1, 0:nx+1], shading='auto',vmin=vmin,vmax=vmax,cmap=cmap)
    #     else:
    #         im = ax.pcolormesh(X[buff:mesh["ngy"]+buff+1, buff:mesh["ngx"]+buff+1], Y[buff:mesh["ngy"]+buff+1, buff:mesh["ngx"]+buff+1], field[buff:mesh["ngy"]+buff+1, buff:mesh["ngx"]+buff+1], shading='auto',vmin=vmin,vmax=vmax,cmap=cmap)
    if bool_buff:
        im = ax.pcolormesh(X[0:ny+1, 0:nx+1], Y[0:ny+1, 0:nx+1], field[0:ny+1, 0:nx+1], shading='auto',vmin=vmin,vmax=vmax,cmap=cmap)
    else:
        im = ax.pcolormesh(X[buff:mesh["ngy"]+buff+1, buff:mesh["ngx"]+buff+1], Y[buff:mesh["ngy"]+buff+1, buff:mesh["ngx"]+buff+1], field[buff:mesh["ngy"]+buff+1, buff:mesh["ngx"]+buff+1], shading='auto',vmin=vmin,vmax=vmax,cmap=cmap)
    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.04)
    cbar.set_label('Electric Potential [V]')
    ax.set_xlabel('x [cm]')
    ax.set_ylabel('y [cm]')
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
    ax.plot(x[0:mesh["ngx"]+mesh["ngb"]+buff+1], field[j, 0:mesh["ngx"]+mesh["ngb"]+buff+1], label = "j = "+str(j))
    ax.plot(x[0:mesh["ngx"]+mesh["ngb"]+buff+1], field[j+1, 0:mesh["ngx"]+mesh["ngb"]+buff+1], label = "j = "+str(j+1))
    ax.plot(x[0:mesh["ngx"]+mesh["ngb"]+buff+1], field[j+2, 0:mesh["ngx"]+mesh["ngb"]+buff+1], label = "j = "+str(j+2))
    ax.plot(x[0:mesh["ngx"]+mesh["ngb"]+buff+1], np.mean(field[:, 0:mesh["ngx"]+mesh["ngb"]+buff+1], axis=0), label = "mean")
    ax.legend(loc=1)
    ax.set_xlabel('x [cm]')
    ax.set_ylabel('valq')
    ax.set_xlim(-(mesh["ngb"]+2)*mesh["dx"], (mesh["ngx"]+mesh["ngb"]+2)*mesh["dx"])
    ax.set_title(title)
    # メッシュの描画
    ax.axvline(0, color='gray', linestyle=':', linewidth=0.5)
    ax.axvline(mesh["ngx"]*mesh["dx"], color='gray', linestyle=':', linewidth=0.5)
    ax.axvline(2.4, color='gray', linestyle=':', linewidth=0.5) ## cathodeline
    fig.tight_layout()
    # fig.savefig(figname)
    # plt.close(fig)
    plt.show()

# プロット関数
def plotField1dy(field, title, figname, type_id, i):
    # 解析領域の ikmin, jkmin を復元
    if type_id == 4:
        buff = mesh["ngb"]+2
    else:
        buff = mesh["ngb"]
    # 定義点を調整
    if type_id == 0 or type_id == 4:
        y = np.linspace(-buff*mesh["dy"], (mesh["ngy"]+mesh["ngb"])*mesh["dy"], mesh["ngy"]+mesh["ngb"]+buff+1)
    elif type_id == 1:
        ## shifted to left
        y = np.linspace(-(buff+0.5)*mesh["dy"], (mesh["ngy"]+mesh["ngb"]-0.5)*mesh["dy"], mesh["ngy"]+mesh["ngb"]+buff+2)
    elif type_id == 2:
        ## shifted to bottom
        y = np.linspace(-buff*mesh["dy"], (mesh["ngy"]+mesh["ngb"])*mesh["dy"], mesh["ngy"]+mesh["ngb"]+buff+1)
    elif type_id == 3:
        ## shifted to left and bottom
        y = np.linspace(-(buff+0.5)*mesh["dy"], (mesh["ngy"]+mesh["ngb"]-0.5)*mesh["dy"], mesh["ngy"]+mesh["ngb"]+buff+2)

    # プロット
    fig, ax = plt.subplots()
    ax.plot(y[0:mesh["ngy"]+mesh["ngb"]+buff+1], field[0:mesh["ngy"]+mesh["ngb"]+buff+1,i  ], label = "i = "+str(i))
    ax.plot(y[0:mesh["ngy"]+mesh["ngb"]+buff+1], field[0:mesh["ngy"]+mesh["ngb"]+buff+1,i+int(i*1/10.0)], label = "i = "+str(i+int(i*1/10.0)))
    ax.plot(y[0:mesh["ngy"]+mesh["ngb"]+buff+1], field[0:mesh["ngy"]+mesh["ngb"]+buff+1,i+int(i*2/10.0)], label = "i = "+str(i+int(i*2/10.0)))
    ax.plot(y[0:mesh["ngy"]+mesh["ngb"]+buff+1], field[0:mesh["ngy"]+mesh["ngb"]+buff+1,i+int(i*3/10.0)], label = "i = "+str(i+int(i*3/10.0)))
    # ax.plot(y[0:mesh["ngy"]+mesh["ngb"]+buff+1], np.mean(field[0:mesh["ngx"]+mesh["ngb"]+buff+1,:], axis=1), label = "mean")
    ax.legend(loc=1)
    ax.set_xlabel('y [cm]')
    ax.set_ylabel('val')
    ax.set_xlim(-(mesh["ngb"]+2)*mesh["dy"], (mesh["ngy"]+mesh["ngb"]+2)*mesh["dy"])
    ax.set_title(title)
    # メッシュの描画
    ax.axvline(0, color='gray', linestyle=':', linewidth=0.5)
    ax.axvline(mesh["ngy"]*mesh["dy"], color='gray', linestyle=':', linewidth=0.5)
    fig.tight_layout()
    # fig.savefig(figname)
    # plt.close(fig)
    plt.show()

def save_field_video(
    bin_pattern="bin/large_Moments_electron_{cycle:08}.bin",
    name="electron",
    unit="[-]",
    type_id=0,
    dtoutput=1,
    Sframe=0,
    Eframe=1,
    outfile="ion_Xe1_n.mp4",
    vmin=None,
    vmax=None,
    fps=5,
    dpi=150
):
    """
    bin_pattern : 読み込みファイル名のフォーマット文字列（cycle を受け取る）
    name        : 変数名
    unit        : 単位
    type_id     : 同じく型 ID
    dtoutput    : ファイル名の cycle 間隔
    Sframe      : 開始フレーム
    Eframe      : 終了フレーム
    outfile     : 出力動画ファイル名
    vmin   : カラーバーの最小値（None の場合は最初のフレームから自動設定）
    vmax   : カラーバーの最大値（None の場合は最初のフレームから自動設定）
    fps         : 動画のフレームレート
    dpi         : 出力解像度
    """
    # time evolution per cycle
    dt_sim = 5e-3

    # 最初のフレーム読み込み
    mesh0, fields0 = getField(bin_pattern.format(cycle=Sframe))
    arr0 = next(arr for nm, tid, arr in fields0 if nm == name and tid == type_id)

    # カラーバー範囲を決定
    cvar_vmin = arr0.min() if vmin is None else vmin
    cvar_vmax = arr0.max() if vmax is None else vmax
    # ゼロ中心ならbwr, それ以外ならplasma
    cmap = 'bwr' if abs(cvar_vmin+cvar_vmax) < 1e-20 else 'plasma'

    ngb, dx, dy = mesh0["ngb"], mesh0["dx"], mesh0["dy"]
    buff = ngb + 1 if type_id == 4 else ngb
    nrows, ncols = arr0.shape

    xmin, xmax = -buff * dx, -buff * dx + dx * (ncols - 1)
    ymin, ymax = -buff * dy, -buff * dy + dy * (nrows - 1)

    fig, ax = plt.subplots()
    im = ax.imshow(
        arr0,
        origin="lower",
        extent=(xmin, xmax, ymin, ymax),
        vmin=cvar_vmin,
        vmax=cvar_vmax,
        interpolation="nearest",
        cmap=cmap,
    )
    ax.set_title(name)
    # メッシュ
    if (nx < 50 or ny < 50):
        for i in range(mesh0["ngx"]+1):
            ax.axvline(i*mesh0["dx"], color='gray', linestyle=':', linewidth=0.5)
        for j in range(mesh0["ngy"]+1):
            ax.axhline(j*mesh0["dy"], color='gray', linestyle=':', linewidth=0.5)

    cbar = fig.colorbar(im, ax=ax, label=name+" "+unit, shrink=0.6)

    # writer の作成
    writer = FFMpegWriter(fps=fps, metadata={"title": name})
    with writer.saving(fig, outfile, dpi=dpi):
        for cycle in range(Sframe, Eframe+dtoutput, dtoutput):
            mesh, fields = getField(bin_pattern.format(cycle=cycle))
            # データ更新
            arr = next(arr for nm, tid, arr in fields if nm == name and tid == type_id)
            im.set_data(arr)
            # カラーバー更新
            cvar_vmin = arr.min() if vmin is None else vmin
            cvar_vmax = arr.max() if vmax is None else vmax
            im.set_clim(cvar_vmin,cvar_vmax)
            # ax.set_title(f'{name}: time = {int(cycle*dt_sim)} [ns]')
            ax.set_title(f'{name}: time = {cycle*dt_sim:.3f} [ns]')
            writer.grab_frame()

    plt.close(fig)
    print(f"Saved video to {outfile}")

def save_field_video_1d(
    bin_pattern="bin/large_Moments_electron_{cycle:08}.bin",
    name="electron_n",
    unit="[-]",
    type_id=0,
    dtoutput=1,
    Sframe=0,
    Eframe=1,
    j=0,
    outfile="electron_n_1d.mp4",
    vmin=None,
    vmax=None,
    fps=10,
    dpi=150,
    show_grid=True,
):
    """
    2Dフィールドの1D断面を動画化
    - bin_pattern: 読み込みファイル名フォーマット（{cycle} 必須）
    - name, type_id: 変数名とタイプID
    - dtoutput: cycle 間隔
    - Sframe, Eframe: 開始と終了 cycle（両端含む）
    - j: 抽出する行の下端インデックス（j, j+1, j+2 を描画）
    - vmin, vmax: y軸の固定範囲（None で自動）
    """
    # 時間刻み（1cycle 当たりの実時間 [ns]）。環境に応じて調整
    dt_sim = 5e-3  # [ns/cycle]

    # 最初のフレーム
    mesh0, fields0 = getField(bin_pattern.format(cycle=Sframe))
    arr0 = next(arr for nm, tid, arr in fields0 if nm == name and tid == type_id)

    # バッファ厚の扱い（plotField1dx と同じ）
    if type_id == 4:
        buff = mesh0["ngb"] + 2
    else:
        buff = mesh0["ngb"]

    # x 定義点（type_id に応じたシフト）
    if type_id in (0, 4):
        x = np.linspace(-buff*mesh0["dx"], (mesh0["ngx"]+mesh0["ngb"])*mesh0["dx"],
                        mesh0["ngx"]+mesh0["ngb"]+buff+1)
    elif type_id == 1:
        x = np.linspace(-(buff+0.5)*mesh0["dx"], (mesh0["ngx"]+mesh0["ngb"]-0.5)*mesh0["dx"],
                        mesh0["ngx"]+mesh0["ngb"]+buff+2)
    elif type_id == 2:
        x = np.linspace(-buff*mesh0["dx"], (mesh0["ngx"]+mesh0["ngb"])*mesh0["dx"],
                        mesh0["ngx"]+mesh0["ngb"]+buff+1)
    elif type_id == 3:
        x = np.linspace(-(buff+0.5)*mesh0["dx"], (mesh0["ngx"]+mesh0["ngb"]-0.5)*mesh0["dx"],
                        mesh0["ngx"]+mesh0["ngb"]+buff+2)
    nxvis = mesh0["ngx"] + mesh0["ngb"] + buff + 1
    xs = x[:nxvis]

    # 初期 y データ
    def slice_rows(A, j0):
        return (A[j0, :nxvis],
                A[j0+1, :nxvis] if j0+1 < A.shape[0] else A[j0, :nxvis],
                A[j0+2, :nxvis] if j0+2 < A.shape[0] else A[j0, :nxvis],
                np.mean(A[:, :nxvis], axis=0))

    y0, y1, y2, ymean = slice_rows(arr0, j)

    # y 範囲
    if vmin is None or vmax is None:
        ymin = np.min([y0.min(), y1.min(), y2.min(), ymean.min()])
        ymax = np.max([y0.max(), y1.max(), y2.max(), ymean.max()])
        pad = 0.05*(ymax - ymin + 1e-30)
        auto_ymin, auto_ymax = ymin - pad, ymax + pad
    ymin = auto_ymin if vmin is None else vmin
    ymax = auto_ymax if vmax is None else vmax

    # Figure 準備
    fig, ax = plt.subplots()
    ln0, = ax.plot(xs, y0, label=f"j={j}")
    ln1, = ax.plot(xs, y1, label=f"j={j+1}")
    ln2, = ax.plot(xs, y2, label=f"j={j+2}")
    lnM, = ax.plot(xs, ymean, label="mean")
    ax.legend(loc=1)
    ax.set_xlabel("x (units)")
    ax.set_ylabel(f"{name} {unit}")
    ax.set_xlim(-(mesh0["ngb"]+2)*mesh0["dx"], (mesh0["ngx"]+mesh0["ngb"]+2)*mesh0["dx"])
    ax.set_ylim(ymin, ymax)
    ax.set_title(f'{name}: time = {Sframe*dt_sim:.3f} [ns]')

    # 参考ライン
    if show_grid:
        ax.axvline(0, color="gray", linestyle=":", linewidth=0.5)
        ax.axvline(mesh0["ngx"]*mesh0["dx"], color="gray", linestyle=":", linewidth=0.5)
        ax.axvline(2.4, color="gray", linestyle=":", linewidth=0.5)  # 任意の目安線
    fig.tight_layout()

    writer = FFMpegWriter(fps=fps, metadata={"title": f"{name} 1D"})
    with writer.saving(fig, outfile, dpi=dpi):
        for cycle in range(Sframe, Eframe+dtoutput, dtoutput):
            mesh, fields = getField(bin_pattern.format(cycle=cycle))
            arr = next(arr for nm, tid, arr in fields if nm == name and tid == type_id)
            y0, y1, y2, ymean = slice_rows(arr, j)
            ln0.set_ydata(y0)
            ln1.set_ydata(y1)
            ln2.set_ydata(y2)
            lnM.set_ydata(ymean)
            # y 範囲
            if vmin is None or vmax is None:
                ymin = np.min([y0.min(), y1.min(), y2.min(), ymean.min()])
                ymax = np.max([y0.max(), y1.max(), y2.max(), ymean.max()])
                pad = 0.05*(ymax - ymin + 1e-30)
                auto_ymin, auto_ymax = ymin - pad, ymax + pad
            ymin = auto_ymin if vmin is None else vmin
            ymax = auto_ymax if vmax is None else vmax
            ax.set_ylim(ymin, ymax)
            ax.set_title(f'{name}: time = {cycle*dt_sim:.3f} [ns]')
            writer.grab_frame()

    plt.close(fig)
    print(f"Saved 1D video to {outfile}")

if __name__ == '__main__':
    args = sys.argv

    # カラーバー調整辞書(get で取り出した場合、無指定の変数は None を渡す)
    dict_cmin = {}
    dict_cmax = {}
    # dict_cmin = {"electron_n":5e9-1e3}
    # dict_cmax = {"electron_n":5e9+1e3}
    # dict_cmin = {"rho":-20, "phi":-100, "Ex":-1e5, "Ey":-1e5, "electron_uy":-1e9}
    # dict_cmax = {"rho":20,  "phi":700,  "Ex":1e5,  "Ey":1e5,  "electron_uy":1e9, "electron_T":60, "DebyeRes":0.2}
    # dict_cmin = {"rho":-1}
    # dict_cmax = {"rho":1}

    ## projectname
    pname = args[1]
    if args[2].isdigit():
        cycle = int(args[2])
    else:
        cycle = None
    mesh, fields = getField("bin/"+pname+f"_EMField_{cycle:08}.bin")

    print(f"Mesh info: ngx={mesh["ngx"]}, ngy={mesh["ngy"]}, ngb={mesh["ngb"]}, dx={mesh["dx"]}, dy={mesh["dy"]}")
    nx = mesh["ngx"]+2*mesh["ngb"]
    ny = mesh["ngy"]+2*mesh["ngb"]

    ## 1dplot position
    ypos_1dxplot = mesh["ngb"]+int(mesh["ngy"]/2)
    xpos_1dyplot = mesh["ngb"]+int(mesh["ngx"]/4*3)

    # plot EMField
    ave_phi = np.zeros(nx+3,ny+3)
    ave_Ex = np.zeros(nx+2,ny+1)
    for name, type_id, arr in fields:
        print(f"{name}: type_id={type_id}, shape={arr.shape}, min={arr.min()}, max={arr.max()}")
        # for idx, val in np.ndenumerate(arr[1,:]):
        #     print(idx[0],",",val)
        # plotField2d(arr, name, f"fig/{name}_wBuff.png" , type_id , True,dict_cmin.get(name),dict_cmax.get(name))
        plotField2d(arr, name, f"fig/{name}.png" , type_id, True,dict_cmin.get(name),dict_cmax.get(name))
        plotField1dx(arr, name, f"fig/{name}_1d_min.png" , type_id , ypos_1dxplot)
        plotField1dy(arr, name, f"fig/{name}_1d_min.png" , type_id , xpos_1dyplot)
    #     if('phi' in name): ave_phi += arr
    #     if('Ex' in name): ave_Ex += arr
    # ave_phi /= fields.shape[0]
    # ave_Ex /= fields.shape[0]
    # plotField1dx(ave_phi, name, f"fig/{name}_1d_min.png" , 4 , ypos_1dxplot)
    # plotField1dx(ave_Ex, name, f"fig/{name}_1d_min.png" , 1 , ypos_1dxplot)
    
    # # electron
    mesh, fields = getField("bin/"+pname+f"_Moments_electron_{cycle:08}.bin")
    Pe = np.zeros((ny+1,nx+1))
    for name, type_id, arr in fields:
        ## 圧力テンソルはスキップ
        if('electron_P' not in name):
            print(f"{name}: type_id={type_id}, shape={arr.shape}, min={arr.min()}, max={arr.max()}")
            # for idx, val in np.ndenumerate(arr[3,:]):
            #     print(idx[0],",",val)
            # plotField2d(arr, name, f"fig/{name}_wBuff.png" , type_id , True)
            plotField2d(arr, name, f"fig/{name}.png" , type_id, True,dict_cmin.get(name),dict_cmax.get(name))
            plotField1dx(arr, name, f"fig/{name}_1d_min.png" , type_id , ypos_1dxplot)
            plotField1dy(arr, name, f"fig/{name}_1d_min.png" , type_id , xpos_1dyplot)
        ## 温度計算
        if(name == "electron_n"): ne = arr # 1/cm3
        if(name == "electron_Pxx"): Pe += arr
        if(name == "electron_Pyy"): Pe += arr
        if(name == "electron_Pzz"): Pe += arr
    # temperature
    name = "electron_T"
    kb = 1.380649e-23 # [J/K]
    KtoeV = 1.0/11604
    type_id = 0
    arr = np.where(ne > 0, Pe/(3*kb*ne*1e6)*KtoeV, 0) ## ne > 0 の要素だけ計算
    print(f"{name}: min={arr.min()}, max={arr.max()}")
    plotField2d(arr, name, f"fig/{name}.png" , type_id, True,dict_cmin.get(name),dict_cmax.get(name))
    plotField1dx(arr, name, f"fig/{name}_1d_min.png" , type_id , ypos_1dxplot)
    # Debye length
    dx = 0.005e-2 # [m]
    name = "DebyeRes"
    type_id = 0
    arr = np.zeros_like(Pe)
    np.sqrt(eps0*Pe/(3*(ne*1e6)**2*ec**2)/dx, out=arr, where=(Pe > 0)&(ne > 0)) ## where を満たす要素だけ sqrt する
    print(f"{name}: min={Pe.min()}, max={Pe.max()}")
    print(f"{name}: min={arr.min()}, max={arr.max()}")
    plotField2d(arr, name, f"fig/{name}.png" , type_id, True,dict_cmin.get(name),dict_cmax.get(name))
    plotField1dx(arr, name, f"fig/{name}_1d_min.png" , type_id , ypos_1dxplot)
    plotField1dy(arr, name, f"fig/{name}_1d_min.png" , type_id , xpos_1dyplot)

    # # ion
    mesh, fields = getField("bin/"+pname+f"_Moments_ion_Xe1_{cycle:08}.bin")
    for name, type_id, arr in fields:
        if('ion_Xe1_P' not in name):
            print(f"{name}: type_id={type_id}, shape={arr.shape}, min={arr.min()}, max={arr.max()}")
            # plotField2d(arr, name, f"fig/{name}_wBuff.png" , type_id , True)
            plotField2d(arr, name, f"fig/{name}.png" , type_id, True,dict_cmin.get(name),dict_cmax.get(name))
            plotField1dx(arr, name, f"fig/{name}_1d_min.png" , type_id , ypos_1dxplot)
            plotField1dy(arr, name, f"fig/{name}_1d_min.png" , type_id , xpos_1dyplot)

    # 動画保存(開始フレーム、最終フレームを数字で受け取り動画保存)
    if (len(args)==5):
        start = int(args[3])
        dt = int(args[4])
        end = int(args[2])
        ## 2d plot
        # save_field_video("bin/"+pname+"_Moments_electron_{cycle:08}.bin" ,"electron_n"   ,"[1/cm3]"  ,0,dt,start,end ,"fig/"+pname+"_from"+args[3]+"_to"+args[2]+"_electron_n.mp4"  ,vmin=dict_cmin.get("electron_n")  ,vmax=dict_cmax.get("electron_n") ,fps=15)
        # save_field_video("bin/"+pname+"_Moments_electron_{cycle:08}.bin" ,"electron_ux"  ,"[cm/s]"   ,0,dt,start,end ,"fig/"+pname+"_from"+args[3]+"_to"+args[2]+"_electron_ux.mp4" ,vmin=dict_cmin.get("electron_ux") ,vmax=dict_cmax.get("electron_ux"),fps=15)
        # save_field_video("bin/"+pname+"_Moments_electron_{cycle:08}.bin" ,"electron_uy"  ,"[cm/s]"   ,0,dt,start,end ,"fig/"+pname+"_from"+args[3]+"_to"+args[2]+"_electron_uy.mp4" ,vmin=dict_cmin.get("electron_uy") ,vmax=dict_cmax.get("electron_uy"),fps=15)
        # save_field_video("bin/"+pname+"_Moments_ion_Xe1_{cycle:08}.bin"  ,"ion_Xe1_n"    ,"[1/cm3]"  ,0,dt,start,end ,"fig/"+pname+"_from"+args[3]+"_to"+args[2]+"_ion_Xe1_n.mp4"   ,vmin=dict_cmin.get("ion_Xe1_n")   ,vmax=dict_cmax.get("ion_Xe1_n")  ,fps=15)
        # save_field_video("bin/"+pname+"_Moments_ion_Xe1_{cycle:08}.bin"  ,"ion_Xe1_ux"   ,"[cm/s]"   ,0,dt,start,end ,"fig/"+pname+"_from"+args[3]+"_to"+args[2]+"_ion_Xe1_ux.mp4"  ,vmin=dict_cmin.get("ion_Xe1_ux")  ,vmax=dict_cmax.get("ion_Xe1_ux") ,fps=15)
        # save_field_video("bin/"+pname+"_Moments_ion_Xe1_{cycle:08}.bin"  ,"ion_Xe1_uy"   ,"[cm/s]"   ,0,dt,start,end ,"fig/"+pname+"_from"+args[3]+"_to"+args[2]+"_ion_Xe1_uy.mp4"  ,vmin=dict_cmin.get("ion_Xe1_uy")  ,vmax=dict_cmax.get("ion_Xe1_uy") ,fps=15)
        save_field_video("bin/"+pname+"_EMField_{cycle:08}.bin"          ,"rho"          ,"[esu/m3]" ,0,dt,start,end ,"fig/"+pname+"_from"+args[3]+"_to"+args[2]+"_rho.mp4"         ,vmin=dict_cmin.get("rho")         ,vmax=dict_cmax.get("rho")        ,fps=30)
        save_field_video("bin/"+pname+"_EMField_{cycle:08}.bin"          ,"Ex"           ,"[V/m]"    ,1,dt,start,end ,"fig/"+pname+"_from"+args[3]+"_to"+args[2]+"_Ex.mp4"          ,vmin=dict_cmin.get("Ex")          ,vmax=dict_cmax.get("Ex")         ,fps=30)
        save_field_video("bin/"+pname+"_EMField_{cycle:08}.bin"          ,"Ey"           ,"[V/m]"    ,2,dt,start,end ,"fig/"+pname+"_from"+args[3]+"_to"+args[2]+"_Ey.mp4"          ,vmin=dict_cmin.get("Ey")          ,vmax=dict_cmax.get("Ey")         ,fps=30)
        ## 1d plot
        save_field_video_1d("bin/"+pname+"_Moments_electron_{cycle:08}.bin" ,"electron_n"   ,"[1/cm3]"  ,0,dt,start,end, ypos_1dxplot,"fig/"+pname+"_from"+args[3]+"_to"+args[2]+"_electron_n_1d.mp4"  ,vmin=dict_cmin.get("electron_n")  ,vmax=dict_cmax.get("electron_n") ,fps=30)
        save_field_video_1d("bin/"+pname+"_Moments_electron_{cycle:08}.bin" ,"electron_ux"  ,"[cm/s]"   ,0,dt,start,end, ypos_1dxplot,"fig/"+pname+"_from"+args[3]+"_to"+args[2]+"_electron_ux_1d.mp4" ,vmin=dict_cmin.get("electron_ux") ,vmax=dict_cmax.get("electron_ux"),fps=30)
        # save_field_video_1d("bin/"+pname+"_Moments_electron_{cycle:08}.bin" ,"electron_uy"  ,"[cm/s]"   ,0,dt,start,end, ypos_1dxplot,"fig/"+pname+"_from"+args[3]+"_to"+args[2]+"_electron_uy_1d.mp4" ,vmin=dict_cmin.get("electron_uy") ,vmax=dict_cmax.get("electron_uy"),fps=30)
        # save_field_video_1d("bin/"+pname+"_Moments_ion_Xe1_{cycle:08}.bin"  ,"ion_Xe1_n"    ,"[1/cm3]"  ,0,dt,start,end, ypos_1dxplot,"fig/"+pname+"_from"+args[3]+"_to"+args[2]+"_ion_Xe1_n_1d.mp4"   ,vmin=dict_cmin.get("ion_Xe1_n")   ,vmax=dict_cmax.get("ion_Xe1_n")  ,fps=30)
        # save_field_video_1d("bin/"+pname+"_Moments_ion_Xe1_{cycle:08}.bin"  ,"ion_Xe1_ux"   ,"[cm/s]"   ,0,dt,start,end, ypos_1dxplot,"fig/"+pname+"_from"+args[3]+"_to"+args[2]+"_ion_Xe1_ux_1d.mp4"  ,vmin=dict_cmin.get("ion_Xe1_ux")  ,vmax=dict_cmax.get("ion_Xe1_ux") ,fps=30)
        # save_field_video_1d("bin/"+pname+"_Moments_ion_Xe1_{cycle:08}.bin"  ,"ion_Xe1_uy"   ,"[cm/s]"   ,0,dt,start,end, ypos_1dxplot,"fig/"+pname+"_from"+args[3]+"_to"+args[2]+"_ion_Xe1_uy_1d.mp4"  ,vmin=dict_cmin.get("ion_Xe1_uy")  ,vmax=dict_cmax.get("ion_Xe1_uy") ,fps=30)
        save_field_video_1d("bin/"+pname+"_EMField_{cycle:08}.bin"          ,"rho"          ,"[esu/m3]" ,0,dt,start,end, ypos_1dxplot,"fig/"+pname+"_from"+args[3]+"_to"+args[2]+"_rho_1d.mp4"         ,vmin=dict_cmin.get("rho")         ,vmax=dict_cmax.get("rho")        ,fps=30)
        save_field_video_1d("bin/"+pname+"_EMField_{cycle:08}.bin"          ,"Ex"           ,"[V/m]"    ,1,dt,start,end, ypos_1dxplot,"fig/"+pname+"_from"+args[3]+"_to"+args[2]+"_Ex_1d.mp4"          ,vmin=dict_cmin.get("Ex")          ,vmax=dict_cmax.get("Ex")         ,fps=30)
        # save_field_video_1d("bin/"+pname+"_EMField_{cycle:08}.bin"          ,"Ey"           ,"[V/m]"    ,2,dt,start,end, ypos_1dxplot,"fig/"+pname+"_from"+args[3]+"_to"+args[2]+"_Ey_1d.mp4"          ,vmin=dict_cmin.get("Ey")          ,vmax=dict_cmax.get("Ey")         ,fps=30)
        save_field_video_1d("bin/"+pname+"_EMField_{cycle:08}.bin"          ,"phi"          ,"[V]"      ,4,dt,start,end, ypos_1dxplot,"fig/"+pname+"_from"+args[3]+"_to"+args[2]+"_phi_1d.mp4"         ,vmin=dict_cmin.get("phi")         ,vmax=dict_cmax.get("phi")        ,fps=30)
