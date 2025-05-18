import struct
import math
import numpy as np
from typing import Sequence

def write_1d_field(
    file_path: str,
    n: int,
    d: float,
    data: Sequence[float]
) -> None:
    """
    EMField.load1dField: が読み込めるバイナリを生成する。

    Parameters
    ----------
    file_path : str
        出力先ファイルパス
    n : int
        配列サイズ n（実データは後で n+1 要素）
    d : float
        ヘッダの浮動小数点値
    data : Sequence[float]
        長さ n+1 のデータシーケンス
    """
    if len(data) != n + 1:
        raise ValueError(f"data length must be n+1 ({n+1}), but got {len(data)}")

    # '<' はリトルエンディアン。'i' は 4 バイト int、'f' は 4 バイト float。
    with open(file_path, 'wb') as f:
        # ヘッダ部分
        f.write(struct.pack('<i', n))
        f.write(struct.pack('<f', d))
        # データ部分
        f.write(struct.pack(f'<{n+1}f', *data))



if __name__ == "__main__":
    dx = 1e-2
    ngx = 250
    ngb = 2
    nx = ngx+2*ngb

    ## set param
    bm_1    = 6e-3
    bm_max  = 10e-3
    bm_2    = 1e-3
    x_1     = ngb*dx
    x_max   = ngb*dx + 0.75
    x_2     = (ngx+ngb)*dx
    sgm1    = 0.625
    sgm2    = 0.625
    # bm_wnum = 8e-3
    # bm_wamp = 0.2e-3
    aa1 = (bm_max-bm_1)/(1-math.exp(-0.5*(x_max-x_1)**2/sgm1**2))
    bb1 = bm_max - aa1
    aa2 = (bm_max-bm_2)/(1-math.exp(-0.5*(x_max-x_2)**2/sgm2**2))
    bb2 = bm_max - aa2
    ## create mesh
    i = np.arange(nx+1)
    bmx = (i + 0.5) * dx
    bs1d = np.where(
        bmx < x_max,
        aa1 * np.exp(-0.5 * (bmx - x_max)**2 / (sgm1**2)) + bb1,
        aa2 * np.exp(-0.5 * (bmx - x_max)**2 / (sgm2**2)) + bb2
    )
    write_1d_field("data/Bz.bin", nx, dx, bs1d)
    ## output zeros
    write_1d_field("data/Bx.bin", nx, dx, np.zeros(nx+1))
    write_1d_field("data/By.bin", nx, dx, np.zeros(nx+1))
    print("done")