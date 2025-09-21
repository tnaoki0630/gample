#!/usr/bin/env python3
import numpy as np
import pandas as pd
import sys
import matplotlib.pyplot as plt
from typing import Optional, Tuple, Dict

# ========= 計算関数 =========
def compute_fft_from_phi(phi: np.ndarray, dt: float) -> Dict[str, np.ndarray]:
    """phi(t) を受け取り rFFT を返す"""
    N = phi.size
    t = np.arange(1, N + 1) * dt  # 2行目→時刻1の仕様
    Phi = np.fft.rfft(phi)
    freq = np.fft.rfftfreq(N, d=dt)
    amp = np.abs(Phi)
    phase = np.angle(Phi)
    return {"t": t, "phi": phi, "freq": freq, "amp": amp, "phase": phase}

# ========= プロット関数 =========
def plot_time_series(
    t: np.ndarray,
    phi: np.ndarray,
    savepath: str = "plot_time_series.png",
    xlim: Optional[Tuple[float, float]] = None,
    ylim: Optional[Tuple[float, float]] = None,
    title: str = "Time Series: phi(t)"
) -> None:
    plt.figure()
    plt.plot(t, phi)
    plt.xlabel("time")
    plt.ylabel("phi")
    plt.title(title)
    if xlim is not None:
        plt.xlim(*xlim)
    if ylim is not None:
        plt.ylim(*ylim)
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(savepath, dpi=200)
    plt.close()

def plot_spectrum(
    freq: np.ndarray,
    amp: np.ndarray,
    savepath: str = "plot_spectrum.png",
    xlim: Optional[Tuple[float, float]] = None,
    ylim: Optional[Tuple[float, float]] = None,
    title: str = "Amplitude Spectrum of phi"
) -> None:
    plt.figure()
    plt.plot(freq, amp)
    plt.xlabel("frequency")
    plt.ylabel("|Phi(f)|")
    plt.title(title)
    if xlim is not None:
        plt.xlim(*xlim)
    if ylim is not None:
        plt.ylim(*ylim)
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(savepath, dpi=200)
    plt.close()

# ========= I/O ユーティリティ =========
def load_phi_from_csv(csv_path: str, phi_col: str = "phi") -> np.ndarray:
    df = pd.read_csv(csv_path)
    if phi_col not in df.columns:
        raise ValueError(f"CSVに '{phi_col}' 列がありません。列名: {df.columns.tolist()}")
    return df[phi_col].to_numpy()

def save_outputs(data: Dict[str, np.ndarray],
                 ts_csv: str = "time_series_used.csv",
                 spec_csv: str = "spectrum_phi.csv") -> None:
    pd.DataFrame({"t": data["t"], "phi": data["phi"]}).to_csv(ts_csv, index=False)
    pd.DataFrame({"freq": data["freq"], "amp": data["amp"], "phase": data["phase"]}).to_csv(spec_csv, index=False)

# ========= main =========
def main():
    # 入力設定
    args = sys.argv
    if (len(args)<3):
        raise ValueError("Number of args is not matched.\n usage: (arg1)project name, (arg2)value name, (arg3,optional)exec fft flag")
    pname = args[1] # ファイル名
    valName = args[2] # 変数名
    filo_ts = pname+"_used_"+valName+".csv"
    filo_sp = pname+"_spectrum_"+valName+".csv"
    bool_fft = (args[3].lower() == "true") if len(args) >= 4 else True # fft 実行フラグ（デフォルトTrue, True 以外の場合はFalse）
    dt = 5e-13  # [s]

    # 描画範囲（必要なければ None）
    t_xlim = None           # 例: (0, 100)
    t_ylim = None           # 例: (-1, 1)
    f_xlim = (0,1e10)       # 例: (0, 50)
    f_ylim = None           # 例: (0, None)

    # 読み込み
    phi = load_phi_from_csv(pname+".csv", phi_col=valName)

    if bool_fft:
        # fft 計算
        data = compute_fft_from_phi(phi, dt)

        # 結果保存
        save_outputs(data, filo_ts, filo_sp)
    else:
        # 既存ファイルから読み込み
        ts_df = pd.read_csv(filo_ts)
        sp_df = pd.read_csv(filo_sp)
        data = {
            "t": ts_df["t"].to_numpy(),
            "phi": ts_df["phi"].to_numpy(),
            "freq": sp_df["freq"].to_numpy(),
            "amp": sp_df["amp"].to_numpy(),
            "phase": sp_df["phase"].to_numpy(),
        }

    # プロット
    plot_time_series(data["t"], data["phi"], savepath=f"fig/{pname}_time_series_{valName}.png", xlim=t_xlim, ylim=t_ylim)
    plot_spectrum(data["freq"], data["amp"], savepath=f"fig/{pname}_spectrum_{valName}.png", xlim=f_xlim, ylim=f_ylim)

    # 情報表示
    print(f"N = {phi.size}")
    print(f"dt = {dt}")
    print(f"time range: {data['t'][0]} .. {data['t'][-1]}")
    print(f"freq range : {data['freq'][1]} .. {data['freq'][-1]}")
    print(f"出力: {filo_ts}, {filo_sp}, fig/{pname}_time_series_{valName}.png, fig/{pname}_spectrum_{valName}.png")

if __name__ == "__main__":
    main()
