import numpy as np
import glob
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def getParticleTrajectory(fileName):
    # バイナリファイルのレコード形式（time:int, x, y, vx, vy, vz:double）
    record_dtype = np.dtype([
        ("time", np.int32),
        ("x", np.float32),
        ("y", np.float32),
        ("vx", np.float32),
        ("vy", np.float32),
        ("vz", np.float32)
    ])

    # ファイル名のパターンに合わせて全粒子のファイルを取得
    files = sorted(glob.glob(fileName))
    if not files:
        raise FileNotFoundError(fileName)

    # 各粒子のデータを読み込み（各ファイルは1粒子分の時系列データ）
    particle_data = [np.fromfile(f, dtype=record_dtype) for f in files]

    num_particles = len(particle_data)
    num_frames = particle_data[0].shape[0]
    for d in particle_data:
        if d.shape[0] != num_frames:
            raise ValueError("粒子間でフレーム数が一致していません。")

    # 各粒子の (x, y) をベクトル化して frame_positions を作成（shape: (num_frames, num_particles, 2)）
    frame_positions = np.stack([np.column_stack((pd["x"], pd["y"])) for pd in particle_data], axis=1)
    # 時刻情報は最初の粒子のデータから取得
    times = particle_data[0]["time"]

    return times, particle_data

if __name__ == "__main__":
    # バイナリデータ抽出
    fileName = "bin/proj1_PhaseSpace_electron_*.bin"
    times, particleData = getParticleTrajectory(fileName)
    fileName = "bin/proj1_restart_PhaseSpace_electron_*.bin"
    times_restart, particleData_restart = getParticleTrajectory(fileName)
    
    # numpy 配列として取得
    timeSeriesX = np.array([pd["x"] for pd in particleData]).T
    timeSeriesY = np.array([pd["y"] for pd in particleData]).T
    timeSeriesX_restart = np.array([pd["x"] for pd in particleData_restart]).T
    timeSeriesY_restart = np.array([pd["y"] for pd in particleData_restart]).T
    
    # プロットの設定
    fig, ax = plt.subplots()
    delay = 29
    scat = ax.scatter(timeSeriesX[delay,:], timeSeriesY[delay,:], c="r", alpha=0.3, label="org")
    scat2 = ax.scatter(timeSeriesX_restart[0,:], timeSeriesY_restart[0,:], c="b", alpha=0.3, label="restart")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(f"Time step: {times[0]}")
    ax.set_xlim(0, 2.5)
    ax.set_ylim(0, 1)
    ax.set_aspect('equal', adjustable='box')
    ax.legend(loc=4)
    # plt.show()

    def update(frame):
        # フレームごとに粒子の位置を更新
        # note:: np.c_() : shape (N,), shape (N,) -> shape (N,2)
        scat.set_offsets(np.c_[timeSeriesX[frame+delay,:], timeSeriesY[frame+delay,:]])
        scat2.set_offsets(np.c_[timeSeriesX_restart[frame,:], timeSeriesY_restart[frame,:]])
        # print(timeSeriesX[frame-1+delay,:])
        # print(timeSeriesX_restart[frame-1,:])
        ax.set_title(f"Time step: {times_restart[frame]}")
        return scat, scat2

    # アニメーション作成
    ani = animation.FuncAnimation(fig, update, frames=timeSeriesX_restart.shape[0], interval=100, blit=True)

    # 動画として保存（ffmpeg のインストールが必要）
    ani.save("particle_trajectories.mp4", writer="ffmpeg", fps=5)
