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
    fileName = "bin/PhaseSpace_restart_electron_*.bin"
    times, particleData = getParticleTrajectory(fileName)
    fileName = "bin/PhaseSpace_restart2_electron_*.bin"
    times_restart, particleData_restart = getParticleTrajectory(fileName)

    timeSeriesX = np.array([pd["x"] for pd in particleData]).T
    timeSeriesX_restart = np.array([pd["x"] for pd in particleData_restart]).T
    print(f"times: {times, times_restart}")
    print(f"timeSeriesX: {timeSeriesX}")
    print(f"timeSeriesX_restart: {timeSeriesX_restart}")
    print(f"diff: {timeSeriesX - timeSeriesX_restart}")