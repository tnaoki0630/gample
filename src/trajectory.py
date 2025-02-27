import numpy as np
import glob
import matplotlib.pyplot as plt
import matplotlib.animation as animation

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
files = sorted(glob.glob("bin/PhaseSpace_electron_*.bin"))
if not files:
    raise FileNotFoundError("PhaseSpace_electron_*.bin が見つかりません。")

# 各粒子のデータをリストに読み込み（各ファイルは1粒子分の時系列データ）
particle_data = []
for f in files:
    data = np.fromfile(f, dtype=record_dtype)
    particle_data.append(data)

num_particles = len(particle_data)
# 各粒子の記録数（全粒子で同じ時刻ステップであることを前提）
num_frames = particle_data[0].shape[0]
for d in particle_data:
    if d.shape[0] != num_frames:
        raise ValueError("粒子間でフレーム数が一致していません。")

# 各フレームごとに全粒子の位置 (x, y) をまとめる配列を作成
# frame_positions: shape = (num_frames, num_particles, 2)
frame_positions = np.empty((num_frames, num_particles, 2))
# 各フレームの時刻情報（全粒子で同一のはず）
times = np.empty(num_frames, dtype=np.int32)
for j in range(num_frames):
    times[j] = particle_data[0]["time"][j]
    for i in range(num_particles):
        frame_positions[j, i, 0] = particle_data[i]["x"][j]
        frame_positions[j, i, 1] = particle_data[i]["y"][j]

# プロットの設定
fig, ax = plt.subplots()
scat = ax.scatter(frame_positions[0, :, 0], frame_positions[0, :, 1])
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_title(f"Time step: {times[0]}")

# 軸の範囲はシミュレーション領域に合わせて調整してください
ax.set_xlim(0, 2)
ax.set_ylim(0, 2)

def update(frame):
    # フレームごとに粒子の位置を更新
    scat.set_offsets(frame_positions[frame])
    print(frame_positions[frame])
    ax.set_title(f"Time step: {times[frame]}")
    return scat,

# アニメーション作成（intervalは各フレームの表示間隔[ms]、fpsは動画保存時のフレームレート）
ani = animation.FuncAnimation(fig, update, frames=num_frames, interval=100, blit=True)

# 動画として保存（ffmpegがインストールされている必要があります）
ani.save("particle_trajectories.mp4", writer="ffmpeg", fps=10)


