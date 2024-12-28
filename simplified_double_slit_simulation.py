"""
This code was refactored by Claude.
このコードはClaudeによってリファクタリングされました。
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from lib.QMLib import QMLib

def main():
    # Grid settings / グリッド設定
    dx = 0.1    # Grid spacing in x / x方向の格子間隔
    dy = 0.1    # Grid spacing in y / y方向の格子間隔
    nx = 800    # Number of x points / x方向の点数
    ny = 600    # Number of y points / y方向の点数
    
    # Initialize quantum simulation / 量子シミュレーションの初期化
    qm = QMLib()
    qm.set_dx(dx)
    qm.set_dy(dy)
    qm.set_nx(nx)
    qm.set_ny(ny)
    qm.init_2d()
    
    # Get grid parameters / グリッドパラメータの取得
    dx = qm.get_dx()  # x方向の間隔
    dy = qm.get_dy()  # y方向の間隔
    
    # Initialize arrays / 配列の初期化
    field_x = np.zeros(nx)     # x座標
    field_y = np.zeros(ny)     # y座標
    psi2d_real = np.zeros((nx, ny))  # 波動関数実部
    psi2d_imag = np.zeros((nx, ny))  # 波動関数虚部
    v2d_real = np.zeros((nx, ny))    # ポテンシャル実部
    v2d_imag = np.zeros((nx, ny))    # ポテンシャル虚部

    # Set initial wave packet / 初期波束の設定
    qm.put_gaussian(
        x=-15.0,   # Initial x position / 初期x位置
        y=0.0,     # Initial y position / 初期y位置
        kx=6.0,    # Wave number in x / x方向の波数
        ky=0.0,    # Wave number in y / y方向の波数
        c=0.25     # Wave packet width / 波束の幅
    )

    # Set double-slit potential barrier / 二重スリットポテンシャルの設定
    slit_w = 1.0   # Slit thickness / スリットの厚さ
    slit_a = 0.5   # Slit width / スリットの幅
    slit_d = 2.0   # Distance between slits / スリット間距離
    slit_p = 50.0  # Potential strength / ポテンシャルの強さ
    
    # Upper barrier / 上部の障壁
    qm.set_v2d_box(-slit_w/2.0, -dy*ny/2.0, slit_w/2.0, -slit_d/2.0-slit_a, slit_p, 0.0)
    # Middle barrier / 中央の障壁
    qm.set_v2d_box(-slit_w/2.0, -slit_d/2.0, slit_w/2.0, slit_d/2.0, slit_p, 0.0)
    # Lower barrier / 下部の障壁
    qm.set_v2d_box(-slit_w/2.0, slit_d/2.0+slit_a, slit_w/2.0, dy*ny/2.0, slit_p, 0.0)

    # Get field data / 場のデータ取得
    qm.get_field_x(field_x)
    qm.get_field_y(field_y)
    qm.get_v2d(v2d_real, v2d_imag)

    # Create potential mask / ポテンシャルマスクの作成
    potential_abs = np.sqrt(v2d_real**2 + v2d_imag**2)
    potential_mask = np.where(potential_abs > 0, 1.0, np.nan)

    # Set up plot / プロットのセットアップ
    fig, ax = plt.subplots(figsize=(10, 8))

    # Plot potential barrier / ポテンシャル障壁の表示
    ax.imshow(
        potential_mask.T,
        extent=(field_x.min(), field_x.max(), field_y.min(), field_y.max()),
        origin="lower",
        cmap="gray",
        aspect="equal",
        alpha=1.0
    )

    # Initialize probability density plot / 確率密度プロットの初期化
    heatmap = ax.imshow(
        np.zeros((nx, ny)).T,
        extent=(field_x.min(), field_x.max(), field_y.min(), field_y.max()),
        origin="lower",
        cmap="viridis",
        aspect="equal",
        alpha=0.8,
        vmin=0,
        vmax=0.01
    )

    # Add labels and colorbar / ラベルとカラーバーの追加
    ax.set_title("Probability Density |Ψ|^2")
    ax.set_xlabel("Position x")
    ax.set_ylabel("Position y")
    plt.colorbar(heatmap)

    def update(frame):
        """Animation update function / アニメーション更新関数"""
        # Time evolution / 時間発展
        qm.step_2d(10)
        
        # Get wave function and calculate probability density
        # 波動関数を取得し確率密度を計算
        qm.get_psi2d(psi2d_real, psi2d_imag)
        probability_density = psi2d_real**2 + psi2d_imag**2

        # Update plot / プロット更新
        heatmap.set_data(probability_density.T)
        return heatmap,

    # Create and show animation / アニメーション作成と表示
    anim = FuncAnimation(
        fig,
        update,
        frames=1500,
        interval=50,
        blit=False
    )
    plt.show()

if __name__ == "__main__":
    main()
