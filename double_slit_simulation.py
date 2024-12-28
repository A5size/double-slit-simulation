"""
This code was refactored by Claude.
このコードはClaudeによってリファクタリングされました。
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.animation import FuncAnimation
from mpl_toolkits.axes_grid1 import make_axes_locatable
from dataclasses import dataclass
from typing import Tuple, Optional
from tqdm import tqdm
from lib.QMLib import QMLib

@dataclass
class SimulationConfig:
    """
    Configuration parameters for the quantum simulation
    量子シミュレーションの設定パラメータ
    """
    # Grid settings / グリッド設定
    dx: float = 0.1  # Spatial step in x direction / x方向の空間ステップ
    dy: float = 0.1  # Spatial step in y direction / y方向の空間ステップ
    nx: int = 800    # Number of grid points in x / x方向のグリッド点数
    ny: int = 600    # Number of grid points in y / y方向のグリッド点数
    margin_nx: int = 300  # Margin for x visualization / x方向の可視化マージン
    margin_ny: int = 300  # Margin for y visualization / y方向の可視化マージン
    
    # Gaussian wave packet parameters / ガウス波束のパラメータ
    gauss_x: float = -15.0  # Initial x position / 初期x位置
    gauss_y: float = 0.0    # Initial y position / 初期y位置
    gauss_kx: float = 6.0   # Wave number in x / x方向の波数
    gauss_ky: float = 0.0   # Wave number in y / y方向の波数
    gauss_c: float = 0.25   # Wave packet width parameter / 波束の幅パラメータ
    
    # Double slit parameters / 二重スリットのパラメータ
    slit_a: float = 0.5   # Slit width / スリットの幅
    slit_d: float = 2.0   # Distance between slits / スリット間距離
    slit_l: float = 10.0  # Distance to screen / スクリーンまでの距離
    slit_w: float = 1.0   # Slit thickness / スリットの厚さ
    slit_p: float = 50.0  # Potential strength / ポテンシャルの強さ
    
    # Slice coordinates / スライス座標
    x_slice_coord_1: float = 0.0    # First slice position / 1番目のスライス位置
    x_slice_coord_2: float = 10.0   # Second slice position (screen) / 2番目のスライス位置（スクリーン）
    
    # Animation settings / アニメーション設定
    replot_interval: int = 10      # Steps between plots / プロット間のステップ数
    total_frames: int = 150       # Total animation frames / アニメーションの総フレーム数
    animation_interval: int = 50    # Interval between frames (ms) / フレーム間隔（ミリ秒）

class DoubleSlitSimulation:
    """
    Handles the quantum mechanical simulation of the double-slit experiment
    二重スリット実験の量子力学シミュレーションを処理するクラス
    """
    def __init__(self, config: SimulationConfig):
        self.config = config
        self.qm = QMLib()
        self.initialize_simulation()
        
    def initialize_simulation(self):
        """
        Initialize all components of the simulation
        シミュレーションの全コンポーネントを初期化
        """
        self._setup_grid()
        self._initialize_arrays()
        self._setup_wave_function()
        self._setup_potential()
        
    def _setup_grid(self):
        """
        Set up the computational grid
        計算グリッドのセットアップ
        """
        self.qm.set_dx(self.config.dx)
        self.qm.set_dy(self.config.dy)
        self.qm.set_nx(self.config.nx)
        self.qm.set_ny(self.config.ny)
        self.qm.init_2d()
        
        # Get actual values (might be adjusted by the Fortran code)
        # 実際の値を取得（Fortranコードによって調整される可能性がある）
        self.dx = self.qm.get_dx()
        self.dy = self.qm.get_dy()
        
    def _initialize_arrays(self):
        """
        Initialize arrays for fields and wave functions
        場と波動関数の配列を初期化
        """
        self.field_x = np.zeros(self.config.nx)
        self.field_y = np.zeros(self.config.ny)
        self.psi2d_real = np.zeros((self.config.nx, self.config.ny))
        self.psi2d_imag = np.zeros((self.config.nx, self.config.ny))
        self.v2d_real = np.zeros((self.config.nx, self.config.ny))
        self.v2d_imag = np.zeros((self.config.nx, self.config.ny))
        
    def _setup_wave_function(self):
        """
        Initialize the Gaussian wave packet
        ガウス波束を初期化
        """
        self.qm.put_gaussian(
            self.config.gauss_x, self.config.gauss_y,
            self.config.gauss_kx, self.config.gauss_ky,
            self.config.gauss_c
        )
        
    def _setup_potential(self):
        """
        Set up the double-slit potential barrier
        二重スリットのポテンシャル障壁をセットアップ
        """
        # Upper barrier / 上部の障壁
        self.qm.set_v2d_box(
            -self.config.slit_w/2.0, -self.config.dy*self.config.ny/2.0,
            self.config.slit_w/2.0, -self.config.slit_d/2.0-self.config.slit_a,
            self.config.slit_p, 0.0
        )
        # Middle barrier / 中央の障壁
        self.qm.set_v2d_box(
            -self.config.slit_w/2.0, -self.config.slit_d/2.0,
            self.config.slit_w/2.0, self.config.slit_d/2.0,
            self.config.slit_p, 0.0
        )
        # Lower barrier / 下部の障壁
        self.qm.set_v2d_box(
            -self.config.slit_w/2.0, self.config.slit_d/2.0+self.config.slit_a,
            self.config.slit_w/2.0, self.config.dy*self.config.ny/2.0,
            self.config.slit_p, 0.0
        )
        
    def get_field_data(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
        Get field and potential data for visualization
        可視化のための場とポテンシャルデータを取得
        """
        self.qm.get_field_x(self.field_x)
        self.qm.get_field_y(self.field_y)
        self.qm.get_v2d(self.v2d_real, self.v2d_imag)
        return self.field_x, self.field_y, self.v2d_real, self.v2d_imag

class Visualization:
    """
    Handles all visualization aspects of the simulation
    シミュレーションの可視化を処理するクラス
    """
    def __init__(self, simulation: DoubleSlitSimulation):
        self.sim = simulation
        self.config = simulation.config
        self._setup_figure()
        self._setup_progress()
        
    @staticmethod
    def calculate_theory_intensity(y_array: np.ndarray, a: float, d: float, L: float,
                                 wavelength: float, I0: float = 1.0) -> np.ndarray:
        """
        Calculate theoretical intensity for Fraunhofer diffraction
        フラウンホーファー回折の理論強度を計算

        Parameters:
            y_array: Observation points / 観測点
            a: Slit width / スリット幅
            d: Slit separation / スリット間隔
            L: Distance to screen / スクリーンまでの距離
            wavelength: Wavelength / 波長
            I0: Initial intensity / 初期強度
        """
        alpha = np.pi * a * y_array / (wavelength * L)
        beta = np.pi * d * y_array / (wavelength * L)
        envelope = np.sinc(alpha/np.pi)**2  # Single slit diffraction / 単スリット回折
        interference = np.cos(beta)**2      # Double slit interference / 二重スリット干渉
        return I0 * envelope * interference
    
    def _setup_progress(self):
        """
        Set up progress bar and statistics tracking
        進捗バーと統計追跡のセットアップ
        """
        self.pbar = tqdm(
            total=self.config.total_frames,
            desc="Simulation Progress",
            bar_format="{desc}: {percentage:3.0f}%|{bar}| {n_fmt}/{total_fmt} " +
                      "[Time: {elapsed}<{remaining}] " +
                      "{postfix}",
            ncols=100
        )
        self.pbar.postfix = {
            "t": 0.0,    # Simulation time / シミュレーション時間
            "P": 1.0     # Total probability / 全確率
        }

    def _calculate_total_probability(self) -> float:
        """
        Calculate total probability (should be conserved)
        全確率を計算（保存されるべき量）
        """
        probability_density = self.sim.psi2d_real**2 + self.sim.psi2d_imag**2
        return np.sum(probability_density) * self.sim.dx * self.sim.dy
    
    def _setup_figure(self):
        """
        Set up the matplotlib figure and axes
        Matplotlibの図とグラフ軸をセットアップ
        """
        self.fig = plt.figure(figsize=(12, 6))
        gs = GridSpec(2, 2, figure=self.fig, width_ratios=[3, 2], height_ratios=[1, 1])
        
        # Get field data / 場のデータを取得
        field_x, field_y, v2d_real, v2d_imag = self.sim.get_field_data()
        
        # Calculate slice indices / スライスのインデックスを計算
        self.x_slice_1 = np.argmin(np.abs(field_x - self.config.x_slice_coord_1))
        self.x_slice_2 = np.argmin(np.abs(field_x - self.config.x_slice_coord_2))
        
        # Set up main heatmap and cross-section plots
        # ヒートマップと断面プロットをセットアップ
        self._setup_heatmap(gs, field_x, field_y, v2d_real, v2d_imag)
        self._setup_cross_sections(gs, field_x, field_y)
        
        plt.tight_layout()
        
    def _setup_heatmap(self, gs: GridSpec, field_x: np.ndarray, field_y: np.ndarray,
                      v2d_real: np.ndarray, v2d_imag: np.ndarray):
        """
        Set up the probability density heatmap
        確率密度のヒートマップをセットアップ
        """
        self.ax_map = self.fig.add_subplot(gs[:, 0])
        
        # Create potential mask / ポテンシャルマスクを作成
        potential_abs = np.sqrt(v2d_real**2 + v2d_imag**2)
        potential_mask = np.where(potential_abs > 0, 1.0, np.nan)
        
        # Calculate plot extents / プロット範囲を計算
        margin_nx, margin_ny = self.config.margin_nx, self.config.margin_ny
        extent = (
            field_x.min()+self.sim.dx*margin_nx/2.0,
            field_x.max()-self.sim.dx*margin_nx/2.0,
            field_y.min()+self.sim.dy*margin_ny/2.0,
            field_y.max()-self.sim.dy*margin_ny/2.0
        )
        
        # Plot potential barrier / ポテンシャル障壁をプロット
        self.ax_map.imshow(
            potential_mask[margin_nx//2:self.config.nx-margin_nx//2,
                         margin_ny//2:self.config.ny-margin_ny//2].T,
            extent=extent,
            origin="lower",
            cmap="gray",
            aspect="equal",
            alpha=1.0
        )
        
        # Initialize probability density heatmap / 確率密度ヒートマップを初期化
        self.heatmap = self.ax_map.imshow(
            np.zeros((self.config.nx-margin_nx, self.config.ny-margin_ny)).T,
            extent=extent,
            origin="lower",
            cmap="viridis",
            aspect="equal",
            alpha=0.8,
            vmin=0,
            vmax=0.01
        )
        
        # Add colorbar / カラーバーを追加
        divider = make_axes_locatable(self.ax_map)
        cax = divider.append_axes("bottom", size="5%", pad=0.5)
        self.fig.colorbar(self.heatmap, cax=cax, orientation="horizontal")
        
        # Add slice lines / スライス線を追加
        self.line_map_1 = self.ax_map.axvline(x=field_x[self.x_slice_1],
                                             color="#ff6347",
                                             label="x_slice_1",
                                             alpha=0.5)
        self.line_map_2 = self.ax_map.axvline(x=field_x[self.x_slice_2],
                                             color="#4169e1",
                                             label="x_slice_2",
                                             alpha=0.5)
        
        self.ax_map.set_title("Probability Density |Ψ|^2")
        self.ax_map.set_xlabel("x")
        self.ax_map.set_ylabel("y")
        
    def _setup_cross_sections(self, gs: GridSpec, field_x: np.ndarray, field_y: np.ndarray):
        """
        Set up the cross-sectional plots showing probability density at specific x positions
        特定のx位置での確率密度を示す断面プロットをセットアップ
        """
        margin_ny = self.config.margin_ny
        y_screen = field_y[margin_ny//2:self.config.ny-margin_ny//2]
        
        # Calculate theoretical intensity / 理論強度を計算
        theory_intensity = self.calculate_theory_intensity(
            y_array=y_screen,
            a=self.config.slit_a,
            d=self.config.slit_d,
            L=self.config.slit_l,
            wavelength=2.0*np.pi/self.config.gauss_kx,
            I0=0.005
        )
        
        # First cross-section (at slit position) / 1番目の断面（スリット位置）
        self.ax_cross1 = self.fig.add_subplot(gs[0, 1])
        self.line1, = self.ax_cross1.plot(
            field_y[margin_ny//2:self.config.ny-margin_ny//2],
            np.zeros(self.config.ny-margin_ny),
            label=f"x = {field_x[self.x_slice_1]:.2f}",
            color="#ff6347"
        )
        self._setup_cross_section_plot(self.ax_cross1, field_x[self.x_slice_1])
        
        # Second cross-section (at screen position) / 2番目の断面（スクリーン位置）
        self.ax_cross2 = self.fig.add_subplot(gs[1, 1])
        self.line2, = self.ax_cross2.plot(
            field_y[margin_ny//2:self.config.ny-margin_ny//2],
            np.zeros(self.config.ny-margin_ny),
            label=f"x = {field_x[self.x_slice_2]:.2f}",
            color="#4169e1"
        )
        self.line_theory, = self.ax_cross2.plot(y_screen, theory_intensity, "k--", label="Theory")
        self._setup_cross_section_plot(self.ax_cross2, field_x[self.x_slice_2])
        
    def _setup_cross_section_plot(self, ax: plt.Axes, x_coord: float):
        """
        Configure a single cross-section plot
        1つの断面プロットを設定
        """
        ax.set_title(f"Cross-section at x = {x_coord:.2f}")
        ax.set_xlabel("y")
        ax.set_ylabel("Probability Density")
        ax.set_ylim([0, 0.008])
        ax.legend()
        
    def update(self, frame):
        """
        Update function for animation
        アニメーションの更新関数
        """
        # Step simulation forward / シミュレーションを進める
        self.sim.qm.step_2d(self.config.replot_interval)
        current_time = self.sim.qm.get_time()
        
        # Get updated wave function / 波動関数を更新
        self.sim.qm.get_psi2d(self.sim.psi2d_real, self.sim.psi2d_imag)
        probability_density = self.sim.psi2d_real**2 + self.sim.psi2d_imag**2
        
        # Calculate total probability / 全確率を計算
        total_prob = self._calculate_total_probability()
        
        # Update progress bar / 進捗バーを更新
        self.pbar.update(1)
        self.pbar.set_postfix({
            "t": f"{current_time:.3f}",
            "P": f"{total_prob:.6f}"
        })
        
        margin_nx, margin_ny = self.config.margin_nx, self.config.margin_ny
        
        # Update probability density heatmap / 確率密度ヒートマップを更新
        self.heatmap.set_data(
            probability_density[margin_nx//2:self.config.nx-margin_nx//2,
                              margin_ny//2:self.config.ny-margin_ny//2].T
        )
        
        # Update cross-sectional plots / 断面プロットを更新
        cross_section_1 = probability_density[self.x_slice_1,
                                           margin_ny//2:self.config.ny-margin_ny//2]
        cross_section_2 = probability_density[self.x_slice_2,
                                           margin_ny//2:self.config.ny-margin_ny//2]
        self.line1.set_ydata(cross_section_1)
        self.line2.set_ydata(cross_section_2)
        
        return self.heatmap, self.line1, self.line2, self.line_map_1, self.line_map_2
    
    def animate(self):
        """
        Create and display the animation
        アニメーションを作成して表示
        """
        self.anim = FuncAnimation(
            self.fig,
            self.update,
            frames=self.config.total_frames,
            interval=self.config.animation_interval,
            blit=False
        )
        try:
            plt.show()
        finally:
            self.pbar.close()

def main():
    """
    Main function to run the simulation
    シミュレーションを実行するメイン関数
    """
    config = SimulationConfig()
    simulation = DoubleSlitSimulation(config)
    viz = Visualization(simulation)
    viz.animate()

if __name__ == "__main__":
    main()
