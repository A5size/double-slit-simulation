# Double Slit Quantum Simulation

![Simulation Image](fig/double_slit.png)

---

## Overview (English)

This project demonstrates a **double-slit quantum simulation** using Fortran for high-performance numerical calculations and Python (via `ctypes`) for visualization. A 2D Schrödinger equation is solved with a finite-difference + 4th-order Runge–Kutta (RK4) method. The resulting wavefunction (probability density) is animated using Matplotlib in Python.

### Features

- **Fortran Code (qm.f90)**  
  Implements the numerical solution of the 2D time-dependent Schrödinger equation:
  \[
    i \frac{\partial \Psi}{\partial t} = -\frac{1}{2} \nabla^2 \Psi + V \Psi
  \]
  - Initializes Gaussian wave packets  
  - Sets up double-slit potential barriers  
  - Performs time evolution via RK4

- **Python Wrapper (QMLib.py)**  
  Provides a convenient Python interface to call the Fortran routines using `ctypes`.

- **Simulation Scripts**
  - **double_slit_simulation.py**  
    Offers advanced visualization of the quantum system (e.g., slicing through the wavefunction, comparing to analytical Fraunhofer diffraction, etc.).
  - **simplified_double_slit_simulation.py**  
    A minimalistic script that focuses on displaying a simple 2D heatmap animation of the wavefunction.

### Requirements

- A Fortran compiler (e.g., gfortran)
- Python 3.x
  - numpy
  - matplotlib
  - tqdm (used in the more detailed simulation script)
- Linux or macOS preferred (Windows works with WSL or an appropriate Fortran setup)

### Setup Instructions

1. **Build the Fortran Library**  
   In the `lib` directory, compile the Fortran code to create `libfort.so` (shared library):

   ```
   cd lib
   ./compile.sh
   ```

   or run a similar command manually:

   ```
   gfortran -shared -fPIC qm.f90 -o libfort.so
   ```

2. **Install Python Dependencies**

   ```
   pip install numpy matplotlib tqdm
   ```

3. **Run the Simulation**
   - For full-featured visualization:

     ```
     python double_slit_simulation.py
     ```

   - For a simplified animation:

     ```
     python simplified_double_slit_simulation.py
     ```

### Usage Highlights

- **double_slit_simulation.py**  
  - Uses a `SimulationConfig` dataclass for easy adjustment of parameters (grid size, Gaussian wave packet parameters, slit width, slit separation, etc.).
  - Displays heatmaps, cross-sectional plots, and a theoretical Fraunhofer diffraction curve for comparison.
  - Animates the time evolution of the wave packet passing through the double-slit.

- **simplified_double_slit_simulation.py**  
  - A concise script demonstrating only the essential steps:
    1. Initialize the Fortran solver  
    2. Set up a Gaussian wave packet  
    3. Define the double-slit potential  
    4. Animate the probability density in 2D

### License

This project is released under the [MIT License](https://opensource.org/licenses/MIT). Feel free to use or modify the code for your own purposes.

---

## 概要 (日本語)

このプロジェクトは、Fortranによる高速な数値計算と、`ctypes` を介してPythonから呼び出した結果をMatplotlibで可視化することで、**二重スリットの量子シミュレーション**を行うデモコードです。2次元シュレディンガー方程式を有限差分 + 4次のルンゲ＝クッタ法で解き、その確率密度(|Ψ|²)をアニメーションとして表示します。

### 特徴

- **Fortranコード (qm.f90)**  
  2次元シュレディンガー方程式  
  \[
    i \frac{\partial \Psi}{\partial t} = -\frac{1}{2} \nabla^2 \Psi + V \Psi
  \]
  を数値的に解く。  
  - ガウス波束の初期化  
  - 二重スリットポテンシャルの設定  
  - ルンゲ＝クッタ法(RK4)を用いた時間発展

- **Pythonラッパー (QMLib.py)**  
  `ctypes` を使い、FortranのサブルーチンをPython側で簡単に呼び出すためのインターフェース。

- **シミュレーションスクリプト**
  - **double_slit_simulation.py**  
    詳細な可視化（断面プロット、フラウンホーファー回折の理論比較など）を行うためのサンプル。
  - **simplified_double_slit_simulation.py**  
    最小限の設定とアニメーション表示のみを行うシンプルなサンプル。

### 必要環境

- Fortranコンパイラ (例: gfortran)
- Python 3系
  - numpy
  - matplotlib
  - tqdm (詳細版シミュレーションで使用)
- LinuxまたはmacOS推奨 (Windowsの場合、WSLや適切なFortranの環境が必要)

### セットアップ手順

1. **Fortranライブラリのビルド**  
   `lib` ディレクトリ内で共有ライブラリ `libfort.so` を生成します:

   ```
   cd lib
   ./compile.sh
   ```

   または下記のように直接コンパイルしてください:

   ```
   gfortran -shared -fPIC qm.f90 -o libfort.so
   ```

2. **Python依存関係のインストール**

   ```
   pip install numpy matplotlib tqdm
   ```

3. **シミュレーションの実行**
   - 詳細な可視化:

     ```
     python double_slit_simulation.py
     ```

   - シンプルなアニメーション:

     ```
     python simplified_double_slit_simulation.py
     ```

### 使い方のポイント

- **double_slit_simulation.py**  
  - `SimulationConfig` (dataclass) でグリッドサイズやガウス波束のパラメータ、スリットの幅・間隔などを一括管理。  
  - ヒートマップや断面プロット、フラウンホーファー回折の理論曲線を同時に可視化。  
  - 時間発展した波束が二重スリットを通過し、干渉パターンが形成される様子をアニメーションで観察可能。

- **simplified_double_slit_simulation.py**  
  - Fortranライブラリを読み込み、ガウス波束と二重スリットポテンシャルを設定し、確率密度をアニメ表示する最小限のサンプル。
  - コードの学習や動作確認に適した簡易版。

### ライセンス

本プロジェクトは [MIT License](https://opensource.org/licenses/MIT) のもとで公開されています。自由に利用・改変いただけます。
