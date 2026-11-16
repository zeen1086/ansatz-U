# 冻结 V1 → V_func_tmp → cal_RVR_split → cal_Hsk

这是第30轮局域核、v/B/完整N读出与反馈模型的正式接口。模型包采用 `LANIO2_LOCAL_VFUNC_V2`，原生投影和组装采用 `LANIO2_LRVR_V2`。本目录只提供当前 V1 的 `local/` 程序；旧五记录和 scalar8/B40 RVR 接口已退役，依赖与历史入口见 [目录说明](../README.md)。

当前本机正式模型入口：

`/Users/lijianyu/Desktop/doctor/ansatz_workflow_20261030/current/models/V_func_tmp`

导出包是自包含目录：原始 `C_v_local.npy`、`C_B_local.npy`、`local_grid.npz`、模型参数及六个冻结推理源文件，另加 `V_FUNC.json` 和 `native_geometry.bin`。不需要原训练目录、旧四Ni采样缓存或运行时DFT占据。参数和空间核逐字节复制，未重拟合；包移动后仍可调用。

## 1. 一次导出和编译

在 `lanio2_ansatzU_fortran` 下，使用主机已有Python环境。25_li为 `/home/25_li/ED/.conda/bin/python3`；本机为 `/Users/lijianyu/.local/share/venvs/codex-shared/bin/python`。

```sh
"$ANSATZ_PYTHON" local_vfunc.py export \
  --source /absolute/frozen_V1/result/final \
  --code /absolute/frozen_V1/model_code \
  --structure /absolute/original_L02_structure.npz \
  --scalar-offset-eV=-7.187822007469299 \
  --output /absolute/export/V_func_tmp

make local-rvr BUILD=/absolute/external/build
"$ANSATZ_PYTHON" local_vfunc.py verify /absolute/export/V_func_tmp
```

标量常数来自当前V1已约定的Néel规范 `−0.528295173387505 Ry`，固定用于所有纹理。导出几何为原始四Ni结构；`structure.npz` 包含 `cell, positions, species, ni, spins`，Å单位、原子序数、零基Ni索引与单位方向。导出、预测与原生程序均拒绝覆盖已有输出。

macOS默认gfortran/Accelerate。25_li的Intel编译示例：

```sh
make local-rvr BUILD=/absolute/external/build FC=ifort \
  FFLAGS='-O2 -check bounds -free -extend-source -assume byterecl' LIBS='-mkl=sequential'
```

## 2. 绑定纯H0、共同Wannier基底并准备投影

已有DFT-only `primitive.bin` 可直接导出原生 `.au` 背景：

```sh
"$ANSATZ_PYTHON" local_vfunc.py export-background \
  --source /absolute/DFT_only/primitive.bin --output /absolute/basis/pure_H0.au
```

导出沿用原来的WS变换与 `p(q)=raw_projector*dis*U`，舍弃输入中的标量势和N0，不读取旧拟合势；输出包括 `.au.json` 来源/单位记录。原接口需要30空间轨道、原始8×8×8 q网格。其他来源不可仅改名为纯H0。

Wannier输入沿用原始 `wan_1 ... wan_30`：little-endian complex128流文件、Fortran排列 `(n_translation, n_grid)`，原始 `+iqR` 翻译约定；投影器自行使用逆平移。网格必须与导出包一致。权重文件首行翻译数，随后每行 `Rx Ry Rz weight`，显式记录保留支持和边界权重。**不自动裁切支持，不自动给重复边界加权，不重归一化Wannier。**

```sh
"$ANSATZ_PYTHON" local_vfunc.py prepare-projection \
  --model /absolute/export/V_func_tmp \
  --wan-dir /absolute/common_gauge/wan \
  --translations /absolute/weighted_translations.dat \
  --background /absolute/basis/pure_H0.au \
  --grid-weight 1 \
  --kz 0 0.125 0.25 0.375 0.5 0.625 0.75 0.875 \
  --output /absolute/RVR
```

`grid-weight=1`适用于既有离散归一化来源；其他基底需使用其实际积分权重。准备器校验字节数、记录模型/Wannier/背景身份，并估算投影文件大小；不执行投影。`PROJECTION.json`与`projection.nml`一起保留。来源绑定防止误混文件；共同gauge的物理来源仍需由实际数据生产记录核验。

## 3. 原生 cal_RVR_split 和 combine

```sh
/absolute/build/local/cal_RVR_split /absolute/RVR/projection.nml 0
/absolute/build/local/cal_RVR_split /absolute/RVR/projection.nml 1
# 对空间轨道2..30执行同一命令；可按轨道行并行。
/absolute/build/local/combine /absolute/RVR/projection.nml
```

行0直接读取 `V_func_tmp` 的局域空间核，只做一次三线性插值、紧支撑窗口与z周期求和，写共享 `spatial_fields`。行1..30读取缓存与共同Wannier，执行Fortran/BLAS投影；第三个可选参数指定一个1基kz索引，支持按轨道/kz分片。合并器要求所有分片完整且模型、基底、平移和权重一致，输出 `B_general`。

投影结果保存所有 scalar62+bias1、B62空间通道。z方向周期固定为1，按显式kz投影；x/y保持真实原胞平移，所以同一个B_general可用于不同L与自旋纹理。新kz要补投影；不对未投影kz默默插值。

标量去均值通过同一Wannier支持上的重叠矩阵处理，不假设截断基底的重叠恰为单位阵。不同真实平移折叠到同一超胞位点时仍逐项累加。

## 4. 一条命令从结构/方向到原生 cal_Hsk

```sh
"$ANSATZ_PYTHON" local_vfunc.py cal_Hsk \
  --model /absolute/export/V_func_tmp \
  --structure /absolute/L04_structure.npz \
  --rvr /absolute/RVR/B_general \
  --background /absolute/basis/pure_H0.au \
  --kpoints /absolute/kpoints.dat \
  --binary /absolute/build/local/cal_Hsk \
  --derivative-sites 1 2 6 --step 0.001 \
  --components --output /absolute/H_L04_case01
```

读取器支持固定原胞几何的 L×L×1，L=1..100，30空间轨道/原胞，显式自旋后维数60L²。结构变化需要重建相应空间投影，运行时会拒绝与包几何不匹配的原子结构。Ni输入自动重排为 `x*L+y+1`，排列保存在输出中；2×2对应旧标签5、7、6、8。

Python只执行冻结的原始特征、非线性读出和反馈，输出 `value.bin` 及选定Ni的 `dsite_XXXX.bin`。每个角度导数均对收敛后的完整模型作五点差分，所有Ni的p和完整N响应均保留；导数站点是1基原生Ni编号，方向为绕全局y轴旋转。省略 `--derivative-sites` 只输出H。

随后调用原生Fortran `cal_Hsk`：

`H=H0+Hv+HB+Σ A_i†Ueff(I/2−N_i)A_i`。

H0、v、B、HU各计一次；导数模式不加入H0、偏置常数或I/2项。全部矩阵为complex128/eV，不强制对称化。k点文件首行点数，随后每行 `kx ky kz weight`，权重和为1。Hubbard顺序为 `[up d1..5,down d1..5]`，Hamiltonian为 `(cell,orbital,up/down)`。

每个构型输出：

- `prediction.npz`：p、完整N、VU、dp/dN、Ni排序。
- `value.bin`、`dsite_XXXX.bin`：原生模型输入，含模型ID，防止混用RVR。
- `value.k000001.h` 等：逐k总H；导数对应 `dsite_XXXX.k000001.h`。
- 开启 `--components` 时另存 `.h0/.hv/.hb/.hu`。
- `.index`：维数、k点/权重、标量均值及Hermiticity；`HAMILTONIAN.json`仅在全部调用成功后写出。

原生矩阵读取示例：`np.fromfile(path,dtype='<c16').reshape(nd,nd,order='F')`。也可单独用 `local_vfunc.py evaluate` 生成系数，再用原生 `cal_Hsk run.nml` 调用；公共入口额外核验基底与背景来源。

## 验证边界

`tests/test_local_vfunc.py` 的原生链测试覆盖L=1/2/3、非零kx/ky/kz、完整复数N、标量均值、三势分量、H及dH，并与独立直接实空间积分比较。`tests/test_v1_runtime_bridge.py` 以真实V1反馈和明确的测试投影核核验整个调用器；它不是材料Hamiltonian验收。

真实V1导出、预测和原生局域核求值的核对记录在项目 `current/validation/V1_native_interface_20261114_r01/`。2026-11-15 的真实共同 Wannier 完整 RVR 及 Néel/FM 的 2×2 H/能带回归已有通过记录；同一已知标量规范下最大差约 5.2e−12 eV。这是模型实现的一致性回归，不是独立 DFT 物理精度验证。批量生产与 2026-11-16 完整能量入口见 [production](../../7_spin/production/README.md)。
