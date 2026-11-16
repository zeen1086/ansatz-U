# 当前 V1 的可复用 RVR

本目录只保留第30轮冻结 V1 所需的空间核投影与背景支持。四个源码集中在本目录，操作与输出说明统一如下。

```text
6_RVR/
├── README.md
├── local_vfunc.f90
├── local_rvr.f90
├── rvr.f90
└── frozen_io.f90
```

`rvr.f90` 是统一命令入口，提供 `prepare/project/combine`；其余三个文件是模块。Makefile 输出 `local/rvr` 和下游 `local/cal_Hsk`。旧 `local/cal_RVR_split`、`local/combine` 命令已移除，作业脚本须采用下文的新命令。

| 文件 | 作用与直接依赖 |
|---|---|
| `local_vfunc.f90` | 读取局域空间核和几何，插值并作周期求和；依赖 `common/numerics.f90`、`4_potential/fit/joint_npy.f90` |
| `local_rvr.f90` | Wannier 投影、重叠矩阵、分片校验与合并、RVR 块读取；依赖 `local_vfunc` |
| `rvr.f90` | 统一准备、投影和合并命令，支持补齐缺失分片；依赖 `local_rvr` |
| `frozen_io.f90` | 读取纯 H0/projector `.au`，执行 Bloch 插值、超胞背景和站点 projector 组装；依赖 `numerics`，由 `7_spin/local/cal_Hsk.f90` 调用 |

上游 `local_vfunc.py` 负责模型导出、纯 H0 导出、投影准备和冻结 V1 反馈求值。下游 `7_spin/local/cal_Hsk.f90` 读取 `local_rvr` 与 `frozen_io`，组装 H0/v/B/HU；`7_spin/production` 在该原生入口上计算批量能量与转矩。`joint_npy` 只是共用的 NPY 读取器，保留它不表示仍有 `joint` RVR 依赖。

一次投影包含 scalar62+bias1 与 B62 共 125 个空间核，并保存同一 Wannier 支持上的重叠矩阵。同一模型、几何、基底、支持和已投影 kz 下，RVR 可复用于不同 L 和自旋纹理；新 kz 或改变几何需要相应新投影。

```sh
make local-rvr BUILD=/absolute/external/build
```

输出为 `local/rvr`、`local/cal_Hsk`。两个程序为串行/BLAS 程序，投影按轨道/kz 从外部调度并行，不使用 `mpirun` 启动多份写同一分片的进程。Python 构建入口 `build.py` 同样只链接当前模块。

新计算默认精简输出：合并成功并检查最终文件后，`rvr combine` 清除空间场缓存及轨道/kz 分片；`local_vfunc.py cal_Hsk` 校验全部结果后只保留矩阵和 `HAMILTONIAN.json`。两者添加 `--keep-intermediates` 可保留调试文件。失败时不清理；独立 `evaluate` 的系数和生产续算检查点仍保留。详细规则见下文。

`LANIO2_LOCAL_VFUNC_V2`、`LANIO2_LRVR_V2` 和纯 H0 的 `LANIO2_AU_V1` 布局未变。背景读取器保留旧 KB/N0 槽位的长度、有限值和 Hermiticity 校验，读取后丢弃；当前 B 来自 V1 RVR，N 来自完整反馈模型。`background_projectors`、`project_at` 从旧 `local_ansatz` 原样迁入，周期相位、Ni/轨道/自旋顺序不变。

## 执行流程与文件

不同 L×L×1 超胞共用同一份 `B_general`、冻结 V1 模型包和 `pure_H0.au`；各自重新计算系数、完整 N 和 H。RVR 投影中断后可重跑相同命令，或用 `rvr project projection.nml` 顺序补齐全部缺失分片；单构型 `cal_Hsk` 保留失败现场，已有输出目录不能直接覆盖重跑。生产批量驱动另有检查点复用。只有相应阶段成功并完成校验后才清理其临时输出。

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

标量常数来自当前V1已约定的Néel规范 `−0.528295173387505 Ry`，固定用于所有纹理。导出几何为原始四Ni结构；`structure.npz` 包含 `cell, positions, species, ni, spins`，Å单位、原子序数、零基Ni索引与单位方向。导出和预测拒绝覆盖已有输出；RVR 程序会校验并复用已完成文件，未完成的 `.partial` 在重试时重新生成。

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

## 3. 统一 RVR 命令与中断续算

```sh
/absolute/build/local/rvr prepare /absolute/RVR/projection.nml
/absolute/build/local/rvr project /absolute/RVR/projection.nml
/absolute/build/local/rvr combine /absolute/RVR/projection.nml
```

`prepare` 读取局域空间核，做三线性插值、紧支撑窗口与 z 周期求和，生成共享 `spatial_fields`。`project` 读取空间场和 Wannier，执行 Fortran/BLAS 投影；不指定轨道时按顺序检查全部轨道/kz，只补算缺失分片。也可指定任务以便外部并行调度：

```sh
# 第 1 个轨道的全部 kz
/absolute/build/local/rvr project /absolute/RVR/projection.nml 1
# 第 1 个轨道、第 3 个 kz（均为 1 基编号）
/absolute/build/local/rvr project /absolute/RVR/projection.nml 1 3
# 保留分片以便调试或重新合并
/absolute/build/local/rvr combine /absolute/RVR/projection.nml --keep-intermediates
```

每个缓存、分片或最终文件只允许一个写入进程；可并行投影不同轨道/kz，所有投影进程结束后才执行合并。不要让全量 `project` 与覆盖相同分片的任务同时运行。

恢复规则：

- `spatial_fields`、`B_general_i_kzk` 和 `B_general` 都先写入同目录的 `.partial`；关闭成功并通过对应格式/长度检查后再重命名发布。生成的空间场和分片包含有限值校验。
- 重新执行 `prepare` 会校验并复用完整缓存。`project` 校验正式分片的模型/基底、轨道/kz、平移、权重、长度及全部数据有限性，合格文件不重写。没有正式文件时重新计算，并替换该任务遗留的 `.partial`。
- 正式文件若截断、身份不符或含非有限值，程序报错并保留文件；先将异常文件归档到其他位置，再重新运行。旧版遗留的不完整正式文件也按此处理，不会被当成成功检查点。
- `combine` 只读取空间缓存头、空间域及末尾的标量核均值，跳过整个网格场，不分配完整 `fields(ng,nc,nd)`。分片逐块读取、验证和写入，未完成的合并从头重建 `.partial`，不重做已完成投影。
- 最终 `B_general` 已存在时，三个命令都会核验其模型/基底、几何、尺寸、kz 和完整数据有限性，再复用。若发布后清理被中断，重跑 `combine` 可继续清理，不依赖已删除的空间缓存。

复用期间须保持模型、Wannier、平移/权重及准备记录不变；变更输入应重新执行 `prepare-projection`，生成新的来源绑定和输出目录。检查过程需要读取数据，不是只检查文件名，也不是物理精度验收。

默认精简输出：最终 `B_general` 发布成功后，`combine` 清除该任务的空间缓存、分片及遗留 `.partial`，目录中的其他文件不受影响。由 `prepare-projection` 创建的目录最终保留：

```text
B_general          # 可复用矩阵元
PROJECTION.json    # 来源绑定及准备记录
projection.nml     # 投影输入，便于追溯和重建
kz.dat             # 投影 kz 列表
```

30 个轨道、8 个 kz 时，完成后的文件数仍为 4 个。`--keep-intermediates` 保留调试文件；计算或校验失败时不进入清理。`PROJECTION.json` 仍是准备记录，成功合并输出 `B_GENERAL_COMPLETE`。完成后的文件布局和数值数据格式与之前兼容。

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

`local_vfunc.py cal_Hsk` 默认完成后只保留：

- `value.k000001.h` 等：逐k总H；要求导数时另有 `dsite_XXXX.k000001.h`。
- 开启 `--components` 时另存 `.h0/.hv/.hb/.hu`；默认不生成分量文件。
- `HAMILTONIAN.json`：预测和组装完成记录、输入路径、模型/基底身份及所有原生索引。`native_indices` 按 `value/dsite_XXXX` 保存原 `.index` 的全部文本行，含维数、k点/权重、标量均值及Hermiticity。

只有全部原生调用成功，且索引及所有预期矩阵的字节数通过检查后，才删除 `value.bin/dsite_XXXX.bin`、对应 `.nml/.index` 和中间的 `PREDICTION.json`。失败时保留现场。添加 `--keep-intermediates` 保留这些文件，并生成诊断用 `prediction.npz`（p、完整N、VU、dp/dN、Ni排序）；不改变数值计算。

单独运行 `local_vfunc.py evaluate` 时，`value.bin/dsite_XXXX.bin` 和 `PREDICTION.json` 是后续组装、转矩与断点续算的输入，因此始终保留。只有额外的 `prediction.npz` 改为 `--keep-intermediates` 时输出。生产批量驱动使用这个独立入口，其必要检查点不受单构型 `cal_Hsk` 清理影响。

原生矩阵读取示例：`np.fromfile(path,dtype='<c16').reshape(nd,nd,order='F')`。也可单独用 `local_vfunc.py evaluate` 生成系数，再用原生 `cal_Hsk run.nml` 调用；公共入口额外核验基底与背景来源。

## 验证边界

`tests/test_local_vfunc.py` 的原生链测试覆盖L=1/2/3、非零kx/ky/kz、完整复数N、标量均值、三势分量、H及dH，并与独立直接实空间积分比较。`tests/test_v1_runtime_bridge.py` 以真实V1反馈和明确的测试投影核核验整个调用器；它不是材料Hamiltonian验收。

真实V1导出、预测和原生局域核求值的核对记录在项目 `current/validation/V1_native_interface_20261114_r01/`。2026-11-15 的真实共同 Wannier 完整 RVR 及 Néel/FM 的 2×2 H/能带回归已有通过记录；同一已知标量规范下最大差约 5.2e−12 eV。这是模型实现的一致性回归，不是独立 DFT 物理精度验证。批量生产与 2026-11-16 完整能量入口见 [production](../7_spin/production/README.md)。

## 历史接口与验证记录

2026-11-16 删除的接口：

- 根目录三分量 `cal_RVR_split/combine` 及 `_lanio2` 版本、`rearrange_HR` 和 `speedup`。
- scalar8/B40 的 `joint_rvr` 和 `joint/` 投影/合并入口。
- `tensor_rvr`、旧 `export_frozen.py`，以及其专用 `7_spin/local_ansatz`、`run_ansatzU`、`joint/cal_Hsk`、`run_spin.sh` 与运行模板。

旧源代码及说明仍可从 [清理前提交 3e1307a](https://github.com/zeen1086/ansatz-U/tree/3e1307ac721e91a66eb3851dae2df6bc4c9c99e1) 找回。拟合、DFT、Wannier 与读取已有历史矩阵的独立 `7_spin` 工具保留；当前分支不再生成旧 RVR。

2026-11-16 清理后验证：

| 检查 | 结果 |
|---|---|
| 外部空目录 `make all`，启用 `-fcheck=all`；`build.py` 单独构建三个 V1 程序 | 通过 |
| `check.py` | 23 项通过，1 项因本机缺少 MPI 跳过；删除的两个旧 RVR 测试不再计入 |
| `test_local_vfunc.py`，L=1/2/3、非零 k、完整复数 N 和分量/导数 | 120 个矩阵对照最大差 1.36e−16 eV；40 个批量/单点对照差为 0 |
| `test_v1_runtime_bridge.py`，真实 V1 反馈与已知测试投影核 | 最大差 9.16e−16 eV；15 个矩阵输出与清理前版本逐字节一致 |
| `test_production_adjoint.py`，含非局域 projector、Bloch 相位和批量导数 | 28 项通过，最大转矩差 2.43e−17 eV/rad |
| 已删除路径引用、文档相对链接、背景/projector 函数迁移 | 代码/配置无残留引用，链接有效；两个迁移函数与原文逐字节一致 |

小型测试覆盖实现一致性，不代替完整材料投影、HPC 生产运行或独立 DFT 物理验证。本次未重新运行生产计算。

同日默认输出精简回归：离线测试 24 项通过、1 项 MPI 跳过；原生合并的精简/调试结果 SHA256 一致，截断分片导致失败时全部中间文件保留；真实 V1 的 15 个 H/dH 分量在两种模式下逐字节一致，索引完整归入 JSON。原生进程失败或后续导数矩阵截断时，已成功的前序结果与所有预测输入均保留。生产批量检查点测试继续通过。

同日统一入口与续算验证：`rvr prepare/project/combine` 覆盖遗留 `.partial` 重算、缺失分片补齐、正式分片不重写、错误身份/截断/非有限值拒绝，以及最终发布后中断清理的恢复。合并在已投影的测试缓存网格区注入 NaN 后仍得到相同结果，证明该阶段只依赖元数据和均值；`prepare` 则拒绝该损坏缓存。与改动前原生程序独立投影得到的 `B_general` SHA256 一致；120 个 H/dH 分量及 40 个批量/单点对照通过。该测试模拟中断留下的文件状态，未进行集群断电或文件系统故障试验。
