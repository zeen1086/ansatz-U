# V1 批量 RVR 与完整能量生产入口

本目录整理了 2026-11-15 的实际批量计算驱动和 2026-11-16 的完整能量实现。
计算核沿用生产源码；路径由外部 `ANSATZ_RUN_CONFIG` 指定。
`PROVENANCE.json` 记录导入来源和源文件摘要。V1 参数、空间核、Wannier、RVR、密度及输出均为外部数据。

| 文件 | 职责 |
|---|---|
| `run_material.py` | 原生 H 批量组装、检查点、公共化学势、带能 HF 转矩与壳层 J |
| `torque_contract.py` | 同 kz 最多 16 个条目共享一次 RVR 读取，收缩所有站点的完整 dp/dN |
| `kz_production.py` | 无损节点本地 kz 切片、显式谱/转矩批处理阶段 |
| `complete_energy.py` | 八个 kz worker，完整试探密度能量与 J12；支持此次计算的 L=4/6/8 |
| `density_map.py`、`complete_readout.py` | WS 波函数、密度/势伴随、重叠矩阵、固定电子数响应 |
| `energy_torque.py` | 完整泛函及 Dudarev 全复数 N 的基本能量/导数 |
| `native_large.f90`、`native_input.py` | 已有 QE 环境中的 Hartree/XC/core/离子势求值，无 SCF 迭代 |
| `reference_projection.py` | 与原四胞 Galerkin 实现对照的独立投影检查 |
| `summarize_energy.py`、`check_large_native.py` | J12/能量曲率汇总，以及真实铺展密度的原生回归 |

## 输入与构建

使用已有 Python 3.11+、NumPy、SciPy 环境。原生局域程序需要 Fortran 与 BLAS/LAPACK。
完整能量还需要与原始标签相符、已编译的 QE 7.6 环境、赝势和参考密度。

```sh
make local-rvr BUILD=/absolute/build
# Intel/OpenMP 生产构建示例：
make local-rvr BUILD=/absolute/intel-build FC=ifort \
  FFLAGS='-O2 -check bounds -free -extend-source -assume byterecl -qopenmp' \
  LIBS='-mkl=parallel'
```

V1 模型导出见 [local RVR 接口](../../6_RVR/README.md)。本仓库的
`4_potential/local/` 保存逐字节复制的六个冻结推理源，可作为导出命令的 `--code`。
导出后的 `interface/V_func_tmp` 含这些源及全部冻结模型参数，运行不依赖训练目录。

独立 `local_vfunc.py evaluate` 默认只生成生产所需的 `value.bin/dsite_XXXX.bin` 和 `PREDICTION.json`；重复存储这些数组的 `prediction.npz` 仅在传入 `--keep-intermediates` 时生成。批量驱动仍保留实际续算所需的系数、索引和完成标记，不使用单构型 `cal_Hsk` 命令完成后的临时文件清理。

复制 [配置模板](../../examples/7_spin/production.json) 到外部运行目录，替换所有绝对路径：

```sh
export ANSATZ_RUN_CONFIG=/absolute/run/production.json
```

配置中的路径及数据约定：

- `python`：执行主机已有解释器；`native_build`：含 `local/cal_Hsk` 的构建根目录。
- `interface`：含 `V_func_tmp/V_FUNC.json` 的目录；模型及参数 ID 必须与模板中的冻结 V1 一致。
- `input_root`：含 `inputs/Lxx/structure.npz`、`inputs/Lxx/CASE.json`。
  CASE 沿用生产配对、角度和站点定义；超胞及 Ni 顺序见 local 接口。
- `references`：共同 Wannier 的 `data.npz`、`k_0001.npz ...` 和 `wr_0001.npy ...`。
- `work`：已有材料投影根目录，含 `RVR/B_general`、`RVR/PROJECTION.json`、
  `basis/pure_H0.au`、`basis/bloch.npy` 和通过的 `regression/RESULT.json`。
- `rvr_slices`：可选的 `{"0":"/node-local/z0/B_general", ...}`。缓存路径在对应 worker
  节点上可见；必须与完整 RVR 的模型、基底、尺寸及 kz 匹配。未指定时从完整源截取。
- `old_J`：可选旧四 Ni 对照文件，文件不存在时不进行该比较。

完整能量的 `complete_energy` 配置另外需要：

- `work`：新的完整能量输出目录；每次更换模型、输入、角步或网格使用新目录。
- `reference_run`：原 2×2 参考包，含 `results/BANDS.npz`、`results/FIELDS.npz`、
  `inputs/reference_actual/DENSITY.bin`。原生铺展回归还读取
  `results/native/reference_actual_ENERGY.dat` 与 `reference_actual_TRUE.bin`。
- `pseudo_dir`：相同来源的 `La.upf`、`Ni.upf`、`O.upf`。
- `native_binary`：下面编译得到的原生 QE 适配器。
- `sizes`：此次明确计算的尺寸列表，仅接受 4、6、8；不自动增加尺寸或采样。
- `reuse_cases`：可选的 `{"4":"/existing/L04/mesh2_h0.001", ...}`，读取其
  `case_00 ... case_10/coeff`。已有本征态若存在并通过残差检查，也可复用。
  省略时在新输出目录生成这些构型的完整 V1 系数。

本实现仍固定 LaNiO2 的原胞 40×40×36 网格、30 空间轨道、完整 8³ WS 数据、
Ueff=1 eV、活性电子数 21/Ni 和全电子价电子数 41/Ni；不是任意材料入口。
原生 QE 适配器的构建在已有 QE 的 `PP/src` 下执行，产物放在独立目录：

```sh
mkdir -p /absolute/qe_adapter_build
make -C /absolute/qe/PP/src \
  -f /absolute/repo/7_spin/production/Makefile.native \
  TASK_DIR=/absolute/qe_adapter_build native
```

该 make 配方沿用生产 Intel 编译器的 `-heap-arrays -check bounds`；需使用相同的
QE `make.inc` 和静态库，不能与其他版本的模块文件混编。

## 批量带能读出

在已分配的 PBS 计算任务中运行。`run_material.py` 沿用原来的阶段接口，
按 kz 分组组装 H、求全部 k 点的公共化学势，再计算 HF 转矩：

```sh
/absolute/python 7_spin/production/run_material.py spectrum --L 4 --case 1 --mesh 2
/absolute/python 7_spin/production/run_material.py torque --L 4 --case 1 --mesh 2
/absolute/python 7_spin/production/run_material.py collect --L 4 --mesh 2
```

`collect` 需要 CASE 中全部规定构型；默认不会把单构型输出标为完整 J。
生产计算要求已有材料回归通过。`V1_K_WORKERS` 控制同 kz 批量原生组装的线程数，
`K_WORKERS` 控制批量收缩的并行度；数学库线程数按实际节点配置。

跨构型的本地切片计算使用显式阶段，任务文件为
[批量示例](../../examples/7_spin/batch_tasks.json) 中的 `[L, mesh, step_rad, case]` 列表：

```sh
/absolute/python 7_spin/production/kz_production.py stage --z 0 --rvr /local/z0.bin
/absolute/python 7_spin/production/kz_production.py prepare --z 0 --rvr /local/z0.bin --tasks /run/tasks.json
/absolute/python 7_spin/production/kz_production.py spectrum --z 0 --rvr /local/z0.bin --tasks /run/tasks.json
```

`prepare` 只由一个任务执行；`spectrum` 对 z=0..7 各执行一次。
所有谱完成后，一个任务执行 `spectrum-finish` 求公共化学势，随后 z=0..7
分别执行 `torque`，全部结束后一个任务执行 `torque-finish`。
这些阶段不提交、取消或释放 PBS 作业；历史迁移作业号没有进入新入口。

## 完整能量读出

这是不同于带能的读出。在固定参考补空间密度下，忽略构型无关常数：

```text
F = Fband - Tr(P H_local+U) + EH[rho_active+rho_bg]
    + Exc[rho_active+rho_bg] - Exc_core
    + integral((n_active+n_bg) vion) + EU[N_active+N_bg]
```

令 `A=H_local+U−W†(V_Hxc+vion)W−P_U†VU[Ntotal]P_U`，完整转矩为
`Tr(A dP/dtheta)`。固定电子数的 Gaussian 响应包含全部谱的化学势导数；
响应的伴随矩阵再进入同一 RVR 收缩，保留所有 v/B/N 反馈导数。
非原始采样点按 `T=S^(-1/2)`、`THT`、`Q=TU` 一致计算密度与响应。

八个 PBS worker 各绑定一个 kz，推荐沿用每节点 16 核及节点本地切片：

```sh
# 每个计算节点的 z 依次为 0..7，所有任务读取同一个外部配置。
/absolute/python 7_spin/production/complete_energy.py 0
# 全部任务完成后：
/absolute/python 7_spin/production/summarize_energy.py
```

固定 2×2×8 k 网格、Gaussian 0.01 Ry、内层五点导数步长 0.001 rad、
外层 ±0.25°/±0.5°/±1°，只计算最近邻 J12 及反向配对检查。
输出保留完整密度、谱、原生能量、逐 k 转矩及逐方向 J。
不同尺寸或更细采样的物理收敛不由完成标记自动保证。

## 验证与边界

`tests/test_batch_checkpoint.py` 检查断点复用、截断输出拒绝及完整批次发布。
`tests/test_production_energy.py` 检查密度/势伴随、FFT gauge、固定电子数响应、
非正交全复数 N 的完整能量导数、参考投影、密度记录和无损 kz 切片。
`tests/test_local_vfunc.py` 提供原生 H/dH 与独立实空间积分及批量/单点对照。
`tests/test_production_adjoint.py` 进一步对照显式 dH，并核对不同条目的导数排列。

2026-11-15 的生产记录中，Néel/FM 同规范 H/能带回归最大差约 5.2e−12 eV；
测试 RVR 上批量伴随与显式 dH 最大差约 2.7e−17 eV。
11 月 16 日完整能量 4×4/6×6 已有计算结果，但尺寸与有限角度检查尚不支持
收敛的独立 DFT 交换参数结论。这些历史记录与本次移植测试分别报告。

密度、伪势、模型及背景必须沿用一致来源；本次没有重拟合、重新 Wannier 化或新 SCF。

本次移植验证（2026-11-16）：

| 检查 | 本次结果 |
|---|---|
| `make local-rvr`，gfortran/Accelerate | PASS |
| `check.py` 全仓离线检查 | 24 通过，1 项 MPI 测试因缺少 MPI 编译器/mpirun 跳过；随后新增外部配置路径测试单独通过 |
| 原生 L=1/2/3，120 项 H/dH 分量对照 | 最大误差 1.353e−16 eV |
| 40 项批量/单点对照 | 最大差 0 |
| 真实 V1 反馈＋测试 RVR | H 最大差 9.156e−16 eV，dH 最大差 1.003e−18 eV |
| 28 项伴随/显式 dH 及混合导数检查 | 最大差 2.289e−17 eV |
| 外部配置、密度/势、固定电子数、完整泛函、记录和切片 | PASS |
| QE 适配器重新编译与材料八节点全流程 | NOT_RUN；`native_large.f90` 与已生产使用的源码逐字节一致 |

首次全仓检查发现切片模块缺少 `os` 导入，补齐后重跑通过。
这里的软件检查不替代大规模生产复跑或独立 DFT 物理验证。
