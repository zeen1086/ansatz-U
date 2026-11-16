# V1 批量能量与转矩使用说明

## 输入与构建

使用已有 Python 3.11+、NumPy、SciPy 环境。原生局域程序需要 Fortran 与 BLAS/LAPACK。
完整能量还需要与参考数据相符、已编译的 QE 7.6 环境、赝势和参考密度。

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
- `sizes`：计算尺寸列表，仅接受 4、6、8；不自动增加尺寸或采样。
- `reuse_cases`：可选的 `{"4":"/existing/L04/mesh2_h0.001", ...}`，读取其
  `case_00 ... case_10/coeff`。已有本征态若存在并通过残差检查，也可复用。
  省略时在新输出目录生成这些构型的完整 V1 系数。

本实现固定 LaNiO2 的原胞 40×40×36 网格、30 空间轨道、完整 8³ WS 数据、
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

## 批量带自由能读出

在已分配的 PBS 计算任务中运行。`run_material.py` 提供以下阶段接口，
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
这些阶段需在已有 PBS 计算任务中执行，不负责提交或取消作业。

## 完整能量读出

能量定义与单位见 [7_spin 使用说明](../README.md)。这是不同于带自由能的读出。在固定参考补空间密度下，忽略构型无关常数：

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

主要结果：带自由能保存在各构型的 `TORQUE.json/F_eV`；完整自由能保存在
`Lxx/case_XX/ENERGY.json`，J12 保存在 `Lxx/RESULT.json`。
运行汇总后生成 `J12_by_size.csv`、各尺寸的 `energy_components.csv`
和 `ENERGY_TORQUE_CHECK.json`。后者比较单站点能量曲率与转矩导数，不是 J12。

带自由能分阶段入口可重跑相同命令复用检查点。完整能量入口包含按文件存在性跳过的阶段，
不应把残缺文件当作有效检查点；异常退出后先检查失败阶段及 `FAILURE.json`，
不要直接复用损坏目录。更换输入、模型、网格或步长时使用新输出目录。
