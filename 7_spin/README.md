# V1 Hamiltonian 与能量使用说明

## 目录与数据放置

所有 Python、Fortran 源码和构建文件直接放在 `7_spin/`，运行数据分为：

```text
7_spin/
├── *.py、*.f90、Makefile.native、README.md
├── shared/                  # 同一模型/几何/基底下，不随 spin 改变
│   ├── model/V_func_tmp/
│   ├── RVR/                 # B_general、PROJECTION.json 等
│   ├── basis/               # pure_H0.au、bloch.npy 等
│   ├── wannier/references/
│   ├── complete_reference/、pseudopotentials/
│   └── complete_energy/     # 固定补空间密度、重叠检查、切片来源记录
└── spins/                   # 构型输入、结果及任务状态
    ├── inputs/Lxx/          # structure.npz、CASE.json
    ├── H_case/              # 单构型 H/dH
    ├── regression/          # 指定参考自旋构型的 H/能带检查
    ├── J/Lxx/.../case_XX/   # 系数、本征态、化学势、带自由能/转矩
    └── complete_energy/     # 构型密度、完整能量、转矩、J 和检查点
```

这两个数据目录受 Git 忽略。可将它们链接到外部存储，或在配置中指定外部绝对路径；不要把大型计算数据推送到 GitHub。更换模型、几何、基底或参考密度时使用新的一组目录，不能仅因文件名相同而复用。

节点本地的波函数、重叠变换和 kz 切片仍放在任务临时缓存，避免共享存储重复 I/O；它们属于共用数据，不是构型结果。结构文件含自旋方向，因此整份文件放入 `spins/`。

## 单构型 Hamiltonian

先按 [RVR 使用说明](../6_RVR/README.md) 准备模型、`B_general` 和 `pure_H0.au`，在仓库根目录编译并运行：

```sh
make local-rvr BUILD=/absolute/build
/absolute/python local_vfunc.py cal_Hsk \
  --model /absolute/model/V_func_tmp \
  --structure /absolute/L04_structure.npz \
  --rvr /absolute/RVR/B_general \
  --background /absolute/basis/pure_H0.au \
  --kpoints /absolute/kpoints.dat \
  --binary /absolute/build/local/cal_Hsk \
  --output /absolute/repo/7_spin/spins/H_case
```

`cal_Hsk.f90` 组装 `H=H0+Hv+HB+HU`。默认输出逐 k 的 H 和 `HAMILTONIAN.json`；添加 `--derivative-sites 1 2 --step 0.001` 计算角度导数，`--components` 输出分量，`--keep-intermediates` 保留调试文件。使用新的输出目录。

## 批量能量、转矩与 J

按 下文批量运行说明 填写外部配置，并在 PBS 计算任务中运行。选择读出时区分：

| 入口 | 能量与输出 |
|---|---|
| `run_material.py`、`kz_production.py` | Gaussian 带自由能 `Fband`；`TORQUE.json` 中的 `F_eV`，对应转矩 `−dFband/dθ` |
| `complete_energy.py` | 固定补空间密度下的完整试探密度自由能；`ENERGY.json` 中的 `energy`，对应转矩 `−dF/dθ` |
| `summarize_energy.py` | 汇总 J12，并比较完整能量与转矩的单站点曲率；该曲率不是 J12 |

两种自由能均使用固定电子数、全 k 公共化学势和 Gaussian 展宽 `σ=0.01 Ry`：

```text
f = erfc((ε−μ)/σ)/2
Fband = Σk,n wk [ε f − σ exp(−((ε−μ)/σ)²)/(2√π)]
F = Fband − Σk wk Tr[P(H−H0)]
    + EH[ρtotal] + Exc[ρtotal] − Exc_core
    + ∫ n_total vion + Ueff/2 Tr[Ntotal−Ntotal²]
```

`P` 是活性空间密度矩阵，`ρtotal/Ntotal` 加入固定参考补空间。输入势扣除项包含 v/B/HU；Hubbard 泛函只加一次。完整能量省略构型无关参考常数，只比较相同尺寸、参考数据和能量规范下的构型差。能量为每个超胞的 eV，转矩为 eV/rad，J 为 meV；两种读出的能量和 J 不应混用。

完整能量采用 Wannier 重叠度量；带自由能入口直接对角化 H，要求输入基底正交归一。非正交插值点不能仅凭相同展宽参数将两种入口的带项直接等同。

## 独立能带工具

`band.f90` 读取已有 `tot_HR/hr_*` 和平移表，仍用于 Wannier 流程检查，不读取 V1 的逐 k `.h` 文件。使用匹配的 `parameter`（TOML）、`plt_path.in` 和已有 MPI 环境：

```sh
ANSATZ_BUILD_DIR=/absolute/band-build /absolute/python build.py \
  7_spin/band.f90 --fc mpifort
cd /absolute/band_run
mpirun -np 2 /absolute/band-build/bin/7_spin/band.x
```

在运行目录准备 `tot_HR/` 和 `plt_path.in`，后者可按本目录模板填写；平移表位置由 `super_h_dir` 指定。输出 `plt_bnd_<title>`。

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

V1 模型导出见 [local RVR 接口](../6_RVR/README.md)。本仓库的
`4_potential/local/` 保存逐字节复制的六个冻结推理源，可作为导出命令的 `--code`。
导出后的 `interface/V_func_tmp` 含这些源及全部冻结模型参数，运行不依赖训练目录。

独立 `local_vfunc.py evaluate` 默认只生成生产所需的 `value.bin/dsite_XXXX.bin` 和 `PREDICTION.json`；重复存储这些数组的 `prediction.npz` 仅在传入 `--keep-intermediates` 时生成。批量驱动仍保留实际续算所需的系数、索引和完成标记，不使用单构型 `cal_Hsk` 命令完成后的临时文件清理。

复制 [配置模板](../examples/7_spin/production.json) 到外部运行目录，替换所有绝对路径：

```sh
export ANSATZ_RUN_CONFIG=/absolute/run/production.json
```

配置中的路径及数据约定：

- `python`：执行主机已有解释器；`native_build`：含 `local/cal_Hsk` 的构建根目录。
- `interface`：含 `V_func_tmp/V_FUNC.json` 的目录；模型及参数 ID 必须与模板中的冻结 V1 一致。
- `input_root`：含 `inputs/Lxx/structure.npz`、`inputs/Lxx/CASE.json`。
  CASE 沿用生产配对、角度和站点定义；超胞及 Ni 顺序见 local 接口。
- `references`：共同 Wannier 的 `data.npz`、`k_0001.npz ...` 和 `wr_0001.npy ...`。
- `work`：共用数据目录 `shared/`，含 `RVR/B_general`、`RVR/PROJECTION.json`、
  `basis/pure_H0.au` 和 `basis/bloch.npy`。
- `rvr_slices`：可选的 `{"0":"/node-local/z0/B_general", ...}`。缓存路径在对应 worker
  节点上可见；必须与完整 RVR 的模型、基底、尺寸及 kz 匹配。未指定时从完整源截取。
- `spin_work`：构型结果根目录 `spins/`，保存 `J/` 和 `regression/`，计算前需有通过的 `regression/RESULT.json`。
- `old_J`：可选旧四 Ni 对照文件，文件不存在时不进行该比较。

完整能量的 `complete_energy` 配置另外需要：

- `work`：随构型变化的完整能量输出目录 `spins/complete_energy/`；每次更换模型、输入、角步或网格使用新目录。
- `shared_work`：同次计算共用的参考密度、重叠检查等输出目录 `shared/complete_energy/`；更换模型或参考输入时也须使用新目录。
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
  -f /absolute/repo/7_spin/Makefile.native \
  TASK_DIR=/absolute/qe_adapter_build native
```

该 make 配方沿用生产 Intel 编译器的 `-heap-arrays -check bounds`；需使用相同的
QE `make.inc` 和静态库，不能与其他版本的模块文件混编。

## 批量带自由能读出

在已分配的 PBS 计算任务中运行。`run_material.py` 提供以下阶段接口，
按 kz 分组组装 H、求全部 k 点的公共化学势，再计算 HF 转矩：

```sh
/absolute/python 7_spin/run_material.py spectrum --L 4 --case 1 --mesh 2
/absolute/python 7_spin/run_material.py torque --L 4 --case 1 --mesh 2
/absolute/python 7_spin/run_material.py collect --L 4 --mesh 2
```

`collect` 需要 CASE 中全部规定构型；默认不会把单构型输出标为完整 J。
生产计算要求已有材料回归通过。`V1_K_WORKERS` 控制同 kz 批量原生组装的线程数，
`K_WORKERS` 控制批量收缩的并行度；数学库线程数按实际节点配置。

跨构型的本地切片计算使用显式阶段，任务文件为
[批量示例](../examples/7_spin/batch_tasks.json) 中的 `[L, mesh, step_rad, case]` 列表：

```sh
/absolute/python 7_spin/kz_production.py stage --z 0 --rvr /local/z0.bin
/absolute/python 7_spin/kz_production.py prepare --z 0 --rvr /local/z0.bin --tasks /run/tasks.json
/absolute/python 7_spin/kz_production.py spectrum --z 0 --rvr /local/z0.bin --tasks /run/tasks.json
```

`prepare` 只由一个任务执行；`spectrum` 对 z=0..7 各执行一次。
所有谱完成后，一个任务执行 `spectrum-finish` 求公共化学势，随后 z=0..7
分别执行 `torque`，全部结束后一个任务执行 `torque-finish`。
这些阶段需在已有 PBS 计算任务中执行，不负责提交或取消作业。

## 完整能量读出

能量定义见上文。完整能量在非原始采样点使用重叠变换，按固定电子数计算密度响应。

八个 PBS worker 各绑定一个 kz，推荐沿用每节点 16 核及节点本地切片：

```sh
# 每个计算节点的 z 依次为 0..7，所有任务读取同一个外部配置。
/absolute/python 7_spin/complete_energy.py 0
# 全部任务完成后：
/absolute/python 7_spin/summarize_energy.py
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

未设置 `spin_work` 的旧配置沿用 `work` 保存构型输出；未设置 `complete_energy.shared_work` 时沿用完整能量的 `work`。已有任务可保留旧配置续算；新任务使用模板中的分离路径。
