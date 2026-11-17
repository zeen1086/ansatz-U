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
│   └── complete_energy/     # 固定补空间密度、持久波函数/算符缓存
└── spins/                   # 构型输入、结果及任务状态
    ├── inputs/Lxx/          # structure.npz、CASE.json
    ├── H_case/              # 单构型 H/dH
    ├── regression/          # 指定参考自旋构型的 H/能带检查
    ├── J/Lxx/.../case_XX/   # 系数、本征态、化学势、带自由能/转矩
    └── complete_energy/     # 构型密度、完整能量、转矩、J 和检查点
```

这两个数据目录受 Git 忽略。可将它们链接到外部存储，或在配置中指定外部绝对路径；不要把大型计算数据推送到 GitHub。更换模型、几何、基底或参考密度时使用新的一组目录，不能仅因文件名相同而复用。

完整能量流程的波函数、重叠矩阵及变换、H0 和 projector 持久保存在 `shared/complete_energy/cache/`，不同 spin 共用。结构文件含自旋方向，因此整份文件放入 `spins/`。

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

按 下文批量运行说明 填写外部配置，在本地或 PBS 计算任务中运行。选择读出时区分：

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

原生组装将固定标量零点写成 `cI`。完整能量入口在折叠基底中替换为 `cS`，再用 `S^(-1/2)` 求解并保存修正后的 H 供转矩计算使用。用此版本复算旧结果时，配置新的 `complete_energy.work` 和 `complete_energy.shared_work`。原模型系数仍可通过 `reuse_cases` 提供。

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
- `interface`：含 `V_func_tmp/V_FUNC.json` 的目录；将模板中全零占位的 `model_id` 替换为 `V_FUNC.json` 的实际 ID。
- `input_root`：含 `inputs/Lxx/structure.npz`、`inputs/Lxx/CASE.json`。
  CASE 沿用生产配对、角度和站点定义；超胞及 Ni 顺序见 local 接口。
- `references`：共同 Wannier 的 `data.npz`、`k_0001.npz ...` 和 `wr_0001.npy ...`。
- `work`：共用数据目录 `shared/`，含 `RVR/B_general`、`RVR/PROJECTION.json`、
  `basis/pure_H0.au` 和 `basis/bloch.npy`。
- `rvr_slices`：可选的 `{"0":"/node-local/z0/B_general", ...}`。缓存路径在对应 worker
  节点上可见；必须与完整 RVR 的模型、基底、尺寸及 kz 匹配。未指定时从完整源截取。
- `spin_work`：构型结果根目录 `spins/`，保存 `J/` 和 `regression/`，`regression/RESULT.json` 是可选诊断报告，不作为计算前置条件。
- `old_J`：可选旧四 Ni 对照文件，文件不存在时不进行该比较。

完整能量的 `complete_energy` 配置另外需要：

- `work`：随构型变化的完整能量输出目录 `spins/complete_energy/`；每次更换模型、输入、角步或网格使用新目录。
- `shared_work`：共用的参考密度、波函数、重叠变换及算符缓存目录 `shared/complete_energy/`；更换模型或参考输入时也须使用新目录。
- `reference_run`：原 2×2 参考包，含 `results/BANDS.npz`、`results/FIELDS.npz`、
  `inputs/reference_actual/DENSITY.bin`。原生铺展回归还读取
  `results/native/reference_actual_ENERGY.dat` 与 `reference_actual_TRUE.bin`。
- `pseudo_dir`：相同来源的 `La.upf`、`Ni.upf`、`O.upf`。
- `native_binary`：下面编译得到的原生 QE 适配器。
- `sizes`：计算尺寸列表，仅接受 4、6、8；不自动增加尺寸或采样。
- `reuse_cases`：可选的 `{"4":"/existing/L04/mesh2_h0.001", ...}`，读取其
  `case_00 ... case_10/coeff`。这些系数必须带有本版本生成的 `CHECKPOINT.json`，
  其结构/spin、模型、求导位点、步长及文件校验须与当前计算匹配；外部目录只读。
  已有本征态若完整且通过残差检查，也可复用，否则重新对角化。
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

可在本地或已分配的 PBS 计算任务中运行。`run_material.py` 提供以下阶段接口，
按 kz 分组组装 H、求全部 k 点的公共化学势，再计算 HF 转矩：

```sh
/absolute/python 7_spin/run_material.py spectrum --L 4 --case 1 --mesh 2
/absolute/python 7_spin/run_material.py torque --L 4 --case 1 --mesh 2
/absolute/python 7_spin/run_material.py collect --L 4 --mesh 2
```

`collect` 需要 CASE 中全部规定构型；默认不会把单构型输出标为完整 J。
材料回归缺失、未通过或报告不匹配只输出 `QUALITY WARNING`，计算继续。`V1_K_WORKERS` 控制同 kz 批量原生组装的线程数，
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
这些阶段不负责提交或取消作业；未设置 PBS_JOBID 或 ANSATZ_RUN_CONFIG 时只打印 NOTE，后者缺失时读取仓库配置模板。模板中的实际数据路径仍需可访问。

## 完整能量读出

能量定义见上文。完整能量在非原始采样点使用重叠变换，按固定电子数计算密度响应。

八个 worker 各绑定一个 kz，推荐沿用每节点 16 核及节点本地切片：

```sh
# 可选：在统一启动脚本中指定同一个编号，传给全部八个 worker。
export ANSATZ_RUN_ID="$(/absolute/python -c 'import uuid; print(uuid.uuid4().hex)')"
# 每个计算节点的 z 依次为 0..7，使用相同 ANSATZ_RUN_CONFIG 和 ANSATZ_RUN_ID。
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

### 中断后重跑

1. 确认上一轮全部 worker 已停止，保持输入、配置、模型与可执行文件不变。
2. 带自由能流程重跑对应阶段命令。完整能量未设置 `ANSATZ_RUN_ID` 时，按 `PBS_JOBID`（本地按父进程 PID）在 `run_ids/<key>` 加锁生成并复用编号。同一启动器的八个 worker 使用同一编号；不同 PBS 作业或不同父进程共同运行时需显式指定相同编号。
3. 已绑定且完整的检查点直接复用。损坏/缺失的谱、密度或逐 k 转矩会重算；未完成的系数、Hamiltonian 或 QE 目录留为 `*.incomplete-*` 后重新生成。
4. 查看 `attempts/<ANSATZ_RUN_ID>/FAILURE.json` 和作业日志。失败标记只打印 NOTE，等待没有硬超时；同编号启动标记可覆盖，活动 worker 仍由 flock 防并发。需要隔离旧 barrier 时显式指定新编号；重用编号不会自动清除旧 barrier。

结果文件先写临时文件，再原子替换。QE 阶段使用全新临时目录，输出验证后发布；若原生阶段已完成而能量 JSON 未写完，重跑时直接补写摘要。带自由能流程只在整个构型的 `TORQUE.json` 发布后清理本征态等临时数据。

`INPUTS.json` 绑定结构/spin、CASE、模型、基底、网格/步长及相关代码；系数和原生输出另外记录输入与文件校验。缺少 `INPUTS.json` 的已有目录通常会打印 NOTE 并写入当前身份。完整能量目录中已有构型结果时，须有 `INPUTS.json` 才能续跑。已有 manifest 与当前身份不同仍报错，需要使用相应输入或新的目录。文件身份采用大小/时间戳，结构/spin 使用实际数组描述，不记录或校验文件内容摘要。旧哈希身份与新记录不同也不会被静默覆盖。

失败留存的 `*.incomplete-*`、进程被强制终止留下的隐藏临时目录，不是可复用检查点。确认没有活动写入者后可自行清理，正常重跑无需先删除它们。运行过程中保持输入和缓存不可修改。

未设置 `spin_work` 的配置仍沿用 `work` 保存构型输出；完整能量请显式设置独立的 `shared_work` 和 `work`，分别保存共用缓存与本轮结果。

## 小型检查与原生回归

```sh
/absolute/python 7_spin/complete_readout.py
/absolute/python -B check.py
```

`complete_readout.py` 包含能量、Hubbard 项和固定电子数响应的有限差分检查；
`density_map.py` 保留独立参考投影算法，由完整能量初始化及离线测试调用。
在配置好参考密度、QE 适配器并分配 PBS 计算任务后，可运行原生铺展密度检查：

```sh
/absolute/python 7_spin/check_large_native.py
```

结果写入 `complete_energy.shared_work/native_check/`，目录须不存在。
这些检查不替代材料收敛和独立 DFT 精度验证。

## 完整能量的共享缓存

先重新编译 `make local-rvr BUILD=/absolute/build`，在配置中指定
`complete_energy.shared_work`（推荐 `7_spin/shared/complete_energy`）。运行
`complete_energy.py 0` 至 `complete_energy.py 7` 时自动准备并复用缓存，不需单独执行缓存命令。

```text
shared/complete_energy/
├── INPUTS.json                 # 输入内容及计算代码的身份记录
├── reference/                  # 固定参考补空间
└── cache/
    ├── Lxx/kXXXX/
    │   ├── wave.npy            # 插值波函数
    │   ├── overlap.npy         # S
    │   ├── transform.npy       # S^(-1/2)
    │   ├── metric.json
    │   └── CACHE.json
    ├── Lxx/operators_<k>/
    │   ├── h0.npy              # 折叠基底 H0
    │   ├── projectors.npy      # 折叠基底 projector
    │   ├── native_projectors.npy
    │   ├── native.bin          # 原生 H 组装用 H0/projector
    │   └── CACHE.json
    └── z*/RVR.bin              # 未指定外部切片时生成
```

首次使用按 L、k 构建；后续 spin 的谱计算复用重叠变换，密度与响应读出复用
H0/projector，原生 H 组装和转矩收缩也读取对应基底的缓存。H、模型系数、
本征态、化学势、密度、能量和转矩仍按 spin 重算并保存到 `spins/`。
单构型 `local_vfunc.py cal_Hsk` 和独立带自由能入口保留原来的直接组装方式。

输入身份使用文件大小、mtime/ctime 及已有 model/basis ID，`cache_identity()` 不再读取整份输入计算哈希。共享缓存复用检查身份、文件集合与字节数，不保留内容摘要字段、空字段或摘要清单。模型包 `files` 记录文件长度，`verify` 检查文件存在和长度；新模型/基底导出使用随机的 64 字符 ID，原生二进制头部长度不变。

新建 `.npy` 缓存不超过 256 MiB 时检查全部行，超过时只检查首行。输入与缓存应在运行期间保持不变；这些检查不能发现所有等长内容损坏。

缓存按条目加文件锁，写到临时目录后整体原子发布。中断留下的
`*.partial-*` 目录不会被读取，重试会重新生成未发布的条目；确认没有写入
任务后可删除残留临时目录。正式缓存缺文件、字节数不符或身份不符时会报错，
不会静默覆盖。损坏条目先归档，再在相同输入下重建；输入或源码改变时用新的
`shared_work` 和构型输出目录。无 `INPUTS.json` 的旧参考/缓存目录会打印 NOTE 并登记当前身份；已有身份不匹配仍报错。

缓存不会在计算完成后自动删除。32 个 k 点的波函数本身约占
`32 × L² × 57600 × 30 × 16` 字节：L=4/6/8 分别约 14.2/31.9/56.6 GB，
另需重叠和算符矩阵空间。构建波函数前检查当前可用空间；并发任务仍需预留
总量。只更换 spin 时复用同一套缓存，无需为每个构型复制。

## 运行检查

材料回归、有限的电子数偏差和较小的正定重叠本征值只给出质量提示，不阻止运行。电子数偏差写入结果文件或运行日志；运行完成不代表这些诊断通过。

文件缺失或损坏、维度/格式错误、模型与缓存身份不匹配、非有限数值、非正定重叠矩阵仍属于功能错误。遇到这些错误须修复输入或重新生成对应结果，不能把错误缓存作为有效检查点继续使用。

Hamiltonian 的 Hermiticity 残差超过 1e-5、参考投影残差超过 1e-6 只给诊断提示；伴随/直接导数残差超过 1e-6 仍报错。

重新导出模型或重建投影会生成新的 ID；更新配置并使用对应的新 RVR，不要混用不同导出的文件。旧格式元数据可读取，但不会在新生成的记录中复制内容摘要。
