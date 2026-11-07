# ansatz+U — LaNiO₂ Fortran workflow

这是独立的 ansatz+U 研究源码项目，GitHub 仓库名为 `ansatz-U`。`physical_release=false`。原始源码、冻结模型和历史结果均保留在目录外；不在本目录运行计算或保存编译产物。本目录没有自动提交作业或启动训练的入口。

生产数值均由 Fortran 执行：空间插值/RVR、H_sk、Hubbard 组装、完整复数 Hamiltonian 对角化、固定电子数 FD 求 EF、B/U 转矩和单 Ni 倾斜 Φ。Python 仅用于一次性存储格式转换与目录外的独立核对。

## 目录与阶段

按 [zeen1086/ansatz](https://github.com/zeen1086/ansatz) 的编号和阶段边界整理；参考树为 `05fa36cbf6fb7740bd55293466e79fa13af545c7`。这里只迁移已有 LaNiO₂ 代码，不引入参考仓库的 CrI₃ 模型、参数系统或 Python 运行器。

```text
lanio2_ansatzU_fortran/
├── 0_QE/                 # 外部 QE 程序来源说明
├── 1_scf/                # 已有训练场和保存态来源说明
├── 2_band/               # 已有磁性 DFT bands 来源说明
├── 3_band_no_mag/        # 已有非磁背景来源说明
├── 4_potential/fit/      # 对角/全分量拟合、V_func_tmp 重启
├── 5_wannier/            # 外部同 gauge Wannier 数据契约
├── 6_RVR/                # RVR、冻结数据读取/转换、背景 HR 重排
├── 7_spin/               # H_sk、ansatz+U、能带、转矩和 Phi
│   └── AFM/              # 保留的旧 AFM 源码入口
├── common/               # real64、BLAS/LAPACK、FD 数目求解
├── examples/             # THE.in 与 4_potential/、7_spin/ 文本模板
├── Makefile
├── PROVENANCE.md
└── VALIDATION.md
```

`0_QE` 至 `3_band_no_mag` 以及 `5_wannier` 仅记录外部依赖，不包含新 DFT/Wannier 生产程序。实际源码从 `4_potential` 开始。原 `src/` 层已移除；所有编译目标名、数值模型、输入格式和外部模型路径保留。完整旧目录小型快照和移动表保存在维护者的外部验证归档中，不纳入 Git 仓库。

## 构建和最小运行

需要 gfortran、双精度 BLAS/LAPACK。默认 macOS 链接 Accelerate，Linux 链接 `-llapack -lblas`。没有例行安装依赖。

```sh
make -C /path/to/lanio2_ansatzU_fortran BUILD=/absolute/external/build
# 需要重新投影或旧 MPI 能带/能量程序时另行编译：
make -C /path/to/lanio2_ansatzU_fortran mpi BUILD=/absolute/external/build MPIFC=mpif90
```

ISSP 已有 MKL 的链接参数如下，使用执行主机上的编译器和 MPI：

```sh
make -C /path/to/lanio2_ansatzU_fortran BUILD=/absolute/external/build \
 LIBS='-L/opt/intel/oneapi/mkl/latest/lib/intel64 -Wl,-rpath,/opt/intel/oneapi/mkl/latest/lib/intel64 -lmkl_rt -lpthread -lm -ldl'
```

先将 `examples/7_spin/run.nml`、`texture.dat`、`kpoints.dat`、`readouts.dat` 复制到独立运行目录，并填写显式数据路径。所有路径相对**调用进程的工作目录**解析；建议使用绝对路径。输出父目录须已存在，输出必须使用未占用的前缀。

```sh
/absolute/external/build/ansatzU /absolute/external/run/run.nml
```

`mode='torque', delta=0` 读取任意给定纹理，不附加扰动。`mode='table3'` 在输入纹理及指定 Ni 沿其给定轴旋转 `delta` 后分别计算，输出差分 Φ；示例为 Néel 单 Ni 绕 +y 倾斜 1°。`mode='torque', delta!=0` 可用于显式旋转和有限差分。生产接口不限于该单倾斜纹理。

输入文件：

- `texture.dat`：L² 行，逐行 `ex ey ez axis_x axis_y axis_z`；Ni 顺序 x 外层、y 内层，y 最快。方向和转矩轴都要求单位范数，程序不做隐式归一化。
- `kpoints.dat`：首行点数，后续 `Kx Ky Kz weight`。K 是超胞面内/primitive z 的约化坐标，权重非负且总和为 1；支持非均匀权重。能带路径也可用此输入，但路径点上的 EF 不是三维布里渊区积分 EF。
- `readouts.dat`：首行读出数，后续零基 Ni 坐标 `x y`。轴取该 Ni 的输入轴。所有读出使用同一个完整 H 的本征态和 FD 占据。
- `run.nml`：`ueff` 默认 1 eV、`sigma` 默认 0.05 eV、`electrons_per_ni` 默认 21；均实际参与计算。电子数为 `electrons_per_ni*L*L`，没有读取后覆盖为常数。旧 THE.in 的 46 不是这个活性电子数。

输出 `.p1.bands`（以及 table3 的 `.p2.bands`）、`.result`、`.phi` 均为文本双精度；`dump_matrices=true` 额外输出小型逐分量矩阵，为外部核对而设，维数限 600。矩阵流为 little-endian complex128、列优先。已有输出会报错，不覆盖。

## 冻结数据：只转换一次

完整模型在外部显式提供，目录内不含模型、波函数或大投影缓存。冻结源 NPZ 的身份为：

`6fe8cb35d6271638fea0fe9227456a4f1b04910fb8cf181fe30cc212420c2ba0`

```sh
/path/to/existing/python 6_RVR/export_frozen.py /external/local_model.npz /external/local_model.au
```

工具核验冻结 SHA、维度、复数精度和写后逐元素读回；若目标存在则停止。生成外部 `.au.json` 格式证书。Fortran 此后直接读取 `.au`，每次计算不调用 Python。ISSP 使用已有 `cri3-edge66` Python 环境进行转换，无需安装包。完整数据在 ISSP 原位处理。

固定二进制 V1：16 字节 magic `LANIO2_AU_V1` 后补四空格，8 个 little-endian int32 为 `0x01020304,729,30,5,10,810,1,0`（最后两项为 eV、背景不含 HU）。随后按 Fortran 列优先连续写入：

| 数组 | 形状 | 类型/含义 |
|---|---|---|
| rr | (3,729) | int32，WS [-4,4]³，原插值权重已含在系数中 |
| hc | (30,30,729) | complex128，H0+固定标量 Hv，eV |
| p | (5,30,729) | complex128，原 Hubbard projector 映射，无量纲 |
| KB | (810,810) | complex128，27 平移 × 30 轨道，局域 B 空间核，eV |
| N0 | (10,10) | complex128，无量纲完整参考占据 |

复数为相邻 real64 实/虚部，总长度严格核验。projector 占据顺序为 `[up d1..d5, down d1..d5]`；Hamiltonian 顺序为 `(cell,orbital,up/down)`，自旋交错。KB 平移按 x,y,z（z 最快）各 -1..1，轨道最快。原 Wannier gauge/projector 全部保留。

**H0 与 Hv 的独立数值分解不在所给冻结 NPZ 中**，因此本接口明确报告 `H0+Hv` 合并背景，不伪造两者的单独输出。上游原算法中它们分别来自 nonmagnetic Wannier 背景及固定标量势投影。此冻结背景无 HU；B/HU 分别组装，Hubbard 只加一次。若输入另一个已含 HU 的背景，必须先在其来源处完成有据可查的扣除，不能将其标作本格式的无 HU 背景。

## 原 Fortran 阶段与全分量张量接口

保留主线 `fit → V_func_tmp → RVR → H_sk → 合并背景 → 对角化/读出`。冻结局域核已完成 fit/空间投影，日常运行从已冻结 RVR 进入，尺寸变化不重拟合、不重 Wannier 化。

`fit_new` 使用 `examples/4_potential/fit.nml`。`diagonal_iters` 指定先只更新 B_general 的自旋分量对角项，`full_iters` 指定随后更新全部九个分量；`lambda_DM`、迭代数必须由用户显式给定，模板为不执行训练的占位值。阶段切换不清空非对角项。可设 `diagonal_iters=0, full_iters>0`，从已有 `restart_file` 直接继续全部分量更新。输出为五个原顺序 sequential unformatted 记录：VA_ave、V0、Vn、B_general、JH；`restart_file=''` 为新初始化，输出必须是新路径。标量/旧模型不在本任务重训。

原拟合是逐样本 B 梯度更新；保留原算法，不擅自填历史迭代数、学习率或训练结果。有限核查发现归档 CrI₃ `fit.f90` 中存在取消对角限制的全分量版本；没有找到可确认身份的 LaNiO₂ 后续训练源码/参数。来源表记录此边界。

旧阶段使用显式 THE.in 参数；其内部仍按原来的阶段工作目录布局写入**外部运行树**，不应在源码目录启动：

```text
/external/run/THE.in
/external/run/6_RVR/tmp/
/external/run/6_RVR/super_H/
/external/run/7_spin/data/
/external/run/7_spin/tot_HR/
```

THE.in 的 `func_dir/nscf_dir/THE_dir` 是显式外部目录，保留末尾 `/`。拟合输入为 `fields_dir/spinN/V_HYC_dp` 和显式 spin_info。RVR 的网格参考另在 `func_dir/reference_V_HYC_dp` 提供，不再偷偷读相邻旧 QE 目录。只需要投影时使用原 Wannier `make_wan/wan_data/wan_*` 与 `cal_TF_*`；已有合格核不重做此步。

```sh
/external/build/fit_new /external/run/fit.nml                 # 仅在明确要训练时执行
cd /external/run/6_RVR
mpirun -np 2 /external/build/cal_RVR_split /external/run/THE.in 1
# 按需要对每个 spinor-Wannier column 2..60 运行；至少两个 MPI ranks。
/external/build/rvr_combine /external/run/THE.in
/external/build/rearrange_HR /external/run/THE.in
cd /external/run/7_spin
/external/build/tensor_Hsk /external/run/THE.in /external/run/expanded_texture.dat
/external/build/combine /external/run/THE.in
```

`tensor_Hsk_AFM` 保留同样的全分量接口，两者现在均从显式 texture 读取，不生成内置 FM/Néel，不隐含乘 0.372/0.4。旧纹理文件每行三项，顺序 x,y,z（z 最快），需覆盖 `super_H/ir_vec_HYC` 所定义的展开区域；这些数值**直接是原张量系数约定下的自旋分量**。`examples/7_spin/legacy_texture.dat` 仅展示幅值 0.372 约定；0.4 也是另一个历史约定，不能混用。

要把原张量 RVR 与新的 HU/FD/转矩连起来，使用 `examples/7_spin/tensor_run.nml`：`tensor_path` 指向原 `rvr_combine` 输出的 full-component `B_general`，`tensor_unit` 必须显式声明（原 Fortran 势/RVR 为 Ry，适配器按统一常量转 eV；已是 eV 的核声明 eV），`tensor_range` 是原对称投影范围，`tensor_spin_scale` 必须显式给出（例如对应核的 0.372）。它**替换**局域 KB 的 HB，不叠加第二个 B；共用同一外部 hc/p/N0、同一 gauge 和同一 H 的本征态。输入方向仍为单位向量，幅值仅乘一次。此适配器支持当前 LaNiO₂ 的 `n_func=3,nwan=30,nspin=1`，全部张量分量均保留。`B_general` 必须是 little-endian complex128、单条记录、两端4字节正 int32长度标记的格式；程序核验总长度及两端标记，超过单记录2 GiB的布局不接受。其它 nwan/n_dep 的旧模型仍使用原阶段接口，不能标作已验证的该冻结模型。

`combine` 的旧输出已改为 complex128 stream，H_sk 的 Ry 按 `13.605693122994 eV/Ry` 转换。新的 `band`/`calc_tot_ene` 读同一精度并调用 ZHEEV；MPI 数据为 DOUBLE_PRECISION/DOUBLE_COMPLEX。旧单精度文件不能直接混入新读入链。`calc_tot_ene` 保留原 sharp-filling 带能算法，THE.in 的 `n_electron` 为总活性电子数且真正使用；其结果不是完整 DFT 总能。当前 FD 主线使用 `ansatzU`。

## 模型和转矩定义

`H = (H0+Hv) + HB + Σ A_R† Ueff(I10/2−N_R) A_R`，所有 H、B、U、sigma 为 eV。N 保留完整复数矩阵及轨道/自旋非对角元，不裁剪本征值。

初始方向 e 对应 `D=Rz(phi) Ry(theta)`，随后扰动用 `exp[−i delta(axis·sigma)/2] D` 左乘；因此方向和占据使用同一旋转。`N_R=D N0 D†`，导数为 `[G,N_R]`，`G=−i(axis·sigma)/2 ⊗ I5`。这是声明的刚性局域模型，不将旧四自旋非刚性 N 关系外推到任意超胞。保持无 SOC、原 projector/Wannier gauge、局域 B 的 [-1,1]³ 范围与 scalar/projector 的原 WS 插值。

`∂HU/∂theta = −Ueff A†[G,N]A`；局域 B 的导数由 `axis × e` 得到，张量 B 则对原三个系数逐项求导。`tau_B/U=−Σ w f_FD <∂HB/HU>`，总转矩为两者之和，单位 eV/rad。Phi 为 `[tau(delta)−tau(0)]/delta`，输出 meV/rad²。表3历史表头单位与此导数定义的差异不用于重标度。两种读出是同一 H 的导数贡献，不能称为 U 开/关效应。

局域 KB 已含 B 幅值，不能额外乘 0.372、0.4 或 S²。程序的带能求和是 `Σ w f E`；FD 转矩对应固定电子数的独立粒子 Helmholtz 自由能导数，不是有限温度带能和的全导数，也不是 DFT 总能量。

使用及验证结论见 `VALIDATION.md`，来源与功能对应见 `PROVENANCE.md`。模型、测试代码与数值输出保存在维护者的外部验证归档中；本仓库不分发大模型或历史作业数据。
