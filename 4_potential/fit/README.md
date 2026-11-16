# 四自旋三势拟合

生产链仅用 Fortran、shell 和 BLAS/LAPACK：

```sh
make joint BUILD=/absolute/build
4_potential/fit/run_joint.sh /absolute/states.tsv \
  examples/4_potential/joint_fit.nml /absolute/new_run /absolute/build
# 先在当前目录保存仓库绝对路径，再进入外部运行目录求值：
repo="$PWD"
cd /absolute/new_run
/absolute/build/eval_joint "$repo/examples/4_potential/joint_eval.nml" \
  V_func_tmp "$repo/examples/4_potential/joint_texture.dat"
```

`RUN`、`BUILD` 必须在源码树外，运行目录须不存在。已有输出不覆盖。
没有 Python 生产依赖、自动安装、DFT/Wannier 重算或作业提交。

## 输入与读取

`states.tsv` 每行五个 TAB 分隔字段：`name weight anchor unit archive_path`。
路径相对清单文件或为绝对路径；`unit=Ry/eV`，权重非负。
`anchor=1` 是精确约束；`weight=0,anchor=0` 是只报告残差的留出态。
示例在 `examples/4_potential/joint_archives.tsv`。

shell 用 `unzip -p` 按固定成员名读取已有 DFT 标签 NPZ，不计算或转换数值：

| 成员 | dtype / shape / 约定 |
|---|---|
| `grid.npy` | int32/int64，`(3,)`；展平实空间网格 x 最快 |
| `atom_ids.npy` | int32/int64，`(4,)`；默认 Ni5,6,7,8 |
| `v_scalar_Ry.npy`, `B_material_Ry.npy` | float64，`(ngrid,)`、`(ngrid,3)` |
| eV 标签的对应成员 | `v_scalar_eV.npy`, `B_material_eV.npy` |
| `moments_muB.npy` | float64，`(4,3)`；显式按配置归一化方向，原幅值写入报告 |
| `N.npy` | complex128，`(4,10,10)`；每Ni `[up d1..5,down d1..5]` |
| 可选 `V_U_eV.npy` | 同 N；校验 `Ueff*(I/2-N)`，不叠加进标量/B标签 |

Fortran 原生解析 NPY v1/v2/v3、C/Fortran 顺序和 little-endian 数据；对象、单精度、截断、尺寸/原子顺序错配、非有限值及非 Hermitian N 均报错。
归一化后使用固定长度的自旋方向；DFT幅值已包含在标签和拟合系数中，不再乘0.372、0.4或S²，也未引入独立的纵向幅值自由度。
也可直接调用 `fit_joint fit.nml dataset.manifest [V_func_tmp [fit_report.txt]]`。
省略后两个参数时，在调用目录新建 `V_func_tmp` 和 `fit_report.txt`；也可显式指定外部路径。
manifest 首行态数，每态两行：`name weight anchor unit`，以及显式解包目录；目录内为 `grid/atom_ids/scalar/B/moments/N.npy` 和可选 `VU.npy`。
shell 直接写入 `dataset.manifest`，不再生成 `rows.manifest`；既有运行结果不改动。

## 模型与拟合目标

- 标量：原 scalar8 特征，即常数、六个 `ei·ej-1` 和 `alpha*(1-alpha)`。
- B：原12组基础矢量和28组非线性残差矢量，三分量同时拟合。
- N：原 `scalar/chi/vector/cross` 特征与实对称/虚反对称 d 轨道矩阵基，重构完整复数10×10矩阵；四个Ni共同依赖全部四个输入自旋，允许占据本征值变化。
- 沿用原64个方形四Ni晶体作用、Ni顺序和QE实d轨道约定；全局自旋 SO(3)/SU(2) 和时间反演性质由特征保持。
- 网格 `nx=ny` 可被4整除、`nz` 为偶数，Ni位于原四分之一晶胞位置；不自动将四输入非刚性N推广到任意大超胞。

目标为显式态权重下的原始标量场、B场和N矩阵三个最小二乘块。
它采用旧 Python 的特征族，**不是旧投影算符联合损失或其冻结训练结果的复现**。
SVD 直接求解，硬锚点用零空间消元；报告有效相对截断、秩、奇异值及逐态残差。
`fit_report.txt` 集中保留拟合参数、态权重与锚点、秩/奇异值、逐态残差和占据检查。
未识别方向取声明的 SVD 最小范数，不静默增加正则项或重调权重。
`remove_scalar_mean` 默认关闭；开启时每态被移除的均值明确记入报告，不能与另一势零点直接比较。
原 `fit_new` 的对角/全分量迭代和重启接口保留；新 SVD 入口不使用迭代次数或学习率。

新版 `V_func_tmp`（原名 `model.joint`）保留新 schema：16字节 `LANIO2_JOINT_V1 `、12个int32、两个float64、两个int32标志，接着列优先float64数组 `(8,ngrid)`、`(40,ngrid)`、`(8,15)`、`(4,10)`、`(12,15)`、`(6,10)`。
标量/B系数为eV，N系数无量纲；训练占据越界时保留 gate=0 候选和报告，求值器拒绝使用。N不裁剪。
读取仍检查格式头、版本、字节序、尺寸和完整字节数；原 `model.joint` 内容无需转换，只改变默认文件名。

## 求值与角度导数

输入四行 `ex ey ez axis_x axis_y axis_z`，按 DFT Ni 顺序；方向与轴均为单位向量。
`derivative_site=1..4` 选择扰动输入。导数按 `dei/dtheta=axis×ei` 解析求出，包含标量/B及所有受影响Ni的N。
`ueff` 从求值输入实际读取，`VU=ueff*(I/2-N)`、`dVU=-ueff*dN`；改变此参数并不重新求解DFT密度。

省略求值的第四参数时，在调用目录输出下列文件及 `eval_report.txt`。保留原二进制布局，不拆分记录：

| 原名 | 新名 | 单位 |
|---|---|---|
| `fields_Ry.bin` | `V_fit` | Ry，scalar/Bx/By/Bz |
| `dfields_Ry_per_rad.bin` | `dV_fit` | Ry/rad，同上 |
| `N_VU_eV.bin` | `N_VU_fit` | N无量纲，VU为eV |
| `dN_dVU_eV_per_rad.bin` | `dN_dVU_fit` | dN为1/rad，dVU为eV/rad |

仍可传入第四参数 `output_prefix`，例如 `/absolute/new_run/probe`，生成 `probe.V_fit` 等四个文件及 `probe.txt`。
场文件首记录为网格，第二记录为 `(ngrid,4)` float64，网格x最快；N文件首记录为原子顺序及Ueff，导数文件首记录为原子顺序及扰动Ni。
N/VU及其导数各保留两个 `(10,10,4)` complex128记录，包含全部轨道/自旋非对角元。

本入口完成 DFT标签→三势拟合→势和解析导数求值。
新版 `V_func_tmp` **不兼容**旧五记录同名文件或局域 `.au` 格式，不能作为旧 `fit_new` 的重启输入。
该 scalar8/B40 拟合与网格求值入口继续保留；其专用 `joint` RVR 和 `ansatzU` 后端已于 2026-11-16 退役。当前 [V1 RVR](../../6_RVR/README.md)读取冻结局域核与反馈模型，不接收此处生成的参数。旧后端可从 [清理前版本](https://github.com/zeen1086/ansatz-U/tree/3e1307ac721e91a66eb3851dae2df6bc4c9c99e1/6_RVR/joint) 查阅。
