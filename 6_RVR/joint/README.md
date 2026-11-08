# 三势 RVR、Hamiltonian 与转矩

`V_func_tmp → cal_RVR_split → combine → B_general → cal_Hsk → 能带/EF/转矩/Φ` 全部由 Fortran 执行。三势入口放在 `joint/` 子目录，沿用阶段文件名；原五记录模型、MPI RVR 和张量入口保留。

```sh
make joint-rvr BUILD=/absolute/external/build
mkdir -p /absolute/new_run/6_RVR/tmp
# 编辑 examples/6_RVR/joint_projection.nml，填写真实输入和新的 rvr_dir。
/absolute/external/build/joint/cal_RVR_split /absolute/run/projection.nml
/absolute/external/build/joint/combine /absolute/run/projection.nml
# 编辑 examples/7_spin/joint_run.nml，再运行任意四自旋纹理：
/absolute/external/build/ansatzU /absolute/run/joint_run.nml
```

省略投影的第二参数时一次处理全部120个空间轨道；也可传入1..120，生成单行 `tmp/B_general_i`。`combine` 检查行号、完整性、输入绑定和 Hermiticity 后新建 `B_general`。已有文件拒绝覆盖。构建和运行目录放在源码树外；没有自动拟合、DFT/Wannier计算或作业提交。

## 背景与输入

本入口明确采用标量势差方案：

`H = (H0 + Hv_ref) + [Hv(e) - Hv_ref] + HB(e) + Σ_i A_i† Ueff(I/2-N_i(e)) A_i`。

- `frozen_path` 为已有 `.au` 冻结背景：30个primitive空间轨道、原 `rr/hc_R/p_R`；格式必须声明无HU。沿用 `hc_R` 的固定标量背景及原 Hubbard projector，旧 KB/N0 不参与新三势的求值。
- `reference_potential_path` 必须是**实际生成该 hc_R 时使用的固定标量场**，不能用新拟合在参考纹理上的预测代替。读取原 `V_HYC_dp` 格式：little-endian顺序记录，首记录3个int32网格，第二记录 `(ngrid,4)` float64，单位Ry；仅使用第一列标量，转为eV投影。它可能来自旧拟合，不能仅凭文件名换成某个DFT标签。
- 必须声明 `background='hc_plus_fixed_Hv_no_HU'`、`reference_scalar_id` 和 `gauge_id`。ID是调用者对参考场/共同gauge来源的明确绑定，程序不能从任意势文件推断其DFT身份；真实输入仍须有匹配的来源记录。
- `model_path` 读取新版 `LANIO2_JOINT_V1` 的 `V_func_tmp`，保留标量8、B40及完整复数N模型；不读旧五记录格式。
- `wan_dir` 中 `wan_1...` 使用原stream complex128空间Wannier布局 `(nR,ngrid,1)`，不带记录标记。原 `make_wan_fun_mpi` 用 `+ikR` 写出这些平移；投影沿用原RVR的逆平移索引。无SOC，哈密顿量自旋次序为每空间轨道up/down；N和projector行次序为 `[up d1..5,down d1..5]`。
- `translations_path` 首行平移数，随后逐行整数 `Rx Ry Rz`，须包含反向平移、不得重复。它定义四Ni胞的Wannier支持；差平移全部保留，不做阈值裁剪。
- `support_path` 非空时严格读取其中的 `cal_TF_i`，每原轨道一行/组nR个0/1；空字符串表示使用列表中全部支持。输入网格x最快，需与势匹配。
- `grid_weight=1` 沿用原已归一化Wannier的网格求和，不额外乘体积、0.372、0.4或S²。若源数据归一化不同，必须显式声明其求积权重并记录依据。

既有冻结局域模型的协议将该参考标量定义为旧 scalar8 在理想Néel态的预测；它不是自动等于原始DFT标量标签。本方案对合并背景加入同一投影支持上的势差，不宣称从hc_R恢复了独立的DFT纯H0；真实输入的参考场身份及支持误差仍须核对。

`wannier_layout='primitive30'` 直接复用30个原primitive Wannier文件。primitive网格为 `(nx/2,ny/2,nz)`，平移由 `primitive_translations_path` 显式给出。程序按原约定把120个四胞轨道重排到Ni5,6,7,8对应的胞 `(0,0),(1,0),(0,1),(1,1)`，无插值、重拟合或重标度；输出平移若截断原支持则报错。`four_cell120` 则直接读取120个同一gauge的四胞Wannier文件。模板平移列表可按实际支持修改；primitive与四胞平移分别以各自晶胞为单位。

## 一次投影与复用

投影分别计算8个标量空间核、40个B空间核和1个实际固定标量参考核。对逻辑Wannier平移a/b，计算 `Σ_g W_a(g)* V_c(g) W_b(g) grid_weight`，按 `R_b-R_a` 累加；任意K点用 `exp(+2πi K·R)` 插值，B的三个Pauli分量在组装时完整加入。

新版 `B_general` 是一个带 `LANIO2_RVR_V1` 头的little-endian容器，包含单位eV的complex128空间核、完整N系数、原背景/projector及来源记录；不兼容旧 `B_general`。模型、背景、参考场、投影控制文件的双模整数指纹用于拒绝混用行文件，不作为安全哈希。Wannier和支持数据须在投影期间保持不变，gauge必须与冻结背景/projector一致。

每套势系数、Wannier/gauge、参考场、求积权重和空间支持只投影一次。完整 `B_general` 可移到其他运行目录；后续改变纹理、K点、Ueff、sigma、电子数或角步只读取它，不再读取势网格、Wannier或旧模型路径。更换这些空间输入需要新的运行目录与投影。

缓存的主要空间为 `16*120²*nR_difference*49` 字节，另有原背景/projector和少量系数；例如125个差平移约1.41GB。该公式不含投影工作数组或求解器工作区，不是作业资源承诺。

## 三项读出

`7_spin/joint/cal_Hsk.f90` 由原 `ansatzU` 程序通过 `joint_rvr_path` 调用。此入口要求 `L=2`，不能同时提供旧 `model_path/tensor_path`。四行纹理及旋转轴均为单位三维向量，按DFT Ni5..8顺序，即x最快；它与旧局域入口的y最快排序不同，不能混用。`perturb_site/readouts` 仍用零基x/y坐标。

`Ueff`、FD宽度和每Ni电子数均实际读取；默认1eV、0.05eV、21，固定电子数目标为84。所有分量使用同一个完整H的本征态和同一EF/FD占据：

`tau_v/B/U = -Σ_kn w_k f_kn <kn|dHv/dHB/dHU|kn>`，`tau_total = tau_v+tau_B+tau_U`。

N和dN由同一四自旋特征解析求出，保留全部复数元且不裁剪本征值；单个输入方向变化时，`dHU=-Ueff Σ_all_four_Ni A_i† dN_i A_i`。参考场固定，标量导数为dHv。这里是同一H的导数分解，不是U开/关效应。

- `.p1/.p2.bands`：完整能带、FD占据和K权重。
- `.p1/.p2.result`：EF、实际电子数、带能和，以及 `tau_v tau_B tau_U tau_total`，单位eV/rad。
- `mode='table3'`：`.phi` 输出 `Phi_v Phi_B Phi_U Phi_total`，定义 `[tau(delta)-tau(0)]/delta*1000`，单位meV/rad²，角度为rad，无S²因子。
- `dump_matrices=.true.`：输出 `.hbg/.hv/.hb/.hu/.h` 及每读出 `.dv/.db/.du`；`.hbg=H0+Hv_ref`，`.hv=Hv-Hv_ref`，全为complex128/eV或eV/rad，不假称纯H0。

带能和 `ΣwfE` 不是完整DFT总能量。本实现保持原四自旋模型的空间和gauge定义，不把非刚性N自动延拓到4/6/8/10等大超胞；真实数据、空间支持和物理精度资格须另行核对。
