# ansatz+U

LaNiO₂ 的 Fortran 势拟合、Wannier/RVR、Hamiltonian、能带和转矩工具。沿用 [ansatz](https://github.com/zeen1086/ansatz) 的阶段划分。

| 目录 | 内容 |
|---|---|
| `0_QE` | QE 7.6 源码修改、安装与编译脚本 |
| `1_scf`—`3_band_no_mag` | 结构、SCF 和 bands 输入模板 |
| `4_potential` | 显式通道映射、势拆分、两阶段拟合与网格转换 |
| `5_wannier` | NSCF/Wannier 模板、波函数与投影格式转换 |
| `6_RVR` | 标准/加速投影、HR 重排、冻结模型读取 |
| `7_spin` | 自旋组装、能带、能量、转矩和 ansatz+U |
| `common`、`examples`、`tests` | 公共数值模块、输入模板、离线测试 |

`*_lanio2.f90` 保留原 LaNiO₂ 阶段接口。`parameter` 用于公共工具和DFT模板；+U运行参数来自 `examples` 中的namelist。二者各有构建入口，互不覆盖。

四自旋的标量、B和非刚性Hubbard三势拟合见 [原生 Fortran/Shell 入口](4_potential/fit/README.md)。`make joint BUILD=/absolute/build` 构建 `fit_joint/eval_joint`；shell直接解包已有DFT标签，由Fortran读取、带权SVD拟合和求解析导数。采用 scalar8、B12+28、完整复数N特征；这是显式原始标签损失，未复现旧Python的投影算符联合损失。新三势参数默认写入 `V_func_tmp`，保留新版格式头及完整系数，**不兼容旧五记录同名文件**或局域 `.au` 模型；新特征到RVR的接口尚未完成。求值默认输出 `V_fit`、`dV_fit`、`N_VU_fit`、`dN_dVU_fit`，单位和布局见上述入口说明。

## LaNiO₂ ansatz+U

需要 gfortran 与 BLAS/LAPACK；MPI 程序另外需要 MPI 编译器。macOS 默认链接 Accelerate，Linux 默认链接 BLAS/LAPACK。

```sh
make BUILD=/absolute/external/build
# 需要原 MPI 阶段时：
make mpi BUILD=/absolute/external/build MPIFC=mpif90
/absolute/external/build/ansatzU /absolute/run/run.nml
```

从 `examples/7_spin/` 复制输入并填写外部路径：

- `run.nml`：Ueff（默认1 eV）、FD宽度（0.05 eV）、每Ni活性电子数（21）、计算模式和输出前缀。
- `texture.dat`：每Ni一行 `ex ey ez axis_x axis_y axis_z`，均为单位向量；按x、y排序，y最快。
- `kpoints.dat`：首行点数，后续 `Kx Ky Kz weight`，权重总和为1。
- `readouts.dat`：首行读出数，后续零基坐标 `x y`。

`mode='torque', delta=0` 计算给定纹理；`mode='table3'` 对指定Ni旋转 `delta` 弧度后取转矩差分。输出能带、EF、B/U与总转矩及Φ，已有文件不覆盖。

冻结模型只转换一次；计算阶段直接由 Fortran 读取：

```sh
/path/to/existing/python 6_RVR/export_frozen.py /external/local_model.npz /external/local_model.au
```

模型包含 `hc_R=H0+Hv`、`KB`、`p_R`、完整复数 `N0` 和 `rr`。二进制为little-endian、列优先、int32/complex128；头为16字节magic及8个int32，数组依次为 `(3,729)`、`(30,30,729)`、`(5,30,729)`、`(810,810)`、`(10,10)`。转换工具提供SHA与读回证书。

采用 `H=(H0+Hv)+HB+Σ A†Ueff(I/2−N)A`，`N=D N0 D†`；B/U导数在同一完整H本征态下读出。能量为eV、角度为rad、转矩为eV/rad、Φ为meV/rad²。带能和为 `ΣwfE`，不等于完整DFT总能。

| 输入 | 计算路径 |
|---|---|
| `ueff` | 同时乘HU及 `−A†(∂N/∂θ)A` |
| `sigma`、`electrons_per_ni` | 固定电子数 `electrons_per_ni*L²` 求EF和FD占据 |
| 纹理、旋转轴、`perturb_site`、`delta` | Rodrigues自旋旋转与同一SU(2)占据旋转；`∂N=[−i(axis·σ)/2,N]` |
| k点及权重、读出列表 | 完整H对角化，`τ=−Σwf⟨∂HB+∂HU⟩`；按列表顺序输出 |
| 张量路径、单位、空间范围、幅值 | 替换HB并同步替换B导数；HU保持相同projector与占据定义 |

局域KB已含B幅值，不额外乘0.372、0.4或S²。接入原三系数张量RVR时用 `tensor_run.nml`，显式填写 `tensor_unit` 与 `tensor_spin_scale`；该HB替换局域HB。RVR记录须为单条complex128记录，带4字节little-endian长度标记。

当前冻结 +U 模型为30轨道、单Ni、L×L×1，N保留全部复数元且不裁剪本征值。`hc_R`为H0与固定Hv合并背景；输入背景不应再次包含同一HU。

原LaNiO₂拟合入口为 `fit_new`，模板在 `examples/4_potential/fit.nml`。`diagonal_iters` 后接 `full_iters`，学习率显式指定；`restart_file`读取原五记录V_func_tmp。原阶段使用 `examples/THE.in`，外部运行树为 `6_RVR/super_H`、`7_spin/data` 和 `7_spin/tot_HR`。

原接口顺序为 `fit_new → cal_RVR_split → rvr_combine → tensor_Hsk → combine → band/calc_tot_ene`；MPI程序的Wannier列索引为第二个参数，`tensor_Hsk`第二参数为显式三分量纹理。张量阶段只支持单Ni的 `n_func=3`，其他通道数明确报错。其纹理是含幅值的系数，+U纹理是单位方向，不能混用。

`rvr_combine`生成的 `B_general` 可直接作为+U的 `tensor_path`。+U内部完成H_sk、H0+Hv合并、HU和导数读出，不再叠加旧 `combine` 的已组装H。`tensor_unit`、`tensor_range`、`tensor_spin_scale`须与原投影匹配；projector、N0和背景仍取同一冻结模型与Wannier gauge。

## 公共 ansatz 流程

Python 3.11+、NumPy、Fortran、BLAS/LAPACK；MPI/Wannier阶段使用相应程序。用已有Python解释器执行以下命令，生成和构建目录均放在仓库外：

```sh
export ANSATZ_BUILD_DIR=/absolute/external/ansatz-build
/path/to/existing/python configure.py --output-dir /absolute/external/case
/path/to/existing/python build.py 4_potential/V_split.f90 4_potential/make_R_vec.f90
```

修改 `parameter` 后重新生成输入并构建。`run.py` 核验构建签名并保留调用目录。QE安装见 `0_QE/README.md`；模板中的赝势名和数值网格须按实际输入填写。结构、投影和晶格模板来自原LaNiO₂输入，不改变冻结模型。

可用 `ANSATZ_PARAMETER=/absolute/external/parameter` 选择外部配置。Wannier运行脚本保留调用目录，先进入外部运行目录再执行；仓库内默认路径为显式占位路径，不读取相邻历史目录。

`V_split`读取 `polar2` 和显式映射文件：首行 `n_spins n_channels n_cells`，随后每个拆分胞一行1基自旋索引，胞按x/y/z排序、z最快；例见 `examples/4_potential/channel_map.dat`。输出保留原 `V1...`，并生成拟合直接读取的 `spinN/V_HYC_dp` 与 `spin_info`。

若标签包含周期端点，拟合后用 `make_R_vec base_V_HYC_dp fitted_extended_model output_model output_R_vec` 去掉端点；五条模型记录完整保留。拟合重启须使用与标签网格一致的模型。原LaNiO₂ RVR按完整网格读取；公共RVR按 `R_vec` 排序，输入模型与Wannier网格必须对应。

`check.py` 执行离线测试；数值结果摘要见 `VALIDATION.md`，来源见 `PROVENANCE.md`。模型、波函数、缓存和运行输出通过外部路径提供。
