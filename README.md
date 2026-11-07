# ansatz+U

Fortran 势拟合、Wannier/RVR、Hamiltonian、能带和转矩工具。包含 [ansatz](https://github.com/zeen1086/ansatz) 的完整公共流程，以及 LaNiO₂ 的局域 ansatz+U 实现。

| 目录 | 内容 |
|---|---|
| `0_QE` | QE 7.6 源码修改、安装与编译脚本 |
| `1_scf`—`3_band_no_mag` | 结构、SCF 和 bands 输入模板 |
| `4_potential` | 势拆分、拟合、对称处理与导出 |
| `5_wannier` | NSCF/Wannier 模板、波函数与投影格式转换 |
| `6_RVR` | 标准/加速投影、HR 重排、冻结模型读取 |
| `7_spin` | 自旋组装、能带、能量、转矩和 ansatz+U |
| `common`、`examples`、`tests` | 公共数值模块、输入模板、离线测试 |

标准文件名沿用上游实现；`*_lanio2.f90` 保留原 LaNiO₂ 阶段接口。上游 `parameter` 是 CrI₃ 示例；LaNiO₂ +U 使用 `examples` 中的 namelist，不读取该示例参数。

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

局域KB已含B幅值，不额外乘0.372、0.4或S²。接入原三系数张量RVR时用 `tensor_run.nml`，显式填写 `tensor_unit` 与 `tensor_spin_scale`；该HB替换局域HB。RVR记录须为单条complex128记录，带4字节little-endian长度标记。

当前冻结 +U 模型为30轨道、单Ni、L×L×1，N保留全部复数元且不裁剪本征值。`hc_R`为H0与固定Hv合并背景；输入背景不应再次包含同一HU。

原LaNiO₂拟合入口为 `fit_new`，模板在 `examples/4_potential/fit.nml`。`diagonal_iters` 后接 `full_iters`，学习率显式指定；`restart_file`读取原五记录V_func_tmp。原阶段使用 `examples/THE.in`，外部运行树为 `6_RVR/super_H`、`7_spin/data` 和 `7_spin/tot_HR`。

## 公共 ansatz 流程

Python 3.11+、NumPy、Fortran、BLAS/LAPACK；MPI/Wannier阶段使用相应程序。用已有Python解释器执行以下命令，生成和构建目录均放在仓库外：

```sh
export ANSATZ_BUILD_DIR=/absolute/external/ansatz-build
/path/to/existing/python configure.py --output-dir /absolute/external/case
/path/to/existing/python build.py 7_spin/cal_Hsk.f90 7_spin/combine.f90 7_spin/cal_torque_hf.f90
```

修改 `parameter` 后重新生成输入并构建。`run.py` 核验构建签名并保留调用目录；原流程、QE安装和Wannier转换命令见 [上游说明](https://github.com/zeen1086/ansatz/blob/main/README.md) 与 `0_QE/README.md`。

`check.py` 执行离线测试；数值结果摘要见 `VALIDATION.md`，来源见 `PROVENANCE.md`。模型、波函数、缓存和运行输出通过外部路径提供。
