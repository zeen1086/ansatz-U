# ansatz+U

LaNiO₂ 的 Fortran 势拟合、Wannier/RVR、Hamiltonian、能带和转矩工具。沿用 [ansatz](https://github.com/zeen1086/ansatz) 的阶段划分。

| 目录 | 内容 |
|---|---|
| `0_QE` | QE 7.6 源码修改、安装与编译脚本 |
| `1_scf`—`3_band_no_mag` | 结构、SCF 和 bands 输入模板 |
| `4_potential` | 势拆分、拟合与网格转换；`local/` 为冻结 V1 推理源 |
| `5_wannier` | NSCF/Wannier 模板、波函数与投影格式转换 |
| `6_RVR` | 当前 V1 局域空间核投影、合并、背景和 projector 读取 |
| `7_spin/local`、`7_spin/production` | V1 Hamiltonian、批量计算及完整能量读出 |
| `common`、`examples`、`tests` | 公共数值模块、输入模板、离线测试 |

## 当前 V1 流程

[RVR 目录与依赖](6_RVR/README.md)只保留第30轮冻结 V1 需要的实现。完整操作见 [V1 → V_func_tmp → 原生 RVR → cal_Hsk](6_RVR/local/README.md)。固定几何 L×L×1 使用全部 V1 通道、完整复数 N 和反馈响应：

```text
local_vfunc.py export / export-background
  → prepare-projection
  → local/cal_RVR_split → local/combine → B_general
  → local_vfunc.py cal_Hsk → local/cal_Hsk
  → 7_spin/production 批量能量、转矩和 J12 读出
```

需要已有 Python 3.11+、NumPy/SciPy、gfortran 和 BLAS/LAPACK；macOS 默认链接 Accelerate。只编译当前 V1：

```sh
make local-rvr BUILD=/absolute/external/build
```

`local_vfunc.py export` 的 `--code` 可指向本仓库 `4_potential/local/`。模型参数、Wannier、RVR 和计算输出均使用外部路径。背景由 `export-background` 生成纯 H0，当前组装为 `H0 + Hv + HB + HU`；`.au` 保持已有二进制布局，旧 KB/N0 占位字段只校验，不参与当前计算。

[7_spin/production](7_spin/production/README.md)提供同 kz 批量组装、RVR 伴随收缩、节点本地切片，以及冻结补空间密度下的完整能量/转矩/J12。使用 [外部 JSON 配置模板](examples/7_spin/production.json)。完整能量读出与带自由能读出有不同定义；实现回归通过不等于独立 DFT 物理验证。

## 已退役接口与其他阶段

2026-11-16 清理了旧三分量 RVR、scalar8/B40 `joint` RVR、`speedup`、旧 HR 重排和刚性 N0/张量 `ansatzU` 运行入口。旧命令和模板不再受本分支支持，源代码可从 [清理前版本 3e1307a](https://github.com/zeen1086/ansatz-U/tree/3e1307ac721e91a66eb3851dae2df6bc4c9c99e1) 查阅；格式、背景定义不同，旧模型不能改名后交给当前 V1。

本次保留拟合阶段：`make joint BUILD=/absolute/build` 构建 `fit_joint/eval_joint`，其 scalar8/B40/复数 N 拟合和求值说明见 [fit](4_potential/fit/README.md)，参数不适用于 V1。`make legacy` 保留 `fit_new` 及 `7_spin` 中读取已有外部历史矩阵的独立工具；旧 RVR 生成链已删除。`make mpi` 仅增加历史 band/能量读出，需已有 MPI 编译器。`make all` 构建 V1 与这些保留的非 MPI 工具。

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

若标签包含周期端点，历史五记录拟合模型可用 `make_R_vec base_V_HYC_dp fitted_extended_model output_model output_R_vec` 去掉端点；五条记录完整保留。拟合重启须使用与标签网格一致的模型。当前 V1 使用自身导出包，空间核与 Wannier 网格、平移及权重的匹配见 `6_RVR/local/README.md`。

`check.py` 执行离线测试；数值结果摘要见 `VALIDATION.md`，来源见 `PROVENANCE.md`。模型、波函数、缓存和运行输出通过外部路径提供。
