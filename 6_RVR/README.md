# 当前 V1 的可复用 RVR

本目录只保留第30轮冻结 V1 所需的空间核投影与背景支持。实际输入、命令和格式说明见 [local/README.md](local/README.md)。

| 文件 | 作用与直接依赖 |
|---|---|
| `local_vfunc.f90` | 读取局域空间核和几何，插值并作周期求和；依赖 `common/numerics.f90`、`4_potential/fit/joint_npy.f90` |
| `local_rvr.f90` | Wannier 投影、重叠矩阵、分片校验与合并、RVR 块读取；依赖 `local_vfunc` |
| `local/cal_RVR_split.f90` | 生成空间场或指定轨道/kz 的投影行；依赖 `local_rvr` |
| `local/combine.f90` | 合并全部投影行；依赖 `local_rvr` |
| `frozen_io.f90` | 读取纯 H0/projector `.au`，执行 Bloch 插值、超胞背景和站点 projector 组装；依赖 `numerics`，由 `7_spin/local/cal_Hsk.f90` 调用 |

上游 `local_vfunc.py` 负责模型导出、纯 H0 导出、投影准备和冻结 V1 反馈求值。下游 `7_spin/local/cal_Hsk.f90` 读取 `local_rvr` 与 `frozen_io`，组装 H0/v/B/HU；`7_spin/production` 在该原生入口上计算批量能量与转矩。`joint_npy` 只是共用的 NPY 读取器，保留它不表示仍有 `joint` RVR 依赖。

一次投影包含 scalar62+bias1 与 B62 共 125 个空间核，并保存同一 Wannier 支持上的重叠矩阵。同一模型、几何、基底、支持和已投影 kz 下，RVR 可复用于不同 L 和自旋纹理；新 kz 或改变几何需要相应新投影。

```sh
make local-rvr BUILD=/absolute/external/build
```

输出为 `local/cal_RVR_split`、`local/combine`、`local/cal_Hsk`。三个程序为串行/BLAS 程序，投影按轨道/kz 从外部调度并行，不使用 `mpirun` 启动多份写同一分片的进程。Python 构建入口 `build.py` 同样只链接当前模块。

新计算默认精简输出：合并成功并检查最终文件后，`combine` 清除空间场缓存及轨道/kz 分片；`local_vfunc.py cal_Hsk` 校验全部结果后只保留矩阵和 `HAMILTONIAN.json`。两者添加 `--keep-intermediates` 可保留调试文件。失败时不清理；独立 `evaluate` 的系数和生产续算检查点仍保留。详见 [文件保留规则](local/README.md)。

`LANIO2_LOCAL_VFUNC_V2`、`LANIO2_LRVR_V2` 和纯 H0 的 `LANIO2_AU_V1` 布局未变。背景读取器保留旧 KB/N0 槽位的长度、有限值和 Hermiticity 校验，读取后丢弃；当前 B 来自 V1 RVR，N 来自完整反馈模型。`background_projectors`、`project_at` 从旧 `local_ansatz` 原样迁入，周期相位、Ni/轨道/自旋顺序不变。

2026-11-16 删除的接口：

- 根目录三分量 `cal_RVR_split/combine` 及 `_lanio2` 版本、`rearrange_HR` 和 `speedup`。
- scalar8/B40 的 `joint_rvr` 和 `joint/` 投影/合并入口。
- `tensor_rvr`、旧 `export_frozen.py`，以及其专用 `7_spin/local_ansatz`、`run_ansatzU`、`joint/cal_Hsk`、`run_spin.sh` 与运行模板。

旧源代码及说明仍可从 [清理前提交 3e1307a](https://github.com/zeen1086/ansatz-U/tree/3e1307ac721e91a66eb3851dae2df6bc4c9c99e1) 找回。拟合、DFT、Wannier 与读取已有历史矩阵的独立 `7_spin` 工具保留；当前分支不再生成旧 RVR。

2026-11-16 清理后验证：

| 检查 | 结果 |
|---|---|
| 外部空目录 `make all`，启用 `-fcheck=all`；`build.py` 单独构建三个 V1 程序 | 通过 |
| `check.py` | 23 项通过，1 项因本机缺少 MPI 跳过；删除的两个旧 RVR 测试不再计入 |
| `test_local_vfunc.py`，L=1/2/3、非零 k、完整复数 N 和分量/导数 | 120 个矩阵对照最大差 1.36e−16 eV；40 个批量/单点对照差为 0 |
| `test_v1_runtime_bridge.py`，真实 V1 反馈与已知测试投影核 | 最大差 9.16e−16 eV；15 个矩阵输出与清理前版本逐字节一致 |
| `test_production_adjoint.py`，含非局域 projector、Bloch 相位和批量导数 | 28 项通过，最大转矩差 2.43e−17 eV/rad |
| 已删除路径引用、文档相对链接、背景/projector 函数迁移 | 代码/配置无残留引用，链接有效；两个迁移函数与原文逐字节一致 |

小型测试覆盖实现一致性，不代替完整材料投影、HPC 生产运行或独立 DFT 物理验证。本次未重新运行生产计算。

同日默认输出精简回归：离线测试 24 项通过、1 项 MPI 跳过；原生合并的精简/调试结果 SHA256 一致，截断分片导致失败时全部中间文件保留；真实 V1 的 15 个 H/dH 分量在两种模式下逐字节一致，索引完整归入 JSON。原生进程失败或后续导数矩阵截断时，已成功的前序结果与所有预测输入均保留。生产批量检查点测试继续通过。
