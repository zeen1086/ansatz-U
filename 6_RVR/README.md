# V1 RVR 使用说明

在仓库根目录执行以下命令。将示例中的绝对路径替换为实际路径，`ANSATZ_PYTHON` 指向已安装 NumPy、SciPy 的 Python 环境。模型、构建和计算输出放在仓库外。

## 1. 编译与模型准备

```sh
export ANSATZ_PYTHON=/absolute/existing/python
make local-rvr BUILD=/absolute/build
```

生成 `/absolute/build/local/rvr` 和 `/absolute/build/local/cal_Hsk`。

已有导出模型和纯 H0 背景时可跳过下面两步；否则执行：

```sh
"$ANSATZ_PYTHON" local_vfunc.py export \
  --source /absolute/frozen_V1/result/final \
  --code /absolute/frozen_V1/model_code \
  --structure /absolute/L02_structure.npz \
  --scalar-offset-eV=-7.187822007469299 \
  --output /absolute/model/V_func_tmp

"$ANSATZ_PYTHON" local_vfunc.py export-background \
  --source /absolute/DFT_only/primitive.bin \
  --output /absolute/basis/pure_H0.au
```

上述标量常数用于当前冻结 V1。导出结构为原始四 Ni 几何；结构 NPZ 包含 `cell, positions, species, ni, spins`，分别使用 Å、原子序数、零基 Ni 索引和单位自旋方向。背景输入须为兼容的 DFT-only `primitive.bin`，与投影 Wannier 使用共同基底。

## 2. 准备投影输入

需要 `wan_1 ... wan_30`：little-endian complex128 流文件，按 Fortran 顺序排列为 `(n_translation, n_grid)`，网格与模型一致。`weighted_translations.dat` 首行为平移数，随后每行 `Rx Ry Rz weight`。

```sh
"$ANSATZ_PYTHON" local_vfunc.py prepare-projection \
  --model /absolute/model/V_func_tmp \
  --wan-dir /absolute/wan \
  --translations /absolute/weighted_translations.dat \
  --background /absolute/basis/pure_H0.au \
  --grid-weight 1 \
  --kz 0 0.125 0.25 0.375 0.5 0.625 0.75 0.875 \
  --output /absolute/RVR
```

`--grid-weight 1` 适用于既有离散归一化来源；其他来源须填写实际积分权重。程序不会自动裁切 Wannier 支持、补边界权重或重新归一化。输出目录须为新目录。

## 3. 计算、续算与输出

```sh
/absolute/build/local/rvr prepare /absolute/RVR/projection.nml
/absolute/build/local/rvr project /absolute/RVR/projection.nml
/absolute/build/local/rvr combine /absolute/RVR/projection.nml
```

`prepare` 生成共享空间场缓存；`project` 计算轨道/kz 分片；`combine` 合并为可复用的 `B_general`。可指定分片任务，轨道和 kz 均为 1 基编号：

```sh
# 第 1 个轨道的全部 kz
/absolute/build/local/rvr project /absolute/RVR/projection.nml 1
# 第 1 个轨道的第 3 个 kz
/absolute/build/local/rvr project /absolute/RVR/projection.nml 1 3
# 调试时保留中间文件
/absolute/build/local/rvr combine /absolute/RVR/projection.nml --keep-intermediates
```

可以并行运行不同轨道/kz 任务，但每个文件只允许一个写入进程；全部投影完成后再合并。不要用 `mpirun` 启动多份写同一分片的进程。

中断后重跑相同命令：空间缓存检查头部与均值后复用，分片检查头部/字节数并抽检首尾块后复用，缺失分片补算，未完成的 `.partial` 从头生成。检测到格式/字节数错误、首尾块非有限值或身份不符时会报错；先将异常文件移到其他目录，再重试。最终文件已发布但清理中断时，重跑 `combine` 可继续清理。

续算期间保持模型、几何、Wannier、平移/权重及准备记录不变；输入变化时使用新的输出目录重新准备投影。

默认仅在合并成功后清理缓存和分片，失败时保留现场。最终目录包含：

```text
B_general          # 可复用矩阵元
PROJECTION.json    # 输入来源与准备记录
projection.nml     # 投影配置
kz.dat             # 已投影的 kz
```

`PROJECTION.json` 不是完成标记；成功合并会输出 `B_GENERAL_COMPLETE`。

## 4. 不同超胞与自旋构型复用

固定原胞几何下，不同 L×L×1 超胞和自旋纹理共用同一个 `B_general`、模型包和 `pure_H0.au`，分别计算各自的 Hamiltonian。kx/ky 可改变，kz 必须在已投影列表内；更换模型、几何或 Wannier 基底需要重新投影。

```sh
"$ANSATZ_PYTHON" local_vfunc.py cal_Hsk \
  --model /absolute/model/V_func_tmp \
  --structure /absolute/L04_structure.npz \
  --rvr /absolute/RVR/B_general \
  --background /absolute/basis/pure_H0.au \
  --kpoints /absolute/kpoints.dat \
  --binary /absolute/build/local/cal_Hsk \
  --derivative-sites 1 2 --step 0.001 \
  --output /absolute/H_case
```

k 点文件首行为点数，随后每行 `kx ky kz weight`，权重和为 1。导数站点采用原生 Ni 编号 `x*L+y+1`，从 1 开始，表示绕全局 y 轴的角度导数；省略 `--derivative-sites` 则只计算 H。

默认输出逐 k 的 `value.k000001.h`、请求的 `dsite_XXXX.k000001.h` 和 `HAMILTONIAN.json`。矩阵为 little-endian complex128、Fortran 顺序、单位 eV，维数为 `60L²`；索引与 k 点信息保存在 JSON 的 `native_indices` 中。

添加 `--components` 输出 H0/Hv/HB/HU 分量；添加 `--keep-intermediates` 保留系数、原生输入和诊断文件。成功校验全部结果后才清理中间文件，失败时保留现场。单构型 `cal_Hsk` 不自动续算已有输出目录。

批量能量、转矩与生产检查点用法见 [production 使用说明](../7_spin/README.md)。

已完成的 `B_general` 仅检查头部、维度、几何和字节数，不再次全量扫描；分片首尾抽检不能检测所有中间块损坏。空间缓存写完不再整文件回读。
