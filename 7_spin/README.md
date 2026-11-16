# V1 Hamiltonian 与能量使用说明

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
  --output /absolute/H_case
```

`local/cal_Hsk.f90` 组装 `H=H0+Hv+HB+HU`。默认输出逐 k 的 H 和 `HAMILTONIAN.json`；添加 `--derivative-sites 1 2 --step 0.001` 计算角度导数，`--components` 输出分量，`--keep-intermediates` 保留调试文件。使用新的输出目录。

## 批量能量、转矩与 J

按 [production 使用说明](production/README.md) 填写外部配置，并在 PBS 计算任务中运行。选择读出时区分：

| 入口 | 能量与输出 |
|---|---|
| `production/run_material.py`、`production/kz_production.py` | Gaussian 带自由能 `Fband`；`TORQUE.json` 中的 `F_eV`，对应转矩 `−dFband/dθ` |
| `production/complete_energy.py` | 固定补空间密度下的完整试探密度自由能；`ENERGY.json` 中的 `energy`，对应转矩 `−dF/dθ` |
| `production/summarize_energy.py` | 汇总 J12，并比较完整能量与转矩的单站点曲率；该曲率不是 J12 |

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
