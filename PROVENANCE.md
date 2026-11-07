# 代码来源与对应

9个原始文件先从维护者的 LaNiO₂ 工作树 `source/` 原样复制，再仅修改副本。`source/` 表示外部来源，不是本仓库中的缺失目录。原始文件的 SHA256 记录如下；原数据和历史结果没有写入本仓库。

| 原始来源 | 本目录文件 | 合并/修复内容 |
|---|---|---|
| source/fit/fit_new.f90 | 4_potential/fit/fit_new.f90、4_potential/fit/fit_mask.f90 | 原逐样本 B 更新、五记录 V_func_tmp 重启；显式对角/全分量两阶段、显式迭代/学习率/输入输出 |
| source/1_make_RVR/cal_RVR_split.f90 | 6_RVR/cal_RVR_split.f90 | 保留全部 B 张量分量；运行时列索引替代模板 HYC_change；显式 THE.in/网格来源、原 MPI 投影 |
| source/1_make_RVR/combine.f90 | 6_RVR/combine.f90 | 保留原列拼接和 B_general 布局；显式 THE.in |
| source/9_super_wan/rearrange_HR_new.f90 | 6_RVR/rearrange_HR.f90 | 原背景超胞重排；显式输入、双精度构建 |
| source/9_super_wan/H_sk/cal_Hsk_new.f90 | 7_spin/cal_Hsk.f90 | 原张量 RVR→H_sk 累加及原空间展开；显式任意纹理，移除隐含幅值与小矩阵元丢弃 |
| source/9_super_wan_AFM/H_sk/cal_Hsk_new.f90 | 7_spin/AFM/cal_Hsk.f90 | 同上；保留旧 AFM 入口但不再内置纹理 |
| source/9_super_wan/combine.f90 | 7_spin/combine.f90 | 原背景+H_sk；complex128 输出、声明 Ry/eV 换算 |
| source/9_super_wan/band.f90 | 7_spin/band.f90 | 原能带路径，全部数组/MPI/ZHEEV 双精度 |
| source/9_super_wan/calc_tot_ene.f90 | 7_spin/calc_tot_ene.f90 | 原 sharp-filling 带能，读取的总活性电子数生效；双精度 MPI/ZHEEV |
| local_model.py 的 interp/keff/assemble | 6_RVR/frozen_io.f90、7_spin/local_ansatz.f90 | 原冻结插值、空间核、projector；广义逐 Ni 纹理、完整复数 N、显式 Ueff |
| 原全分量 B_general 投影布局 | 6_RVR/tensor_rvr.f90 | 原 RVR 与新增 HU/FD/转矩的直接接口；显式 Ry/eV 转换，不将张量能力缩减为 b(r)e |
| local_model.py 的 derivatives/filling；table3_f4.py、extend.py | common/numerics.f90、7_spin/run_ansatzU.f90 | 三维旋转/交换子、B/U 同 H 读出、固定电子数 FD、能带、转矩/Φ；网格全部由输入给定 |
| qualify.py 的冻结结果 | 6_RVR/export_frozen.py（一次性格式转换） | 复用已资格化的模型；不重投影；无损读回，生产不调用 Python |

所列 Python 参考均来自 `current/runs/seki_table3_local_ansatzU_20261107_r01/code/`，extend.py 来自对应 `..._L8L10_20261107_r01/code/`。已阅读两个 run 的 PROTOCOL、LOCAL_MODEL_FROZEN 和报告；资格证据复用，数值移植另外按双精度门槛核对。

源码身份的有限核查：当前 LaNiO₂ `source/fit/fit_new.f90` 含对角更新限制；归档 `archive/cri3_original/source/pbs/4_potential_2x2/fit.f90` 含已注释的同一限制，而 fit0.f90 保留限制。后者证明归档中存在全分量更新源码，不能据此认定其为 LaNiO₂ 后续版本或推断历史训练次数/学习率。本任务显式保留用户要求的两个阶段；没有重训。

目录结构参考 [zeen1086/ansatz](https://github.com/zeen1086/ansatz)，main 树 `05fa36cbf6fb7740bd55293466e79fa13af545c7`，README blob `2bc0c0878f25abf2c0e1f243f827e6e975b9156b`。采用其 0_QE 至 7_spin 的阶段编号、拟合子目录、RVR 背景重排与 spin/data/tot_HR 位置；没有导入其 CrI₃ 数值模型、源码或 Python 配置/运行系统。移动记录、完整逐行补丁和旧源码快照由维护者保存在仓库外。

原始来源身份（未修改的输入源码）：

| 原文件 | SHA256 |
|---|---|
| source/fit/fit_new.f90 | `1ef3f9e1471d7942ad47a945e973f1e6755d96db4ebaf6c751ab591acb5f6a8b` |
| source/1_make_RVR/cal_RVR_split.f90 | `ebec429a295ecddc46b6cefc5a4bb8c546bcd93906f08e1a573adf750e429bf5` |
| source/1_make_RVR/combine.f90 | `30a7a21ccf963056eeca3c7f183ee7c7a0e8d1271f0eb8d2e9848585a49fff1b` |
| source/9_super_wan/calc_tot_ene.f90 | `fe3e27ff780e689bb0cd3ebe5fe3921f15ec0a4a5e9fa144595f3ebec9f892de` |
| source/9_super_wan/rearrange_HR_new.f90 | `590b1ae20899e351b2d98898aaf3af54a07df9414708947240334c72959dc58d` |
| source/9_super_wan/combine.f90 | `3ecc2a2898b3c238ae38ca419d914e7602d7d6fb8aa5be5a37094768af2cae67` |
| source/9_super_wan/band.f90 | `0d815e1162c5f47a004f91d1a4f040919b978b43854539067052f9514e2662cb` |
| source/9_super_wan/H_sk/cal_Hsk_new.f90 | `49f7009e7ccd38bc713cfdd57e0618f181f18082adc7624e868bc3d6e6f95698` |
| source/9_super_wan_AFM/H_sk/cal_Hsk_new.f90 | `b8fca8dd347c0d474d7e89418ac4328a36d6e004f174b65889c1afe30933f01a` |
