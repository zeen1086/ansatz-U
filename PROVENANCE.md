# 来源

2026-11-16：`6_RVR` 已收敛到当前 V1。下列旧三分量、scalar8/B40 和刚性 N0 入口的来源记录对应 [清理前版本 3e1307a](https://github.com/zeen1086/ansatz-U/tree/3e1307ac721e91a66eb3851dae2df6bc4c9c99e1)，不是当前可执行接口。当前文件和依赖见 [6_RVR](6_RVR/README.md)；`background_projectors`、`project_at` 从旧 `local_ansatz` 原样移入 `frozen_io`，`.au` 格式和校验保持兼容。

公共工具来自 [zeen1086/ansatz](https://github.com/zeen1086/ansatz)，commit `05fa36cbf6fb7740bd55293466e79fa13af545c7`。保留QE修改、Wannier、RVR及构建工具；材料专用拟合、对称投影和自旋通道分支已删除，输入模板改用原LaNiO₂结构与投影。QE许可证保留在 `0_QE/License`。

LaNiO₂ 分支来源对应：

| 原工作树文件 | 历史导入位置 |
|---|---|
| `fit/fit_new.f90` | `4_potential/fit/fit_new.f90`、`fit_mask.f90` |
| `1_make_RVR/cal_RVR_split.f90` | `6_RVR/cal_RVR_split_lanio2.f90` |
| `1_make_RVR/combine.f90` | `6_RVR/combine_lanio2.f90` |
| `9_super_wan/rearrange_HR_new.f90` | `6_RVR/rearrange_HR_lanio2.f90` |
| `9_super_wan/H_sk/cal_Hsk_new.f90` | `7_spin/cal_Hsk_lanio2.f90` |
| `9_super_wan_AFM/H_sk/cal_Hsk_new.f90` | `7_spin/AFM/cal_Hsk_lanio2.f90` |
| `9_super_wan/combine.f90` | `7_spin/combine_lanio2.f90` |
| `9_super_wan/band.f90` | `7_spin/band_lanio2.f90` |
| `9_super_wan/calc_tot_ene.f90` | `7_spin/calc_tot_ene_lanio2.f90` |

新增 `frozen_io`、`tensor_rvr`、`local_ansatz`、`numerics`、`run_ansatzU` 移植并扩展了既有LaNiO₂局域Python参考。冻结模型源SHA256：`6fe8cb35d6271638fea0fe9227456a4f1b04910fb8cf181fe30cc212420c2ba0`。完整差异、原源码身份和数值记录保存在外部归档。

四自旋三势入口 `joint_features/joint_model` 对应既有 `covariant_basis.py`、`b_features.py` 中的 scalar8、B12+28、复数N及64个晶体作用。`joint_npy` 原生读取既有DFT标签，`joint_svd/fit_joint` 为新带权SVD目标，`eval_joint` 输出三势及解析角度导数，`run_joint.sh` 仅解包、构建和启动；生产链不调用Python。该目标按用户选择实现，不继承旧算符损失的数值资格或训练结果。

三势参数默认名由 `model.joint` 改为 `V_func_tmp`，求值文件沿用 `V_fit/dV_fit/N_VU_fit/dN_dVU_fit` 命名；新版格式与旧五记录同名模型不兼容。

三势后续阶段沿用原名称，在 `joint/` 子目录区分新格式：

| 原阶段/共用函数 | 新入口或扩展 | 功能 |
|---|---|---|
| `cal_RVR_split` | `6_RVR/joint/cal_RVR_split.f90`、`joint_rvr.f90` | 沿用逆平移Wannier索引及网格收缩，读取新版V_func_tmp，投影标量8/B40/固定标量参考核 |
| `combine` | `6_RVR/joint/combine.f90` | 验证分行绑定和Hermiticity，生成新版B_general并保留原背景/projector与N系数 |
| `local_ansatz::assemble` | 共用 `background_projectors` | 提取原H0+固定Hv与p_R的插值/周期相位，无数值定义改动；原局域入口仍调用同一函数 |
| `cal_Hsk` | `7_spin/joint/cal_Hsk.f90` | 组装标量势差、B及全四Ni HU；同一H下读出三项转矩、能带/EF及Φ |
| `run_ansatzU` | 新增显式 `joint_rvr_path` | 选择四自旋原2×2胞后端；拒绝重复背景/B输入或隐式大胞N延拓 |

primitive30读取映射直接遵循 `5_wannier/make_wan/make_wan_fun_mpi.f90` 的 `+ikR` 输出和原RVR的逆平移约定；不重新Wannier化、不增加幅值因子。新阶段的生产数值均为Fortran/BLAS/LAPACK，未改原拟合特征或权重。
