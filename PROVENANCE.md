# 来源

公共工具来自 [zeen1086/ansatz](https://github.com/zeen1086/ansatz)，commit `05fa36cbf6fb7740bd55293466e79fa13af545c7`。保留QE修改、Wannier、RVR及构建工具；材料专用拟合、对称投影和自旋通道分支已删除，输入模板改用原LaNiO₂结构与投影。QE许可证保留在 `0_QE/License`。

LaNiO₂ 分支来源对应：

| 原工作树文件 | 当前文件 |
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
