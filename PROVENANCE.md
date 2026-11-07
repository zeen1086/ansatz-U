# 来源

公共流程来自 [zeen1086/ansatz](https://github.com/zeen1086/ansatz)，commit `05fa36cbf6fb7740bd55293466e79fa13af545c7`。已复制其QE修改、输入模板、拟合、Wannier、RVR、spin、构建配置和测试；根README与忽略规则作合并。QE许可证保留在 `0_QE/License`。

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
