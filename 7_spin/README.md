# 7_spin：Hamiltonian 与读出

run_ansatzU.f90/local_ansatz.f90 执行冻结模型或全分量张量 RVR 的 +U、FD/EF、能带和 B/U 转矩/Phi；所有数值输入和输出路径显式给定。模板在 ../examples/7_spin/。

原阶段链的 cal_Hsk.f90、combine.f90、band.f90 和 calc_tot_ene.f90 从外部运行树的 7_spin 工作目录调用；参考仓库同样将 Hsk_* 放在 data/，合并结果放在 tot_HR/。背景通过 THE_dir 明确读取 6_RVR/super_H，RVR 读取 6_RVR/B_general。AFM/cal_Hsk.f90 保留旧源码入口，也从同一 7_spin 工作目录调用。编译和数值输出仍须在交付目录之外。
