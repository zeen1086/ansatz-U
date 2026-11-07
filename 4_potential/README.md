# 4_potential：势拟合

fit/fit_new.f90 保留原 B_general 拟合及五记录 V_func_tmp 重启；fit_mask.f90 实现对角阶段后进入全分量阶段。输入模板为 ../examples/4_potential/fit.nml。需要训练时显式设置迭代数和学习率；整理目录不启动训练。
