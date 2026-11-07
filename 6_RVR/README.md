# 6_RVR：投影核与背景重排

cal_RVR_split.f90/combine.f90 保留原全分量张量投影和列合并；rearrange_HR.f90 在外部运行树的 6_RVR/super_H 写入重排背景。frozen_io.f90 读取已冻结局域模型，tensor_rvr.f90 将原 B_general 接入 +U，必须声明 Ry/eV 单位。

export_frozen.py 只用于一次性无损存储转换；日常计算由 Fortran 直接读取外部 .au 文件。使用已资格化冻结模型时不重建投影。运行命令与格式见根 README.md。
