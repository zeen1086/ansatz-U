# 5_wannier：外部同 gauge 数据

对应 Wannier/实空间转换阶段。本交付复用已有 Wannier gauge、原 Hubbard projector 和冻结空间核；没有复制波函数、完整 H/W 或投影缓存。原 RVR 分支通过显式 nscf_dir 引用 make_wan/wan_data/wan_*、cal_TF_*、ir_vec_HYC、H_deg_HYC 与 hr_HYC_*，这些记录必须属于同一来源。完整冻结模型运行无需重新执行此阶段。
