"""
PIES 调度模型参数
"""
# ==========MTG 相关==========
# MTG 数量
n_G = 10
# MTG 最大允许爬坡速率
delt_P_MTG_U_max = 0.5
# MTG 最大允许滑坡速率
delt_P_MTG_D_max = 0.5
# MTG 最大储能
P_MTG_i_max = 10
# ==========ESS 相关==========
# ESS 最大充电量
P_ESS_ch_max = 100
# ESS 最大放电量
P_ESS_dc_max = 100
# ESS 充电效率
e_ch = 0.8
# ESS 放电效率
e_dc = 0.8
# ESS 最低储能量
C_ESS_min = 0
# ESS 最高储能量
C_ESS_max = 100
# ==========电网相关==========
# t 时 PIES 从电网买电量上界
P_b_gr_max = 100
# t 时 PIES 向电网卖电量上界
P_s_gr_max = 100
# ==========其他相关==========
# 调度周期
T = 24
# 对电负载的灵活 TL 补偿价格
rho_TL = 1.5
# 对电负载的灵活 IL 补偿价格
rho_IL = 1.5
# 对中断冷负载的灵活冷补偿价格
rho_c = 1.5
P_TL_in_min = 0
P_TL_in_max = 100
P_TL_out_min = 0
P_TL_out_max = 100

# IL 中断负载下界
P_IL_t_min = 0
# IL 终端负载上界
P_IL_t_max = 100

# PV cell processing 产线数
n_c = 10
n_cp = 10
# PV module encapsulation 产线数
n_m = 10
n_mp = 10

# 建筑数
W = 4

# T_m
T_m_min = 0
T_m_max = 23
# T_c
T_c_min = 0
T_c_max = 23

# 场景数
n_s = 3*5

F_m = 10
F_c = 10
F_cp = 10
G_cp = 10
F_mp = 10
G_mp = 10
