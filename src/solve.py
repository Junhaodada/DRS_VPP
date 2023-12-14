"""
PIES 主问题与子问题数学规划求解
"""
from params import *
from docplex.mp.model import Model


def solve_milp(type, gama, gap):
    model = Model('solve')
    # ======================✅决策变量=====================
    # ======================一阶段决策变量=====================
    # 第一阶段，开发组件的运行状态策略，组件包括 MTG、ESS 和工厂中的机器
    # 控制 MTG 开关状态
    eta_MTG_i_t = {(i, t): model.binary_var(name=f'eta_MTG_{i}_{t}') for i in range(n_G) for t in range(T)}
    # 控制 MTG 开关状态
    y_i_t = {(i, t): model.binary_var(name=f'y_{i}_{t}') for i in range(n_G) for t in range(T)}
    z_i_t = {(i, t): model.binary_var(name=f'z_{i}_{t}') for i in range(n_G) for t in range(T)}
    # 控制 ESS 充电状态
    eta_ch_t = {t: model.binary_var(name=f'eta_ch_{t}') for t in range(T)}
    # 控制 ESS 放电状态
    eta_dc_t = {t: model.binary_var(name=f'eta_dc_{t}') for t in range(T)}

    # 控制生产线运行状态-PV module encapsulation 产线
    v_m_t = {(m, t): model.binary_var(name=f'v_{m}_{t}') for m in range(n_m) for t in range(T)}
    v_mp_t = {(mp, t): model.binary_var(name=f'v_{mp}_{t}') for mp in range(n_mp) for t in range(T)}

    # 控制生产线运行状态-PV cell processing 产线
    v_c_t = {(c, t): model.binary_var(name=f'v_{c}_{t}') for c in range(n_c) for t in range(T)}
    v_cp_t = {(cp, t): model.binary_var(name=f'v_{cp}_{t}') for cp in range(n_cp) for t in range(T)}

    # 0-1 整数变量
    U = [
        y_i_t,
        z_i_t,
        eta_dc_t,
        eta_ch_t,
        v_m_t,
        v_mp_t,
        v_c_t,
        v_cp_t
    ]

    # ======================二阶段决策变量=====================
    # 第二阶段，在揭示光伏和冷负荷需求的不确定性后，制定调度方案，包括 MTG 输出、从主网购电量、卖电量等
    # t 时第 i 个 MTG 的生产需要的能量
    P_MTG_i_t = {(i, t): model.continuous_var(name=f'P_MTG_{i}_{t}') for i in range(n_G) for t in range(T)}
    # t 时第 i 个 MTG 的维护需要的能量
    u_MTG_i_t = {(i, t): 0 for i in range(n_G) for t in range(T)}
    # t 时 ESS 充电量
    P_ESS_ch_t = {t: model.continuous_var(name=f'P_ESS_ch_{t}') for t in range(T)}
    # t 时 ESS 放电量
    P_ESS_dc_t = {t: model.continuous_var(name=f'P_ESS_dc_{t}') for t in range(T)}
    # t 时 PIES 从电网买电量
    P_b_gr_t = {t: model.continuous_var(lb=0, ub=P_b_gr_max, name=f'P_b_gr_{t}') for t in range(T)}
    # t 时 PIES 向电网卖电量
    P_s_gr_t = {t: model.continuous_var(lb=0, ub=P_s_gr_max, name=f'P_s_gr_{t}') for t in range(T)}
    # t 时转移进的输入负载
    P_TL_in_t = {t: model.continuous_var(name=f'P_TL_in_{t}') for t in range(T)}
    # t 时转移出的输出负载
    P_TL_out_t = {t: model.continuous_var(name=f'P_TL_out_{t}') for t in range(T)}
    # t 时中断负载
    P_IL_t = {t: model.continuous_var(lb=P_IL_t_min, ub=P_IL_t_max, name=f'P_IL_{t}') for t in range(T)}
    # t 时注入建筑 w 的冷空气量
    C_IL_c_w_t ={(w, t): model.continuous_var(name=f'C_IL_c_{w}_{t}') for w in range(W) for t in range(T)}
    # t 是光伏输出电量
    P_PV_t = {t: model.continuous_var(name=f'P_PV_{t}') for t in range(T)}
    # 连续变量
    V = [
        P_MTG_i_t,
        P_ESS_ch_t,
        P_ESS_dc_t,
        P_PV_t,
        P_b_gr_t,
        P_s_gr_t,
        P_TL_in_t,
        P_TL_out_t,
        P_IL_t,
        C_IL_c_w_t
    ]

    # ======================✅目标函数=====================
    # ======================计算一阶段目标函数=====================
    # MTG 的开关成本
    U_i = [0 for _ in range(n_G)]
    D_i = [0 for _ in range(n_G)]
    # 工厂机器的生产维护成本
    C_M_t = 0
    # 一阶段目标函数
    # 总成本=所有的 MTG 的开关成本 + 工厂机器的维护成本
    obj_first_stage = model.sum(
        C_M_t + model.sum(U_i[i] * y_i_t[i, t] + D_i[i] * z_i_t[i, t] for i in range(n_G)) for t in range(T))
    # =====================计算 MTG 操作成本=====================
    # 主系数
    a_i = [0 for i in range(n_G)]
    # 常系数
    b_i = [0 for i in range(n_G)]
    # MTG 操作成本
    C_MTG = model.sum(model.sum(a_i[i] * P_MTG_i_t[i, t] + b_i * u_MTG_i_t[i, t] for i in range(n_G)) for t in range(T))
    # =====================计算 ESS 操作成本=====================

    # ESS 充/放电单价
    epsilon_ESS = 0

    # ESS 操作成本
    C_ESS = model.sum(epsilon_ESS * (P_ESS_ch_t[t] + P_ESS_dc_t[t]) for t in range(T))

    # =====================计算 CO2 排放惩罚成本=====================
    # 光伏削减惩罚成本
    epsilon_CO2 = 0
    # MTG 生产每单位电量的 COR2 排放量
    k_MTG = 0
    # 电网生产每单位电量的 COR2 排放量
    k_gr = 0
    # CO2 排放惩罚成本
    C_CO2 = model.sum(
        epsilon_CO2 * model.sum(k_MTG * model.sum(P_MTG_i_t[i, t] for i in range(n_G)) + k_gr * P_b_gr_t[t]) for t in
        range(T))
    # =====================计算 DR 补偿成本=====================
    # DR 补偿成本
    C_DR = model.sum(rho_TL * P_TL_in_t[t] + rho_IL * P_IL_t[t] + rho_c * C_IL_c_w_t[w,t] for w in range(W) for t in range(T))
    # =====================计算太阳能削减成本=====================
    # !
    C_PV_loss_t = {t: 0 for t in range(T)}
    # 太阳能削减成本
    C_loss = model.sum(C_PV_loss_t[t] for t in range(T))
    # =====================计算电网购电成本=====================
    # PIES 向电网购电单价
    epsilon_b_gr = 0
    # 电网购电成本
    C_b_gr = model.sum(epsilon_b_gr * P_b_gr_t[t] for t in range(T))
    # =====================计算电网卖电收入=====================
    # PIES 向电网卖电单价
    epsilon_s_gr = 0
    # 电网卖电收入
    C_s_gr = model.sum(epsilon_s_gr * P_s_gr_t[t] for t in range(T))
    # ======================计算二阶段目标函数=====================
    # 二阶段目标函数
    obj_second_stage = C_MTG + C_ESS + C_CO2 + C_DR + C_loss + C_b_gr - C_s_gr

    # ======================✅模型约束=====================
    # for i in range(n_G):
    #     for t in range(T):
    #         model.add_constraint(y_i_t[i, t] + z_i_t[i, t] <= 1)
    # 需求响应 - 电负载模型
    P_MTG_i_max = [0 for _ in range(n_G)]
    P_PV_ge_t = {t: 0 for t in range(T)}
    P_load_t = {t: 0 for t in range(T)}
    C_c_load_t = {t: 0 for t in range(T)}
    # 约束 (34)-(45)
    for i in range(n_G):
        for t in range(T):
            model.add_constraint(P_MTG_i_t[i, t] <= eta_MTG_i_t[i, t] * P_MTG_i_max)
            model.add_constraint(y_i_t[i, t] + z_i_t[i, t] <= 1)
            # !
            model.add_constraint(
                eta_MTG_i_t[i, t] * P_MTG_i_t[i, t + 1] - eta_MTG_i_t[i, t + 1] * P_MTG_i_t[i, t] <= delt_P_MTG_U_max
            )
            # !
            model.add_constraint(
                eta_MTG_i_t[i, t] * P_MTG_i_t[i, t - 1] - eta_MTG_i_t[i, t + 1] * P_MTG_i_t[i, t] <= delt_P_MTG_D_max
            )
    for t in range(T):
        model.add_constraint(
            P_ESS_ch_t[t] <= eta_ch_t[t] * P_ESS_ch_max
        )
        model.add_constraint(
            P_ESS_dc_t[t] <= eta_dc_t[t] * P_ESS_dc_max
        )
        model.add_constraint(
            eta_ch_t[t] + eta_dc_t[t] <= 1
        )
        model.add_constraint(
            C_ESS_0 + e_ch * model.sum(P_ESS_ch_t[tt] for tt in range(T)) - (1 / e_dc) * model.sum(
                P_ESS_dc_t[tt] for tt in range(T)) >= C_ESS_min
        )
        model.add_constraint(
            C_ESS_0 + e_ch * model.sum(P_ESS_ch_t[tt] for tt in range(T)) - (1 / e_dc) * model.sum(
                P_ESS_dc_t[tt] for tt in range(T)) <= C_ESS_max
        )
        model.add_constraint(
            C_ESS_0 == C_ESS_end
        )
        model.add_constraint(
            P_PV_ge_t[t] + P_b_gr_t[t] + P_ESS_dc_t[t] + model.sum(P_MTG_i_t[i, t] for i in range(n_G)) == P_load_t[t] +
            C_c_load_t[t] + P_ESS_ch_t[t]
        )
    # 需求响应约束
    # 电负载约束 (1)-(6)
    # TLs-转移负载，ILs-中断负载
    # 控制 TL in
    eta_TL_in_t = {t: model.binary_var(name=f'eta_TL_in_{t}') for t in range(T)}
    # 控制 TL out
    eta_TL_out_t = {t: model.binary_var(name=f'eta_TL_out_{t}') for t in range(T)}
    for t in range(T):
        # 约束 (1)
        model.add_constraint(P_TL_in_t[t] * eta_TL_in_t[t] >= P_TL_in_min)
        model.add_constraint(P_TL_in_t[t] * eta_TL_in_t[t] >= P_TL_in_min)
        # 约束 (2)
        model.add_constraint(P_TL_out_t[t] * eta_TL_out_t[t] <= P_TL_out_max)
        model.add_constraint(P_TL_out_t[t] * eta_TL_out_t[t] <= P_TL_out_max)
        # 约束 (3)
        model.add_constraint(eta_TL_in_t[t] + eta_TL_out_t[t] <= 1)
        # 约束 (5) P_IL_t 定义时添加上下界约束
        # 约束 (6)
        model.add_constraint(P_load_t[t] == P_load_t[0] + P_TL_in_t[t] - P_TL_out_t[t] - P_IL_t[t])
    # 约束 (4)
    model.add_constraint(model.sum(P_TL_in_t[t] - P_TL_out_t[t] for t in range(T)) == 0)

    # 生成维护约束
