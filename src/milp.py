"""
PIES 主问题与子问题数学规划求解
"""
from params import *
import numpy as np
from docplex.mp.model import Model
from docplex.mp.solution import SolveSolution


def solve_mp(L_bar=0.5, epsilon_mp=0.01) -> SolveSolution:
    """
    求解主问题

    参数：
        L_bar: Au≥L_bar 右端向量对应的值
        epsilon_mp: 求解的容差
    返回：
        主问题的最优解
    """
    model = Model('mp')
    model.parameters.mip.tolerances.mipgap = epsilon_mp
    # =========================================================================================
    # ✅一阶段决策变量 U (46)
    # =========================================================================================
    # 控制 MTG 开关状态
    y_i_t = {(i, t): model.binary_var(name=f'y_{i}_{t}')
             for i in range(n_G) for t in range(T)}
    z_i_t = {(i, t): model.binary_var(name=f'z_{i}_{t}')
             for i in range(n_G) for t in range(T)}
    # 控制 ESS 放电状态
    eta_dc_t = {t: model.binary_var(name=f'eta_dc_{t}') for t in range(T)}
    # 控制 ESS 充电状态
    eta_ch_t = {t: model.binary_var(name=f'eta_ch_{t}') for t in range(T)}
    # 控制生产线运行状态-PV module encapsulation 产线
    v_m_t = {(m, t): model.binary_var(name=f'v_m{m}_{t}')
             for m in range(n_m) for t in range(T)}
    v_mp_t = {(mp, t): model.binary_var(name=f'v_mp{mp}_{t}')
              for mp in range(n_mp) for t in range(T)}
    # 控制生产线运行状态-PV cell processing 产线
    v_c_t = {(c, t): model.binary_var(name=f'v_c{c}_{t}')
             for c in range(n_c) for t in range(T)}
    v_cp_t = {(cp, t): model.binary_var(name=f'v_cp{cp}_{t}')
              for cp in range(n_cp) for t in range(T)}
    # 控制 MTG 开关状态
    eta_MTG_i_t = {(i, t): model.binary_var(name=f'eta_MTG_{i}_{t}')
                   for i in range(n_G) for t in range(T)}

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
    # =========================================================================================
    # ✅二阶段决策变量 V
    # =========================================================================================
    # t 时第 i 个 MTG 的生产需要的能量
    P_MTG_i_t = {(s, i, t): model.continuous_var(name=f'P_MTG_{s}_{i}_{t}')
                 for s in range(n_s) for i in range(n_G) for t in range(T)}
    # t 时 ESS 充电量
    P_ESS_ch_t = {(s, t): model.continuous_var(
        name=f'P_ESS_ch_{s}_{t}') for s in range(n_s) for t in range(T)}
    # t 时 ESS 放电量
    P_ESS_dc_t = {(s, t): model.continuous_var(
        name=f'P_ESS_dc_{s}_{t}') for s in range(n_s) for t in range(T)}
    # t 时 PIES 从电网买电量
    P_b_gr_t = {(s, t): model.continuous_var(
        lb=0, ub=P_b_gr_max, name=f'P_b_gr_{s}_{t}') for s in range(n_s) for t in range(T)}
    # t 时 PIES 向电网卖电量
    P_s_gr_t = {(s, t): model.continuous_var(
        lb=0, ub=P_s_gr_max, name=f'P_s_gr_{s}_{t}') for s in range(n_s) for t in range(T)}
    # t 是光伏输出电量
    P_PV_t = {(s, t): model.continuous_var(name=f'P_PV_{s}_{t}')
              for s in range(n_s) for t in range(T)}
    # t 时转移进的输入负载
    P_TL_in_t = {(s, t): model.continuous_var(
        name=f'P_TL_in_{s}_{t}') for s in range(n_s) for t in range(T)}
    # t 时转移出的输出负载
    P_TL_out_t = {(s, t): model.continuous_var(
        name=f'P_TL_out_{s}_{t}') for s in range(n_s) for t in range(T)}
    # t 时中断负载
    P_IL_t = {(s, t): model.continuous_var(
        lb=P_IL_t_min, ub=P_IL_t_max, name=f'P_IL_{s}_{t}') for s in range(n_s) for t in range(T)}
    # t 时注入建筑 w 的冷空气量
    C_IL_c_w_t = {(s, w, t): model.continuous_var(name=f'C_IL_c_{s}_{w}_{t}')
                  for s in range(n_s) for w in range(W) for t in range(T)}
    # t 时第 i 个 MTG 的维护需要的能量
    u_MTG_i_t = {(s, i, t): 1 for s in range(n_s)
                 for i in range(n_G) for t in range(T)}
    # 连续变量
    V_vars = [
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
    gama = model.continuous_var(name='gama')
    # =========================================================================================
    # ✅一阶段目标函数 min:Au+γ (32)
    # =========================================================================================
    # MTG 的开关成本
    U_i = [1 for _ in range(n_G)]
    D_i = [1 for _ in range(n_G)]
    # 工厂机器的生产维护成本 (25)
    # C_c C_m C_cp C_mp 4条产线
    C_c = 1
    C_m = 1
    C_cp = 1
    C_mp = 1
    C_M_t = {t: model.sum(v_c_t[c, t]*C_c for c in range(n_c)) + model.sum(v_m_t[m, t]*C_m for m in range(n_m)) +
             model.sum(v_cp_t[cp, t]*C_cp for cp in range(n_cp)) +
             model.sum(v_c_t[mp, t]*C_mp for mp in range(n_mp)) for t in range(T)}
    # 一阶段目标函数 Au+γ
    # 总成本=所有的 MTG 的开关成本 + 工厂机器的维护成本
    obj_first_stage = model.sum(
        C_M_t[t] + model.sum(U_i[i] * y_i_t[i, t] + D_i[i] * z_i_t[i, t] for i in range(n_G)) for t in range(T)) + gama
    # =========================================================================================
    # ✅二阶段目标函数 ∑{B·vs}
    # =========================================================================================
    # 主系数
    a_i = [1 for i in range(n_G)]
    # 常系数
    b_i = [1 for i in range(n_G)]
    # MTG 操作成本
    C_MTG = model.sum(
        model.sum(
            model.sum(
                a_i[i] * P_MTG_i_t[s, i, t] + b_i[i] * u_MTG_i_t[s, i, t] for i in range(n_G)
            ) for t in range(T)
        ) for s in range(n_s)
    )
    # 计算 ESS 操作成本
    # ESS 充/放电单价
    epsilon_ESS = 0.1
    # ESS 操作成本
    C_ESS = model.sum(model.sum(
        epsilon_ESS * (P_ESS_ch_t[s, t] + P_ESS_dc_t[s, t]) for t in range(T)) for s in range(n_s))
    # 计算 CO2 排放惩罚成本
    # 光伏削减惩罚成本
    epsilon_CO2 = 0.1
    # MTG 生产每单位电量的 COR2 排放量
    k_MTG = 1
    # 电网生产每单位电量的 COR2 排放量
    k_gr = 1
    # CO2 排放惩罚成本
    C_CO2 = model.sum(model.sum(
        epsilon_CO2 * model.sum(k_MTG * model.sum(P_MTG_i_t[s, i, t] for i in range(n_G)) + k_gr * P_b_gr_t[s, t]) for t in
        range(T)) for s in range(n_s))
    # 计算 DR 补偿成本
    # DR 补偿成本
    C_DR = model.sum(model.sum(
        rho_TL * P_TL_in_t[s, t] + rho_IL * P_IL_t[s, t] + rho_c * C_IL_c_w_t[s, w, t] for w in range(W) for t in range(T)) for s in range(n_s))
    # 计算太阳能削减成本
    C_PV_loss_t = {(s, t): 0 for s in range(n_s) for t in range(T)}
    # 太阳能削减成本
    C_loss = model.sum(model.sum(C_PV_loss_t[s, t]
                       for t in range(T)) for s in range(n_s))
    # 计算电网购电成本
    # PIES 向电网购电单价
    epsilon_b_gr = 0.1
    # 电网购电成本
    C_b_gr = model.sum(
        model.sum(epsilon_b_gr * P_b_gr_t[s, t] for t in range(T)) for s in range(n_s))
    # 计算电网卖电收入
    # PIES 向电网卖电单价
    epsilon_s_gr = 0.1
    # 电网卖电收入
    C_s_gr = model.sum(
        model.sum(epsilon_s_gr * P_s_gr_t[s, t] for t in range(T)) for s in range(n_s))
    # 计算二阶段目标函数
    # 二阶段目标函数
    # ∑Bv
    p_s = 0.2
    obj_second_stage = C_MTG + C_ESS + C_CO2 + C_DR + C_loss + C_b_gr - C_s_gr
    # =========================================================================================
    # ✅约束
    # =========================================================================================
    # todo: γ≥∑pBv
    model.add_constraint(
        gama >= p_s * obj_second_stage
    )
    # =========================================================================================
    # todo: Au+γ≥L_bar
    model.add_constraint(
        obj_first_stage + gama >= L_bar
    )
    # =========================================================================================
    # todo: Cu=c
    # 无此类约束
    # =========================================================================================
    # todo: Du≤d
    # 约束 (35)
    for i in range(n_G):
        for t in range(T):
            model.add_constraint(
                y_i_t[i, t] + z_i_t[i, t] <= 1
            )
    # 约束 (40)
    for t in range(T):
        model.add_constraint(
            eta_ch_t[t] + eta_dc_t[t] <= 1
        )
    # (17)~(20)
    # 约束 (17)
    h = 0
    for m in range(n_m):
        model.add_constraint(
            model.sum(v_m_t[m, t] for t in range(h, T)) >= T_m_min
        )
        model.add_constraint(
            model.sum(v_m_t[m, t] for t in range(h, T)) <= T_m_max
        )
    for c in range(n_c):
        model.add_constraint(
            model.sum(v_c_t[c, t] for t in range(h, T)) >= T_c_min
        )
        model.add_constraint(
            model.sum(v_c_t[c, t] for t in range(h, T)) <= T_c_max
        )

    # 约束 (18)
    a = 0
    b = 0
    for h in range(T):
        model.add_constraint(model.sum(model.sum(v_cp_t[cp, t]*F_cp*G_cp for cp in range(1, n_cp))
                             for t in range(1, h+a)) <= model.sum(model.sum(v_c_t[c, t]*F_c for c in range(1, n_c))
                             for t in range(1, h)))
    # 约束 (19)
    for h in range(T):
        model.add_constraint(model.sum(model.sum(v_m_t[m, t]*F_m for m in range(1, n_m))
                             for t in range(1, h+a)) <= model.sum(model.sum(v_cp_t[cp, t]*F_cp*G_cp for cp in range(1, n_cp))
                             for t in range(1, h)))
    # 约束 (20)
    for h in range(T):
        model.add_constraint(model.sum(model.sum(v_mp_t[mp, t]*F_mp*G_mp for mp in range(1, n_mp))
                             for t in range(1, h+1)) <= model.sum(model.sum(v_m_t[m, t]*F_m for m in range(1, n_m))
                             for t in range(1, h)))
    # =========================================================================================
    # todo: Ev=e
    # 约束 (42)
    for s in range(n_s):
        C_ESS_0 = epsilon_ESS * (P_ESS_ch_t[s, 0] + P_ESS_dc_t[s, 0])
        C_ESS_end = epsilon_ESS * (P_ESS_ch_t[s, T-1] + P_ESS_dc_t[s, T-1])
        model.add_constraint(
            C_ESS_0 == C_ESS_end
        )
    # =========================================================================================
    # todo: Fv≤f
    # 约束 (34)
    for s in range(n_s):
        for i in range(n_G):
            for t in range(T):
                model.add_constraint(
                    P_MTG_i_t[s, i, t] <= eta_MTG_i_t[i, t] * P_MTG_i_max)

    # 约束 (36)
    for i in range(n_G):
        for t in range(T-1):
            for s in range(n_s):
                model.add_constraint(
                    eta_MTG_i_t[i, t] * P_MTG_i_t[s, i, t + 1] -
                    eta_MTG_i_t[i, t + 1] *
                    P_MTG_i_t[s, i, t] <= delt_P_MTG_U_max
                )
    # 约束 (37)
    for i in range(n_G):
        for t in range(1, T-1):
            for s in range(n_s):
                model.add_constraint(
                    eta_MTG_i_t[i, t] * P_MTG_i_t[s, i, t - 1] -
                    eta_MTG_i_t[i, t + 1] *
                    P_MTG_i_t[s, i, t] <= delt_P_MTG_D_max
                )
    # 约束 (44)
    # 见 P_b_gr_t 定义
    # 约束 (45)
    # 见 P_s_gr_t 定义
    # =========================================================================================
    # todo: Gu+Hv≤g
    # 约束 (38)
    for s in range(n_s):
        for t in range(T):
            model.add_constraint(
                P_ESS_ch_t[s, t] <= eta_ch_t[t] * P_ESS_ch_max
            )
    # 约束 (39)
    for s in range(n_s):
        for t in range(T):
            model.add_constraint(
                P_ESS_dc_t[s, t] <= eta_dc_t[t] * P_ESS_dc_max
            )
    # 约束 (43)
    P_PV_ge_t = {t: 10 for t in range(T)}
    # 约束(15)
    # 冷负载(输入场景)
    C_c_load_t = {t: 10 for t in range(T)}
    # 约束(6)
    # 电负载
    P_load_t = {t: 10 for t in range(T)}
    # 电负载约束 (1)~(6)
    # TLs-转移负载，ILs-中断负载
    # 控制 TL in
    eta_TL_in_t = {t: model.binary_var(
        name=f'eta_TL_in_{t}') for t in range(T)}
    # 控制 TL out
    eta_TL_out_t = {t: model.binary_var(
        name=f'eta_TL_out_{t}') for t in range(T)}
    for s in range(n_s):
        for t in range(T):
            # 约束 (1)
            model.add_constraint(
                P_TL_in_t[s, t] * eta_TL_in_t[t] >= P_TL_in_min)
            model.add_constraint(
                P_TL_in_t[s, t] * eta_TL_in_t[t] >= P_TL_in_min)
            # 约束 (2)
            model.add_constraint(
                P_TL_out_t[s, t] * eta_TL_out_t[t] <= P_TL_out_max)
            model.add_constraint(
                P_TL_out_t[s, t] * eta_TL_out_t[t] <= P_TL_out_max)
            # 约束 (3)
            model.add_constraint(eta_TL_in_t[t] + eta_TL_out_t[t] <= 1)
            # 约束 (5) P_IL_t 定义时添加上下界约束
            # 约束 (6)
            model.add_constraint(
                P_load_t[t] == P_load_t[0] + P_TL_in_t[s, t] - P_TL_out_t[s, t] - P_IL_t[s, t])
        # 约束 (4)
        model.add_constraint(
            model.sum(P_TL_in_t[s, t] - P_TL_out_t[s, t] for t in range(T)) == 0)

    for s in range(n_s):
        for t in range(T):
            model.add_constraint(
                P_PV_ge_t[t] + P_b_gr_t[s, t] + P_ESS_dc_t[s, t] + model.sum(
                    P_MTG_i_t[s, i, t] for i in range(n_G)) == P_load_t[t] + C_c_load_t[t] + P_ESS_ch_t[s, t]
            )
    # =========================================================================================
    # 极小化目标
    # =========================================================================================
    model.minimize(obj_first_stage+gama)
    solution: SolveSolution = model.solution
    if solution:
        # print(model.solve_details())
        return solution
    else:
        print('mp问题无解')
        exit(-1)


def solve_sp(u: SolveSolution) -> SolveSolution:
    """
    求解子问题

    参数：
        u: 主问题 u 的最优解
    返回：
        子问题的最优解
    """
    model = Model('sp')
    # =========================================================================================
    # ✅一阶段决策变量 U
    # =========================================================================================
    # 控制 MTG 开关状态
    y_i_t = {(i, t): u[f'y_{i}_{t}']
             for i in range(n_G) for t in range(T)}
    z_i_t = {(i, t): u[f'z_{i}_{t}']
             for i in range(n_G) for t in range(T)}
    # 控制 ESS 放电状态
    eta_dc_t = {t: u[f'eta_dc_{t}'] for t in range(T)}
    # 控制 ESS 充电状态
    eta_ch_t = {t: u[f'eta_ch_{t}'] for t in range(T)}
    # 控制生产线运行状态-PV module encapsulation 产线
    v_m_t = {(m, t): u[f'v_m{m}_{t}']
             for m in range(n_m) for t in range(T)}
    v_mp_t = {(mp, t): u[f'v_mp{mp}_{t}']
              for mp in range(n_mp) for t in range(T)}
    # 控制生产线运行状态-PV cell processing 产线
    v_c_t = {(c, t): u[f'v_c{c}_{t}']
             for c in range(n_c) for t in range(T)}
    v_cp_t = {(cp, t): u[f'v_cp{cp}_{t}']
              for cp in range(n_cp) for t in range(T)}
    # 控制 MTG 开关状态
    eta_MTG_i_t = {(i, t): u[f'eta_MTG_{i}_{t}']
                   for i in range(n_G) for t in range(T)}

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
    # =========================================================================================
    # ✅二阶段决策变量 V
    # =========================================================================================
    # t 时第 i 个 MTG 的生产需要的能量
    P_MTG_i_t = {(s, i, t): model.continuous_var(name=f'P_MTG_{s}_{i}_{t}')
                 for s in range(n_s) for i in range(n_G) for t in range(T)}
    # t 时 ESS 充电量
    P_ESS_ch_t = {(s, t): model.continuous_var(
        name=f'P_ESS_ch_{s}_{t}') for s in range(n_s) for t in range(T)}
    # t 时 ESS 放电量
    P_ESS_dc_t = {(s, t): model.continuous_var(
        name=f'P_ESS_dc_{s}_{t}') for s in range(n_s) for t in range(T)}
    # t 时 PIES 从电网买电量
    P_b_gr_t = {(s, t): model.continuous_var(
        lb=0, ub=P_b_gr_max, name=f'P_b_gr_{s}_{t}') for s in range(n_s) for t in range(T)}
    # t 时 PIES 向电网卖电量
    P_s_gr_t = {(s, t): model.continuous_var(
        lb=0, ub=P_s_gr_max, name=f'P_s_gr_{s}_{t}') for s in range(n_s) for t in range(T)}
    # t 是光伏输出电量
    P_PV_t = {(s, t): model.continuous_var(name=f'P_PV_{s}_{t}')
              for s in range(n_s) for t in range(T)}
    # t 时转移进的输入负载
    P_TL_in_t = {(s, t): model.continuous_var(
        name=f'P_TL_in_{s}_{t}') for s in range(n_s) for t in range(T)}
    # t 时转移出的输出负载
    P_TL_out_t = {(s, t): model.continuous_var(
        name=f'P_TL_out_{s}_{t}') for s in range(n_s) for t in range(T)}
    # t 时中断负载
    P_IL_t = {(s, t): model.continuous_var(
        lb=P_IL_t_min, ub=P_IL_t_max, name=f'P_IL_{s}_{t}') for s in range(n_s) for t in range(T)}
    # t 时注入建筑 w 的冷空气量
    C_IL_c_w_t = {(s, w, t): model.continuous_var(name=f'C_IL_c_{s}_{w}_{t}')
                  for s in range(n_s) for w in range(W) for t in range(T)}
    # t 时第 i 个 MTG 的维护需要的能量
    u_MTG_i_t = {(s, i, t): 1 for s in range(n_s)
                 for i in range(n_G) for t in range(T)}
    # 连续变量
    V_vars = [
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
    # =========================================================================================
    # ✅二阶段目标函数 ∑{B·vs}
    # =========================================================================================
    # 主系数
    a_i = [1 for i in range(n_G)]
    # 常系数
    b_i = [1 for i in range(n_G)]
    # MTG 操作成本
    C_MTG = model.sum(
        model.sum(
            model.sum(
                a_i[i] * P_MTG_i_t[s, i, t] + b_i[i] * u_MTG_i_t[s, i, t] for i in range(n_G)
            ) for t in range(T)
        ) for s in range(n_s)
    )
    # 计算 ESS 操作成本
    # ESS 充/放电单价
    epsilon_ESS = 0.1
    # ESS 操作成本
    C_ESS = model.sum(model.sum(
        epsilon_ESS * (P_ESS_ch_t[s, t] + P_ESS_dc_t[s, t]) for t in range(T)) for s in range(n_s))
    # 计算 CO2 排放惩罚成本
    # 光伏削减惩罚成本
    epsilon_CO2 = 0.1
    # MTG 生产每单位电量的 COR2 排放量
    k_MTG = 1
    # 电网生产每单位电量的 COR2 排放量
    k_gr = 1
    # CO2 排放惩罚成本
    C_CO2 = model.sum(model.sum(
        epsilon_CO2 * model.sum(k_MTG * model.sum(P_MTG_i_t[s, i, t] for i in range(n_G)) + k_gr * P_b_gr_t[s, t]) for t in
        range(T)) for s in range(n_s))
    # 计算 DR 补偿成本
    # DR 补偿成本
    C_DR = model.sum(model.sum(
        rho_TL * P_TL_in_t[s, t] + rho_IL * P_IL_t[s, t] + rho_c * C_IL_c_w_t[s, w, t] for w in range(W) for t in range(T)) for s in range(n_s))
    # 计算太阳能削减成本
    C_PV_loss_t = {(s, t): 0 for s in range(n_s) for t in range(T)}
    # 太阳能削减成本
    C_loss = model.sum(model.sum(C_PV_loss_t[s, t]
                       for t in range(T)) for s in range(n_s))
    # 计算电网购电成本
    # PIES 向电网购电单价
    epsilon_b_gr = 0.1
    # 电网购电成本
    C_b_gr = model.sum(
        model.sum(epsilon_b_gr * P_b_gr_t[s, t] for t in range(T)) for s in range(n_s))
    # 计算电网卖电收入
    # PIES 向电网卖电单价
    epsilon_s_gr = 0.1
    # 电网卖电收入
    C_s_gr = model.sum(
        model.sum(epsilon_s_gr * P_s_gr_t[s, t] for t in range(T)) for s in range(n_s))
    # 计算二阶段目标函数
    # 二阶段目标函数
    # ∑Bv
    p_s = 0.2
    obj_second_stage = C_MTG + C_ESS + C_CO2 + C_DR + C_loss + C_b_gr - C_s_gr
    # =========================================================================================
    # ✅约束
    # =========================================================================================
    # todo: Ev=e
    # 约束 (42)
    for s in range(n_s):
        C_ESS_0 = epsilon_ESS * (P_ESS_ch_t[s, 0] + P_ESS_dc_t[s, 0])
        C_ESS_end = epsilon_ESS * (P_ESS_ch_t[s, T-1] + P_ESS_dc_t[s, T-1])
        model.add_constraint(
            C_ESS_0 == C_ESS_end
        )
    # =========================================================================================
    # todo: Fv≤f
    # 约束 (34)
    for s in range(n_s):
        for i in range(n_G):
            for t in range(T):
                model.add_constraint(
                    P_MTG_i_t[s, i, t] <= eta_MTG_i_t[i, t] * P_MTG_i_max)

    # 约束 (36)
    for i in range(n_G):
        for t in range(T-1):
            for s in range(n_s):
                model.add_constraint(
                    eta_MTG_i_t[i, t] * P_MTG_i_t[s, i, t + 1] -
                    eta_MTG_i_t[i, t + 1] *
                    P_MTG_i_t[s, i, t] <= delt_P_MTG_U_max
                )
    # 约束 (37)
    for i in range(n_G):
        for t in range(1, T-1):
            for s in range(n_s):
                model.add_constraint(
                    eta_MTG_i_t[i, t] * P_MTG_i_t[s, i, t - 1] -
                    eta_MTG_i_t[i, t + 1] *
                    P_MTG_i_t[s, i, t] <= delt_P_MTG_D_max
                )
    # 约束 (44)
    # 见 P_b_gr_t 定义
    # 约束 (45)
    # 见 P_s_gr_t 定义
    # =========================================================================================
    # todo: Gu+Hv≤g
    # 约束 (38)
    for s in range(n_s):
        for t in range(T):
            model.add_constraint(
                P_ESS_ch_t[s, t] <= eta_ch_t[t] * P_ESS_ch_max
            )
    # 约束 (39)
    for s in range(n_s):
        for t in range(T):
            model.add_constraint(
                P_ESS_dc_t[s, t] <= eta_dc_t[t] * P_ESS_dc_max
            )
    # 约束 (43)
    P_PV_ge_t = {t: 10 for t in range(T)}
    # 约束(15)
    # 冷负载(输入场景)
    C_c_load_t = {t: 10 for t in range(T)}
    # 约束(6)
    # 电负载
    P_load_t = {t: 10 for t in range(T)}
    # 电负载约束 (1)~(6)
    # TLs-转移负载，ILs-中断负载
    # 控制 TL in
    eta_TL_in_t = {t: model.binary_var(
        name=f'eta_TL_in_{t}') for t in range(T)}
    # 控制 TL out
    eta_TL_out_t = {t: model.binary_var(
        name=f'eta_TL_out_{t}') for t in range(T)}
    for s in range(n_s):
        for t in range(T):
            # 约束 (1)
            model.add_constraint(
                P_TL_in_t[s, t] * eta_TL_in_t[t] >= P_TL_in_min)
            model.add_constraint(
                P_TL_in_t[s, t] * eta_TL_in_t[t] >= P_TL_in_min)
            # 约束 (2)
            model.add_constraint(
                P_TL_out_t[s, t] * eta_TL_out_t[t] <= P_TL_out_max)
            model.add_constraint(
                P_TL_out_t[s, t] * eta_TL_out_t[t] <= P_TL_out_max)
            # 约束 (3)
            model.add_constraint(eta_TL_in_t[t] + eta_TL_out_t[t] <= 1)
            # 约束 (5) P_IL_t 定义时添加上下界约束
            # 约束 (6)
            model.add_constraint(
                P_load_t[t] == P_load_t[0] + P_TL_in_t[s, t] - P_TL_out_t[s, t] - P_IL_t[s, t])
        # 约束 (4)
        model.add_constraint(
            model.sum(P_TL_in_t[s, t] - P_TL_out_t[s, t] for t in range(T)) == 0)

    for s in range(n_s):
        for t in range(T):
            model.add_constraint(
                P_PV_ge_t[t] + P_b_gr_t[s, t] + P_ESS_dc_t[s, t] + model.sum(
                    P_MTG_i_t[s, i, t] for i in range(n_G)) == P_load_t[t] + C_c_load_t[t] + P_ESS_ch_t[s, t]
            )
    # =========================================================================================
    # 极大化目标
    # =========================================================================================
    #! 求和和概率怎么理解
    model.maximize(model.sum(p_s*model.min(obj_second_stage)))
    solution: SolveSolution = model.solution
    if solution:
        # print(model.solve_details())
        return solution
    else:
        print('sp问题无解')
        exit(-1)


if __name__ == "__main__":
    print(solve_mp())
    print(solve_sp(solve_mp()))
