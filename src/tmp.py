"""
临时存放
"""
from params import *
import numpy as np
from docplex.mp.model import Model
from docplex.mp.solution import SolveSolution


def other_constraint():
    pass
    # # =========================================================================================
    # # 其他约束：与 P_TL、P_IL、C_IL_c_t_w 相关的约束，论文未标出
    # # =========================================================================================
    # # todo: 需求响应约束:(1) 电负载约束 (1)~(6)
    # # TLs-转移负载，ILs-中断负载
    # # 控制 TL in
    # eta_TL_in_t = {t: model.binary_var(
    #     name=f'eta_TL_in_{t}') for t in range(T)}
    # # 控制 TL out
    # eta_TL_out_t = {t: model.binary_var(
    #     name=f'eta_TL_out_{t}') for t in range(T)}
    # for s in range(n_s):
    #     for t in range(T):
    #         # 约束 (1)
    #         model.add_constraint(
    #             P_TL_in_t[s, t] * eta_TL_in_t[t] >= P_TL_in_min)
    #         model.add_constraint(
    #             P_TL_in_t[s, t] * eta_TL_in_t[t] >= P_TL_in_min)
    #         # 约束 (2)
    #         model.add_constraint(
    #             P_TL_out_t[s, t] * eta_TL_out_t[t] <= P_TL_out_max)
    #         model.add_constraint(
    #             P_TL_out_t[s, t] * eta_TL_out_t[t] <= P_TL_out_max)
    #         # 约束 (3)
    #         model.add_constraint(eta_TL_in_t[t] + eta_TL_out_t[t] <= 1)
    #         # 约束 (5) P_IL_t 定义时添加上下界约束
    #         # 约束 (6)
    #         model.add_constraint(
    #             P_load_t[t] == P_load_t[0] + P_TL_in_t[s, t] - P_TL_out_t[s, t] - P_IL_t[s, t])
    #     # 约束 (4)
    #     model.add_constraint(
    #         model.sum(P_TL_in_t[s, t] - P_TL_out_t[s, t] for t in range(T)) == 0)
    # # =========================================================================================
    # # todo: 需求响应约束:(2) 冷负载约束 (7)~(15)
    # T_in_c_t_w = {(w, t): model.continuous_var()
    #               for w in range(W) for t in range(T)}
    # T_out_c_t_w = {(w, t): model.continuous_var()
    #                for w in range(W) for t in range(T)}
    # C_c_t_w = {(w, t): model.continuous_var()
    #            for w in range(W) for t in range(T)}
    # tol = 1
    # R = 1
    # KF = 1
    # C_air = 1
    # d_air = 1
    # V = 1
    # PMV = model.continuous_var()
    # T_in_a_w = {w: model.continuous_var() for w in range(W)}
    # H = 1
    # I_cl = 1
    # for w in range(W):
    #     for t in range(1, T):
    #         # 约束 (7)
    #         model.add_constraint(T_in_c_t_w[w, t] == T_in_c_t_w[w, t-1]*np.exp(-1/tol)+(
    #             T_out_c_t_w[w, t-1]-C_c_t_w[w, t-1]*R)*(1-np.exp(-1/tol)))
    #         # 约束 (8)
    #         # model.add_constraint(C_c_t_w[w, t-1] == T_out_c_t_w[w, t-1] - np.exp(-1/tol)*(
    #         #     T_in_c_t_w[w, t]-T_in_c_t_w[w, t-1])/(1-np.exp(-1/tol)))/R
    #         # # 约束 (9)
    #         model.add_constraint(C_c_t_w[w, t-1] == ((T_out_c_t_w[w, t]-T_in_c_t_w[w, t])+KF*(
    #             T_out_c_t_w[w, t]-T_in_c_t_w[w, t-1])/C_air*d_air*V)/(1/KF)+(1/(C_air*d_air*V)))
    #         # 约束 (10)
    #         model.add_constraint(PMV == 2.34 - 3.76 *
    #                              (T_in_a_w[w]-T_in_c_t_w[w, t])/(H*(I_cl+0.1)))
    # # 约束 (11)
    # for t in range(T):
    #     if 7 <= t <= 18:
    #         model.add_constraint(PMV == 0.5)
    #     else:
    #         model.add_constraint(PMV == 0.9)
    # # =========================================================================================
    # # todo: 生成维护约束 (16)~(25)
    # # 约束 (16)
    # F_bar_m = model.continuous_var()
    # for h in range(T):
    #     model.add_constraint(F_bar_m == model.sum(
    #         model.sum(v_m_t[m, t]*F_m for m in range(1, n_m)) for t in range(h, T)))
    # # (17)~(20) 已添加
    # # 约束 (21)
    # # 见 S_t_c 定义
    # S_c_min = 0
    # S_c_max = 10
    # S_t_c = {t: model.continuous_var(lb=S_c_min, ub=S_c_max) for t in range(T)}
    # # 约束 (22)
    # # 见 S_t_m 定义
    # S_m_min = 0
    # S_m_max = 10
    # S_t_m = {t: model.continuous_var(lb=S_m_min, ub=S_m_max) for t in range(T)}
    # for t in range(1, T):
    #     # 约束 (23)
    #     model.add_constraint(S_t_c[t] == S_t_c[t-1] + model.sum(v_c_t[c, t]*F_c for c in range(
    #         n_c) - model.sum(v_cp_t[cp, t]*F_cp*G_cp for cp in range(n_cp))))
    #     # 约束 (24)
    #     model.add_constraint(S_t_m[t] == S_t_m[t-1] + model.sum(v_mp_t[mp, t]*F_mp *
    #                          G_mp for mp in range(n_mp)) - model.sum(v_m_t[m, t]*F_m for m in range(n_m)))
    # # =========================================================================================
    # # todo: 光伏发电成本约束 (26)~(30)
    # # code
    # # =========================================================================================
    # # todo: 其他约束 (34)-(45)
    # P_MTG_i_max = [0 for _ in range(n_G)]
    # P_PV_ge_t = {t: 0 for t in range(T)}
    # P_load_t = {t: 0 for t in range(T)}
    # C_c_load_t = {t: 0 for t in range(T)}
    # for s in range(n_s):
    #     for i in range(n_G):
    #         for t in range(1, T-1):
    #             model.add_constraint(
    #                 P_MTG_i_t[s, i, t] <= eta_MTG_i_t[i, t] * P_MTG_i_max)
    #             model.add_constraint(y_i_t[i, t] + z_i_t[i, t] <= 1)
    #             # !
    #             model.add_constraint(
    #                 eta_MTG_i_t[i, t] * P_MTG_i_t[s, i, t + 1] - eta_MTG_i_t[i,
    #                                                                          t + 1] * P_MTG_i_t[s, i, t] <= delt_P_MTG_U_max
    #             )
    #             # !
    #             model.add_constraint(
    #                 eta_MTG_i_t[i, t] * P_MTG_i_t[s, i, t - 1] - eta_MTG_i_t[i,
    #                                                                          t + 1] * P_MTG_i_t[s, i, t] <= delt_P_MTG_D_max
    #             )
    # for s in range(n_s):
    #     for t in range(T):
    #         model.add_constraint(
    #             P_ESS_ch_t[s, t] <= eta_ch_t[t] * P_ESS_ch_max
    #         )
    #         model.add_constraint(
    #             P_ESS_dc_t[s, t] <= eta_dc_t[t] * P_ESS_dc_max
    #         )
    #         model.add_constraint(
    #             eta_ch_t[t] + eta_dc_t[t] <= 1
    #         )
    #         model.add_constraint(
    #             C_ESS_0 + e_ch * model.sum(P_ESS_ch_t[s, tt] for tt in range(T)) - (1 / e_dc) * model.sum(
    #                 P_ESS_dc_t[s, tt] for tt in range(T)) >= C_ESS_min
    #         )
    #         model.add_constraint(
    #             C_ESS_0 + e_ch * model.sum(P_ESS_ch_t[s, tt] for tt in range(T)) - (1 / e_dc) * model.sum(
    #                 P_ESS_dc_t[s, tt] for tt in range(T)) <= C_ESS_max
    #         )
    #         model.add_constraint(
    #             C_ESS_0 == C_ESS_end
    #         )
    #         model.add_constraint(
    #             P_PV_ge_t[t] + P_b_gr_t[s, t] + P_ESS_dc_t[s, t] + model.sum(P_MTG_i_t[s, i, t] for i in range(n_G)) == P_load_t[s, t] +
    #             C_c_load_t[s, t] + P_ESS_ch_t[s, t]
    #         )
