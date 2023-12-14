"""
IC&CG 算法
"""
import numpy as np
from milp import solve_mp, solve_sp


def iccg():
    """iccg 算法"""
    # Initialization
    L_bar = 0
    U_bar = np.inf
    epsilon = np.random.uniform()
    epsilon_tilde = np.random.random() * epsilon / (1 + epsilon)
    sigma = np.random.random()
    epsilon_MP = []
    epsilon_MP.append(np.random.random())
    LB = []
    UB = []
    U = []
    k = 0
    l = 0
    while True:
        # MP solving
        solution_mp = solve_mp(L_bar, epsilon_MP[k])
        # 修改
        LB_k = solution_mp.objective_value-0.2
        LB.append(LB_k)
        UB_k = solution_mp.objective_value+solution_mp['gama']
        UB.append(UB_k)
        if LB_k > L_bar:
            l = k
        L_bar = UB[k]
        # SP solving
        solution_sp = solve_sp(u=solution_mp)
        U_bar = min(U_bar, solution_mp.objective_value +
                    solution_sp.objective_value)  # 修改
        if (U_bar - LB[l]) / U_bar >= epsilon:
            if (U_bar - UB[k]) / U_bar < epsilon_tilde:
                k = l
                L_bar = LB[l]
                for i in range(l, len(epsilon_MP)):
                    epsilon_MP[i] = sigma*epsilon_MP[i]
            else:
                # ! Ω的作用
                k += 1
            continue
        else:
            break
    return solution_mp.solve_details(), solution_sp.solve_details()
