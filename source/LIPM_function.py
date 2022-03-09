import math
import numpy as np
from datetime import datetime
from pathlib import Path
import time
import matplotlib.pyplot as plt
import os

def modifiable_x_OSG_RampDSP_RampSSP(dt, wn_T, A, a, T, k, BoolPlot):
    vec_sum = []
    if T <= 0:
        print("input error: T must be greater than zero")
    else:
        N = len(A)
        size = (1, N)
        d1 = np.zeros(size, np.float64)
        d2 = np.zeros(size, np.float64)
        d1[0][0] = 0
        d2[0][N-1] = 0

        for i in range(N-1):
            d2[0][i] = (A[i + 1] - A[i]) / 2
            d1[0][i+1] = d2[0][i]

        t_length = math.floor(N*T/dt)
        totalT = t_length * dt

        for i in range(0, N):
            xn_vn_pn = x_OSCG_ZMPBC2_RampDSP_RampSSP(dt, wn_T, A[i], a[i], d1[0][i], d2[0][i], k, T, 0, False)

            for j in range(0, len(xn_vn_pn)):
                vec_three = []
                [xn, vn, pn] = xn_vn_pn[j]
                vec_three.append(xn)
                vec_three.append(vn)
                vec_three.append(pn)

                vec_sum.append(vec_three)

    return vec_sum

def modifiable_y_OSG_RampDSP_SinSSP(dt, wn_T, B, b, T, k, BoolFinalZero, BoolPlot):
    # ok
    vec_sum = []
    if T <= 0:
        print("input error: T must be greater than zero")
    else:
        N = len(B)
        size = (1, N)
        d1 = np.zeros(size, np.float64)
        d2 = np.zeros(size, np.float64)
        d1[0][0] = 0
        if BoolFinalZero:
            d2[0][N-1] = 0
        else:
            d2[0][N-1] = (B[N-1] + B[N-2]) / 2
        
        for i in range(0, N-1):
            d2[0][i] = (B[i+1] + B[i])/2
            d1[0][i+1] = d2[0][i]
        t_length = math.floor(N*T/dt)
        totalT = t_length * dt

        for i in range(0, N):
            yn_vn_qn = y_OSCG_ZMPBC2_RampDSP_SinSSP(dt, wn_T, B[i], b[i], d1[0][i], d2[0][i], k, T, 0, False)
            for j in range(len(yn_vn_qn)):
                vec_three = []

                [yn, vn, qn] = yn_vn_qn[j]

                vec_three.append(yn)
                vec_three.append(vn)
                vec_three.append(qn)

                vec_sum.append(vec_three)
    return vec_sum

def x_OSCG_ZMPBC2_RampDSP_RampSSP(dt, wn_T, A, a, d1, d2, k, T, Bias_x, BoolPlot):
    # ok
    vec_sum = []
    g = 9.8
    sum_wn_T = sum(wn_T)

    wn = sum_wn_T/ len(wn_T)

    vec_wn_temp = []
    vec_wn_temp.append(wn)
    vec_T_temp = []
    vec_T_temp.append(T)

    if k>0:
        # run
        t1 = k * T
        t2 = (1 - k)*T
        L1 = A - d1
        K1 = (d1 - a) / (k*T)
        L2 = A-a-2*a*t1/(1-2*k)/T
        K2 = 2 * a / ((1 - (2 * k)) * T)
        L3 = A + a - ((d2 - a)*t2) / (k * T)
        K3 = (d2 - a) / (k * T)

        P1S_four = StepZMP2CoM(vec_T_temp, vec_wn_temp, L1, 0) #P1S,X1S,V1S,A1S
        P1R_four = RampZMP2CoM(vec_T_temp, vec_wn_temp, K1, 0) #P1R,X1R,V1R,A1R
        P2S_four = StepZMP2CoM(vec_T_temp, vec_wn_temp, L2 - L1, t1) #P2S,X2S,V2S,A2S
        P2R_four = RampZMP2CoM(vec_T_temp, vec_wn_temp, K2 - K1, t1) #P2R,X2R,V2R,A2R
        P3S_four = StepZMP2CoM(vec_T_temp, vec_wn_temp, L3 - L2, t2) #P3S,X3S,V3S,A3S
        P3R_four = RampZMP2CoM(vec_T_temp, vec_wn_temp, K3 - K2, t2) #P3R,X3R,V3R,A3R

        X1S = P1S_four[0][1]
        X1R = P1R_four[0][1]
        X2S = P2S_four[0][1]
        X2R = P2R_four[0][1]
        X3S = P3S_four[0][1]
        X3R = P3R_four[0][1]

        X1 = X1S + X1R
        X2 = X2S + X2R
        X3 = X3S + X3R

        C = math.cosh(wn*T)
        S = math.sinh(wn*T)
        x0 = A - d1
        xf = A + d2
        x_v0 = (xf - x0 * C - X1 - X2 - X3)*wn / S

    else:
        alpha = 2 * math.pi / T
        L1 = A - a
        K1 = 2 * a / T
        P1S_four = StepZMP2CoM(vec_T_temp, vec_wn_temp, L1, 0) #P1S,X1S,V1S,A1S  P1S[0][1]
        P1R_four = RampZMP2CoM(vec_T_temp, vec_wn_temp, K1, 0) #P1R,X1R,V1R,A1R
        X1S = P1S_four[0][1]
        X1R = P1R_four[0][1]
        X1 = X1S + X1R
        X2 = 0
        X3 = 0
        C = math.cosh(wn*T)
        S = math.sinh(wn*T)
        x0 = A - d1
        xf = A + d2
        x_v0 = (xf - x0 * C - X1 - X2 - X3)*wn / S

    t = np.arange(dt, T+dt, dt)
    ## debug here
    if k > 0:
        P1S_four = StepZMP2CoM(t, wn_T, L1, 0) #$P1S,X1S,V1S,A1S
        P1R_four = RampZMP2CoM(t, wn_T, K1, 0) #P1R,X1R,V1R,A1R
        P2S_four = StepZMP2CoM(t, wn_T, L2 - L1, t1) #P2S,X2S,V2S,A2S
        P2R_four = RampZMP2CoM(t, wn_T, K2 - K1, t1) #P2R,X2R,V2R,A2R
        P3S_four = StepZMP2CoM(t, wn_T, L3 - L2, t2) #P3S,X3S,V3S,A3S
        P3R_four = RampZMP2CoM(t, wn_T, K3 - K2, t2) #P3R,X3R,V3R,A3R

        for i in range(len(t)):
            # [P1S,X1S,V1S,A1S] 
            [_, X1S, V1S, A1S] = P1S_four[i]
            # [P1R,X1R,V1R,A1R]
            [_, X1R, V1R, A1R] = P1R_four[i]
            # [P2S,X2S,V2S,A2S]
            [_, X2S, V2S, A2S] = P2S_four[i]
            # [P2R,X2R,V2R,A2R]
            [_, X2R, V2R, A2R] = P2R_four[i]
            # [P3S,X3S,V3S,A3S]
            [_, X3S, V3S, A3S] = P3S_four[i]
            # [P3R,X3R,V3R,A3R]
            [_, X3R, V3R, A3R] = P3R_four[i]	
            vec_three_num = []
            x_temp = x0 * math.cosh(wn_T[i] * t[i]) + x_v0 / wn_T[i] * math.sinh(wn_T[i] * t[i]) + X1S + X1R + X2S + X2R + X3S + X3R
            v_temp = x0 * wn_T[i] * math.sinh(wn_T[i] * t[i]) + x_v0 * math.cosh(wn_T[i] * t[i]) + V1S + V1R + V2S + V2R + V3S + V3R
            a_temp = x0 * (pow(wn_T[i], 2))*math.cosh(wn_T[i] * t[i]) + x_v0 * wn_T[i] * math.sinh(wn_T[i] * t[i]) + A1S + A1R + A2S + A2R + A3S + A3R
            p_temp = x_temp - a_temp / pow(wn_T[i], 2)
            x_temp = x_temp + Bias_x
            p_temp = p_temp + Bias_x

            vec_three_num.append(x_temp)
            vec_three_num.append(v_temp)
            vec_three_num.append(p_temp)
            vec_sum.append(vec_three_num)

    else:
        P1S_four = StepZMP2CoM(t, wn_T, L1, 0) #P1S,X1S,V1S,A1S 
        P1R_four = RampZMP2CoM(t, wn_T, K1, 0) #P1R,X1R,V1R,A1R			
        for i in range(len(t)):
            # [P1S,X1S,V1S,A1S] 
            _, X1S, V1S, A1S = P1S_four
            # [P1R,X1R,V1R,A1R]
            _, X1R, V1R, A1R = P1R_four

            vec_three_num = []
            x_temp = x0 * math.cosh(wn_T[i] * t[i]) + x_v0 / wn_T[i] * math.sinh(wn_T[i] * t[i]) + X1S + X1R
            v_temp = x0 * wn_T[i] * math.sinh(wn_T[i] * t[i]) + x_v0 * math.cosh(wn_T[i] * t[i]) + V1S + V1R
            a_temp = x0 * (pow(wn_T[i], 2))*math.cosh(wn_T[i] * t[i]) + x_v0 * wn_T[i] * math.sinh(wn_T[i] * t[i]) + A1S + A1R
            p_temp = x_temp - a_temp / pow(wn_T[i], 2)
            x_temp = x_temp + Bias_x
            p_temp = p_temp + Bias_x

            vec_three_num.append(x_temp)
            vec_three_num.append(v_temp)
            vec_three_num.append(p_temp)
            vec_sum.append(vec_three_num)

    return vec_sum

def y_OSCG_ZMPBC2_RampDSP_SinSSP(dt, wn_T, B, b, y0, yf, k, T, Bias_y, BoolPlot):
    # ok
    vec_sum = []
    g = 9.8
    wn_T_sum = sum(wn_T)
    wn = wn_T_sum/len(wn_T)
    vec_T_temp = []
    vec_T_temp.append(T)
    vec_wn_temp = []
    vec_wn_temp.append(wn)

    if k>0:
        t1 = k * T
        t2 = (1 - k)*T
        L1 = y0
        K1 = (B - y0) / k / T
        L2 = B
        M2 = b
        alpha = math.pi / (1 - 2 * k) / T
        T_d = alpha * t1
        L3 = B - (yf - B)*t2 / k / T
        K3 = (yf - B) / k / T

        # [Q1S,Y1S,V1S,A1S]
        Q1S_four = StepZMP2CoM(vec_T_temp, vec_wn_temp, L1, 0)
        # [Q1R,Y1R,V1R,A1R]
        Q1R_four = RampZMP2CoM(vec_T_temp, vec_wn_temp, K1, 0)
        # [Q2S,Y2S,V2S,A2S]
        Q2S_four = StepZMP2CoM(vec_T_temp, vec_wn_temp, L2 - L1, t1)
        # [Q2R,Y2R,V2R,A2R]
        Q2R_four = RampZMP2CoM(vec_T_temp, vec_wn_temp, -K1, t1)
        # [Q2Sin,Y2Sin,V2Sin,A2Sin]
        Q2Sin_four = SinTdZMP2CoM(vec_T_temp, vec_wn_temp, M2, alpha, T_d, t1)
        # [Q3S,Y3S,V3S,A3S]
        Q3S_four = StepZMP2CoM(vec_T_temp, vec_wn_temp, L3 - L2, t2)
        # [Q3R,Y3R,V3R,A3R]
        Q3R_four = RampZMP2CoM(vec_T_temp, vec_wn_temp, K3, t2)
        # [Q3Sin,Y3Sin,V3Sin,A3Sin]
        Q3Sin_four = SinTdZMP2CoM(vec_T_temp, vec_wn_temp, -M2, alpha, T_d, t2)

        Y1S = Q1S_four[0][1]
        Y1R = Q1R_four[0][1]

        Y2S = Q2S_four[0][1]
        Y2R = Q2R_four[0][1]
        Y2Sin = Q2Sin_four[0][1]

        Y3S = Q3S_four[0][1]
        Y3R = Q3R_four[0][1]
        Y3Sin = Q3Sin_four[0][1]

        Y1 = Y1S + Y1R
        Y2 = Y2S + Y2R + Y2Sin
        Y3 = Y3S + Y3R + Y3Sin

        C = math.cosh(wn*T)
        S = math.sinh(wn*T)

        y_v0 = (yf - y0 * C - Y1 - Y2 - Y3)*wn / S
    else:
        b = 0
        K1 = 2 * b / T
        L1 = (B - b)
        C = math.cosh(wn*T)
        S = math.sinh(wn*T)
        y_v0 = (yf - y0 * C - K1 * (T - S / wn) - L1 * (1 - C))*wn / S
    t = np.arange(dt, T+dt, dt)	
    if k>0:
        ## [Q1S,Y1S,V1S,A1S]
        Q1S_four = StepZMP2CoM(t, wn_T, L1, 0)
        ## [Q1R,Y1R,V1R,A1R]
        Q1R_four = RampZMP2CoM(t, wn_T, K1, 0)
        ## [Q2S,Y2S,V2S,A2S]
        Q2S_four = StepZMP2CoM(t, wn_T, L2 - L1, t1)
        ## [Q2R,Y2R,V2R,A2R]
        Q2R_four = RampZMP2CoM(t, wn_T, -K1, t1)
        ##----------------------------------------------------
        ## [Q2Sin,Y2Sin,V2Sin,A2Sin]
        Q2Sin_four = SinTdZMP2CoM(t, wn_T, M2, alpha, T_d, t1)
        ## [Q3S,Y3S,V3S,A3S]
        Q3S_four = StepZMP2CoM(t, wn_T, L3 - L2, t2)
        ## [Q3R,Y3R,V3R,A3R]
        Q3R_four = RampZMP2CoM(t, wn_T, K3, t2)
        ## [Q3Sin,Y3Sin,V3Sin,A3Sin]
        Q3Sin_four = SinTdZMP2CoM(t, wn_T, -M2, alpha, T_d, t2)		

        for i in range(len(t)):
            # [Q1S,Y1S,V1S,A1S]
            [_, Y1S, V1S, A1S] = Q1S_four[i]
            # [Q1R,Y1R,V1R,A1R]
            [_, Y1R, V1R, A1R] = Q1R_four[i]
            # [Q2S,Y2S,V2S,A2S]
            [_, Y2S, V2S, A2S] = Q2S_four[i]
            # [Q2R,Y2R,V2R,A2R]
            [_, Y2R, V2R, A2R] = Q2R_four[i]
            #-----------------------------------
            # [Q2Sin,Y2Sin,V2Sin,A2Sin]
            [_, Y2Sin, V2Sin, A2Sin] = Q2Sin_four[i]
            # [Q3S,Y3S,V3S,A3S]
            [_, Y3S, V3S, A3S] = Q3S_four[i]
            # [Q3R,Y3R,V3R,A3R]
            [_, Y3R, V3R, A3R] = Q3R_four[i]
            # [Q3Sin,Y3Sin,V3Sin,A3Sin]
            [_, Y3Sin, V3Sin, A3Sin] = Q3Sin_four[i]

            vec_three_num = []
            y_temp = y0 * math.cosh(wn_T[i] * t[i]) + y_v0 / wn_T[i] * math.sinh(wn_T[i] * t[i])\
                + Y1S + Y1R\
                + Y2S + Y2R + Y2Sin\
                + Y3S + Y3R + Y3Sin	
            v_temp = y0 * wn_T[i] * math.sinh(wn_T[i] * t[i]) + y_v0 * math.cosh(wn_T[i] * t[i])\
                + V1S + V1R\
                + V2S + V2R + V2Sin\
                + V3S + V3R + V3Sin		
            a_temp = y0 * (pow(wn_T[i], 2))*math.cosh(wn_T[i] * t[i]) + y_v0 * wn_T[i] * math.sinh(wn_T[i] * t[i])\
                + A1S + A1R\
                + A2S + A2R + A2Sin\
                + A3S + A3R + A3Sin

            q_temp = y_temp - a_temp / pow(wn_T[i], 2)
            y_temp = y_temp + Bias_y
            q_temp = q_temp + Bias_y

            vec_three_num.append(y_temp)
            vec_three_num.append(v_temp)
            vec_three_num.append(q_temp)
            vec_sum.append(vec_three_num)
    else:
        for i in range(0, len(t)):
            vec_three_num = []
            y_temp = y0 * math.cosh(wn_T[i] * t[i]) + y_v0 / wn_T[i] * math.sinh(wn_T[i] * t[i]) + K1 * (t[i] - 1 / wn_T[i] * math.sinh(wn_T[i] * t[i])) + L1 * (1 - math.cosh(wn_T[i] * t[i]))
            v_temp = y0 * wn_T[i] * math.sinh(wn_T[i] * t[i]) + y_v0 * math.cosh(wn_T[i] * t[i]) + K1 * (1 - math.cosh(wn_T[i] * t[i])) - L1 * wn_T[i] * math.sinh(wn_T[i] * t[i])
            a_temp = y0 * (pow(wn_T[i], 2))*math.cosh(wn_T[i] * t[i]) + y_v0 * wn_T[i] * math.sinh(wn_T[i] * t[i]) - K1 * wn_T[i] * math.sinh(wn_T[i] * t[i]) - L1 * (pow(wn_T[i], 2))*math.cosh(wn_T[i] * t[i])

            q_temp = y_temp - a_temp / pow(wn_T[i], 2)
            y_temp = y_temp + Bias_y
            q_temp = q_temp + Bias_y

            vec_three_num.append(y_temp)
            vec_three_num.append(v_temp)
            vec_three_num.append(q_temp)
            vec_sum.append(vec_three_num)
    return vec_sum


def StepZMP2CoM(t, wn, A_s, T_s):
# ok
    vec_sum = []
    for i in range(0, len(t)):
        vec_four_num = []
        p_temp = A_s * Heaviside(t[i] - T_s)
        x_temp = A_s * (1 - math.cosh(wn[i] * (t[i] - T_s))) * Heaviside(t[i] - T_s)
        v_temp = -A_s * wn[i] * math.sinh(wn[i] *(t[i] - T_s)) * Heaviside(t[i] - T_s)
        a_temp = -A_s * ((wn[i]**2))* math.cosh(wn[i] * (t[i] - T_s))* Heaviside(t[i] - T_s)
        vec_four_num.append(p_temp)
        vec_four_num.append(x_temp)
        vec_four_num.append(v_temp)
        vec_four_num.append(a_temp)
        vec_sum.append(vec_four_num)
    return vec_sum

def RampZMP2CoM(t, wn, A_r, T_r):
# ok
    vec_sum = []
    for i in range(len(t)):
        vec_four_num = []
        p_temp = A_r * t[i] * Heaviside(t[i] - T_r)
        x_temp = A_r * (t[i] - T_r * math.cosh(wn[i] * (t[i] - T_r)) - 1 / wn[i] * math.sinh(wn[i] * (t[i] - T_r)))*Heaviside(t[i] - T_r)
        v_temp = A_r * (1 - T_r * wn[i] * math.sinh(wn[i] * (t[i] - T_r)) - math.cosh(wn[i] * (t[i] - T_r)))*Heaviside(t[i] - T_r)
        a_temp = A_r * (-T_r * (pow(wn[i], 2))*math.cosh(wn[i] * (t[i] - T_r)) - wn[i] * math.sinh(wn[i] * (t[i] - T_r)))*Heaviside(t[i] - T_r)
        
        vec_four_num.append(p_temp)
        vec_four_num.append(x_temp)
        vec_four_num.append(v_temp)
        vec_four_num.append(a_temp)

        vec_sum.append(vec_four_num)
    return vec_sum

def SinTdZMP2CoM(t, wn, A_beta, alpha, T_d, T_beta):
    # ok
    # print(f"T_beta = {T_beta}")
    vec_sum = []
    # [ps, xs, vs, as]
    Ps_four = SinZMP2CoM(t, wn, A_beta, alpha, T_beta)
    # [pc, xc, vc, ac] 
    Pc_four = CosZMP2CoM(t, wn, A_beta, alpha, T_beta)

    for i in range(0, len(t)):
        vec_four_num = []

        # [ps, xs, vs, as]
        p_s = Ps_four[i][0]
        x_s = Ps_four[i][1]
        v_s = Ps_four[i][2]
        a_s = Ps_four[i][3]
        # [pc, xc, vc, ac]
        p_c = Pc_four[i][0]
        x_c = Pc_four[i][1]
        v_c = Pc_four[i][2]
        a_c = Pc_four[i][3]

        p_temp = p_s * math.cos(T_d) - p_c * math.sin(T_d)
        x_temp = x_s * math.cos(T_d) - x_c * math.sin(T_d)
        v_temp = v_s * math.cos(T_d) - v_c * math.sin(T_d)
        a_temp = a_s * math.cos(T_d) - a_c * math.sin(T_d)

        vec_four_num.append(p_temp)
        vec_four_num.append(x_temp)
        vec_four_num.append(v_temp)
        vec_four_num.append(a_temp)

        vec_sum.append(vec_four_num)
    return vec_sum

def SinZMP2CoM(t, wn, A_beta, alpha, T_beta):
    # ok
    vec_sum = []
    phi = A_beta * math.sin(alpha*T_beta)
    psi = alpha * A_beta*math.cos(alpha*T_beta)
    for i in range(len(t)):
        vec_four_num = []
        p_temp = A_beta * math.sin(alpha*t[i])*Heaviside(t[i] - T_beta)

        # x_temp =  1 / (pow(alpha,2)+pow(wn[i],2)) * (phi * pow(wn[i],2) * cos( alpha * (t[i]-T_beta) ) + psi / alpha * pow(wn[i],2) * sin( alpha * (t[i]-T_beta) ) - phi * ( pow(wn[i],2) ) * cosh( wn[i] * (t[i]-T_beta) ) - psi * wn[i] * sinh( wn[i] * (t[i]-T_beta) )) * Heaviside( t[i]-T_beta );
        x_temp = 1 / (pow(alpha, 2) + pow(wn[i], 2)) * (phi * pow(wn[i], 2) * math.cos(alpha * (t[i] - T_beta))\
            + psi / alpha * pow(wn[i], 2) * math.sin(alpha * (t[i] - T_beta))\
            - phi * (pow(wn[i], 2)) * math.cosh(wn[i] * (t[i] - T_beta))\
            - psi * wn[i] * math.sinh(wn[i] * (t[i] - T_beta))) * Heaviside(t[i] - T_beta)

        v_temp = 1 / (pow(alpha, 2) + pow(wn[i], 2))*(-alpha * phi*(pow(wn[i], 2))*math.sin(alpha*(t[i] - T_beta)) +\
            psi * (pow(wn[i], 2))*math.cos(alpha*(t[i] - T_beta)) -\
            phi * (pow(wn[i], 3))*math.sinh(wn[i] * (t[i] - T_beta)) -\
            psi * (pow(wn[i], 2))*math.cosh(wn[i] * (t[i] - T_beta)))*Heaviside(t[i] - T_beta)

        a_temp = 1 / (pow(alpha, 2) + pow(wn[i], 2))*(-phi * (pow(alpha, 2))*(pow(wn[i], 2))*math.cos(alpha*(t[i] - T_beta)) -\
            psi * alpha*(pow(wn[i], 2))*math.sin(alpha*(t[i] - T_beta)) -\
            phi * (pow(wn[i], 4))*math.cosh(wn[i] * (t[i] - T_beta)) -\
            psi * (pow(wn[i], 3))*math.sinh(wn[i] * (t[i] - T_beta)))*Heaviside(t[i] - T_beta)
        vec_four_num.append(p_temp)
        vec_four_num.append(x_temp)
        vec_four_num.append(v_temp)
        vec_four_num.append(a_temp)
        vec_sum.append(vec_four_num)
    return vec_sum

def CosZMP2CoM(t, wn, A_beta, alpha, T_beta):
    # ok
    vec_sum = []
    phi = A_beta * math.cos(alpha*T_beta)
    psi = -alpha * A_beta*math.sin(alpha*T_beta)	
    for i in range(len(t)):
        vec_four_num = []
        p_temp = A_beta * math.cos(alpha*t[i])*Heaviside(t[i] - T_beta)
        x_temp = 1 / (pow(alpha, 2) + pow(wn[i], 2)) * (phi * pow(wn[i], 2) * math.cos(alpha * (t[i] - T_beta))\
            + psi / alpha * pow(wn[i], 2) * math.sin(alpha * (t[i] - T_beta))\
            - phi * (pow(wn[i], 2)) * math.cosh(wn[i] * (t[i] - T_beta))\
            - psi * wn[i] * math.sinh(wn[i] * (t[i] - T_beta))) * Heaviside(t[i] - T_beta)
        v_temp = 1 / (pow(alpha, 2) + pow(wn[i], 2))*(-alpha * phi*(pow(wn[i], 2))*math.sin(alpha*(t[i] - T_beta)) +\
            psi * (pow(wn[i], 2))*math.cos(alpha*(t[i] - T_beta)) -\
            phi * (pow(wn[i], 3))*math.sinh(wn[i] * (t[i] - T_beta)) -\
            psi * (pow(wn[i], 2))*math.cosh(wn[i] * (t[i] - T_beta)))*Heaviside(t[i] - T_beta)
        a_temp = 1 / (pow(alpha, 2) + pow(wn[i], 2))*(-phi * (pow(alpha, 2))*(pow(wn[i], 2))*math.cos(alpha*(t[i] - T_beta)) -\
            psi * alpha*(pow(wn[i], 2))*math.sin(alpha*(t[i] - T_beta)) -\
            phi * (pow(wn[i], 4))*math.cosh(wn[i] * (t[i] - T_beta)) -\
            psi * (pow(wn[i], 3))*math.sinh(wn[i] * (t[i] - T_beta)))*Heaviside(t[i] - T_beta)
        
        vec_four_num.append(p_temp)
        vec_four_num.append(x_temp)
        vec_four_num.append(v_temp)
        vec_four_num.append(a_temp)
        vec_sum.append(vec_four_num)
    return vec_sum

def Heaviside(x):
    # ok
    if (x > 0.0000001):
        heavisde_x = 1
    elif ((x <= 0.0000001) and (x >= -0.0000001)):
        heavisde_x = 0.5
    else:
        heavisde_x = 0
    return heavisde_x	