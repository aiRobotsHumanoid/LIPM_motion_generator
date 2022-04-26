import os
import math
import time
import numpy as np
from datetime import datetime
from pathlib import Path
import matplotlib.pyplot as plt
from source.LIPM_function import Heaviside
from source.InverseKinematics import InvK

def mod(x, y):
    if x*y < 0:
        result = x % y + y
    else:
        result = x % y
    return result

def StepSize2StrideLength(S):
    vec_sum = np.array([])
    vec_L = np.array([])
    LR = np.array([])
    LL = np.array([])
    L = np.array([])
    N = len(S)

    vec_L = np.append(vec_L, S[0])
    for i in range(1, N):
        L_temp = S[i] + S[i-1]
        vec_L = np.append(vec_L, L_temp)
    
    for k in range(0, N):
        # vec_three = np.array([])
        LR = np.append(LR, vec_L[k] * mod(k + 1, 2))
        LL = np.append(LL, vec_L[k] * abs(mod(k + 1, 2) -1))
        L = np.append(L, vec_L[k])
    return [LR, LL, L]

def StrideLength2ZMPAmplitude(L):
    vec_A = np.array([])
    N = len(L)
    size = 1, N+2
    mat_A = np.zeros(size, np.float64)
    L_sum = sum(L)
    mat_A[0][N+1] = L_sum
    for i in range(N-1):
        mat_A[0][i+2] = L[i] + mat_A[0][i]
    for k in range(1, N+1):
        vec_A = np.append(vec_A, mat_A[0][k])
    return vec_A

def Vec_Arr_Sam(vec_size, inside_num):
    x = np.array([])
    for i in range(1, vec_size+1):
        x = np.append(x, inside_num)
    return x

def Completed_R_generation(dt, Bias_R, T, k, A, boolPlot):
    # ok
    R = []
    R_temp = []
    row_k = len(k)
    row_A = len(A)
    N = row_k
    col_k = len(k[0])
    col_A = len(A[0])

    if (row_k != row_A) or (col_k != 3) or (col_A !=5):
        print("input error: k is a Nx3 vector; A is a Nx5 vector")
        print("col_k", col_k)
        print("col_A", col_A)
    else:
        for n in range(N):
            R_temp = OneStep_R_generation(dt, Bias_R, T, k[n], A[n], False)
            for i in range(0, len(R_temp)):
                R.append(R_temp[i])
    return R

def OneStep_R_generation(dt, Bias_R, T, k, A, boolPlot):
    # ok
    R = []
    col_k = len(k)
    col_A = len(A)

    if (col_k != 3) or (col_A != 5):
        print("input error: k is a 1x3 vector; A is a 1x5 vector")
    else:
        t = np.arange(dt, T+dt, dt)
        t_k = []
        N = len(t)
        size = 1, N
        mat_R = np.zeros(size, np.float64)

        t_k.append(0*T)
        t_k.append(k[0] * T)
        t_k.append(k[1] * T)
        t_k.append(k[2] * T)
        t_k.append(1*T)

        for n in range(4):
            dA = A[n + 1] - A[n]
            dT = t_k[n + 1] - t_k[n]
            if dT != 0:
                for i in range(len(t)):
                    theta = (t[i] - t_k[n]) * 2 * math.pi / dT
                    K = dA / 2 / math.pi * (theta - math.sin(theta))
                    mat_R[0][i] = mat_R[0][i] + + (A[n] + K) * \
                        (Heaviside(t[i] - t_k[n]) - Heaviside(t[i] - t_k[n + 1]))
        mat_R[0][N-1] = A[len(A) - 1]
        for j in range(N):
            R.append(mat_R[0][j] + Bias_R)
    return R

def modifiable_foot_generation(dt, Amp, T, T_DSP, BoolPlot):
    # ok
    vec_sum = []
    a = len(Amp)
    N = len(Amp[0])

    if a == 3:
        t1 = T_DSP / 2
        t2 = T - T_DSP / 2
        T_SSP = T - T_DSP
        biasx = 0
        biasy = 0
        t_index = 1
        n_index = 1	
        t_real = 0
        size_t_temp = N * T / dt

        for k in range(1, int(size_t_temp)+1):
            t = t_index * dt
            Ampx = Amp[0][n_index - 1]
            Ampy = Amp[1][n_index - 1]
            Ampz = Amp[2][n_index - 1]
            temp_x =  (biasx + Ampx / 2 / math.pi * (2 * math.pi*(t - t1) / T_SSP - math.sin(2 * math.pi*(t - t1) / T_SSP)) * (Heaviside(t - t1 - dt) - Heaviside(t - t2))\
                + Ampx * Heaviside(t - t2))
            temp_y = biasy + Ampy / 2 / math.pi * (2 * math.pi*(t - t1) / T_SSP - math.sin(2 * math.pi*(t - t1) / T_SSP)) * (Heaviside(t - t1 - dt) - Heaviside(t - t2))\
                + Ampy * Heaviside(t - t2)
            temp_z = Ampz / 2 * (1 - math.cos(2 * math.pi*(t - t1) / T_SSP)) * (Heaviside(t - t1 - dt) - Heaviside(t - t2)) 
            t_index = t_index + 1
            t_real = t_real + dt

            vec_sum.append([temp_x, temp_y, temp_z])

            if t_real > n_index*T:
                biasx = biasx + Amp[0][n_index - 1]
                biasy = biasy + Amp[1][n_index - 1]
                t_index = 1	
                n_index = n_index + 1
                if n_index > N: 
                    break
    return vec_sum

def OutputMotion(initR, initL, Rup, T, sampleT, dt, hip, Legs, thetaR_length,\
    PRx, PRy, PRz, PLx, PLy, PLz, Lean_angleR, Lean_angleL, DesiredTheta_R1, DesiredTheta_L1, Turn, index_acc, index_dec, foot_height):
    motor_job_data = 0
    hip_r = hip[0]
    hip_l = hip[1]
    samplek = sampleT / dt
    d2 = 43.42 / 1000  # motor35/45 ＆ 地面的距離
    ####################################
    thetaR_leg = []
    thetaL_leg = []

    """ 
        @brief  Motion Planning : Calculate Inverse Kinetic for all motor

        @param  thetaR_leg - planned motion for right leg
        @param  thetaL_leg - planned motion for right leg
    """
    for i in range(1, thetaR_length+1):
        index_i = int(i * samplek)
        
        P_R = [PRx[index_i - 1], PRy[index_i - 1], PRz[index_i - 1]]
        R_R = [[ math.cos(Lean_angleR[index_i - 1]),	0,	 math.sin(Lean_angleR[index_i - 1]) ],
               [				0,			        	1,							     	  0 ],
               [ -math.sin(Lean_angleR[index_i - 1]),	0,	 math.cos(Lean_angleR[index_i - 1]) ] ]
        P_L = [PLx[index_i - 1], PLy[index_i - 1], PLz[index_i - 1] ]
        R_L = [[ math.cos(Lean_angleL[index_i - 1]),	0,	 math.sin(Lean_angleL[index_i - 1]) ],
               [				0,				        1,							       	 0  ],
               [ -math.sin(Lean_angleL[index_i - 1]),	0,	 math.cos(Lean_angleL[index_i - 1]) ]]
        thetaR_leg.append(InvK(d2, Legs, R_R, P_R))
        thetaL_leg.append(InvK(d2, Legs, R_L, P_L))

        if Turn == True:
            thetaR_leg[i-1][0] = DesiredTheta_R1[index_i-1]
            thetaL_leg[i-1][0] = DesiredTheta_L1[index_i-1]
        # thetaR_leg[i-1][0] = DesiredTheta_R1[index_i-1]
        # thetaL_leg[i-1][0] = DesiredTheta_L1[index_i-1]

        if index_i < 1:
            current_step = 0
        else:
            current_step = math.floor((index_i - 1)*dt / T)
        # print(f"hip_r={hip_r} hip_l={hip_l}")=
        if (mod(current_step, 2) == 0):
            if Rup:
                thetaL_leg[i - 1][1] = hip_l * thetaL_leg[i - 1][1]
            else:		
                thetaR_leg[i - 1][1] = hip_r * thetaR_leg[i - 1][1]
        else:
            if Rup:
                thetaR_leg[i - 1][1] = hip_r * thetaR_leg[i - 1][1]
            else:	
                thetaL_leg[i - 1][1] = hip_l * thetaL_leg[i - 1][1]													


    """
        @brief  Calculate motor scaled and offset in Sim and Real

        @param  ScaledthetaX_leg_test - motor data for simulation
        @param  Motor_ScaledthetaX_leg_test - motor data for real world
    """
    size = 6, thetaR_length
    ScaledthetaR_leg_test = np.zeros(size, np.float64)
    ScaledthetaL_leg_test = np.zeros(size, np.float64)
    Motor_ScaledthetaR_leg_test = np.zeros(size, np.float64)
    Motor_ScaledthetaL_leg_test = np.zeros(size, np.float64)

    # Joint angle: radian -> motor resolution
    # For Sim scaled (rad to rad)
    MScaleR_Leg = [1, 1, 1, 1, 1, 1]
    MScaleL_Leg = [1, 1, 1, 1, 1, 1]
    # For Real scaled (rad to HaiTai/Dyn)
    Motor_MScaleR_Leg = [1, 180/math.pi, 1, 1, 1, 180/math.pi]
    Motor_MScaleL_Leg = [1, 180/math.pi, 1, 1, 1, 180/math.pi]

    # motor angle -> motor position
    Motor_MotorOffsetR = []
    Motor_MotorOffsetL = []
    MotorOffsetR = []
    MotorOffsetL = []	

    # motor rotation +/-
    Motor_Dir_R = [1, 1,  1,-1, -1, 1]
    Motor_Dir_L = [1, 1, -1, 1, -1, 1]
    for j in range(6):
        """ Calculate motor initial offset  """
        """ @param  Offset - Like bias, how far are the planned motion(in time 0) and the Initial pose  """
        # For Sim
        MotorOffsetR.append(initR[j] - Motor_Dir_R[j] * thetaR_leg[0][j])
        MotorOffsetL.append(initL[j] - Motor_Dir_L[j] * thetaL_leg[0][j])
        # For Real
        Motor_MotorOffsetR.append(initR[j] - Motor_Dir_R[j] * thetaR_leg[0][j])
        Motor_MotorOffsetL.append(initL[j] - Motor_Dir_L[j] * thetaL_leg[0][j])
    
    for i in range(thetaR_length):
        """ Calculate the motor data with scaled and offset """
        for j in range(6):
            # For Sim   
            ScaledthetaR_leg_test[j][i] = MotorOffsetR[j] + thetaR_leg[i][j] * MScaleR_Leg[j] * Motor_Dir_R[j]
            ScaledthetaL_leg_test[j][i] = MotorOffsetL[j] + thetaL_leg[i][j] * MScaleL_Leg[j] * Motor_Dir_L[j]
            # For Real
            Motor_ScaledthetaR_leg_test[j][i] = (Motor_MotorOffsetR[j] + thetaR_leg[i][j] * Motor_Dir_R[j]) * Motor_MScaleR_Leg[j]
            Motor_ScaledthetaL_leg_test[j][i] = (Motor_MotorOffsetL[j] + thetaL_leg[i][j] * Motor_Dir_L[j]) * Motor_MScaleL_Leg[j]

            # if i == 0: # 應該要跟 init pose(rad/s) 相同 => 修改 Motor_Dir_R / Motor_Dir_L
            #     print(f"{Motor_ScaledthetaR_leg_test[j][i]/Motor_MScaleR_Leg[j]} = ({Motor_MotorOffsetR[j]} + {thetaR_leg[i][j]} * {Motor_Dir_R[j]})")
            #     print(f"{Motor_ScaledthetaL_leg_test[j][i]/Motor_MScaleL_Leg[j]} = ({Motor_MotorOffsetL[j]} + {thetaL_leg[i][j]} * {Motor_Dir_L[j]})")
            if j == 0:
                Motor_ScaledthetaR_leg_test[j][i] = thetaR_leg[i][j] * Motor_MScaleR_Leg[j] * Motor_Dir_R[j]
                Motor_ScaledthetaL_leg_test[j][i] = thetaL_leg[i][j] * Motor_MScaleL_Leg[j] * Motor_Dir_L[j] 
                # print(f"{Motor_ScaledthetaR_leg_test[j][i]} = {thetaR_leg[i][j]} * {Motor_MScaleR_Leg[j]} * {Motor_Dir_R[j]}")
    
    """ FORCE the robot back to initial pose in the end """
    for i in range(6):
        # For Sim
        ScaledthetaR_leg_test[i][thetaR_length-1] = ScaledthetaR_leg_test[i][0]
        ScaledthetaL_leg_test[i][thetaR_length-1] = ScaledthetaL_leg_test[i][0]
        # For Real
        Motor_ScaledthetaR_leg_test[i][thetaR_length-1] = Motor_ScaledthetaR_leg_test[i][0]
        Motor_ScaledthetaL_leg_test[i][thetaR_length-1] = Motor_ScaledthetaL_leg_test[i][0]

    """ For Sim, let the robot be stable after walking """
    # add more stable step 
    # for i in range(5):
    #     ScaledthetaR_leg_test = np.append(ScaledthetaR_leg_test, ScaledthetaR_leg_test[i][0])
    #     ScaledthetaL_leg_test = np.append(ScaledthetaL_leg_test, ScaledthetaL_leg_test[i][0])
        # Motor_ScaledthetaR_leg_test = np.append(Motor_ScaledthetaR_leg_test, Motor_ScaledthetaR_leg_test[i][0])
        # Motor_ScaledthetaL_leg_test = np.append(Motor_ScaledthetaL_leg_test, Motor_ScaledthetaL_leg_test[i][0])
    

    """ Collect motor data """
    # For Sim
    motor_job_data =  np.concatenate((ScaledthetaR_leg_test, ScaledthetaL_leg_test), axis=0)
    motor_job_data = np.round(motor_job_data.T, 3)
    for i in range(30):
        motor_job_data = np.append(motor_job_data, [motor_job_data[-1]], axis=0)

    # For Real
    Motor_test_data = np.concatenate((Motor_ScaledthetaR_leg_test, Motor_ScaledthetaL_leg_test), axis=0)
    Motor_test_data = np.round(Motor_test_data.T, 3)
    
    # print(f'motor_job_data = {len(motor_job_data)}')
    return motor_job_data , Motor_test_data

