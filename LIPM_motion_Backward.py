#coding=UTF-8

import math
import os
import csv
import numpy as np
from datetime import datetime
from pathlib import Path
import time
from source.command_generator import Cmd_Generator

write_data = True
motion_file = "./source/motordata/B_data.csv"    # LIPM給模擬環境的動作檔 
w_file = "./Motion_cmd/Motion_Backward.csv" # motor的刻度檔案


class LIPM_motion:
    def __init__(self):
        """ LIPM參數設定, 1右 2左 """
        # 大腿擺幅, |0.05| < B < |0.1|
        self.B1 = -105 * 0.001 # (+)往內畫 (-)往外畫 -0.04 
        self.B2 = 90 * 0.001  # (+)往外畫 (-)往內畫

        # Hip側開
        self.h1 = 0.2   # 右腳hip (+)抬左腳更低, 重心偏右 (-)抬左腳更高  -0.8
        self.h2 = 0.18    # 左腳hip (+)抬右腳更低, 重心偏左 (-)抬右腳更高   0.9

        # 膝蓋彎曲程度, 越大抬越高 (H > 0.05 調幅比較明顯) 
        self.H1 = 0.056   # 0.06   
        self.H2 = 0.054 

        # 腳底板傾斜角度, (+)腳尖往下 (-)腳尖往上
        self.CR1 = 0.75 #-0.75 -1.5
        self.CR2 = 0.5
        self.CL1 = -0.5  # 1.5
        self.CL2 = 0
        
        # Lean angle sequence [[step1], [step2]] 
                            #    1      2       3       4       5                   
        # self.A_Lean_angle_R = [[ 7.0,   12.0,    14.0,    12.0,     11.0], \
        #                        [ 12,    9.0,    6,    4,     4] ] # (-)腳尖往下 (+)腳尖往上
        self.A_Lean_angle_R = [[ 3,   8.5,    12,    9.5,     8.5], \
                               [ 6.3,    5.3,    4,    3.7,     3.3] ] # (-)腳尖往下 (+)腳尖往上
        self.A_Lean_angle_L = [[  0,    0.0,   -1.5,    -2.5,     -3], \
                               [  -4,    -3,  -2,   0.0,      -0.1]] # (+)腳尖往下 (-)腳尖往上

        # 向前的跨步距離 (+)forward (-)backward
        self.S1 = 0.08   # 0.1
        self.S2 = 0.0 

        # 側移的跨步距離 (側移時用)
        self.LyR = [0, 0]
        self.LyL = [0, 0]

        # 起步的前後偏移量，(+)腳尖往下 (-)腳尖往上 (單位:degree)
        self.L = 0 
        self.LeanAngleInit_R = self.L 
        self.LeanAngleInit_L = self.L

        # 調整Hip yaw角度, (+)往左 (-)往右
        # self.turn_angle = 0
        self.Turn = False
        self.A_DesiredTheta_R = [[  0,   0,   0,   0,   0],\
                                 [  0,   0,   0,   0,   0]]
        self.A_DesiredTheta_L = [[  0,   0,   0,   0,   0],\
                                 [  0,   0,   0,   0,   0]]

        """ 機器人質心高度設定 """
        self.foot_height = 0.59        # foot_height: 行走時髖關節的高度，影響機器人行走時蹲的幅度 
        self.zc = 0.55 #0.45           # zc: 質心高度，影響LIPM產生的質心軌跡 
        self.CoMx_bias = 0.            # CoMx_bias: 質心前後偏移量

        """ 動作時間設定 """
        self.footstep = 3
        self.k_DSP = 0.4419             # 值越小，抬腳的時間越長
        self.T = 2                      # T: 一步的週期時間(sec)

        self.stepH = [self.H1, self.H2] # stepH: 跨步的高度
        self.S = [self.S1, self.S2]     # S: 前進、後退的距離
        self.B = [self.B1, self.B2]
        self.b = [-0.005, 0.005]
        self.hip = [self.h1, self.h2]
        self.DesiredTheta_R = [[0, 0, 0, 0, 0], [0, 0, 0, 0, 0]]
        self.DesiredTheta_L = [[0, 0, 0, 0, 0], [0, 0, 0, 0, 0]]

        self.Rup = True                 # Rup: 是否先抬右腳
        self.Forward = False # True
        self.Shift = False
        self.NoPeriod_acc = 0
        self.NoPeriod_dec = 0

        """ 機器人基礎姿態設定 """
        # Set Initial Pose
        self.initR = [0,  0.0,   0.4, -0.8,  0.4,  0.0] 
        self.initL = [0,  0.0,   0.4, -0.8,  0.4,  0.0]
        # Legs : Legs link lenghth (mm)
        self.Legs = np.array([102.5, 170, 175, 70.33, 34.7])/1000

    def GenerateMotionData(self):
        LR = []
        LL = []
        L = []
        wn_T = []
        wn_T_D = []
        wn_T_S = []
        # A_DesiredTheta_R = np.ones([2,5])
        # A_DesiredTheta_L = np.ones([2,5])
        A_DesiredTheta_R = [[1, 1, 1, 1, 1], [1, 1, 1, 1, 1]]
        A_DesiredTheta_L = [[1, 1, 1, 1, 1], [1, 1, 1, 1, 1]]
        vector_S = self.S[:self.footstep]
        N = len(vector_S)
        N_Trans = 0

        T = self.T
        k_DSP = self.k_DSP
        zc = self.zc
        #--------------------------------------
        LIPMtype = 'S'
        k = 0.1          									
        dt = 0.01
        t_T = np.arange(dt, T+dt, dt)
        g = 9.81
        sampleT = 0.05
        Deg2Rad = math.pi / 180
        T_DSP = k_DSP * T
        #--------------------------------------
        CoMx = []
        CoMy = []
        vx = []
        vy = []
        ZMPx = []
        ZMPy = []
        k_lean_col = []
        k_desired_col = []
        k_Lean_angle = []
        k_DesiredTheta = []
        hip_r = self.hip[0]
        hip_l = self.hip[1]
        foot_height_R = self.foot_height
        foot_height_L = self.foot_height
        ########################################
        m = 2703.39 / 1000 											# Unit Kg
        M = 4500 / 1000 											# Unit Kg
        lm = 180/1000#166.23 / 1000 								# CoM height of Upperbody:  m, CoM height of roboy: 
        ########################################
        wn = math.sqrt(g/zc)
        HR = [self.stepH[0], 0]
        HL = [0, self.stepH[1]]
        

        if (LIPMtype == 'D'):
            wn_SSP = math.sqrt((M*g + 0.5*m*g) / (M*zc + 0.25*m*zc + M * lm))
            wn_DSP = math.sqrt(g / zc)
            # wn_DSP
            for i in range(len(t_T)):
                temp_wn_T_D = wn_DSP * ((t_T[i] <= T_DSP / 2) | (t_T[i] >  T - T_DSP / 2)) + wn_SSP * ((t_T[i] >  T_DSP / 2) & (t_T[i] <= T - T_DSP / 2))
                wn_T_D = np.append(wn_T_D, temp_wn_T_D)	
            wn_T = wn_T_D
        else:
            for i in range(len(t_T)):
                wn_T_S = np.append(wn_T_S, wn)
            wn_T = wn_T_S

        #########################################################################
        ##                           Gait Generation                           ##
        #########################################################################
        from source.motion_function import StepSize2StrideLength, StrideLength2ZMPAmplitude, Vec_Arr_Sam, \
            Completed_R_generation, modifiable_foot_generation, OutputMotion
        from source.LIPM_function import modifiable_x_OSG_RampDSP_RampSSP, modifiable_y_OSG_RampDSP_SinSSP

        [LR, LL, L] = StepSize2StrideLength(vector_S)
        A = StrideLength2ZMPAmplitude(L)

        a = Vec_Arr_Sam(len(A), 0.005)	

        for i in range(len(A)):
            a[i] = (self.S[i]!=0) * a[i]

        if (self.Forward == False):
            A = -A
            a = -a

        CoMx_vx_ZMPx = modifiable_x_OSG_RampDSP_RampSSP(dt, wn_T, A, a, T, k, 1)
        CoMy_vy_ZMPy = modifiable_y_OSG_RampDSP_SinSSP(dt, wn_T, self.B, self.b, T, k, (~self.Shift), True)	

        for i in range(len(CoMx_vx_ZMPx)):
            if i == 0:
                CoMx.append(CoMx_vx_ZMPx[i][0])
            else:
                CoMx.append(CoMx_vx_ZMPx[i][0] + self.CoMx_bias)
            CoMy.append(CoMy_vy_ZMPy[i][0])
            # vx.append(CoMx_vx_ZMPx[i][1])
            # vy.append(CoMy_vy_ZMPy[i][1])
            # ZMPx.append(CoMx_vx_ZMPx[i][2])
            # ZMPy.append(CoMy_vy_ZMPy[i][2])

        k_lean_col = [k_DSP/2, 0.5, 1-(k_DSP/2)]
        k_desired_col = [k_DSP/2, 0.5, 1-(k_DSP/2)]

        for i in range(N):
            k_Lean_angle.append(k_lean_col)
            k_DesiredTheta.append(k_desired_col)
        # did not check k_DesiredTheta
        Lean_angleR = Completed_R_generation(dt, self.LeanAngleInit_R, T, k_Lean_angle, self.A_Lean_angle_R, True)
        Lean_angleL = Completed_R_generation(dt, self.LeanAngleInit_L, T, k_Lean_angle, self.A_Lean_angle_L, True)

        Lean_angleR = np.array(Lean_angleR) * Deg2Rad
        Lean_angleL = np.array(Lean_angleL) * Deg2Rad

        # ## turn angle ##	
        DesiredTheta_R1 = []
        DesiredTheta_L1 = []        
        
        if self.Turn == True:
            DesiredTheta_R1 = Completed_R_generation(dt, 0, T, k_DesiredTheta, self.A_DesiredTheta_R, True)
            DesiredTheta_L1 = Completed_R_generation(dt, 0, T, k_DesiredTheta, self.A_DesiredTheta_L, True)
            
            DesiredTheta_R1 = np.array(DesiredTheta_R1) * Deg2Rad
            DesiredTheta_L1 = np.array(DesiredTheta_L1) * Deg2Rad

        if self.Forward == False:
            LR = -np.array(LR)
            LL = -np.array(LL)

        LHR = []
        LHL = []
        LR_LyR_HR = []
        LL_LyL_HL = []
        LR_LyR_HR.extend([LR, self.LyR, HR])
        LL_LyL_HL.extend([LL, self.LyL, HL])
        
        if (self.Rup):
            LHR = LR_LyR_HR
            LHL = LL_LyL_HL
        else:
            LHL = LR_LyR_HR
            LHR = LL_LyL_HL

        FootR = modifiable_foot_generation(dt, LHR, T, T_DSP, False) 
        FootL = modifiable_foot_generation(dt, LHL, T, T_DSP, False)

        t = np.arange(dt, N*(T)+dt, dt)	

        Zc_sinusoid = []
                
        for i in range(len(t)):
            temp_Zc = 1.2 / 1000 * -math.cos(2 * math.pi*t[i] / T)
            Zc_sinusoid.append(temp_Zc)
        PRx = []
        PRy = []
        PRz = []
        PLx = []
        PLy = []
        PLz = []

        for i in range(len(FootR)):
            PRx.append(FootR[i][0] - CoMx[i])
            PRy.append(FootR[i][1] - CoMy[i])
            PRz.append(FootR[i][2] - (Zc_sinusoid[i] + foot_height_R))

            PLx.append(FootL[i][0] - CoMx[i])
            PLy.append(FootL[i][1] - CoMy[i])
            PLz.append(FootL[i][2] - (Zc_sinusoid[i] + foot_height_L))

        Lean_angleR_ini = 0
        Lean_angleL_ini = 0	
        #####################################
        arr_PR_ini = [0 , 0, -0.5]
        arr_PL_ini = [0 , 0, -0.5]	
        #####################################
        PR_ini = arr_PR_ini
        PL_ini = arr_PL_ini

        N_Total = N
        index_acc = 0
        index_dec = 0	
        if ((self.NoPeriod_acc > 0) and (self.NoPeriod_dec > 0) and (self.NoPeriod_acc < self.NoPeriod_dec)):
            index_acc = round((self.NoPeriod_acc + N_Trans)*T / sampleT )
            index_dec = round((self.NoPeriod_dec + N_Trans)*T / sampleT) 
                        
        # #########################################################################
        # ##                           Output Motion                             ##
        # #########################################################################
        
        thetaR_length = math.floor(N_Total*T / sampleT)	

        # plt.plot(PRz,label='PRz')
        # plt.plot(PLz,label='PLz')
        # plt.ylabel('PRz (m)')
        # plt.show()

        motor_job_data, Motor_test_data = OutputMotion(self.initR, self.initL, self.Rup, T, sampleT, dt, self.hip, self.Legs, thetaR_length,\
        PRx, PRy, PRz, PLx, PLy, PLz, Lean_angleR, Lean_angleL, DesiredTheta_R1, DesiredTheta_L1, self.Turn, index_acc, index_dec, foot_height_R)
        # print("motorjob\n", motor_job_data)

        return motor_job_data, Motor_test_data	


    
def write_file(Motor_test_data):

    folder = os.path.exists('../Motor_Driver/USB2CAN_CppController/motordata')
    if not folder:
        #如果不存在，則建立新目錄
        os.makedirs('../Motor_Driver/USB2CAN_CppController/motordata')
    # Store the new csv
    with open(motion_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerows(Motor_test_data)
        # for row in Motor_test_data:
        #     writer.writerow(row)


LIPM_obj = LIPM_motion()
jointPoses, Motor_data= LIPM_obj.GenerateMotionData() 

if write_data == True:
    write_file(Motor_data)

CMD_obj = Cmd_Generator(motion_file, w_file)
CMD_obj.write_cmd()
