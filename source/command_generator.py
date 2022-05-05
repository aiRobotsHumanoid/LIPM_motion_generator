import numpy as np
import pandas as pd
from math import *
import time
import csv



class Cmd_Generator:
    def __init__(self, motion_file, w_file):

        """ Presetting """

        self.r_file = motion_file
        self.w_file = w_file

        # 用哪些馬達
        self.add_Dynamixel_data = True # True False
        self.add_Haitai_data = [30, 32, 33, 34, 40, 42, 43, 44] 

        if len(self.add_Haitai_data) == 4:
            self.frame_range = 6
        elif len(self.add_Haitai_data) == 8:
            self.frame_range = 12
        # print(f"len of add_Haitai_data = {len(self.add_Haitai_data)}")

        # 每個timestep間隔
        self.time_step = 0.1


        with open(self.r_file, 'r', newline='') as r_f:
            self.datas = pd.read_csv(r_f, header=None)
            self.datas = self.datas.values.tolist()

        self.motor_list = [30, 31, 32, 33, 34, 35, 40, 41, 42, 43, 44, 45]
        self.pre_pos = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        self.mode_list = np.ones([len(self.datas),len(self.datas[0])]) # 0: pos of this step is same as pos of last step
        self.vel_list = np.ones([len(self.datas),len(self.datas[0])])

    def data_scan(self):
        """
            @brief: Scan the datas to record the Haitai mode and Dyn velocity

            @param mode_list - Haitai mode
            @param vel_list  - Dyn velocity
        """
        for frame_idx, frame in enumerate(self.datas):
            for list_idx in range(self.frame_range):
                motordata = frame[list_idx]
                
                ## Motorid ##
                motorid = self.motor_list[list_idx]

                ## Position command ##
                p_cmd = np.round(motordata, 2)
                
                ## Set HT03 mode 0 (situation:first,last, the same command) ##
                if frame_idx==0 or frame_idx==len(self.datas)-1:
                    self.mode_list[frame_idx][list_idx] = 0

                elif p_cmd==self.pre_pos[list_idx]:
                    self.mode_list[frame_idx-1][list_idx] = 0

                ## Calculate Dynamixel velocity ##
                if  self.add_Dynamixel_data == True:
                    if frame_idx == 0:
                        time_cmd = 0.8
                    else:
                        time_cmd = 0.1

                    if motorid in [31, 41]: # Pro100
                        scaled = 2920 / 175.2  # (max/ (deg/sec))
                        self.vel_list[frame_idx][list_idx] = np.round((p_cmd - self.pre_pos[list_idx]) / time_cmd * scaled)

                    elif motorid in [3, 35, 45]: # Pro20+
                        scaled = 2920 / 175.2  
                        self.vel_list[frame_idx][list_idx] = np.round((p_cmd - self.pre_pos[list_idx]) / time_cmd * scaled)

                    elif motorid in [4, 10, 11, 12, 13, 14, 15, 20, 21, 22, 23, 24, 25]: # Xm540
                        scaled = 1023 / 1405.6 
                        self.vel_list[frame_idx][list_idx] = np.round((p_cmd - self.pre_pos[list_idx]) / time_cmd * scaled)
                
                self.pre_pos[list_idx] = p_cmd

        return self.mode_list, self.vel_list

    def write_cmd(self):
        """ 
            Data Format : 
            HaiTai       : [motortype, motorid, mode, p_cmd, time_cmd, t_ff, kp, kd]
            Dynamixel    : [motortype, motorid, p_cmd, vel]
            Wait Command : [motortype, wait_time]

            @param command_list - All motor cmd
        """
        mode_list, vel_list = self.data_scan()

        command_list = []
        for frame_idx, frame in enumerate(self.datas):
            if frame_idx == 0:
                # command = [motortype, motorid, p_cmd, vel]
                command = [5, 3, -0.3, 8.0]
                if self.add_Dynamixel_data == True:
                    command_list.append(command)

            # frame_range = len(frame)-6 # right leg
            for list_idx in range(self.frame_range):
                motordata = frame[list_idx]

                ## Motorid ##
                motorid = self.motor_list[list_idx]

                ## motortype ##
                # Right Leg
                if motorid == 30 or motorid == 32 or motorid == 33 or motorid == 34:
                    motortype = 3
                # Left Leg
                elif motorid == 40 or motorid == 42 or motorid == 43 or motorid == 44:
                    motortype = 4
                # Dynamixel pro
                elif motorid == 3 or motorid == 31 or motorid == 35 or motorid == 41 or motorid == 45:
                    motortype = 5


                ## Velocity (just for Dynamixel) ##
                vel = vel_list[frame_idx][list_idx]
                
                ## Mode (just for haitai) ##
                mode = mode_list[frame_idx][list_idx]

                ## Position command (check +/-) ##
                if motorid == 31 or motorid == 35 or motorid == 41 or motorid == 45: 
                    p_cmd = np.round(-motordata, 2)
                else: p_cmd = np.round(motordata, 2)

                ## Time command (sec) ##
                if frame_idx == 0:
                    time_cmd = 0.8 
                else: time_cmd = self.time_step

                ## Torque command ##
                t_ff = 0

                ## Kp command, if the motion is jumping, we need to set a smaller kp ##
                if frame_idx >= 0 :
                    kp = 100
                
                ## Kd command ##
                kd = 0.2

                """ Write Command """
                if frame_idx==len(self.datas)-1:
                    if motorid == 3 or motorid == 31 or motorid == 35 or motorid == 41 or motorid == 45:
                        command = [motortype, motorid, p_cmd, vel]
                        if self.add_Dynamixel_data == True:
                            command_list.append(command)
                    elif motorid in self.add_Haitai_data:
                        command = [motortype, motorid, mode, p_cmd, time_cmd, t_ff, kp, kd]
                        command_list.append(command)

                elif p_cmd == self.pre_pos[list_idx] and frame_idx>0:
                    continue

                else:
                    if motorid == 31 or motorid == 35 or motorid == 41 or motorid == 45:
                        command = [motortype, motorid, p_cmd, vel]
                        if self.add_Dynamixel_data == True:
                            command_list.append(command)

                    elif motorid in self.add_Haitai_data:
                        command = [motortype, motorid, mode, p_cmd, time_cmd, t_ff, kp, kd]
                        command_list.append(command)
                self.pre_pos[list_idx] = p_cmd
            

            ## Wait Command = [motortype, wait_time] ##
            wait_command = [100, time_cmd*pow(10, 6)]
            command_list.append(wait_command)
        ending_command = [100, 3*pow(10, 6)]  # c++ usleep(1) = 10^-6 sec
        command_list.append(ending_command)


        with open(self.w_file, 'w', newline='') as w_f:
            writer = csv.writer(w_f)
            for row in command_list:
                writer.writerow(row)

        print('Done Command Generating!')
        
        
