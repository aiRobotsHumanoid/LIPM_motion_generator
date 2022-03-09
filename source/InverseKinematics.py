import math
import numpy as np
from datetime import datetime
from pathlib import Path
import time
import matplotlib.pyplot as plt
import os


def InvK(d2, L, R, P):
    """ 
    @ L : represent Legs lenghth
    @ q : q[0]-q[5] represent each joint angle
    """
    q = np.zeros(6)

    [Px, Py, Pz] = P

    [L2, L3, L4, L5, L6] = L

    [r11, r12, r13] = R[0]
    [r21, r22, r23] = R[1]
    [r31, r32, r33] = R[2]

    P05x = Px + L6 * r13
    P05y = Py + L6 * r23
    P05z = Pz + L6 * r33
    
    q[1] = np.arctan(-(Py + L6 * r23) / (L6*r33 + Pz))
    q[5] = np.arcsin(-r23 * np.cos(q[1]) - r33 * np.sin(q[1]))

    H04x = Px + r13 * (L6 + L5 * np.cos(q[5])) + L5 * r12*np.sin(q[5])
    H04y = Py + r23 * (L6 + L5 * np.cos(q[5])) + L5 * r22*np.sin(q[5])
    H04z = Pz + r33 * (L6 + L5 * np.cos(q[5])) + L5 * r32*np.sin(q[5])

    H02x = d2
    H02y = L2 * np.sin(q[1])
    H02z = -L2 * np.cos(q[1])

    H04 = [ H04x, H04y, H04z ]
    H02 = [ H02x, H02y, H02z ]
    L24 = np.sqrt((H04[0] - H02[0])**2 + (H04[1] - H02[1])**2 + (H04[2] - H02[2])**2)
    theta = ( L24**2 - L3**2 - L4**2) / (2 * L3*L4)
    if theta > 1:
        theta = 1
    elif theta < -1:
        theta = -1
    else:
        theta = theta
    q[3] = np.arccos(theta)

    alpha = Py * np.sin(q[1]) - Pz * np.cos(q[1]) - (L6 + L5 * np.cos(q[5])) * (r33*np.cos(q[1]) - r23 * np.sin(q[1])) - L2 - L5 * np.sin(q[5]) * (r32* np.cos(q[1]) - r22 * np.sin(q[1]))
    beta = d2 - Px - r13 * (L6 + L5 * np.cos(q[5])) - L5 * r12 * np.sin(q[5])
    gamma = L3 + L4 * np.cos(q[3])
    phi = L3 * gamma + L4**2 * np.sin(q[3])**2 + L4 * np.cos(q[3]) * gamma

    test_sin = (beta*gamma - alpha * L4 * np.sin(q[3])) / phi
    if (test_sin > 1 or test_sin< -1):
        print("theta 3 error !")

    q[2] = np.arcsin((beta * gamma - alpha * L4 * np.sin(q[3])) / phi)


    A = np.cos(q[5]) * (r33*np.cos(q[1]) - r23 * np.sin(q[1])) + np.sin(q[5]) * (r32*np.cos(q[1]) - r22 * np.sin(q[1]))
    B = r13 * np.cos(q[5]) + r12 * np.sin(q[5])

    q[4] = np.arcsin(B*np.cos(q[2] + q[3]) - A * np.sin(q[2] + q[3]))

    return q