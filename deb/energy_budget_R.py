import numpy as np

from scipy.integrate import solve_ivp
import pandas as pd
from .dydt_KDEB import dydt_KDEB

import matplotlib.pyplot as plt
import time as TIME

import sys


def energy_budget(aux, p):
    """This function is called from the sim() function in sim.py. The function
    accepts two arguments, aux and p, and returns a NumPy array containing the
    results of the simulation.
    
    
    Args:
        aux (dict): Dictionary of auxiliary variables
        p (dict): Dictionary of species-specific parameters

    Returns:
        list: List containing the following values:
            - tEVHRX:
            - Lw: 
            - Ww: 
            - WwR0: 
            - X: 
            - Fish_en:
    """
    i = 0
    IC = aux["IC"]
    tArr = aux["tArr"]
    tc = tArr[0]
    rezTmp = []
    
    t = tArr
    
    R = solve_ivp(
        dydt_KDEB,
        [t[0], t[-1]],
        IC,
        args=[aux, p, IC],
        method="RK45",
        t_eval=t,
        )

    time = R["t"].T
    time = time.reshape((time.shape[0], 1))
    EVHR = R["y"].T
    data = np.hstack((time, EVHR))
    rezTmp.append(data)

    rezTmp = np.concatenate(rezTmp)
    tEVHRX = rezTmp
    
    # Repro
    Repro = tEVHRX[:, 4]
    dateV = aux["dateArr"]
    
    # finding indices of first days of december and last days of march for 
    # each year
    spwn_idx = np.vstack(
        (
            np.where((dateV.month == 12) & (dateV.day == 1))[0],
            np.where((dateV.month == 3) & (dateV.day == 31))[0])
    ).T
    
    spwn_durations = (spwn_idx[:, 1] - spwn_idx[:, 0])+1
    
    nr_reproductions = len(spwn_idx[:, 0])    
    for i in range(1, nr_reproductions):
        daily_repro_output = Repro[spwn_idx[i, 1]] / spwn_durations[i]
        curr_repro_idx = range(spwn_idx[i, 0], spwn_idx[i, 1])
        cum_curr_repro_output = range(1, spwn_durations[i])*daily_repro_output
        Repro[curr_repro_idx] = np.max(0, Repro[curr_repro_idx] - cum_curr_repro_output)
        
        if spwn_idx[i, 1] < len(Repro):
            Repro[spwn_idx[i, 1]:] = Repro[spwn_idx[i, 1]] + Repro[spwn_idx[i-1, 1]:] - cum_curr_repro_output[-1]
        
        if np.abs(Repro[spwn_idx[i, 1]]) > 1e-07:
            pass


    tEVHRX[:, 4] = Repro

    Lw = np.divide((tEVHRX[:, 2]**(1/3)), p["del_M"])
    Ww = aux["w"]*(p["d_V"]*tEVHRX[:, 2] + (p["w_E"]/p["mu_E"])*(tEVHRX[:, 1]+tEVHRX[:, 4]))
    WwR0 = aux["w"]*(p["d_V"]*tEVHRX[:, 2] + (p["w_E"]/p["mu_E"])*(tEVHRX[:, 1]))
    Fish_en = p["d_v"] * tEVHRX[:, 2] + (p["mu_V"]/p["w_V"]) +tEVHRX[:, 1]+tEVHRX[:, 4]- ((p["d_V"] * tEVHRX(1,3) * p["mu_V"] / p["w_V"]) + tEVHRX(1,2) + tEVHRX(1,5))

    X = tEVHRX[:, 5]

    rez = [rezTmp[:, 0], Lw,  Ww, WwR0, X, Fish_en]

    return rez
