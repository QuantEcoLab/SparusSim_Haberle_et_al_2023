import numpy as np
# from spline1 import spline1
from scipy.interpolate import interp1d

def event_function(t, EVHRX, aux, p, EVHR_):
    
    w = aux["w"] * (p["d_V"] * EVHRX[1] * (1 + (p["ome"] * EVHRX[0]/(p["E_m"] * EVHRX[1]))))
    return w - 400
    
def dydt_KDEB(t, EVHRX, aux, p, EVHR_):
    # Modify for new code
    E = EVHRX[0]  # J, reserve energy
    V = EVHRX[1]  # cm^3, structural volume
    E_H = EVHRX[2]  # J , cumulated energy inversted into maturity
    E_R = EVHRX[3]  # J, reproduction buffer
    X = EVHRX[4]  # -, food input

    if E_H <= p["E_Hb"]:
        s_M = 1
    elif E_H > p["E_Hb"] and E_H < p["E_Hj"]:
        p["L_b"] = 0.363 # total length at birth - value from AmP database
        s_M = V**(1/3)/p["L_b"]
    else:
        s_M = p["s_M"]

    if "tT" in aux:
        Tf = interp1d(aux["tT"][:, 0], aux["tT"][:,1], kind="zero")
        T = Tf(t)
    else:
        T = aux["T"] + 273.15

    # Add the new condition for temperature
    if T < (273.15 + 12):
        f = 0.3
    elif "tf" in aux:
        ff = interp1d(aux["tf"][:, 0], aux["tf"][:,1], kind="zero")
        f = ff(t)
    else:
        f = aux["f"]
    
    if p["T_AH"] != None:
        s_H_ratio = (1+np.exp(p["T_AH"]/p["T_H"] - p["T_AH"]/p["T_ref"]))/(1+np.exp(p["T_AH"]/p["T_H"] - p["T_AH"]/T))
        c_T = np.exp(p["T_A"]/p["T_ref"] - p["T_A"]/T) * ((T >= p["T_ref"]) * s_H_ratio + (T < p["T_ref"])) 
    else:
        c_T = np.exp(p["T_A"]/p["T_ref"] - p["T_A"]/T)

    p_AmT = c_T * p["p_Am"] * s_M
    v_T = c_T * p["v"] * s_M
    p_MT = c_T * p["p_M"]
    p_TT = c_T * p["p_T"]
    k_JT = c_T * p["k_J"]
    p_XmT = p_AmT / p["kap_X"]
    k_tM = c_T * p["k_M"]

    if E_H < p["E_Hb"]:    
        pX = 0  # embryo stage -> f=0
    else:
        pX = f * p_XmT * V**(2/3)

    pA = p["kap_X"] * pX
    pM = p_MT * V
    pT = p_TT * V**(2/3)
    pS = pM + pT
    pC = (E/V) * (p["E_G"] * v_T * V**(2/3) + pS ) / (p["kap"] * E/V + p["E_G"]) 
    pJ = k_JT * E_H
   
    missing_pC = pS - (p["kap"]*pC)
    if missing_pC < 0:
        missing_pC = 0

    if (E_R < 0) or (E_H < 0):
        output = [-E, 0, 0, 0, 0]
        return output

    dE = pA-pC
    dV = np.max((0, (p["kap"]*pC -pS)/p["E_G"]))
    if E_H < p["E_Hp"]:
        dH = (1-p["kap"]) * pC - pJ - missing_pC
        dR = 0
    else:
        dH = 0
        dR = (1-p["kap"]) * pC - pJ - missing_pC
        
    dX = pX

    output = [dE, dV, dH, dR, dX]

    return output