from scipy.io import loadmat

import pandas as pd
import numpy as np
from deb.energy_budget_R import energy_budget
import matplotlib.pyplot as plt

def sim(temp):
    """Within the function, the initial weight of the fish is calculated based 
    on the values of aux["init_Lw"], aux["w"], and several parameters from the 
    p dictionary. The initial weight is then used to calculate the initial 
    reserve density aux["IC"], which is a NumPy array containing two values.

    The function then sets aux["init_Ww"] to the product of aux["w"], p["d_V"], 
    and the second value of aux["IC"], multiplied by a scaling factor that 
    depends on several parameters from the p dictionary.

    The function then enters a try block that calls the energy_budget() 
    function from the deb.energy_budget_R library with the arguments aux and p. 
    The output of this function is assigned to a variable simu.

    The function then extracts several values from the simu output and 
    calculates two values FCR_en and FCR using NumPy functions. The function 
    then creates a commented-out block of code that appears to create a 
    visualization of the simulation results using matplotlib.pyplot.

    Args:
        temp (ndarray): temperature timeseries

    Returns:
        dict: dictionary containing timeseries of the following variables:
            - TTM: time to maturity
            - FCR: feed conversion ratio
            - FCR_en: feed conversion ratio (energy)
            - W2Y: weight at 2 years
            - FCR2Y: feed conversion ratio at 2 years
            - FCR2Y_en: feed conversion ratio at 2 years (energy)
            - MTEMP2Y: mean temperature at 2 years
            - NBRD2Y: number of days with temperature >= 28C at 2 years
    """
    aux = {}
    sparus_data = loadmat(
        "data/allpars_Sparus_aurata.mat", simplify_cells=True
        )["allpars_Sparus_aurata"]
    species_name = "Sparus_aurata"

    temp_data = temp

    fArr = [0.9]
    aux["init_Lw"] = 6.5
    aux["w"] = 5

    SSTdata = temp_data.values
    aux["tT"] = np.vstack((np.arange(0, len(SSTdata), 1), SSTdata)).T

    # setting up figures goes here

    for k in range(len(fArr)):
        aux["f"] = fArr[k]
        if "tT" in aux:
            aux["dateArr"] = temp_data.index

            aux["tArr"] = aux["tT"][:, 0]

        p = sparus_data

        aux["IC"] = [
            (p["E_m"]*(aux["init_Lw"]*p["del_M"])**3),
            (aux["init_Lw"]*p["del_M"])**3,
            np.min(
                [p["E_Hp"], 
                p["E_Hp"]*(
                    aux["init_Lw"]/(
                        p["L_p"]/p["del_M"]
                        )
                    )])
            , 0, 0] 
        aux["init_Ww"] = aux["w"] * (
            p["d_V"] * aux["IC"][1] * (
                1 + (p["ome"] * aux["IC"][0]/(p["E_m"] * aux["IC"][1]))))

    try:
        simu = energy_budget(aux, p)
        Lw = simu[1]
        Ww = simu[2]
        WwR0 = simu[3]
        X = simu[4]
        fish_en = simu[5]

        FCR_en = np.divide(X, fish_en)
        FCR = np.divide(X*(p["w_X"]/p["mu_X"]), ((Ww-aux["init_Ww"])/aux["w"]))

        # Uncomment for parameter visualization
        
        # fig, axs = plt.subplots(3, 1)
        # ax_c = axs[2].twinx()
        # ax_d = axs[1].twinx()
        # axs[0].set_title("FCR")
        # axs[0].plot(FCR_en)
        # axs[0].plot(FCR)
        # axs[0].legend(["FCR_en", "FCR"])
        # axs[0].set_ylabel("FCR")
        # axs[0].set_xlabel("Days")
        # axs[1].set_title("Ww")
        # axs[1].plot(Ww)
        # axs[1].plot(WwR0)
        # axs[1].set_xlabel("Days")
        # axs[1].set_ylabel("Ww")
        # ax_d.plot(temp_data.values, "--")
        # ax_d.set_ylabel("temp (K)")
        # axs[1].legend(["Ww", "WwR0"])
        # axs[2].set_title("Lw")
        # ax_c.plot(temp_data.values, "--")
        # ax_c.set_ylabel("temp (K)")
        # axs[2].plot(Lw)
        # axs[2].set_xlabel("Days")
        # axs[2].set_ylabel("Lw")
        # plt.show()
        
        # TTM (400G)
        ttm = np.where(Ww >= 400)[0][0]

        # FCR[TTM]
        FCR_ = FCR[ttm]
        
        # FCR_en[TTM]
        FCR_en_ = FCR_en[ttm]

        # WEIGHT (2y)
        w2y = Ww[2*365]

        # FCR (2y)
        fcr2y = FCR[2*365]
        
        # FCR_en (2y)
        fcr2y_en = FCR_en[2*365]

        # Mean temp (2y)
        mtemp = np.mean(SSTdata[:2*365])

        # NBRD >= 28C (2y)
        nbrd = np.where(SSTdata[:2*365]-273.15 >= 28)[0].shape[0]

        return {
        "TTM": ttm,
        "FCR": FCR_,
        "FCR_en": FCR_en_,
        "W2Y": w2y,
        "FCR2Y": fcr2y,
        "FCR2Y_en": fcr2y_en,
        "MTEMP2Y": mtemp,
        "NBRD2Y": nbrd
        }

    except IndexError:
        print("FAILED")
        return {
        "TTM": np.nan,
        "FCR": np.nan,
        "FCR_en": np.nan,
        "W2Y": np.nan,
        "FCR2Y": np.nan,
        "FCR2Y_en": np.nan,
        "MTEMP2Y": np.nan,
        "NBRD2Y": np.nan
        }
        

        
        
        

        