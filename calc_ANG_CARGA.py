# calc_UEL.py
import numpy as np
import pandas as pd

def calc_ANG(params_UEL, n_pontos=120):
    XqUEL = params_UEL["XqUEL"]
    DESW = params_UEL["DESW"]
    DESW1 = params_UEL["DESW1"]
    DESW2 = params_UEL["DESW2"]
    XnUEL = params_UEL["XnUEL"]
    Vt = params_UEL["Vt"]
    Xd = params_UEL["Xd"]
    DELTA = (params_UEL["DELTA_UEL"])
    DELTA2 = (params_UEL["DELTA_UEL2"])
    Xq2 = params_UEL["XqUEL2"]
    k = params_UEL["k_uel"]

    # angulo de carga
    Vt2 = Vt**2
    P = np.linspace(0, 0.01*n_pontos, n_pontos+1)
    Q = P/(np.tan(DESW*(np.pi/180)))-(Vt**2)/(XqUEL)
    Q_A1 = P/(np.tan(DESW1*(np.pi/180)))-(Vt**2)/(Xq2)
    Q_A2 = P/(np.tan(DESW2*(np.pi/180)))-(Vt**2)/(XqUEL)

    # UEL
    Q1 = P/(np.tan(DELTA*(np.pi/180)))-Vt*k
    Q2 = P/(np.tan(DELTA2*(np.pi/180)))-Vt*k

    df_resultado = pd.DataFrame({
        "P_ANG_CARGA": P,
        "Q_ANG_CARGA": Q,
        "P_UEL": P,
        "Q_UEL": Q1,
        "P_AUX1": P,
        "P_AUX2": P,
        "Q_AUX1": Q_A1,
        "Q_AUX2": Q_A2,
        "Q_UEL2": Q2
    })

    return {
        "Tabela": df_resultado
       }  

