import numpy as np
import pandas as pd

def curvas(params):
    Xd = params["XD"]
    Xq = params["XQ"]
    Vt = params["VTERM"]
    refMin = params["refMin"]
    refMax = params["refMax"]               
    refSCL = params["refSCL"]
    offset = 0.1

    centro = 0.5 * (Vt**2)*(1/Xd + 1/Xq)
    raio = 0.5 * (Vt**2)*(1/Xq - 1/Xd)
    theta = np.linspace(0, 180, 180)

        # Calculo da curva da saliencia polar
    P = raio * np.sin(theta*np.pi/180)
    Q = raio * np.cos(theta*np.pi/180) - centro

    # Calculo da curva nominal
    Q1 = np.cos(theta*np.pi/180)
    P1 = np.sin(theta*np.pi/180)
 
    # Estabilidade teórica
    theta2 = np.linspace(0, 75, 75)
    P2 = 2 * raio * (np.sin(np.deg2rad(theta2))**2) * np.tan(np.deg2rad(theta2))
    Q2 = -(centro + raio - P2/np.tan(np.deg2rad(theta2)))

    # completa até 180 com NaN
    theta2_full = np.linspace(0, 180, 180)
    P2_full = np.full_like(theta2_full, np.nan, dtype=float)
    Q2_full = np.full_like(theta2_full, np.nan, dtype=float)

    P2_full[:len(P2)] = P2
    Q2_full[:len(Q2)] = Q2

    # Estabilidade prática
    P3 = np.zeros_like(theta2)

    for i in range(75):
        if P2[i] - offset < 0:
            P3[i] = 0
        else:
            P3[i] = P2[i] - offset

    Q3 = np.zeros_like(theta2)    
    for i in range(len(theta2)):   
        if (P2[i]**2 - (P2[i] - offset)**2) >= 0:
            Q3[i] = Q2[i] + np.sqrt(P2[i]**2 - (P2[i] - offset)**2)
        else:
            Q3[i] = Q3[i+1]

    P3_full = np.full_like(theta2_full, np.nan, dtype=float)
    Q3_full = np.full_like(theta2_full, np.nan, dtype=float)

    # completar os valores que são 0
    idx = np.where(Q3 != 0)[0]

    if len(idx) > 0: 
        first_nonzero_index = idx[0]
        Q3[:first_nonzero_index] = Q3[first_nonzero_index]

    P3_full[:len(P3)] = P3
    Q3_full[:len(Q3)] = Q3

    # Curva do SCL
    min = 0.0
    Min_Nom = min*Vt
    Min_max = min*refMax
    Min_min = min*refMin

    Q4 = Vt*refSCL*np.cos(theta*np.pi/180)
    P4_final = Vt*refSCL*np.sin(theta*np.pi/180)
            
    '''P4_1 = np.linspace(0, 1.195, int(1.195/0.005) + 1)
    P4_2 = np.linspace(1.195, 0, int(1.195/0.005) + 1)
    P4 = np.concatenate((P4_1, P4_2))

   
    ang_1 = np.arcsin(-P4_1/refSCL/Vt)
    ang_2 = np.arcsin(P4_2/refSCL/Vt)
    ang = np.concatenate((ang_1, ang_2))

    for i in range(1, len(ang)):
        if np.isnan(ang[i]):
            ang[i] = ang[i-1]

    Vterm = refSCL * Vt * np.cos(ang)

    n = len(Vterm)
    meio = n // 2  # ponto médio do vetor

    Q4 = np.zeros(n)
    for i in range(n):
        if abs(Vterm[i]) <= Min_Nom:  # dentro do limite
            if i < meio:  # primeira metade
                Q4[i] = Min_Nom
            else:         # segunda metade
                Q4[i] = -Min_Nom
        else:  # fora do limite
            if i < meio:  # primeira metade
                Q4[i] = Vterm[i]
            else:         # segunda metade
                Q4[i] = -Vterm[i]

    P4_final = np.zeros_like(Vterm)

    for i in range(n):
        if np.abs(Q4[i]) <= Min_Nom:
            P4_final[i] = 100
        else:
            P4_final[i] = P4[i]'''
            
   
    df_resultado = pd.DataFrame({
        "P_nominal": P1,
        "Q_nominal": Q1,
        "P_polar": P,
        "Q_polar": Q,
        "P_estab": P2_full,
        "Q_estab": Q2_full,
        "P_pratica": P3_full,
        "Q_pratica": Q3_full,
    })

    df_resultado1 = pd.DataFrame({
        "P_SCL": P4_final,
        "Q_SCL": Q4
    })

    return {
        "Tabela1": df_resultado,
        "Tabela2": df_resultado1
    }

