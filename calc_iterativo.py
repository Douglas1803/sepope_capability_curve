# curva_ifd.py
import numpy as np
import pandas as pd

def calc_curva_ifd(params, n_pontos=120):
    Xd   = params["XD"]
    Xq   = params["XQ"]     
    Xld  = params["XlD"]
    Xlld = params["XllD"]
    XL   = params["XL"]
    Ra   = params["RA"]
    AG   = params["AG"]
    BG   = params["BG"]
    PNOM = params["PNOM"]
    QNOM = params["QNOM"]
    VNOM = params["VNOM"]  
    Vterm = params["VTERM"] 

    # --- Cálculos iniciais ---
    Snom = PNOM + 1j * QNOM
    SNOM_CONJ = np.conjugate(Snom)

    EQ = VNOM + (Ra + 1j * Xq) * SNOM_CONJ
    ang_EQ = np.angle(EQ)
    ang_SNOM_CONJ = np.pi/2 if SNOM_CONJ == 0 else np.angle(SNOM_CONJ)

    ID = np.abs(SNOM_CONJ) * np.sin(ang_EQ - ang_SNOM_CONJ)
    IQ = np.abs(SNOM_CONJ) * np.cos(ang_EQ - ang_SNOM_CONJ)
    VD = np.abs(VNOM) * np.sin(ang_EQ)
    VQ = np.abs(VNOM) * np.cos(ang_EQ)


    Ellq = VQ + Xld*ID + Ra*IQ
    Elq = Ellq + ID * (Xlld - Xld)
    ET = Elq + (Xd - Xlld)*ID
    EI_EQ = ET + AG * np.exp(BG * (Ellq - 0.8))
    IFD = EI_EQ / (Xd - XL)
    ET_TRAD = np.sqrt(PNOM**2 + (QNOM + VNOM**2 / Xq)**2)

    # --- Vetores de P ---
    arrP_FORM = np.linspace(0, 0.01*n_pontos, n_pontos+1)
    arrQ_TRAD = np.zeros_like(arrP_FORM)
    arrQ_SAT  = np.zeros_like(arrP_FORM)
    arrDQ     = np.zeros_like(arrP_FORM)
    arrERRO   = np.zeros_like(arrP_FORM)

    # --- limites e iteração ---
    limA = 0.2
    limB = 0.1
    flag = 0
    deltQ_iter = limA
    count = 0
    comp_tabDELTA_EI = 99
    fator_intervalo = 0

    while flag == 0:
        i = 0
        deltQ_iter = limA
        arrQ_TRAD[i] = np.sqrt(ET_TRAD**2 - arrP_FORM[i]**2) - (Vterm**2)/Xq
        arrQ_SAT[i] = arrQ_TRAD[i] + deltQ_iter
        tabIT = np.conjugate((arrP_FORM[i] + 1j*arrQ_SAT[i])/Vterm)
        tabEQ = Vterm + (Ra + 1j*Xq)*tabIT
        tabDELTA = np.angle(tabEQ)
        tabTETA  = np.angle(tabIT)
        tabID = np.abs(tabIT)*np.sin(tabDELTA - tabTETA)
        tabIQ = np.abs(tabIT)*np.cos(tabDELTA - tabTETA)
        tabEllQ = Vterm*np.cos(tabDELTA) + Xld*tabID + Ra*tabIQ
        tabElQ = tabEllQ + tabID*(Xlld - Xld)
        tabE = tabElQ + (Xd - Xlld)*tabID
        tabEI = tabE + AG*np.exp(BG*(tabEllQ - 0.8))
        tabDELTA_EI = tabEI - EI_EQ

        if tabDELTA_EI > 0:
            limA -= 0.05
        else: 
            flag = 1
    
    deltQ_iter = limA + 0.06
    flag = 0

    for i in range(n_pontos + 1):
        arrQ_TRAD[i] = np.sqrt(ET_TRAD**2 - arrP_FORM[i]**2)-(Vterm**2)/Xq
        deltQ_iter = deltQ_iter + 0.01
        limA = deltQ_iter 

        while flag == 0:
            arrQ_SAT[i] = arrQ_TRAD[i] + deltQ_iter
            tabIT = np.conjugate((arrP_FORM[i] + 1j*arrQ_SAT[i])/Vterm)
            tabEQ = Vterm + (Ra + 1j*Xq)*tabIT
            tabDELTA = np.angle(tabEQ)
            tabTETA  = np.angle(tabIT)
            tabID = np.abs(tabIT)*np.sin(tabDELTA - tabTETA)
            tabIQ = np.abs(tabIT)*np.cos(tabDELTA - tabTETA)
            tabEllQ = Vterm*np.cos(tabDELTA) + Xld*tabID + Ra*tabIQ
            tabElQ = tabEllQ + tabID*(Xlld - Xld)
            tabE = tabElQ + (Xd - Xlld)*tabID
            tabEI = tabE + AG*np.exp(BG*(tabEllQ - 0.8))
            tabDELTA_EI = tabEI - EI_EQ   

            if np.abs(tabDELTA_EI) < 0.001 + fator_intervalo:
                flag = 1
                arrDQ[i] = deltQ_iter
                arrERRO[i] = tabDELTA_EI
                count = 0
                comp_tabDELTA_EI = 99
            else:
                if np.abs(tabDELTA_EI) < comp_tabDELTA_EI:
                    deltQ_iter = deltQ_iter - 0.001
                    count += 1
                    comp_tabDELTA_EI = np.abs(tabDELTA_EI)
                else:    
                    flag = 1
                    deltQ_iter = deltQ_iter + 0.001
                    arrDQ[i] = deltQ_iter
                    arrQ_SAT[i] = arrQ_TRAD[i] + deltQ_iter
                    tabIT = np.conjugate((arrP_FORM[i] + 1j*arrQ_SAT[i])/Vterm)
                    tabEQ = Vterm + (Ra + 1j*Xq)*tabIT
                    tabDELTA = np.angle(tabEQ)
                    tabTETA  = np.angle(tabIT)
                    tabID = np.abs(tabIT)*np.sin(tabDELTA - tabTETA)
                    tabIQ = np.abs(tabIT)*np.cos(tabDELTA - tabTETA)
                    tabEllQ = Vterm*np.cos(tabDELTA) + Xld*tabID + Ra*tabIQ
                    tabElQ = tabEllQ + tabID*(Xlld - Xld)
                    tabE = tabElQ + (Xd - Xlld)*tabID
                    tabEI = tabE + AG*np.exp(BG*(tabEllQ - 0.8))
                    tabDELTA_EI = tabEI - EI_EQ                   
                    arrERRO[i] = tabDELTA_EI
                    count = 0
                    comp_tabDELTA_EI = 99
            if count > 800:
                deltQ_iter = limA
                count = 0
                fator_intervalo += 0.001        
        flag = 0
        fator_intervalo = 0

    # --- DataFrame final ---
    df_resultado = pd.DataFrame({
        "P": arrP_FORM,
        "Qtrad": arrQ_TRAD,
        "DQ": arrDQ,
        "Qsat": arrQ_SAT,
        "Erro": arrERRO
    })



    return {
        "DELTA": ang_EQ,
        "ET": ET,
        "Elq": Elq,
        "ET_TRAD": ET_TRAD,
        "IFD": IFD,
        "VQ": VQ,
        "IQ": IQ,
        "ID": ID,
        "EI_EQ": EI_EQ, 
        "EQ": EQ,   
        "Tabela": df_resultado
    }

def calc_curva_OELTH(params, n_pontos=120):
    Xd   = params["XD"]
    Xq   = params["XQ"]     
    Xld  = params["XlD"]
    Xlld = params["XllD"]
    XL   = params["XL"]
    Ra   = params["RA"]
    AG   = params["AG"]
    BG   = params["BG"]
    PNOM = params["PNOM"]
    QNOM = params["QNOM"]
    VNOM = params["VNOM"]  
    Vterm = params["VTERM"]    
    FP = params["FP"] 
    avr_OELTH = params["OELTH"]
    avr_OELPK = params["OELPK"]
    avr_MEL = params["MEL"]

    # --- Cálculos iniciais ---
    Snom = PNOM + 1j * QNOM
    SNOM_CONJ = np.conjugate(Snom)

    EQ = VNOM + (Ra + 1j * Xq) * SNOM_CONJ
    ang_EQ = np.angle(EQ)
    ang_SNOM_CONJ = np.pi/2 if SNOM_CONJ == 0 else np.angle(SNOM_CONJ)

    ID = np.abs(SNOM_CONJ) * np.sin(ang_EQ - ang_SNOM_CONJ)
    IQ = np.abs(SNOM_CONJ) * np.cos(ang_EQ - ang_SNOM_CONJ)
    VD = np.abs(VNOM) * np.sin(ang_EQ)
    VQ = np.abs(VNOM) * np.cos(ang_EQ)

    Ellq = VQ + Xld*ID + Ra*IQ
    Elq = Ellq + ID * (Xlld - Xld)
    POT_DQ = (VD + Ra*ID)*ID + (VQ + Ra*IQ)*IQ
    ET = Elq + (Xd - Xlld)*ID
    ET_Orig = ET
    EI_EQ = ET + AG * np.exp(BG * (Ellq - 0.8))
    IFD = EI_EQ / (Xd - XL)
    ET_TRAD = np.sqrt(PNOM**2 + (QNOM + VNOM**2 / Xq)**2)
    EI_EQ_vz = 1 + AG*np.exp(BG*(1 - 0.8))

    # Ajustando as grandezas para os valores de OELTH
    EI_EQ = EI_EQ*avr_OELTH
    ET_TRAD = ET_TRAD*avr_OELTH

    # --- Vetores de P ---
    arrP_FORM = np.linspace(0, 0.01*n_pontos, n_pontos+1)
    arrQ_TRAD = np.zeros_like(arrP_FORM)
    arrQ_SAT  = np.zeros_like(arrP_FORM)
    arrDQ     = np.zeros_like(arrP_FORM)
    arrERRO   = np.zeros_like(arrP_FORM)

    # --- limites e iteração ---
    limA = 0.2
    limB = 0.1
    flag = 0
    deltQ_iter = limA
    count = 0
    comp_tabDELTA_EI = 99
    fator_intervalo = 0

    while flag == 0:
        i = 0
        deltQ_iter = limA
        arrQ_TRAD[i] = np.sqrt(ET_TRAD**2 - arrP_FORM[i]**2) - (Vterm**2)/Xq
        arrQ_SAT[i] = arrQ_TRAD[i] + deltQ_iter
        tabIT = np.conjugate((arrP_FORM[i] + 1j*arrQ_SAT[i])/Vterm)
        tabEQ = Vterm + (Ra + 1j*Xq)*tabIT
        tabDELTA = np.angle(tabEQ)
        tabTETA  = np.angle(tabIT)
        tabID = np.abs(tabIT)*np.sin(tabDELTA - tabTETA)
        tabIQ = np.abs(tabIT)*np.cos(tabDELTA - tabTETA)
        tabEllQ = Vterm*np.cos(tabDELTA) + Xld*tabID + Ra*tabIQ
        tabElQ = tabEllQ + tabID*(Xlld - Xld)
        tabE = tabElQ + (Xd - Xlld)*tabID
        tabEI = tabE + AG*np.exp(BG*(tabEllQ - 0.8))
        tabDELTA_EI = tabEI - EI_EQ

        if tabDELTA_EI > 0:
            limA -= 0.05
        else: 
            flag = 1
    
    deltQ_iter = limA + 0.06
    flag = 0

    for i in range(n_pontos + 1):
        arrQ_TRAD[i] = np.sqrt(ET_TRAD**2 - arrP_FORM[i]**2)-(Vterm**2)/Xq
        deltQ_iter = deltQ_iter + 0.01
        limA = deltQ_iter 

        while flag == 0:
            arrQ_SAT[i] = arrQ_TRAD[i] + deltQ_iter
            tabIT = np.conjugate((arrP_FORM[i] + 1j*arrQ_SAT[i])/Vterm)
            tabEQ = Vterm + (Ra + 1j*Xq)*tabIT
            tabDELTA = np.angle(tabEQ)
            tabTETA  = np.angle(tabIT)
            tabID = np.abs(tabIT)*np.sin(tabDELTA - tabTETA)
            tabIQ = np.abs(tabIT)*np.cos(tabDELTA - tabTETA)
            tabEllQ = Vterm*np.cos(tabDELTA) + Xld*tabID + Ra*tabIQ
            tabElQ = tabEllQ + tabID*(Xlld - Xld)
            tabE = tabElQ + (Xd - Xlld)*tabID
            tabEI = tabE + AG*np.exp(BG*(tabEllQ - 0.8))
            tabDELTA_EI = tabEI - EI_EQ   

            if np.abs(tabDELTA_EI) < 0.001 + fator_intervalo:
                flag = 1
                arrDQ[i] = deltQ_iter
                arrERRO[i] = tabDELTA_EI
                count = 0
                comp_tabDELTA_EI = 99
            else:
                if np.abs(tabDELTA_EI) < comp_tabDELTA_EI:
                    deltQ_iter = deltQ_iter - 0.001
                    count += 1
                    comp_tabDELTA_EI = np.abs(tabDELTA_EI)
                else:    
                    flag = 1
                    deltQ_iter = deltQ_iter + 0.001
                    arrDQ[i] = deltQ_iter
                    arrQ_SAT[i] = arrQ_TRAD[i] + deltQ_iter
                    tabIT = np.conjugate((arrP_FORM[i] + 1j*arrQ_SAT[i])/Vterm)
                    tabEQ = Vterm + (Ra + 1j*Xq)*tabIT
                    tabDELTA = np.angle(tabEQ)
                    tabTETA  = np.angle(tabIT)
                    tabID = np.abs(tabIT)*np.sin(tabDELTA - tabTETA)
                    tabIQ = np.abs(tabIT)*np.cos(tabDELTA - tabTETA)
                    tabEllQ = Vterm*np.cos(tabDELTA) + Xld*tabID + Ra*tabIQ
                    tabElQ = tabEllQ + tabID*(Xlld - Xld)
                    tabE = tabElQ + (Xd - Xlld)*tabID
                    tabEI = tabE + AG*np.exp(BG*(tabEllQ - 0.8))
                    tabDELTA_EI = tabEI - EI_EQ                   
                    arrERRO[i] = tabDELTA_EI
                    count = 0
                    comp_tabDELTA_EI = 99
            if count > 800:
                deltQ_iter = limA
                count = 0
                fator_intervalo += 0.001        
        flag = 0
        fator_intervalo = 0

    # --- DataFrame final ---
    df_resultado_OEL = pd.DataFrame({
        "P": arrP_FORM,
        "Qtrad": arrQ_TRAD,
        "DQ": arrDQ,
        "Qsat": arrQ_SAT,
        "Erro": arrERRO
    })

    return {
        "DELTA": ang_EQ,
        "ET": ET,
        "Elq": Elq,
        "ET_TRAD": ET_TRAD,
        "IFD": IFD,
        "EI_EQ": EI_EQ,
        "EI_EQ_vz": EI_EQ_vz,
        "Tabela": df_resultado_OEL
    }    

def calc_curva_OELPK(params, n_pontos=120):
    Xd   = params["XD"]
    Xq   = params["XQ"]     
    Xld  = params["XlD"]
    Xlld = params["XllD"]
    XL   = params["XL"]
    Ra   = params["RA"]
    AG   = params["AG"]
    BG   = params["BG"]
    PNOM = params["PNOM"]
    QNOM = params["QNOM"]
    VNOM = params["VNOM"]  
    Vterm = params["VTERM"]    
    FP = params["FP"] 
    avr_OELTH = params["OELTH"]
    avr_OELPK = params["OELPK"]
    avr_MEL = params["MEL"]

    # --- Cálculos iniciais ---
    Snom = PNOM + 1j * QNOM
    SNOM_CONJ = np.conjugate(Snom)

    EQ = VNOM + (Ra + 1j * Xq) * SNOM_CONJ
    ang_EQ = np.angle(EQ)
    ang_SNOM_CONJ = np.pi/2 if SNOM_CONJ == 0 else np.angle(SNOM_CONJ)

    ID = np.abs(SNOM_CONJ) * np.sin(ang_EQ - ang_SNOM_CONJ)
    IQ = np.abs(SNOM_CONJ) * np.cos(ang_EQ - ang_SNOM_CONJ)
    VD = np.abs(VNOM) * np.sin(ang_EQ)
    VQ = np.abs(VNOM) * np.cos(ang_EQ)

    Ellq = VQ + Xld*ID + Ra*IQ
    Elq = Ellq + ID * (Xlld - Xld)
    POT_DQ = (VD + Ra*ID)*ID + (VQ + Ra*IQ)*IQ
    ET = Elq + (Xd - Xlld)*ID
    ET_Orig = ET
    EI_EQ = ET + AG * np.exp(BG * (Ellq - 0.8))
    IFD = EI_EQ / (Xd - XL)
    ET_TRAD = np.sqrt(PNOM**2 + (QNOM + VNOM**2 / Xq)**2)
    EI_EQ_vz = 1 + AG*np.exp(BG*(1 - 0.8))

    # Ajustando as grandezas para os valores de OELTH
    ET = ET/avr_OELTH*avr_OELPK
    EI_EQ = EI_EQ*avr_OELPK
    ET_TRAD = ET_TRAD*avr_OELPK

    # --- Vetores de P ---
    arrP_FORM = np.linspace(0, 0.01*n_pontos, n_pontos+1)
    arrQ_TRAD = np.zeros_like(arrP_FORM)
    arrQ_SAT  = np.zeros_like(arrP_FORM)
    arrDQ     = np.zeros_like(arrP_FORM)
    arrERRO   = np.zeros_like(arrP_FORM)

    # --- limites e iteração ---
    limA = 0.2
    limB = 0.1
    flag = 0
    deltQ_iter = limA
    count = 0
    comp_tabDELTA_EI = 99
    fator_intervalo = 0

    while flag == 0:
        i = 0
        deltQ_iter = limA
        arrQ_TRAD[i] = np.sqrt(ET_TRAD**2 - arrP_FORM[i]**2) - (Vterm**2)/Xq
        arrQ_SAT[i] = arrQ_TRAD[i] + deltQ_iter
        tabIT = np.conjugate((arrP_FORM[i] + 1j*arrQ_SAT[i])/Vterm)
        tabEQ = Vterm + (Ra + 1j*Xq)*tabIT
        tabDELTA = np.angle(tabEQ)
        tabTETA  = np.angle(tabIT)
        tabID = np.abs(tabIT)*np.sin(tabDELTA - tabTETA)
        tabIQ = np.abs(tabIT)*np.cos(tabDELTA - tabTETA)
        tabEllQ = Vterm*np.cos(tabDELTA) + Xld*tabID + Ra*tabIQ
        tabElQ = tabEllQ + tabID*(Xlld - Xld)
        tabE = tabElQ + (Xd - Xlld)*tabID
        tabEI = tabE + AG*np.exp(BG*(tabEllQ - 0.8))
        tabDELTA_EI = tabEI - EI_EQ

        if tabDELTA_EI > 0:
            limA -= 0.05
        else: 
            flag = 1
    
    deltQ_iter = limA + 0.06
    flag = 0

    for i in range(n_pontos + 1):
        arrQ_TRAD[i] = np.sqrt(ET_TRAD**2 - arrP_FORM[i]**2)-(Vterm**2)/Xq
        deltQ_iter = deltQ_iter + 0.01
        limA = deltQ_iter 

        while flag == 0:
            arrQ_SAT[i] = arrQ_TRAD[i] + deltQ_iter
            tabIT = np.conjugate((arrP_FORM[i] + 1j*arrQ_SAT[i])/Vterm)
            tabEQ = Vterm + (Ra + 1j*Xq)*tabIT
            tabDELTA = np.angle(tabEQ)
            tabTETA  = np.angle(tabIT)
            tabID = np.abs(tabIT)*np.sin(tabDELTA - tabTETA)
            tabIQ = np.abs(tabIT)*np.cos(tabDELTA - tabTETA)
            tabEllQ = Vterm*np.cos(tabDELTA) + Xld*tabID + Ra*tabIQ
            tabElQ = tabEllQ + tabID*(Xlld - Xld)
            tabE = tabElQ + (Xd - Xlld)*tabID
            tabEI = tabE + AG*np.exp(BG*(tabEllQ - 0.8))
            tabDELTA_EI = tabEI - EI_EQ   

            if np.abs(tabDELTA_EI) < 0.001 + fator_intervalo:
                flag = 1
                arrDQ[i] = deltQ_iter
                arrERRO[i] = tabDELTA_EI
                count = 0
                comp_tabDELTA_EI = 99
            else:
                if np.abs(tabDELTA_EI) < comp_tabDELTA_EI:
                    deltQ_iter = deltQ_iter - 0.001
                    count += 1
                    comp_tabDELTA_EI = np.abs(tabDELTA_EI)
                else:    
                    flag = 1
                    deltQ_iter = deltQ_iter + 0.001
                    arrDQ[i] = deltQ_iter
                    arrQ_SAT[i] = arrQ_TRAD[i] + deltQ_iter
                    tabIT = np.conjugate((arrP_FORM[i] + 1j*arrQ_SAT[i])/Vterm)
                    tabEQ = Vterm + (Ra + 1j*Xq)*tabIT
                    tabDELTA = np.angle(tabEQ)
                    tabTETA  = np.angle(tabIT)
                    tabID = np.abs(tabIT)*np.sin(tabDELTA - tabTETA)
                    tabIQ = np.abs(tabIT)*np.cos(tabDELTA - tabTETA)
                    tabEllQ = Vterm*np.cos(tabDELTA) + Xld*tabID + Ra*tabIQ
                    tabElQ = tabEllQ + tabID*(Xlld - Xld)
                    tabE = tabElQ + (Xd - Xlld)*tabID
                    tabEI = tabE + AG*np.exp(BG*(tabEllQ - 0.8))
                    tabDELTA_EI = tabEI - EI_EQ                   
                    arrERRO[i] = tabDELTA_EI
                    count = 0
                    comp_tabDELTA_EI = 99
            if count > 800:
                deltQ_iter = limA
                count = 0
                fator_intervalo += 0.001        
        flag = 0
        fator_intervalo = 0

    # --- DataFrame final ---
    df_resultado_OELPK = pd.DataFrame({
        "P": arrP_FORM,
        "Qtrad": arrQ_TRAD,
        "DQ": arrDQ,
        "Qsat": arrQ_SAT,
        "Erro": arrERRO
    })

    return {
        "Tabela": df_resultado_OELPK
    }  

def calc_curva_MEL(params, n_pontos=120):
    Xd   = params["XD"]
    Xq   = params["XQ"]     
    Xld  = params["XlD"]
    Xlld = params["XllD"]
    XL   = params["XL"]
    Ra   = params["RA"]
    AG   = params["AG"]
    BG   = params["BG"]
    PNOM = params["PNOM"]
    QNOM = params["QNOM"]
    VNOM = params["VNOM"]  
    Vterm = params["VTERM"]    
    FP = params["FP"] 
    avr_OELTH = params["OELTH"]
    avr_OELPK = params["OELPK"]
    avr_MEL = params["MEL"]

    # --- Cálculos iniciais ---
    Snom = PNOM + 1j * QNOM
    SNOM_CONJ = np.conjugate(Snom)

    EQ = VNOM + (Ra + 1j * Xq) * SNOM_CONJ
    ang_EQ = np.angle(EQ)
    ang_SNOM_CONJ = np.pi/2 if SNOM_CONJ == 0 else np.angle(SNOM_CONJ)

    ID = np.abs(SNOM_CONJ) * np.sin(ang_EQ - ang_SNOM_CONJ)
    IQ = np.abs(SNOM_CONJ) * np.cos(ang_EQ - ang_SNOM_CONJ)
    VD = np.abs(VNOM) * np.sin(ang_EQ)
    VQ = np.abs(VNOM) * np.cos(ang_EQ)

    Ellq = VQ + Xld*ID + Ra*IQ
    Elq = Ellq + ID * (Xlld - Xld)
    POT_DQ = (VD + Ra*ID)*ID + (VQ + Ra*IQ)*IQ
    ET = Elq + (Xd - Xlld)*ID
    ET_Orig = ET
    EI_EQ = ET + AG * np.exp(BG * (Ellq - 0.8))
    IFD = EI_EQ / (Xd - XL)
    ET_TRAD = np.sqrt(PNOM**2 + (QNOM + VNOM**2 / Xq)**2)
    EI_EQ_vz = 1 + AG*np.exp(BG*(1 - 0.8))

    # Ajustando as grandezas para os valores de OELTH
    ET = ET/avr_OELTH*avr_MEL
    EI_EQ = EI_EQ*avr_MEL
    ET_TRAD = ET_TRAD*avr_MEL

    # Calculos específicos do MEL
    MEL_RAIO = np.sqrt((FP**2)+((Vterm**2/Xq) + np.sqrt(1-FP**2))**2)
    MEL_PPARTIDA = -Vterm**2 / Xq
    MEL_D = (Vterm**2)*(Xd - Xq)/(Xd * Xq)

    x1 = (Vterm**2)*(Xd - Xq)/(Xd * Xq)
    x2 = np.sqrt(PNOM**2 + (QNOM + 1/Xq)**2)
    x3 = avr_MEL * (x2 - x1) + x1
    flagPACEL = 0
    comp_tabDELTA_EI = 99

    # Dimensionando inicialmente os vetores
    arrP_FORM = np.linspace(0, 0.01*n_pontos, n_pontos+1)
    arrANG = np.zeros_like(arrP_FORM)
    arrQ_FORM = np.zeros_like(arrP_FORM)
    arrQ_SAT = np.zeros_like(arrP_FORM)
    arrDQ = np.zeros_like(arrP_FORM)
    arrQ_TRAD = np.zeros_like(arrP_FORM)
    arrP_TRAD = np.zeros_like(arrP_FORM)
    arrERRO = np.zeros_like(arrP_FORM)

    # Definindo limites de procura da DELTA_Q
    limA = +0.3
    limB = -0.3
    flag = 0
    deltQ_iter = limA
    fator_intervalo = 0
    count = 0
    i = 0

    try:
        arrANG[i] = np.arcsin(arrP_FORM[i] / x3)
    except Exception:
        # equivalente ao "Exit For"
        pass
    else:
        arrQ_TRAD[i] = -Vterm**2 / Xq + x3 * np.cos(arrANG[i])

    while flag == 0:
        i = 0
        deltQ_iter = limA
        arrQ_TRAD[i] = -(Vterm**2)/Xq + x3 * np.cos(arrANG[i])
        arrQ_SAT[i] = arrQ_TRAD[i] + deltQ_iter
        tabIT = np.conjugate((arrP_FORM[i] + 1j*arrQ_SAT[i])/Vterm)
        tabEQ = Vterm + (Ra + 1j*Xq)*tabIT 
        tabDELTA = np.angle(tabEQ)
        tabTETA  = np.angle(tabIT)
        tabID = np.abs(tabIT)*np.sin(tabDELTA - tabTETA)
        tabIQ = np.abs(tabIT)*np.cos(tabDELTA - tabTETA)
        tabEllQ = np.abs(Vterm)*np.cos(tabDELTA) + Xld*tabID + Ra*tabIQ
        tabElQ = tabEllQ + tabID*(Xlld - Xld)
        tabE = tabElQ + (Xd - Xlld)*tabID
        tabEI = tabE + AG*np.exp(BG*(tabEllQ - 0.8))
        tabDELTA_EI = tabEI - EI_EQ

        if tabDELTA_EI > 0:
            limA -= 0.05
        else: 
            flag = 1
    
    deltQ_iter = limA + 0.06
             
    for i in range(n_pontos + 1):
        try:
            arrANG[i] = np.arcsin(arrP_FORM[i] / x3)
            if arrANG[i] >= np.pi / 2:
                i = i   # não faz nada, só mantém como no VBA
        except Exception:
            # equivalente ao "Exit For"
            break
        else:
            arrQ_TRAD[i] = - Vterm**2  / Xq + x3 * np.cos(arrANG[i])
            deltQ_iter = deltQ_iter + 0.01
            limA = deltQ_iter
            flag = 0

        while flag == 0:
            if flagPACEL == 0:
                arrQ_SAT[i] = arrQ_TRAD[i] + deltQ_iter
                tabIT = np.conjugate((arrP_FORM[i] + 1j*arrQ_SAT[i])/Vterm)
                tabEQ = Vterm + (Ra + 1j*Xq)*tabIT
                tabDELTA = np.angle(tabEQ)
                tabTETA  = np.angle(tabIT)
                tabID = np.abs(tabIT)*np.sin(tabDELTA - tabTETA)
                tabIQ = np.abs(tabIT)*np.cos(tabDELTA - tabTETA)
                tabEllQ = Vterm*np.cos(tabDELTA) + Xld*tabID + Ra*tabIQ
                tabElQ = tabEllQ + tabID*(Xlld - Xld)
                tabE = tabElQ + (Xd - Xlld)*tabID
                tabEI = tabE + AG*np.exp(BG*(tabEllQ - 0.8))
                tabDELTA_EI = tabEI - EI_EQ   
                tabPACEL = tabE*Vterm/Xd*np.cos(tabDELTA) + ((1 / Xq - 1 / Xd) * Vterm ** 2) * np.cos(2 * tabDELTA)

                if tabPACEL < 0:
                    flag = 1
                    flagPACEL = 1
                elif np.abs(tabDELTA_EI) < 0.001 + fator_intervalo:
                    flag = 1
                    arrDQ[i] = deltQ_iter
                    arrERRO[i] = tabDELTA_EI
                    count = 0
                    comp_tabDELTA_EI = 99
                else:
                    if np.abs(tabDELTA_EI) < np.abs(comp_tabDELTA_EI):
                        deltQ_iter = deltQ_iter - 0.001
                        count += 1
                        comp_tabDELTA_EI = np.abs(tabDELTA_EI)
                    else:    
                        flag = 1
                        deltQ_iter = deltQ_iter + 0.001
                        arrDQ[i] = deltQ_iter
                        arrQ_SAT[i] = arrQ_TRAD[i] + deltQ_iter
                        tabIT = np.conjugate((arrP_FORM[i] + 1j*arrQ_SAT[i])/Vterm)
                        tabEQ = Vterm + (Ra + 1j*Xq)*tabIT
                        tabDELTA = np.angle(tabEQ)
                        tabTETA  = np.angle(tabIT)
                        tabID = np.abs(tabIT)*np.sin(tabDELTA - tabTETA)
                        tabIQ = np.abs(tabIT)*np.cos(tabDELTA - tabTETA)
                        tabEllQ = Vterm*np.cos(tabDELTA) + Xld*tabID + Ra*tabIQ
                        tabElQ = tabEllQ + tabID*(Xlld - Xld)
                        tabE = tabElQ + (Xd - Xlld)*tabID
                        tabEI = tabE + AG*np.exp(BG*(tabEllQ - 0.8))
                        tabDELTA_EI = tabEI - EI_EQ                   
                        arrERRO[i] = tabDELTA_EI
                        count = 0
                        comp_tabDELTA_EI = 99
                if count > 800:
                    deltQ_iter = limA
                    count = 0
                    fator_intervalo += 0.001 
            else:
                flag = 1   

        flag = 0
        fator_intervalo = 0   
        if flagPACEL:
            arrQ_SAT[i] = np.nan  # mantém o tamanho do array

    df_resultado_MEL = pd.DataFrame({
        "P": arrP_FORM,
        "Qtrad": arrQ_TRAD,
        "DQ": arrDQ,
        "Qsat": arrQ_SAT,
        "Erro": arrERRO
    })

    return {
        "Tabela": df_resultado_MEL
    }

def calc_curva_MEL2(params, n_pontos=120):
    Xd   = params["XD"]
    Xq   = params["XQ"]     
    Xld  = params["XlD"]
    Xlld = params["XllD"]
    XL   = params["XL"]
    Ra   = params["RA"]
    AG   = params["AG"]
    BG   = params["BG"]
    PNOM = params["PNOM"]
    QNOM = params["QNOM"]
    VNOM = params["VNOM"]  
    Vterm = params["VTERM"]    
    FP = params["FP"] 
    avr_OELTH = params["OELTH"]
    avr_OELPK = params["OELPK"]
    avr_MEL2 = params["MEL2"]

    # --- Cálculos iniciais ---
    Snom = PNOM + 1j * QNOM
    SNOM_CONJ = np.conjugate(Snom)

    EQ = VNOM + (Ra + 1j * Xq) * SNOM_CONJ
    ang_EQ = np.angle(EQ)
    ang_SNOM_CONJ = np.pi/2 if SNOM_CONJ == 0 else np.angle(SNOM_CONJ)

    ID = np.abs(SNOM_CONJ) * np.sin(ang_EQ - ang_SNOM_CONJ)
    IQ = np.abs(SNOM_CONJ) * np.cos(ang_EQ - ang_SNOM_CONJ)
    VD = np.abs(VNOM) * np.sin(ang_EQ)
    VQ = np.abs(VNOM) * np.cos(ang_EQ)

    Ellq = VQ + Xld*ID + Ra*IQ
    Elq = Ellq + ID * (Xlld - Xld)
    POT_DQ = (VD + Ra*ID)*ID + (VQ + Ra*IQ)*IQ
    ET = Elq + (Xd - Xlld)*ID
    ET_Orig = ET
    EI_EQ = ET + AG * np.exp(BG * (Ellq - 0.8))
    IFD = EI_EQ / (Xd - XL)
    ET_TRAD = np.sqrt(PNOM**2 + (QNOM + VNOM**2 / Xq)**2)
    EI_EQ_vz = 1 + AG*np.exp(BG*(1 - 0.8))

    # Ajustando as grandezas para os valores de OELTH
    ET = ET/avr_OELTH*avr_MEL2
    EI_EQ = EI_EQ*avr_MEL2
    ET_TRAD = ET_TRAD*avr_MEL2

    # Calculos específicos do MEL
    MEL_RAIO = np.sqrt((FP**2)+((Vterm**2/Xq) + np.sqrt(1-FP**2))**2)
    MEL_PPARTIDA = -Vterm**2 / Xq
    MEL_D = (Vterm**2)*(Xd - Xq)/(Xd * Xq)

    x1 = (Vterm**2)*(Xd - Xq)/(Xd * Xq)
    x2 = np.sqrt(PNOM**2 + (QNOM + 1/Xq)**2)
    x3 = avr_MEL2 * (x2 - x1) + x1
    flagPACEL = 0
    comp_tabDELTA_EI = 99

    # Dimensionando inicialmente os vetores
    arrP_FORM = np.linspace(0, 0.01*n_pontos, n_pontos+1)
    arrANG = np.zeros_like(arrP_FORM)
    arrQ_FORM = np.zeros_like(arrP_FORM)
    arrQ_SAT = np.zeros_like(arrP_FORM)
    arrDQ = np.zeros_like(arrP_FORM)
    arrQ_TRAD = np.zeros_like(arrP_FORM)
    arrP_TRAD = np.zeros_like(arrP_FORM)
    arrERRO = np.zeros_like(arrP_FORM)

    # Definindo limites de procura da DELTA_Q
    limA = +0.3
    limB = -0.3
    flag = 0
    deltQ_iter = limA
    fator_intervalo = 0
    count = 0
    i = 0

    try:
        arrANG[i] = np.arcsin(arrP_FORM[i] / x3)
    except Exception:
        # equivalente ao "Exit For"
        pass
    else:
        arrQ_TRAD[i] = -Vterm**2 / Xq + x3 * np.cos(arrANG[i])

    while flag == 0:
        i = 0
        deltQ_iter = limA
        arrQ_TRAD[i] = -(Vterm**2)/Xq + x3 * np.cos(arrANG[i])
        arrQ_SAT[i] = arrQ_TRAD[i] + deltQ_iter
        tabIT = np.conjugate((arrP_FORM[i] + 1j*arrQ_SAT[i])/Vterm)
        tabEQ = Vterm + (Ra + 1j*Xq)*tabIT 
        tabDELTA = np.angle(tabEQ)
        tabTETA  = np.angle(tabIT)
        tabID = np.abs(tabIT)*np.sin(tabDELTA - tabTETA)
        tabIQ = np.abs(tabIT)*np.cos(tabDELTA - tabTETA)
        tabEllQ = np.abs(Vterm)*np.cos(tabDELTA) + Xld*tabID + Ra*tabIQ
        tabElQ = tabEllQ + tabID*(Xlld - Xld)
        tabE = tabElQ + (Xd - Xlld)*tabID
        tabEI = tabE + AG*np.exp(BG*(tabEllQ - 0.8))
        tabDELTA_EI = tabEI - EI_EQ

        if tabDELTA_EI > 0:
            limA -= 0.05
        else: 
            flag = 1
    
    deltQ_iter = limA + 0.06
             
    for i in range(n_pontos + 1):
        try:
            arrANG[i] = np.arcsin(arrP_FORM[i] / x3)
            if arrANG[i] >= np.pi / 2:
                i = i   # não faz nada, só mantém como no VBA
        except Exception:
            # equivalente ao "Exit For"
            break
        else:
            arrQ_TRAD[i] = - Vterm**2  / Xq + x3 * np.cos(arrANG[i])
            deltQ_iter = deltQ_iter + 0.01
            limA = deltQ_iter
            flag = 0

        while flag == 0:
            if flagPACEL == 0:
                arrQ_SAT[i] = arrQ_TRAD[i] + deltQ_iter
                tabIT = np.conjugate((arrP_FORM[i] + 1j*arrQ_SAT[i])/Vterm)
                tabEQ = Vterm + (Ra + 1j*Xq)*tabIT
                tabDELTA = np.angle(tabEQ)
                tabTETA  = np.angle(tabIT)
                tabID = np.abs(tabIT)*np.sin(tabDELTA - tabTETA)
                tabIQ = np.abs(tabIT)*np.cos(tabDELTA - tabTETA)
                tabEllQ = Vterm*np.cos(tabDELTA) + Xld*tabID + Ra*tabIQ
                tabElQ = tabEllQ + tabID*(Xlld - Xld)
                tabE = tabElQ + (Xd - Xlld)*tabID
                tabEI = tabE + AG*np.exp(BG*(tabEllQ - 0.8))
                tabDELTA_EI = tabEI - EI_EQ   
                tabPACEL = tabE*Vterm/Xd*np.cos(tabDELTA) + ((1 / Xq - 1 / Xd) * Vterm ** 2) * np.cos(2 * tabDELTA)

                if tabPACEL < 0:
                    flag = 1
                    flagPACEL = 1
                elif np.abs(tabDELTA_EI) < 0.001 + fator_intervalo:
                    flag = 1
                    arrDQ[i] = deltQ_iter
                    arrERRO[i] = tabDELTA_EI
                    count = 0
                    comp_tabDELTA_EI = 99
                else:
                    if np.abs(tabDELTA_EI) < np.abs(comp_tabDELTA_EI):
                        deltQ_iter = deltQ_iter - 0.001
                        count += 1
                        comp_tabDELTA_EI = np.abs(tabDELTA_EI)
                    else:    
                        flag = 1
                        deltQ_iter = deltQ_iter + 0.001
                        arrDQ[i] = deltQ_iter
                        arrQ_SAT[i] = arrQ_TRAD[i] + deltQ_iter
                        tabIT = np.conjugate((arrP_FORM[i] + 1j*arrQ_SAT[i])/Vterm)
                        tabEQ = Vterm + (Ra + 1j*Xq)*tabIT
                        tabDELTA = np.angle(tabEQ)
                        tabTETA  = np.angle(tabIT)
                        tabID = np.abs(tabIT)*np.sin(tabDELTA - tabTETA)
                        tabIQ = np.abs(tabIT)*np.cos(tabDELTA - tabTETA)
                        tabEllQ = Vterm*np.cos(tabDELTA) + Xld*tabID + Ra*tabIQ
                        tabElQ = tabEllQ + tabID*(Xlld - Xld)
                        tabE = tabElQ + (Xd - Xlld)*tabID
                        tabEI = tabE + AG*np.exp(BG*(tabEllQ - 0.8))
                        tabDELTA_EI = tabEI - EI_EQ                   
                        arrERRO[i] = tabDELTA_EI
                        count = 0
                        comp_tabDELTA_EI = 99
                if count > 800:
                    deltQ_iter = limA
                    count = 0
                    fator_intervalo += 0.001 
            else:
                flag = 1   

        flag = 0
        fator_intervalo = 0   
        if flagPACEL:
            arrQ_SAT[i] = np.nan  # mantém o tamanho do array

    df_resultado_MEL2 = pd.DataFrame({
        "P": arrP_FORM,
        "Qtrad": arrQ_TRAD,
        "DQ": arrDQ,
        "Qsat": arrQ_SAT,
        "Erro": arrERRO
    })

    return {
        "Tabela": df_resultado_MEL2
    } 
def calc_curva_MEL3(params, n_pontos=120):
    Xd   = params["XD"]
    Xq   = params["XQ"]     
    Xld  = params["XlD"]
    Xlld = params["XllD"]
    XL   = params["XL"]
    Ra   = params["RA"]
    AG   = params["AG"]
    BG   = params["BG"]
    PNOM = params["PNOM"]
    QNOM = params["QNOM"]
    VNOM = params["VNOM"]  
    Vterm = params["VTERM"]    
    FP = params["FP"] 
    avr_OELTH = params["OELTH"]
    avr_OELPK = params["OELPK"]
    avr_MEL3 = params["MEL3"]

    # --- Cálculos iniciais ---
    Snom = PNOM + 1j * QNOM
    SNOM_CONJ = np.conjugate(Snom)

    EQ = VNOM + (Ra + 1j * Xq) * SNOM_CONJ
    ang_EQ = np.angle(EQ)
    ang_SNOM_CONJ = np.pi/2 if SNOM_CONJ == 0 else np.angle(SNOM_CONJ)

    ID = np.abs(SNOM_CONJ) * np.sin(ang_EQ - ang_SNOM_CONJ)
    IQ = np.abs(SNOM_CONJ) * np.cos(ang_EQ - ang_SNOM_CONJ)
    VD = np.abs(VNOM) * np.sin(ang_EQ)
    VQ = np.abs(VNOM) * np.cos(ang_EQ)

    Ellq = VQ + Xld*ID + Ra*IQ
    Elq = Ellq + ID * (Xlld - Xld)
    POT_DQ = (VD + Ra*ID)*ID + (VQ + Ra*IQ)*IQ
    ET = Elq + (Xd - Xlld)*ID
    ET_Orig = ET
    EI_EQ = ET + AG * np.exp(BG * (Ellq - 0.8))
    IFD = EI_EQ / (Xd - XL)
    ET_TRAD = np.sqrt(PNOM**2 + (QNOM + VNOM**2 / Xq)**2)
    EI_EQ_vz = 1 + AG*np.exp(BG*(1 - 0.8))

    # Ajustando as grandezas para os valores de OELTH
    ET = ET/avr_OELTH*avr_MEL3
    EI_EQ = EI_EQ*avr_MEL3
    ET_TRAD = ET_TRAD*avr_MEL3

    # Calculos específicos do MEL
    MEL_RAIO = np.sqrt((FP**2)+((Vterm**2/Xq) + np.sqrt(1-FP**2))**2)
    MEL_PPARTIDA = -Vterm**2 / Xq
    MEL_D = (Vterm**2)*(Xd - Xq)/(Xd * Xq)

    x1 = (Vterm**2)*(Xd - Xq)/(Xd * Xq)
    x2 = np.sqrt(PNOM**2 + (QNOM + 1/Xq)**2)
    x3 = avr_MEL3 * (x2 - x1) + x1
    flagPACEL = 0
    comp_tabDELTA_EI = 99

    # Dimensionando inicialmente os vetores
    arrP_FORM = np.linspace(0, 0.01*n_pontos, n_pontos+1)
    arrANG = np.zeros_like(arrP_FORM)
    arrQ_FORM = np.zeros_like(arrP_FORM)
    arrQ_SAT = np.zeros_like(arrP_FORM)
    arrDQ = np.zeros_like(arrP_FORM)
    arrQ_TRAD = np.zeros_like(arrP_FORM)
    arrP_TRAD = np.zeros_like(arrP_FORM)
    arrERRO = np.zeros_like(arrP_FORM)

    # Definindo limites de procura da DELTA_Q
    limA = +0.3
    limB = -0.3
    flag = 0
    deltQ_iter = limA
    fator_intervalo = 0
    count = 0
    i = 0

    try:
        arrANG[i] = np.arcsin(arrP_FORM[i] / x3)
    except Exception:
        # equivalente ao "Exit For"
        pass
    else:
        arrQ_TRAD[i] = -Vterm**2 / Xq + x3 * np.cos(arrANG[i])

    while flag == 0:
        i = 0
        deltQ_iter = limA
        arrQ_TRAD[i] = -(Vterm**2)/Xq + x3 * np.cos(arrANG[i])
        arrQ_SAT[i] = arrQ_TRAD[i] + deltQ_iter
        tabIT = np.conjugate((arrP_FORM[i] + 1j*arrQ_SAT[i])/Vterm)
        tabEQ = Vterm + (Ra + 1j*Xq)*tabIT 
        tabDELTA = np.angle(tabEQ)
        tabTETA  = np.angle(tabIT)
        tabID = np.abs(tabIT)*np.sin(tabDELTA - tabTETA)
        tabIQ = np.abs(tabIT)*np.cos(tabDELTA - tabTETA)
        tabEllQ = np.abs(Vterm)*np.cos(tabDELTA) + Xld*tabID + Ra*tabIQ
        tabElQ = tabEllQ + tabID*(Xlld - Xld)
        tabE = tabElQ + (Xd - Xlld)*tabID
        tabEI = tabE + AG*np.exp(BG*(tabEllQ - 0.8))
        tabDELTA_EI = tabEI - EI_EQ

        if tabDELTA_EI > 0:
            limA -= 0.05
        else: 
            flag = 1
    
    deltQ_iter = limA + 0.06
             
    for i in range(n_pontos + 1):
        try:
            arrANG[i] = np.arcsin(arrP_FORM[i] / x3)
            if arrANG[i] >= np.pi / 2:
                i = i   # não faz nada, só mantém como no VBA
        except Exception:
            # equivalente ao "Exit For"
            break
        else:
            arrQ_TRAD[i] = - Vterm**2  / Xq + x3 * np.cos(arrANG[i])
            deltQ_iter = deltQ_iter + 0.01
            limA = deltQ_iter
            flag = 0

        while flag == 0:
            if flagPACEL == 0:
                arrQ_SAT[i] = arrQ_TRAD[i] + deltQ_iter
                tabIT = np.conjugate((arrP_FORM[i] + 1j*arrQ_SAT[i])/Vterm)
                tabEQ = Vterm + (Ra + 1j*Xq)*tabIT
                tabDELTA = np.angle(tabEQ)
                tabTETA  = np.angle(tabIT)
                tabID = np.abs(tabIT)*np.sin(tabDELTA - tabTETA)
                tabIQ = np.abs(tabIT)*np.cos(tabDELTA - tabTETA)
                tabEllQ = Vterm*np.cos(tabDELTA) + Xld*tabID + Ra*tabIQ
                tabElQ = tabEllQ + tabID*(Xlld - Xld)
                tabE = tabElQ + (Xd - Xlld)*tabID
                tabEI = tabE + AG*np.exp(BG*(tabEllQ - 0.8))
                tabDELTA_EI = tabEI - EI_EQ   
                tabPACEL = tabE*Vterm/Xd*np.cos(tabDELTA) + ((1 / Xq - 1 / Xd) * Vterm ** 2) * np.cos(2 * tabDELTA)

                if tabPACEL < 0:
                    flag = 1
                    flagPACEL = 1
                elif np.abs(tabDELTA_EI) < 0.001 + fator_intervalo:
                    flag = 1
                    arrDQ[i] = deltQ_iter
                    arrERRO[i] = tabDELTA_EI
                    count = 0
                    comp_tabDELTA_EI = 99
                else:
                    if np.abs(tabDELTA_EI) < np.abs(comp_tabDELTA_EI):
                        deltQ_iter = deltQ_iter - 0.001
                        count += 1
                        comp_tabDELTA_EI = np.abs(tabDELTA_EI)
                    else:    
                        flag = 1
                        deltQ_iter = deltQ_iter + 0.001
                        arrDQ[i] = deltQ_iter
                        arrQ_SAT[i] = arrQ_TRAD[i] + deltQ_iter
                        tabIT = np.conjugate((arrP_FORM[i] + 1j*arrQ_SAT[i])/Vterm)
                        tabEQ = Vterm + (Ra + 1j*Xq)*tabIT
                        tabDELTA = np.angle(tabEQ)
                        tabTETA  = np.angle(tabIT)
                        tabID = np.abs(tabIT)*np.sin(tabDELTA - tabTETA)
                        tabIQ = np.abs(tabIT)*np.cos(tabDELTA - tabTETA)
                        tabEllQ = Vterm*np.cos(tabDELTA) + Xld*tabID + Ra*tabIQ
                        tabElQ = tabEllQ + tabID*(Xlld - Xld)
                        tabE = tabElQ + (Xd - Xlld)*tabID
                        tabEI = tabE + AG*np.exp(BG*(tabEllQ - 0.8))
                        tabDELTA_EI = tabEI - EI_EQ                   
                        arrERRO[i] = tabDELTA_EI
                        count = 0
                        comp_tabDELTA_EI = 99
                if count > 800:
                    deltQ_iter = limA
                    count = 0
                    fator_intervalo += 0.001 
            else:
                flag = 1   

        flag = 0
        fator_intervalo = 0   
        if flagPACEL:
            arrQ_SAT[i] = np.nan  # mantém o tamanho do array

    df_resultado_MEL3 = pd.DataFrame({
        "P": arrP_FORM,
        "Qtrad": arrQ_TRAD,
        "DQ": arrDQ,
        "Qsat": arrQ_SAT,
        "Erro": arrERRO
    })

    return {
        "Tabela": df_resultado_MEL3
    }      
     
