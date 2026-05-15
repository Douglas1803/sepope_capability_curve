# Dados Genéricos de Gerador Síncrono

import numpy as np
import matplotlib.pyplot as plt

# --------------------------------------------------------------------------
# Dados do Gerador
Ifd_fl = 1000            # Corrente de campo nominal arbitrária (A)
Sn = 100                 # Potência aparente nominal (MVA)
Vn = 13.8                # Tensão nominal (kV)
FP = 0.95                # Fator de potência nominal
Itn = Sn / np.sqrt(3) / Vn
Ifd_nl = 650             # Corrente de campo a vazio (A)

# --------------------------------------------------------------------------
# Parâmetros do Gerador (valores típicos)
Xd = 1.05
Xdl = 0.30
Xdll = 0.20
Xq = 0.70
Xl = 0.15
Ra = 0.003               # 0.3%

# --------------------------------------------------------------------------
# Parâmetros da Saturação (genéricos)
Ag = 0.02
Bg = 8.0
Ifd_nl_pu = Ifd_nl / Ifd_fl
Ifd_fl_pu = 1.0

# --------------------------------------------------------------------------
# Dados Operacionais
Pmin = 20                # MW
Pmax = 100               # MW

# --------------------------------------------------------------------------
# Sistema de Excitação
Refmax = 1.05
Refmin = 0.95
Vt_base = Vn
Ifd_base = Ifd_fl

# --------------------------------------------------------------------------
# Referências dos limitadores

# OEL (sobrecorrente de campo)
I_OELTH = 1050           # 105%
RefOelTh = I_OELTH / Ifd_base
avr_OELTH = RefOelTh

I_OELPK = 1100           # 110%
RefOelPK = I_OELPK / Ifd_base
avr_OELPK = RefOelPK

# MEL (mínima excitação)
I_MEL = 200              # ~20%
RefMel = I_MEL / Ifd_base
avr_MEL = RefMel

I_MEL2 = 250
RefMel2 = I_MEL2 / Ifd_base
avr_MEL2 = RefMel2

I_MEL3 = 300
RefMel3 = I_MEL3 / Ifd_base
avr_MEL3 = RefMel3

# SCL (corrente do estator)
I_SCL = 1.10 * Itn
RefSCL = I_SCL / Itn

# --------------------------------------------------------------------------
# Ângulo de carga
DESW = 30
DESW_A1 = 40
DESW_A2 = 28
Xq_UEB = 0.50
Xq_UEB2 = Xq
k_uel = 0.90

# --------------------------------------------------------------------------
# UEL
Delta_UEL = 65
Delta_UEL2 = 72

# --------------------------------------------------------------------------
# Dados para curva
Vte = Vn
Vt = Vte / Vn
Sb = 100                 # Base da curva (MVA)

Pn = FP
Qn = np.sin(np.arccos(FP))

# --------------------------------------------------------------------------
# Pontos de operação fictícios
P_op = 80
Q_op = -40

P2 = 60
Q2 = -30

P3 = 95
Q3 = -20

# --------------------------------------------------------------------------
# Dicionário principal
params = {
    "XD": Xd,
    "XQ": Xq,
    "XlD": Xdl,
    "XllD": Xdll,
    "XL": Xl,
    "RA": Ra,
    "AG": Ag,
    "BG": Bg,
    "PNOM": Pn,
    "QNOM": Qn,
    "VNOM": Vt,
    "VTERM": Vt,
    "FP": FP,
    "OELTH": avr_OELTH,
    "OELPK": avr_OELPK,
    "MEL": avr_MEL,
    "MEL2": avr_MEL2,
    "MEL3": avr_MEL3,
    "refMax": Refmax,
    "refMin": Refmin,
    "refSCL": RefSCL,
    "Sb": Sb,
    "Sn": Sn,
    "V": Vte,
    "P_op": P_op,
    "Q_op": Q_op,
    "P2": P2,
    "Q2": Q2,
    "P3": P3,
    "Q3": Q3,
    "Pmin": Pmin,
    "Pmax": Pmax
}

# --------------------------------------------------------------------------
# Dicionário UEL
params_UEL = {
    "XqUEL": Xq_UEB,
    "XqUEL2": Xq_UEB2,
    "Vt": Vt,
    "DESW": DESW,
    "XnUEL": 0.00001,
    "DELTA_UEL": Delta_UEL,
    "DELTA_UEL2": Delta_UEL2,
    "Xd": Xd,
    "DESW1": DESW_A1,
    "DESW2": DESW_A2,
    "k_uel": k_uel
}
