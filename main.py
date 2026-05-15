# main.py
from dados import params
from dados import params_UEL
from calc_iterativo import calc_curva_ifd
from calc_iterativo import calc_curva_OELTH
from calc_iterativo import calc_curva_OELPK
from calc_iterativo import calc_curva_MEL
from calc_iterativo import calc_curva_MEL2
from calc_iterativo import calc_curva_MEL3
from calc_ANG_CARGA import calc_ANG
from curvas import curvas
from plotagem import plot_curva_saturada, plot_curva_n_saturada
import pandas as pd
from openpyxl import load_workbook

resultado_ifd = calc_curva_ifd(params)
df = resultado_ifd["Tabela"]

resultado_OELTH = calc_curva_OELTH(params)
df1 = resultado_OELTH["Tabela"]

resultado_OELPK = calc_curva_OELPK(params)
df2 = resultado_OELPK["Tabela"]

resultado_MEL = calc_curva_MEL(params)
df3 = resultado_MEL["Tabela"]

resultado_MEL2 = calc_curva_MEL2(params)
df7 = resultado_MEL2["Tabela"]

resultado_MEL3 = calc_curva_MEL3(params)
df8 = resultado_MEL3["Tabela"]

resultado_ANG = calc_ANG(params_UEL)
df4 = resultado_ANG["Tabela"]

calc_curvas = curvas(params)
df5 = calc_curvas["Tabela1"]
df6 = calc_curvas["Tabela2"]

print("DELTA:", resultado_ifd["DELTA"])
print("ET:", resultado_OELTH["ET"])
print("Elq:", resultado_ifd["Elq"])
print("ET_TRAD:", resultado_ifd["ET_TRAD"])
print("IFD:", resultado_ifd["IFD"])
print("VQ:", resultado_ifd["VQ"])
print("IQ:", resultado_ifd["IQ"])
print("ID:", resultado_ifd["ID"])
print("EI_EQ:", resultado_ifd["EI_EQ"])
print("EQ:", resultado_ifd["EQ"])
#print("EI_EQ_vz:", resultado_OELTH["EI_EQ_vz"])

'''opcao = input("Você quer a curva com saturação (s) ou sem saturação (n)? ")'''

'''if opcao.lower().startswith("s"):
    plot_curva_saturada(df, df1, df2, df3, df4, df5, df6, params["V"], params["Sb"], params["Sn"], params["P_op"], params["Q_op"])
else:
    plot_curva_n_saturada(df, df1, df2, df3, df4, df5, df6, params["V"], params["Sb"], params["Sn"], params["P_op"], params["Q_op"])'''

#plot_curva_saturada(df, df1, df2, df3, df4, df5, df6,  params["V"], params["Sb"], params["Sn"], params["Pmin"], params["Pmax"])
plot_curva_saturada(df1, df3, df4, df5, df6, params["V"], params["Sb"], params["Sn"], params["Pmin"], params["Pmax"], params["P_op"], params["Q_op"])


arquivo_saida = "Curva_resultados.xlsx"

with pd.ExcelWriter(arquivo_saida, engine="openpyxl") as writer:
    df.to_excel(writer, sheet_name="IFD", index=False)
    df1.to_excel(writer, sheet_name="OELTH", index=False)
    df2.to_excel(writer, sheet_name="OELPK", index=False)
    df3.to_excel(writer, sheet_name="MEL", index=False)
    df4[["P_ANG_CARGA", "Q_ANG_CARGA"]].to_excel(writer, sheet_name="ANGULO DE CARGA", index=False)
    df4[["P_UEL", "Q_UEL"]].to_excel(writer, sheet_name="UEL", index=False)
    df5[["P_nominal", "Q_nominal"]].to_excel(writer, sheet_name="NOMINAL", index=False)
    df5[["P_polar", "Q_polar"]].to_excel(writer, sheet_name="SALIENCIA POLAR", index=False)
    df5[["P_estab", "Q_estab"]].to_excel(writer, sheet_name="ESTABILIDADE TEORICA", index=False)
    df5[["P_pratica", "Q_pratica"]].to_excel(writer, sheet_name="ESTABILIDADE PRÁTICA", index=False) 
    df6[["P_SCL", "Q_SCL"]].to_excel(writer, sheet_name="SCL", index=False)  
print(f"Resultados exportados para {arquivo_saida}")

