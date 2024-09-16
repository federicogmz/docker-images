#%%
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


caudales = pd.read_csv('/Users/federicogomez/Documents/Work/2024_pomca_sanbartolo/Susceptibilidad_AvT/Caudales/datos_caudales.csv')

estaciones = caudales.groupby('CodigoEstacion')

def calcular_indice_variabilidad(caudales, porcentajes):
    Qi = caudales[0]
    Qf = caudales[-1]
    Xi = porcentajes[0]
    Xf = porcentajes[-1]
    return (np.log10(Qi) - np.log10(Qf)) / (np.log10(Xi) - np.log10(Xf))

for nombre, grupo in estaciones:
    caudales_ordenados = sorted(grupo['Valor'][grupo['Valor'] > 0], reverse=True)
    n = len(caudales_ordenados)
    
    prob_excedencia = [(i+1)/n for i in range(n)]
    periodos_retorno = [1/p for p in prob_excedencia]
    
    plt.figure(figsize=(10, 6))
    plt.loglog(periodos_retorno, caudales_ordenados)
    plt.title(f'Curva de Duración de Caudales - Estación {nombre}\nÍndice de Variabilidad: {indice*100:.1f}%')
    plt.xlabel('Período de Retorno (años)')
    plt.ylabel('Caudal (m³/s)')
    plt.show()
    
    indice = calcular_indice_variabilidad(caudales_ordenados, periodos_retorno)

#%%

import geopandas as gpd

cuencas = gpd.read_file('/Users/federicogomez/Documents/Work/2024_pomca_sanbartolo/Susceptibilidad_AvT/Cuencas/Cuencas.shp')


def categorizar_densidad_drenaje(area, densidad):
    if area < 15:
        if densidad < 1.50:
            return 1
        elif 1.50 <= densidad < 2.00:
            return 2
        elif 2.00 <= densidad < 2.50:
            return 3
        elif 2.50 <= densidad <= 3.00:
            return 4
        else:
            return 5
    elif 16 <= area <= 50:
        if densidad < 1.20:
            return 1
        elif 1.20 <= densidad < 1.80:
            return 2
        elif 1.80 <= densidad < 2.00:
            return 3
        elif 2.00 <= densidad <= 2.50:
            return 4
        else:
            return 5
    else:  # area > 50
        if densidad < 1.00:
            return 1
        elif 1.00 <= densidad < 1.50:
            return 2
        elif 1.50 <= densidad < 2.00:
            return 3
        elif 2.00 <= densidad <= 2.50:
            return 4
        else:
            return 5

def categorizar_pendiente(area, pendiente):
    if area < 15:
        if pendiente < 20:
            return 1
        elif 20 <= pendiente < 35:
            return 2
        elif 35 <= pendiente < 50:
            return 3
        elif 50 <= pendiente <= 75:
            return 4
        else:
            return 5
    else:  # area >= 50
        if pendiente < 15:
            return 1
        elif 15 <= pendiente < 30:
            return 2
        elif 30 <= pendiente < 45:
            return 3
        elif 45 <= pendiente <= 65:
            return 4
        else:
            return 5

def categorizar_coeficiente_compacidad(coeficiente):
    if coeficiente < 1.125:
        return 5
    elif 1.125 <= coeficiente < 1.250:
        return 4
    elif 1.250 <= coeficiente < 1.375:
        return 3
    elif 1.375 <= coeficiente <= 1.500:
        return 2
    else:
        return 1

cuencas['Cat_DD'] = cuencas.apply(lambda row: categorizar_densidad_drenaje(row['Area (A)'], row['Densidad_d']), axis=1)
cuencas['Cat_P'] = cuencas.apply(lambda row: categorizar_pendiente(row['Area (A)'], row['Pendmean']), axis=1)
cuencas['Cat_CC'] = cuencas['Coef C'].apply(categorizar_coeficiente_compacidad)

cuencas.to_file('/Users/federicogomez/Documents/Work/2024_pomca_sanbartolo/Susceptibilidad_AvT/Cuencas/Cuencas.shp')
