#%%
# SALGAR

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib

matplotlib.rcParams.update({'font.size': 16})

rain = pd.read_csv(r'C:\Users\Federico Gomez\Documents\rain\rainfall_12m.hdr',
    skiprows=5,
    index_col=2, parse_dates=True,
    infer_datetime_format=True,
    usecols = (1,2,3))[80:-160]

rain[' Lluvia'] *= 60
rain['acumulado'] = (rain[' Lluvia'] / 12).cumsum()

Titulo = u'Precipitación'
fig, ax = plt.subplots(1,figsize=(8,6), sharex=True)
# Primer eje
l1 = ax.plot(rain[' Lluvia'], lw=1.5, color='teal', alpha=0.8)
ax.set_ylim(bottom=0)
ax.set_ylabel('Rainfall ($mmh^{-1}$)')
ax.set_xlabel(u'Date')
# Formato de eje X
ax.xaxis.set_minor_locator(mdates.HourLocator(interval = 6))
plt.setp(ax.xaxis.get_majorticklabels(), ha="center", rotation=0 )
ax.xaxis.set_major_locator(mdates.HourLocator(interval = 6))
time_format = mdates.DateFormatter('%d/%m\n%H:%M')
ax.xaxis.set_major_formatter(time_format)
plt.tight_layout()
plt.show()

#%%
# POTRERITO
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import pandas as pd

matplotlib.rcParams.update({'font.size': 16})

serie = pd.read_excel(r'G:\Mi unidad\Laboral-Academico\Maestría\Tesis\Potrerito\Insumos\Hidrografas.xlsx').loc[:,['Time','Clark','Rain']].set_index('Time')
serie.Rain = serie.Rain / 12
serie.index = pd.to_datetime(serie.index, format='%H:%M:%S')

fig, ax = plt.subplots(1,figsize=(10,8), sharex=True)
#ax.set_title(Titulo)
# Primer eje
l1 = ax.plot(serie['Clark'], label=f'Streamflow', lw=1.5, color='olive', alpha=0.8)
ax.set_ylim(bottom=-0.05)
ax.set_ylabel(u'Streamflow $mm^3$')
ax.set_xlabel(u'Time')
ax1 = ax.twinx()
l2 = ax1.plot(serie['Rain'], label=f'Rainfall', lw=1.5, ls='--', color='teal', alpha=0.8)
ax1.set_ylim(bottom=0)
ax1.set_ylabel(u'Rainfall $mm$')
# Formato de eje X
ax.xaxis.set_minor_locator(mdates.HourLocator(interval = 56))
plt.setp(ax.xaxis.get_majorticklabels(), ha="center", rotation=0 )
ax.xaxis.set_major_locator(mdates.MinuteLocator(interval = 56))
time_format = mdates.DateFormatter('%H:%M')
ax.xaxis.set_major_formatter(time_format)
# Leyenda
ls = l1+l2
ax.legend(ls, [l.get_label() for l in ls], loc='upper right', fancybox=True, shadow=True, ncol=1)
plt.tight_layout()
plt.show()

