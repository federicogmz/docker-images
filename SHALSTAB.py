#%%
import geohazards

path = r'/mnt/g/Unidades compartidas/av_torrencial_dabeiba/3.Estudio_Detalle/SHALTAB'
demPath = path + r'/demCuenca.tif'
geoPath = path + r'/geoCuenca.shp'
shalstabPath = path + r'/asd.tif'
criticalRainPath = path + r'/criticalRain.tif'

geoColumns = ['C_kPa', 'Phi_deg', 'Gamma_Nm3', 'Ks_ms', 'hmin', 'hmax']

q=50 # mm/h

exportShalstab = True
exportCriticalRain = False

shalstab, criticalRain = geohazards.SHALSTAB(demPath, geoPath, geoColumns, q, shalstabPath, criticalRainPath, exportShalstab, exportCriticalRain)

shalstab.plot()

print(shalstab.attrs['Reporte'])
