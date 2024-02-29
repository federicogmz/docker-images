#%%

################
### SHALSTAB ###
################

import geohazards

path = r'/mnt/c/Users/Fede/OneDrive - Universidad Nacional de Colombia/Work/2023/Supia/Workingfolder'
demPath = path + r'/DTM.tif'
geoPath = path + r'/Geo/GeoRural.shp'
zsPath = path + r'/zs.tif'
shalstabPath = path + r'/shalstab.tif'
criticalRainPath = path + r'/criticalRain.tif'

geoColumns = ['C', # kPa
              'Phi', # °
              'Gamma', # kN/m3
              'K mh' # m/
              ]

q = 89.82 # mm/h

exportShalstab = False
exportZs = True
exportCriticalRain = False

shalstab = geohazards.SHALSTAB(demPath, geoPath, geoColumns, q, zsPath, shalstabPath, criticalRainPath, exportZs, exportShalstab, exportCriticalRain)

s, cr = shalstab()

print(s.attrs['Reporte'])