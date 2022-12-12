#%%
from geohazards import Catani

path = r''
demPath = f'{path}/dem.tif'
geoPath = f'{path}/geologia.shp'

# Para usar con valores globales de mínimo y máximo espesor
zs = Catani(demPath, '', 0.1, 3)

# Para usar con valores locales de mínimo y máximo espesor según
# los atributos de la geología (geoPath)
hmin = ''
hmax = ''
zs = Catani(demPath, geoPath, hmin, hmax)

# Exportar mapa de espesor del suelo
zs.rio.to_raster(f'{path}/zs.tif')
