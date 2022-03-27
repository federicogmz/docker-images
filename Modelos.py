#%%
from geohazards.geohazards import FS, Catani, Topog
import rasterio
from geohazards import *
import geopandas as gpd

#%%
path='/Volumes/GoogleDrive-111862222221066300789/Shared drives/Proyectos(Fede)/Cartagena_2021'
dem_path = f'{path}/insumos/DEM.tif'
geo_path = f'{path}/v1/UGsv1.shp'

#%%
dem = rasterio.open(dem_path)
plot_rasterio(dem)

#%%
gpd.read_file(geo_path)

#%%
zs = Catani(dem_path, geo_path, hmin, hmax, output='./zs.tif')
zw = Topog(dem_path, f'{path}/v1/Aacu1.tif', geo_path, zs_path, q, ks='ks', output='./zw.tif')
fs = FS(dem_path, geo_path, zw_path, c, phi, gammas, output='./FS.tif')
