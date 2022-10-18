#%%
from geohazards import Shalstab
import geopandas as gpd
import os
import rioxarray
from geocube.api.core import make_geocube

path=r'G:\Unidades compartidas\Proyectos(Fede)\Tarso\Amenaza\TRIGRS'
dem_path = f'{path}/dem.tif'
acc_path = f'{path}/aacum.tif'
geo_path = r'G:\Unidades compartidas\Proyectos(Fede)\Tarso\Amenaza\insumos/UndGeol25.shp'
slope_path = f'{path}/slope.tif'
#%%
q = 71.67

s = Shalstab(dem_path, acc_path, q, stability='./Str5.tif', qcrit='./qcrit.tif', rm_zs=True, rm_qcrit=True)

#%%
os.remove('slope.tif')
os.remove('Str5.tif')
os.remove('qcrit.tif')


#%%
from geohazards import Catani
zs = Catani(dem_path, slope_path, 0.1, 3)


#%%

dem = rioxarray.open_rasterio(dem_path, masked=True)
gdf = gpd.read_file(geo_path)

out_grid = make_geocube(
    vector_data = gdf,
    resolution = dem.rio.resolution(),
    fill = -9999,
    output_crs = dem.rio.crs
)
