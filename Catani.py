#%%
import rasterio
from geohazards import Catani, plot_rasterio, Shalstab
import geopandas as gpd

#%%
path=r'G:\Unidades compartidas\Proyectos(Fede)\Tarso\Amenaza\TRIGRS'
dem_path = f'{path}/dem.tif'
acc_path = f'{path}/aacum.tif'
geo_path = r'G:\Unidades compartidas\Proyectos(Fede)\Tarso\Amenaza\insumos/UndGeol25.shp'

hmin = 'hmin'
hmax = 'hmax'
c = 'C_kPa'
phy = 'Phi_grad'
k = 'k_cmh'
gammas = 'Gamma_knm3'
q = 71.67

#%%
dem = rasterio.open(dem_path)
plot_rasterio(dem)

#%% CATANI
geo_path='/asd'

zs = Catani(dem_path, geo_path, 0.1, 3, f'{path}/zs.tif')

plot_rasterio(zs)
# %% SHALSTAB

Shalstab(dem_path, acc_path, geo_path, hmin, hmax, c, phy, k, gammas, q, stability='./Str5.tif', qcrit='./qcrit.tif', rm_zs=True, rm_qcrit=True)