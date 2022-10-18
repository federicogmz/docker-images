#%%
from geohazards import Catani, plot_rasterio
import rasterio

path=r'G:\Unidades compartidas\Proyectos(Fede)\Tarso\Urbano/'
dem_path = f'{path}dem.tif'
slope_path = f'{path}slope.tif'
dem = rasterio.open(dem_path)
plot_rasterio(dem)

zs = Catani(dem_path, slope_path, 0.1, 2.3)
zs.plot()
zs.rio.to_raster(f'{path}zs.tif')
# zw = Topog(dem_path, f'{path}/v1/Aacu1.tif', geo_path, zs_path, q, ks='ks', output='./zw.tif')
# fs = FS(dem_path, geo_path, zw_path, c, phi, gammas, output='./FS.tif')
    # %%
