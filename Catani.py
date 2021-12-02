#%%
import rasterio
from geohazards import Catani, plot_rasterio
import geopandas as gpd

#%%
path='/Volumes/GoogleDrive/Geteilte Ablagen/ETCR(9)/3.Mov. Masa/11.TRIGRS/Yarí/Mapas tif'
dem_path = f'{path}/DTM_Yari_v2_Fixed_clip_re.tif'
geo_path = f'{path}/Geo_Yari_re.shp'

#%%
dem = rasterio.open(dem_path)
plot_rasterio(dem)

#%%
zs = Catani(dem_path, geo_path, 0.1, 'Hmax', f'{path}/zs.tif')
plot_rasterio(zs)