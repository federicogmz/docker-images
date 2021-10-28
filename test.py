from modelos import Catani, plot_rasterio
path='/Volumes/GoogleDrive/.shortcut-targets-by-id/18bBJu8Sm269Nzq4Won8o62yPgwdYqT6D/Cartografía/Catany'
dem_path = f'{path}/DEM.tif'
geo_path = f'{path}/UGS.shp'
#zs = Catani(dem_path, geo_path, 'Hmin', 'Hmax', './jhnkklkmjk.tif')
#plot_rasterio(zs)
import rasterio
from rasterio.plot import show
show(rasterio.open(dem_path))