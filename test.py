#%%
import rasterio as rio
import numpy as np

dem = rio.open(r'C:\Users\Federico Gomez\Downloads\DEM.tif')

raster = np.loadtxt(r'C:\Users\Federico Gomez\Downloads\rasterrr')

meta = dem.meta.copy()
meta.update(compress='lzw', nodata=-9999)

#Exports raster output file and assign metadata
with rio.open(r'G:\Mi unidad\raster.tif', 'w+', **meta) as out:
    out.write_band(1, raster)

