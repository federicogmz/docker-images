import os
import numpy as np
import matplotlib.pyplot as plt
import rasterio
from rasterio.plot import show
import geopandas as gpd
from osgeo import gdal
import rioxarray
import xarray as xr
from geocube.api.core import make_geocube
import richdem as rd

np.seterr(divide='ignore', invalid='ignore')

def SHALSTAB(demPath, geoPath, geoColumns, q, shalstabPath, criticalRainPath, exportShalstab, exportCriticalRain):

    #Imports rasters and opens them with rasterio
    dem = rioxarray.open_rasterio(demPath, mask_and_scale=True)

    dem_ = rd.LoadGDAL(demPath)
    accum = rd.FlowAccumulation(dem_, method='D8')
    accum = xr.DataArray(accum, coords=[dem.coords['y'],dem.coords['x']])

    gdf = gpd.read_file(geoPath)

    #Rasterize from geology
    c = make_geocube(vector_data=gdf, measurements = [geoColumns[0]], 
                    like=dem, fill = np.nan)[geoColumns[0]]
    p = make_geocube(vector_data=gdf, measurements = [geoColumns[1]], 
                    like=dem, fill = np.nan)[geoColumns[1]]
    p = np.radians(p)
    g = make_geocube(vector_data=gdf, measurements = [geoColumns[2]], 
                    like=dem, fill = np.nan)[geoColumns[2]]
    k = make_geocube(vector_data=gdf, measurements = [geoColumns[3]], 
                    like = dem, fill = np.nan)[geoColumns[3]]

    #Calculates slope
    slope = rd.TerrainAttribute(dem_, attrib='slope_radians')
    slope = xr.DataArray(slope, coords=[dem.coords['y'],dem.coords['x']])
    slope_rad = slope.where(slope!=-9999.)

    #Calculates zs from Catani
    zs = Catani(demPath, geoPath, hmin=geoColumns[4], hmax=geoColumns[5])

    #Unit weight of water
    gammaw = 9.81

    #Left from unconditional conditions
    izq_un = np.tan(slope_rad)

    #Right from unstable condition
    der_unstable = np.tan(p) + (c / (g * zs * np.cos(slope_rad)**2))

    #Right from stable condition
    der_stable = (1 - (gammaw/g)) * np.tan(p) + (c / (g * zs * np.cos(slope_rad)**2))

    #Left from main equation
    izq = accum / dem.rio.resolution()[0]

    #Right from main equation
    der = ((k * (zs * np.cos(slope_rad)) * np.sin(slope_rad)) / (q/1000)) * ((g / gammaw) * (1 - np.tan(slope_rad) / np.tan(p)) + (c / (gammaw * zs * np.cos(slope_rad)**2 * np.tan(p)))) 

    #Shalstab stability categories
    shalstab = np.where(izq_un >= der_unstable, 2, # Unconditionally unstable
                np.where(izq_un < der_stable, 1, # Unconditionally stable
                np.where(izq > der, 3, # Unstable
                np.where(izq <= der, 4, np.nan)))) # Stable

    shalstab = xr.DataArray(shalstab, coords=[dem.coords['y'],dem.coords['x']])
    shalstab.rio.write_nodata(-9999, inplace=True)

    if exportShalstab:
        shalstab.rio.to_raster(shalstabPath)

    #Calculates critical rainfall
    criticalRain = 1000 * (k * zs * np.cos(slope_rad) * np.sin(slope_rad)) * (dem.rio.resolution()[0] / accum) * ((g / gammaw) * (1 - (np.tan(slope_rad) / np.tan(p))) + c / (gammaw * zs * np.cos(slope_rad)**2 * np.tan(p)))
    criticalRain = criticalRain.where(criticalRain>=0)

    criticalRain = np.where(izq_un >= der_unstable, -2, 
                np.where(izq_un < der_stable, -1, criticalRain))
    criticalRain = xr.DataArray(criticalRain, coords=[dem.coords['y'],dem.coords['x']])

    criticalRain.rio.write_nodata(-9999, inplace=True)

    if exportCriticalRain:
        criticalRain.rio.to_raster(criticalRainPath)
    
    totalCeldas = shalstab.to_series().sum()

    try: 
        incondEstables = shalstab.to_series().value_counts()[1]
    except:
        incondEstables = 0
    try:
        incondInestables = shalstab.to_series().value_counts()[2]
    except:
        incondInestables = 0
    try:
        inestables = shalstab.to_series().value_counts()[3]
    except:
        inestables = 0
    try:
        estables = shalstab.to_series().value_counts()[4]
    except:
        estables = 0

    stabilityReport = ''

    stabilityReport += f'Incondicionalmente estables: {incondEstables*100/totalCeldas:.2f}%\n'
    stabilityReport += f'Incondicionalmente inestables: {incondInestables*100/totalCeldas:.2f}%\n'
    stabilityReport += f'Inestables: {inestables*100/totalCeldas:.2f}%\n'
    stabilityReport += f'Estables: {estables*100/totalCeldas:.2f}%'

    shalstab.attrs['Reporte'] = stabilityReport

    return shalstab, criticalRain

def Catani(demPath, geoPath, hmin=0.1, hmax=3.0):

    """
    Function to get soil thickness from model S Catani et. al (2010)
    """

    #Imports dem and opens it with rasterio
    dem = rioxarray.open_rasterio(demPath, mask_and_scale=True)
    dem_ = rd.LoadGDAL(demPath)

    slope = rd.TerrainAttribute(dem_, attrib='slope_radians')
    slope = xr.DataArray(slope, coords=[dem.coords['y'],dem.coords['x']])
    slope_rad = slope.where(slope!=-9999.)

    gdf = gpd.read_file(geoPath)

    if isinstance(hmin, str):
        hmin = make_geocube(vector_data=gdf, measurements = [hmin], 
                        like = dem, fill = np.nan)[hmin]
    else: hmin = hmin

    if isinstance(hmax, str):
        hmax = make_geocube(vector_data=gdf, measurements = [hmax], 
                        like = dem, fill = np.nan)[hmax]
    else: 
        hmax = hmax

    #Calculates variables with Tangent
    tan_slope = np.tan(slope_rad)
    tan_slope_max = np.tan(np.nanmax(slope_rad))
    tan_slope_min = np.tan(np.nanmin(slope_rad))

    #Calculates soil thickness
    catani = hmax * (1 - ( (tan_slope - tan_slope_min) / (tan_slope_max - tan_slope_min) ) * (1 - (hmin / hmax)) )

    catani = xr.DataArray(catani, coords=[dem.coords['y'], dem.coords['x']])

    #Returns soil thickness rasterio file
    return catani

def Topog(dem_path, acc_path, geo_path, zs_path, q, ks='ks', output='./zw.tif'):

    #Imports dem and opens it with rasterio
    dem = rasterio.open(dem_path)
    b = dem.res[0]

    #Calculates slope
    slope = calculate_slope(dem_path)
    slope_array = slope.read(1)

    #Sets no data to np.nan
    slope_array[slope_array == slope.nodata] = np.nan

    #Converts slope from degrees to radians
    slope_rad = np.radians(slope_array)

    acc_ = rasterio.open(acc_path)
    zs = rasterio.open(zs_path)

    #Rasterize ks from geology
    if isinstance(ks,float) or isinstance(ks,int):
        ks = ks
    else:
        ks = rasterize(dem_path, geo_path, attribute='ks').read(1)

    zw_ = q*(acc_/b)*ks*zs*np.cos(slope_rad)*np.sin(slope_rad)
    zw = np.where(zw_>=zs, zs, zw_)

    #Copies metadata from dem
    meta = dem.meta.copy()
    meta.update(compress='lzw', nodata=-9999)

    #Exports raster output file and assign metadata
    with rasterio.open(output, 'w+', **meta) as out:
        out.write_band(1, zw)
    
    #Returns soil thickness rasterio file
    return rasterio.open(output)

def FS(dem_path, geo_path, zw_path, c, phi, gammas, output='./FS.tif'):

    dem = rasterio.open(dem_path)
    zw = rasterio.open(zw_path)

    #Calculates slope
    slope = calculate_slope(dem_path)
    slope_array = slope.read(1)

    #Sets no data to np.nan
    slope_array[slope_array == slope.nodata] = np.nan

    #Converts slope from degrees to radians
    slope_rad = np.radians(slope_array)

    gammaw = 9.81

    C = rasterize(dem_path, geo_path, attribute=c).read(1)
    phi = rasterize(dem_path, geo_path, attribute=phi).read(1)
    gammas = rasterize(dem_path, geo_path, attribute=c).read(1)


    FS = C + (gammas - gammaw) * zw * (np.cos(slope_rad)**2) * np.tan(phi) / gammas * zw * np.sin(slope_rad) * np.cos(slope_rad)

    #Copies metadata from dem
    meta = dem.meta.copy()
    meta.update(compress='lzw', nodata=-9999)

    #Exports raster output file and assign metadata
    with rasterio.open(output, 'w+', **meta) as out:
        out.write_band(1, FS)
    
    #Returns soil thickness rasterio file
    return rasterio.open(output)

def rasterize(raster, polygon, attribute):

    #Imports raster and shapefile
    poly = gpd.read_file(polygon)
    rst = rasterio.open(raster)

    #Copys metadata from raster
    meta = rst.meta.copy()
    meta.update(compress='lzw')

    #Export raster output file and assign metadata
    with rasterio.open('rasterized.tif', 'w+', **meta) as out:
        out_arr = out.read(1)

        #Iters from each polygon and extracts the value of interest
        shapes = ((geom,value) for geom, value in zip(poly.geometry, poly[attribute]))

        #Rasterize each polygon with the value of interest
        burned = features.rasterize(shapes=shapes, fill=-9999, out=out_arr, transform=out.transform)

        #Exports raster output file
        out.write_band(1, burned)

    #Reads back rasterio file
    raster = rasterio.open('rasterized.tif')

    #Removes raster from system
    os.remove('rasterized.tif')

    #Returns rasterio file
    return raster

def calculate_slope(dem_path):

    """
    Function to calculate slope from DEM
    """

    #Imports dem and opens it with rasterio
    dem = rasterio.open(dem_path)
    
    #Calculates slope
    gdal.DEMProcessing('slope.tif', dem_path, 'slope')
    
    #Copys metadata from raster
    meta = dem.meta.copy()
    meta.update(compress='lzw')
    
    #Exports slope output file
    slope = rasterio.open('slope.tif', **meta)
    
    #Removes raster from system
    #os.remove('slope.tif')
    
    #Returns rasterio file
    return slope

def plot_rasterio(rio_dataset):

    """
    Function to plot a rasterio raster with legend
    """

    #Creates plot
    fig, ax = plt.subplots()

    #Reads rasterio raster as np.array
    array = rio_dataset.read(1)

    #Sets no data to np.nan
    array[array == rio_dataset.nodata] = np.nan

    #Adds hidden plot for colorbar
    image_hidden = ax.imshow(array)

    #Adds the rasterio plot
    show(rio_dataset, ax=ax)

    #Adds the colorbar
    fig.colorbar(image_hidden, ax=ax)

    #Returns the plot
    return plt.show()
