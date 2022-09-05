import os
import numpy as np
import matplotlib.pyplot as plt
import rasterio
from rasterio.plot import show
import geopandas as gpd
from osgeo import gdal
import rioxarray
from geocube.api.core import make_geocube

np.seterr(divide='ignore', invalid='ignore')

def Shalstab(dem_path, acc_path, q, geo_path, stability='./stability.tif', qcrit='./qcrit.tif', rm_zs=False, rm_qcrit=False):

    """
    Function to get stability model and critical rainfall using SHALSTAB model from  Montgomery and Dietrich (1984)
    """

    #Imports rasters and opens them with rasterio
    dem = rioxarray.open_rasterio(dem_path, masked=True)
    acc_ = rioxarray.open_rasterio(acc_path, masked=True)

    gdf = gpd.read_file(geo_path)

    #Rasterize from geology
    parameters = make_geocube(
        vector_data=gdf,
        #measurements=["C_kPa"],
        resolution=(-12.5, 12.5),
        fill=-9999,
        output_crs=dem.rio.crs
    )
    
    #Calculates slope
    slope = calculate_slope(dem_path)
    slope_array = slope.read(1)

    #Sets no data to np.nan
    slope_array[slope_array == slope.nodata] = np.nan

    #Converts slope from degrees to radians
    slope_rad = np.radians(slope_array)

    zs = rasterio.open(r'G:\Unidades compartidas\Proyectos(Fede)\Tarso\Amenaza\TRIGRS\zs.asc').read(1)

    #Calculates zs from Catani
    #zs = Catani(dem_path, geo_path, hmin=hmin, hmax=hmax).read(1)
    #if rm_zs:
    #    pass
    #else:
    #    os.remove('./zs.tif')

    #Unit weight of water
    gammaw = 9.81

    #Left from unconditional conditions
    izq_un = np.tan(slope_rad)

    #Right from unstable condition
    der_unstable = np.tan(friction_rad) + (cohesion / (gamma * zs * np.cos(slope_rad)**2))

    #Right from stable condition
    der_stable = (1 - (gammaw/gamma)) * np.tan(friction_rad) + (cohesion / (gamma * zs * np.cos(slope_rad)**2))

    #Left from main equation
    izq = acc / dem.res[0]

    #Right from main equation
    der = ((0.01 * k * (zs * np.cos(slope_rad)) * np.sin(slope_rad)) / (0.001*q)) * ((gamma / gammaw) * (1 - np.tan(slope_rad) / np.tan(friction_rad)) + (cohesion / (gammaw * zs * np.cos(slope_rad)**2 * np.tan(friction_rad)))) 

    #Shalstab stability categories
    shalstab = np.where(izq_un >= der_unstable, 2, np.where(izq_un < der_stable, 1, np.where(izq > der, 3, np.where(izq <= der, 4, 0))))

    #Calculates critical rainfall
    q_crit = (1000 * 0.01 * k * zs * np.cos(slope_rad) * np.sin(slope_rad)) * (dem.res[0] / acc) * ((gamma / gammaw) * (1 - (np.tan(slope_rad) / np.tan(friction_rad))) + cohesion / (gammaw * zs * np.cos(slope_rad)**2 * np.tan(friction_rad)))
    q_crit[q_crit==np.inf] = np.nanmax(q_crit)

    #Adds unconditional categories
    q_crit = np.where(izq_un >= der_unstable, -2, np.where(izq_un < der_stable, -1, q_crit))

    #Copys metadata from dem
    meta = dem.meta.copy()
    meta.update(compress='lzw', nodata=0)

    #Exports raster output file and assign metadata
    with rasterio.open(stability, 'w+', **meta) as out:
        out.write_band(1, shalstab)

    meta.update(compress='lzw', nodata=-9999)

    with rasterio.open(qcrit, 'w+', **meta) as out:
        out.write_band(1, q_crit)

    st = rasterio.open(stability)
    qc = rasterio.open(qcrit)

    if rm_qcrit:
        pass
    else:
        os.remove(qcrit)

    #Returns soil thickness rasterio file
    return st, qc

def Catani(dem_path, slope_path, hmin, hmax):

    """
    Function to get soil thickness from model S Catani et. al (2010)
    """

    #Imports dem and opens it with rasterio
    dem = rioxarray.open_rasterio(dem_path)
    slope = rioxarray.open_rasterio(slope_path)

    #Converts slope from degrees to radians
    slope_rad = np.radians(slope)

    hmax = hmax
    hmin = hmin

    #Calculates variables with Tangent
    tan_slope = np.tan(slope_rad)
    tan_slope_max = np.tan(np.nanmax(slope_rad))
    tan_slope_min = np.tan(np.nanmin(slope_rad))

    #Calculates soil thickness
    catani = hmax * (1 - ( (tan_slope - tan_slope_min) / (tan_slope_max - tan_slope_min) ) * (1 - (hmin / hmax)) )

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
