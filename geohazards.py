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
