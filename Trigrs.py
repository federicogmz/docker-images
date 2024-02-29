#%%
from geohazards import TRIGRS, geohazards

path = r'/mnt/c/Users/Fede/OneDrive - Universidad Nacional de Colombia/Work/2023/Supia/MM_Urbano/Casco Urbano/Insumos'

# geohazards().preprocess_dem(demPath)

t = TRIGRS()
t(dem_path=f'{path}/DTM_CU.tif',
  geo_path=f'{path}/LitologiaCascoU.shp',
  out_path='/home/geohazards/TRIGRS',
  hora=4,
  cri=30.02,
  fosm=True)

#%%
import rioxarray
import numpy as np
import pandas as pd
import xarray as xr
import geopandas as gpd
from matplotlib import colors
import matplotlib.pyplot as plt

out_path='/home/geohazards/TRIGRS'
path = r'/mnt/c/Users/Fede/OneDrive - Universidad Nacional de Colombia/Work/2023/Supia/MM_Urbano/Casco Urbano/Insumos'

gdf = gpd.read_file(f'{path}/LitologiaCascoU.shp')
zones = gdf['Zona'].max()

cohesion = gdf['Cohesion'].values
friccion = gdf['Friccion'].values
gamma = gdf['Gamma'].values
ks = gdf['Ks'].values

variation_C = 40 / 100
C_means = cohesion
C_vectors = {}
C_vars = {}
C_stds = {}
for i, C_mean in enumerate(C_means, start=1):
    C_vectors[i] = [C_mean, C_mean + (C_mean * variation_C), C_mean - (C_mean * variation_C)]
    C_vars[i] = np.var(C_vectors[i])
    C_stds[i] = np.std(C_vectors[i])

variation_phi = 13 / 100
phi_means = friccion
phi_vectors = {}
phi_vars = {}
phi_stds = {}
for i, phi_mean in enumerate(phi_means, start=1):
    phi_vectors[i] = [phi_mean, phi_mean + (phi_mean * variation_phi), phi_mean - (phi_mean * variation_phi)]
    phi_vars[i] = np.var(phi_vectors[i])
    phi_stds[i] = np.std(phi_vectors[i])

variation_uws = 7 / 100
uws_means = gamma
uws_vectors = {}
uws_vars = {}
uws_stds = {}
for i, uws_mean in enumerate(uws_means, start=1):
    uws_vectors[i] = [uws_mean, uws_mean + (uws_mean * variation_uws), uws_mean - (uws_mean * variation_uws)]
    uws_vars[i] = np.var(uws_vectors[i])
    uws_stds[i] = np.std(uws_vectors[i])

variation_ks = 90 / 100
ks_means = ks
ks_vectors = {}
ks_vars = {}
ks_stds = {}
for i, ks_mean in enumerate(ks_means, start=1):
    ks_vectors[i] = [ks_mean, ks_mean + (ks_mean * variation_ks), ks_mean - (ks_mean * variation_ks)]
    ks_vars[i] = np.var(ks_vectors[i])
    ks_stds[i] = np.std(ks_vectors[i])


def plot_percentage_variance(porc):
    ind = np.arange(4)
    porc.plot.bar(legend=False)
    plt.ylabel('Porcentaje de la varianza')
    plt.xticks(ind, ('Cohesión', 'Ángulo de\nfricción', 'Peso\nEspecífico', 'Ks'), rotation=0)
    plt.tight_layout()
    plt.show()

def plot_Ind_conf_hora(Ind_conf_hora):
    c_map = colors.ListedColormap(['red', 'yellow', 'green'])
    bounds = [0, 1.0, 2.5, 10]
    norm = colors.BoundaryNorm(bounds, c_map.N)
    plt.imshow(Ind_conf_hora, cmap=c_map, norm=norm)

def read_result(file_path):
    data = np.genfromtxt(file_path, skip_header=6, delimiter=' ')
    data = np.where(data == -9999, np.nan, data)
    data = np.where(data >= 10, 10, data)
    return data

# Reading results
FS_medio = read_result(f'{out_path}/Resultados/TRfs_min_M_1.txt')
FS_cohe = read_result(f'{out_path}/Resultados/TRfs_min_C_1.txt')
FS_phi = read_result(f'{out_path}/Resultados/TRfs_min_P_1.txt')
FS_uws = read_result(f'{out_path}/Resultados/TRfs_min_G_1.txt')
FS_ks = read_result(f'{out_path}/Resultados/TRfs_min_K_1.txt')

zonas = np.genfromtxt(f'{out_path}/zonas.asc', skip_header=6, delimiter=' ')
zonas = np.where(zonas == -9999, np.nan, zonas)

def dFS(FS, mean_list):
    for i in range(1, zones + 1):
        if i == 1:
            dFS = np.where(zonas == i, FS / (0.1 * mean_list[i - 1]), FS)
        else:
            dFS = np.where(zonas == i, FS / (0.1 * mean_list[i - 1]), dFS)
    return dFS

dFS_c = dFS(FS_cohe, C_means)
dFS_phi = dFS(FS_phi, phi_means)
dFS_uws = dFS(FS_uws, uws_means)
dFS_ks = dFS(FS_ks, ks_means)

def VF(dF, var_list):
    for i in range(1, zones + 1):
        if i == 1:
            vF = np.where(zonas == i, (dF ** (2)) * var_list[i], dF)
        else:
            vF = np.where(zonas == i, (dF ** (2)) * var_list[i], vF)
    return vF

vF_C = VF(dFS_c, C_vars)
vF_phi = VF(dFS_phi, phi_vars)
vF_uws = VF(dFS_uws, uws_vars)
vF_ks = VF(dFS_ks, ks_vars)

std_Fs = np.sqrt(vF_C + vF_phi + vF_uws + vF_ks)

vF = vF_C + vF_phi + vF_uws + vF_ks

Ind_conf = (FS_medio - 1) / std_Fs
Ind_conf = np.where(Ind_conf < 0, -9999, Ind_conf)
Ind_conf = np.where(Ind_conf == np.math.inf, 9999, Ind_conf)

porc = pd.DataFrame([np.nanmean(vF_C / vF), np.nanmean(vF_phi / vF),
                    np.nanmean(vF_uws / vF), np.nanmean(vF_ks / vF)])

plot_percentage_variance(porc)
plot_Ind_conf_hora(Ind_conf)

dem = rioxarray.open_rasterio(f'{path}/DTM_CU.tif', mask_and_scale=True)
Ind_conf = xr.DataArray(Ind_conf, coords=[dem.coords['y'],dem.coords['x']])
Ind_conf.rio.write_nodata(-9999, inplace=True)
Ind_conf.rio.write_crs(dem.rio.crs, inplace=True)
Ind_conf.rio.to_raster(f'{out_path}/Ind_Conf.tif')