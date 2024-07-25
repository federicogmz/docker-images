#%%

#########################
### ESPESOR DEL SUELO ###
#########################

import richdem as rd
import xarray as xr
import numpy as np

dem_path = '/home/fdgmz/Documents/san_bartolo/DEM.tif'
hmin = 0.1
hmax = 3.0

dem = xr.open_dataarray(dem_path, mask_and_scale=True)
dem_ = rd.LoadGDAL(dem_path)
slope = rd.TerrainAttribute(dem_, attrib='slope_radians')
slope_rad = xr.DataArray(slope, coords=[dem.coords['y'],dem.coords['x']])

tan_slope = np.tan(slope_rad)
tan_slope_max = np.tan(np.nanmax(slope_rad))
tan_slope_min = np.tan(np.nanmin(slope_rad))

catani = hmax * (1 - ( (tan_slope - tan_slope_min) / (tan_slope_max - tan_slope_min) ) * (1 - (hmin / hmax)) )

catani = xr.DataArray(catani, coords=[dem.coords['y'], dem.coords['x']])

catani = catani.where(~np.isnan(dem), np.nan).squeeze()

catani.rio.to_raster('/home/fdgmz/Documents/san_bartolo/Susceptibilidad/catani.tif')


#%%
import rioxarray
import numpy as np
import xarray as xr
import pandas as pd
import seaborn as sns
from scipy import stats
import geopandas as gpd
import matplotlib.pyplot as plt
from geocube.api.core import make_geocube

from scipy.stats import boxcox, skew, ttest_ind
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.metrics import confusion_matrix, cohen_kappa_score, accuracy_score
from sklearn.model_selection import RepeatedStratifiedKFold, GridSearchCV, train_test_split

class POMCA():

    def __init__(self, dem_path:str, geo:gpd.GeoDataFrame):
        """
        Inicializa clase padre y crea variables base.

        Args:
            dem_path (str): Ruta al archivo de datos del modelo de elevación digital.
            geo (gpd.GeoDataFrame): GeoDataFrame de las unidades geológicas.
        """
        self.figsize = (10, 8)
        self.dem = rioxarray.open_rasterio(dem_path)

        self.geo = geo.to_crs(self.dem.rio.crs)
        self.extent = [self.geo.bounds.minx.values[0], self.geo.bounds.maxx.values[0],
                    self.geo.bounds.miny.values[0], self.geo.bounds.maxy.values[0]]

    def rasterize(self, gdf, column):

        raster = make_geocube(vector_data=gdf, measurements = [column], like=self.dem, fill = np.nan)[column]
        raster = raster.where(~np.isnan(self.dem))

        return raster


ruta_insumos = '/home/fdgmz/Documents/san_bartolo/Susceptibilidad/Insumos/'

ruta_muestra = '/home/fdgmz/Documents/san_bartolo/Inventario/Muestra.shp'
ruta_dem = ruta_insumos + 'dem.tif'
ruta_ugs = '/home/fdgmz/Documents/san_bartolo/Ugs_preliminar.shp'

ugs = gpd.read_file(ruta_ugs)

#%%

pomca = POMCA(ruta_dem, ugs)


muestra_gdf = gpd.read_file(ruta_muestra)

filtered_points = muestra_gdf[muestra_gdf['id'] == 1]

variables_vector = [
    'ugs', # UGS str
    'gmf', # SUBUNIDAD float->int
    'cobert', # NOMENCLAT float->int
    'usos', # USO_ACT float->int
    'geo', # NOMENCLAT str
]

rasterized_layers = {}
gdf_dict = {}
for var in variables_vector:
    gdf = gpd.read_file(f"{ruta_insumos}{var}.shp")

    joined = gpd.sjoin(filtered_points, gdf, how="inner", op="within")
    
    density = joined.groupby(joined.index_right).size() / gdf.area
    
    max_density = density.max()
    gdf['weight'] = density / max_density

    gdf['weight'] = gdf['weight'].fillna(0)

    rasterized = pomca.rasterize(gdf, 'weight')
    rasterized_layers[var] = rasterized
    gdf_dict[var] = gdf

variables_raster = [
    'dem', 
    'pend', 
    'pends', 
    'rugos',
    'curvar',
    'orient', 
    'insol', 
    'acuenca', 
    'long',
    'distvias', 
    'distdren', 
    'densdren', 
    'distfalla',
    'densfract', 
    'espesor', 
]

rasters = {}
for var in variables_raster:
    rasters[var] = rioxarray.open_rasterio(f"{ruta_insumos}{var}.tif")

rasters.update(rasterized_layers)
 
muestra_gdf = gpd.read_file(ruta_muestra)
muestra_gdf = muestra_gdf.explode(index_parts=False).reset_index(drop=True)

sampled_data = []
for index, point in muestra_gdf.iterrows():
    point_data = {'id': point['id']}
    for var, raster in rasters.items():
        coords = [(point['geometry'].x, point['geometry'].y)]
        sample_value = raster.sel(x=coords[0][0], y=coords[0][1], method="nearest").values
        point_data[var] = sample_value[0] if sample_value.size > 0 else np.nan
    sampled_data.append(point_data)

df = pd.DataFrame(sampled_data)

df.where(df!=-3.4028234663852886e+38, np.nan, inplace=True)
df.where(df!=4294967295, np.nan, inplace=True)
df.where(df!=-999999.0, np.nan, inplace=True)
df.where(df!=-99999.0, np.nan, inplace=True)
df.where(df!=-9999.0, np.nan, inplace=True)

df['long'] = df['long'][df['long']>0]
df['densfract'] = df['densfract'][df['densfract']>0]

for column in df.columns:
    if column not in ['id']:
        col_min = df[column].min()
        col_max = df[column].max()
        col_mean = df[column].mean()
        df[column] = (df[column] - col_min) / (col_max - col_min)
        non_zero_min = df[column][df[column] != 0].min()
        df.loc[df[column]<non_zero_min, column] = non_zero_min
        df[column] = df[column].fillna(non_zero_min)

# for column in df.columns:
#     if column not in ['id']:
#         col_min = df[column].min()
#         col_max = df[column].max()
#         col_mean = df[column].mean()
#         print(f"Columna: {column}, Min: {col_min}, Max: {col_max}, Media: {col_mean}")


#%%

import numpy as np
import xarray as xr
import rioxarray
import matplotlib.pyplot as plt

# Define la función para limpiar y normalizar los datos del raster
def clean_and_normalize_raster(raster, column_name):
    raster = raster.where(raster != -3.4028234663852886e+38, np.nan)
    raster = raster.where(raster != 4294967295, np.nan)
    raster = raster.where(raster != -999999.0, np.nan)
    raster = raster.where(raster != -99999.0, np.nan)
    raster = raster.where(raster != -9999.0, np.nan)

    # Condiciones específicas para 'long' y 'densfract'
    if column_name == 'long':
        raster = raster.where(raster > 0, np.nan)
    elif column_name == 'densfract':
        raster = raster.where(raster > 0, np.nan)
    
    # Normalización entre 0 y 1
    col_min = raster.min().values
    col_max = raster.max().values
    raster = (raster - col_min) / (col_max - col_min)
    
    non_zero_min = raster.where(raster != 0).min().values
    raster = raster.where(raster >= non_zero_min, non_zero_min)
    raster = raster.fillna(non_zero_min)
    
    return raster

# Transformar y exportar todos los rasters
for var in list(rasters.keys())[18:]:
    print(var)
    # Limpiar y normalizar el raster
    cleaned_raster = clean_and_normalize_raster(rasters[var], var)
    
    # Asegurarse de que las dimensiones están en el orden esperado antes de guardar
    cleaned_raster = cleaned_raster.transpose('band', 'y', 'x')
    
    # Guardar el raster
    cleaned_raster.rio.to_raster(f"{ruta_insumos}{var}_norm.tif")






#%%


pd.options.display.float_format = '{:,.2f}'.format


#%%
# Histogramas
folder = '/home/fdgmz/Documents/san_bartolo/Susceptibilidad/Figuras'
for col in df.iloc[:,1:].columns:

    data = df[col]
    sns.displot(data, kde=True, bins=40)
    plt.xlabel(''); plt.ylabel('Conteo'); plt.suptitle(col)
    plt.savefig(f'{folder}/{col}.png')

#%%Test de normalidad

import pandas as pd
from scipy import stats
from scipy.stats import skew, boxcox
import seaborn as sns
import matplotlib.pyplot as plt

normalidad = []

for col in df.iloc[:, 1:].columns:
    data = df[col].dropna()  # Asegúrate de eliminar los NaN antes de las pruebas

    args = stats.norm.fit(data)
    pvalue = stats.kstest(data, 'norm', args=args).pvalue

    if pvalue > 0.05:
        print(f'{col} es normal: p-value={pvalue:.2f}')
    else:
        print(f'{col} no es normal: p-value={pvalue:.2f}')

    if skew(data) > 0.1:
        sesgo = 'positivo'
        transformada = 'lognormal'
        print(f'{col} tiene sesgo positivo. Se transforma a través de lognormal')
        data_transformed, _ = boxcox(data)
        nuevo_p = stats.kstest(data_transformed, 'norm').pvalue

        if nuevo_p > 0.05:
            print(f'{col} transformada es normal: p-value={nuevo_p:.2f}')
        else:
            print(f'{col} transformada no es normal: p-value={nuevo_p:.2f}')
        
        sns.displot(data, kde=True, bins=40)
        plt.title(f'Datos no transformados {col}')
        plt.show()
        
        sns.displot(data_transformed, kde=True, bins=40)
        plt.title(f'Datos transformados (lognormal) {col}')
        plt.show()

    elif skew(data) < -0.1:
        sesgo = 'negativo'
        transformada = 'cubo'
        print(f'{col} tiene sesgo negativo. Se transforma a través del cubo')
        data_transformed = data ** 3
        nuevo_p = stats.kstest(data_transformed, 'norm').pvalue
        
        if nuevo_p > 0.05:
            print(f'{col} transformada es normal: p-value={nuevo_p:.2f}')
        else:
            print(f'{col} transformada no es normal: p-value={nuevo_p:.2f}')
        
        sns.displot(data, kde=True, bins=40)
        plt.title(f'Datos no transformados {col}')
        plt.show()
        
        sns.displot(data_transformed, kde=True, bins=40)
        plt.title(f'Datos transformados (cubo) {col}')
        plt.show()
    else:
        sesgo = 'normal'
        transformada = np.nan
        nuevo_p = np.nan
    
    normalidad.append([col, pvalue, sesgo, transformada, nuevo_p])

normalidad_df = pd.DataFrame(normalidad, columns=['Variable', 'p-value', 'Sesgo', 'Transformada', 'nuevo p-value'])
print(normalidad_df)


#%%
#Matriz de correlación

corr = df.drop(columns=['id']).corr()

#Quitar diagonal superior
mask = np.zeros_like(corr, dtype=bool)
mask[np.triu_indices_from(mask)] = True
corr[mask] = np.nan

plt.figure(figsize=(8,8))
plt.tick_params(
    which='both',
    left=False,
    bottom=False,
    top=False,
    labelbottom=True)

ax = sns.heatmap(
    corr, 
    vmin = -1, vmax = 1, center = 0,
    cmap = sns.diverging_palette(20, 220, n=200, sep=180),
    square = True)

ax.set_xticklabels(
    ax.get_xticklabels(),
    rotation=90)

plt.savefig(f'{folder}/correlacion.png')
# print('Dan con alta correlación positiva UGS-Espesor y DensFract-Geo\nDan con alta correlación negativa Insol-Pend')
#%%
#T-Test

contraste = []
for col in df.iloc[:,1:].columns:
    media_in = np.mean(df[col][df['id']==1])
    std_in = np.std(df[col][df['id']==1])
    media_es = np.mean(df[col][df['id']==0])
    std_es = np.std(df[col][df['id']==0])
    print(f'{col}\nInestables: media:{media_in}, std:{std_in}')
    print(f'{col}\nEstables: media:{media_es}, std:{std_es}')
    pvalue = ttest_ind(df[col][df['id']==1], df[col][df['id']==0]).pvalue
    if pvalue > 0.05: print(f'Los grupos de {col} son iguales: p-value={pvalue:.2f}')
    else: print(f'Los grupos de {col} son diferentes: p-value={pvalue:.2f}')
    contraste.append([col,media_es,std_es,media_in,std_in,pvalue])

contraste = pd.DataFrame(contraste, columns=['Variable','Media ES','Std ES','Media IN','Std IN','p-value'])
contraste.style.apply(lambda x: ["background: red" if v > 0.05 else "" for v in x], subset=['p-value'])

#%%
#Análisis discriminante
def LDA(df):

    print(df.columns)

    X_train, X_test, y_train, y_test = train_test_split(df.iloc[:,1:], df.iloc[:,0:1], test_size=0.2, random_state=1)

    model = LinearDiscriminantAnalysis()
    cv = RepeatedStratifiedKFold(n_splits=10, n_repeats=3, random_state=1)

    grid = dict()
    grid['solver'] = ['svd', 'lsqr', 'eigen']
    search = GridSearchCV(model, grid, scoring='accuracy', cv=cv, n_jobs=-1)
    results = search.fit(X_train, y_train)
    print('Mean Accuracy: %.3f' % results.best_score_)

    model = LinearDiscriminantAnalysis(solver=results.best_params_['solver']).fit(X_train, y_train)
    cm = confusion_matrix(y_train, model.predict(X_train))

    group_names = ['Verdaderos\nNegativos\n','Falsos\nPositivos\n',
                'Falsos\nNegativos\n','Verdaderos\nPositivos\n']
    group_counts = ["{0:0.0f}".format(value) for value in
                    cm.flatten()]
    group_percentages = ["{0:.1%}".format(value) for value in
                        cm.flatten()/np.sum(cm)]
    labels = [f"{v1}\n{v2}\n{v3}" for v1, v2, v3 in
            zip(group_names,group_counts,group_percentages)]
    labels = np.asarray(labels).reshape(2,2)

    ax = sns.heatmap(cm/np.sum(cm), annot=labels, cmap='Blues', fmt='')

    ax.set_title('Matriz de confusión\n');
    ax.set_xlabel('\nPredecido')
    ax.set_ylabel('Observado');

    ## Ticket labels - List must be in alphabetical order
    ax.xaxis.set_ticklabels(['Falso','Verdadero'])
    ax.yaxis.set_ticklabels(['Falso','Verdadero'])

    ## Display the visualization of the Confusion Matrix.
    plt.show()

    acctrain = accuracy_score(y_train, model.predict(X_train))
    cktrain = cohen_kappa_score(y_train, model.predict(X_train))
    print(f'Entrenamiento. Accuracy: {acctrain:.2f}. Kohen-Cappa: {cktrain:.2f}.')

    acctest = accuracy_score(y_test, model.predict(X_test))
    cktest = cohen_kappa_score(y_test, model.predict(X_test))
    print(f'Validación. Accuracy: {acctest:.2f}. Kohen-Cappa: {cktest:.2f}.')

    return model

#%% C1
#Todas las variables
model = LDA(df)

# #%% C2
# #Sacando INSOL por correlación con pendiente
# df_ = df.drop(columns=['INSOL'])
# model = LDA(df_)

# #%% C3
# #Sacando INSOL por correlación con pendiente y ESPESOR Y DENSFRACT por alta correlación y ser derivadas de UGS y GEO
# df__ = df.drop(columns=['INSOL','ESPESOR','DENSFRACT'])
# model = LDA(df__)

# #%% C4
# #Variables con medias estadísticamente diferentes y y ESPESOR Y DENSFRACT por alta correlación y ser derivadas de UGS y GEO
# df___ = df.drop(columns=['DEM','RUGOS','CURVAR','ORIENT','INSOL',
#                     'ACUENCA','LONG','DISTDREN','ESPESOR','DENSFRACT'])
# model = LDA(df___)

# #%%


# #%%
# #Encuentra todas las combinaciones posibles y analiza LDA
# def LDA(df):
#     X_train, X_test, y_train, y_test = train_test_split(df.iloc[:,1:], df.iloc[:,0:1], test_size=0.2, random_state=1)
#     model = LinearDiscriminantAnalysis(solver='svd').fit(X_train, y_train)
#     acctrain = accuracy_score(y_train, model.predict(X_train))
#     acctest = accuracy_score(y_test, model.predict(X_test))
    
#     return model, acctrain, acctest

# from itertools import combinations

# sample_list = df.columns[1:]
# list_combinations = list()
# for n in range(2, len(sample_list) + 1):
#     list_combinations += list(combinations(sample_list, n))

# # combinaciones = []
# # for i in list_combinations:
# #     df_ = df[['id']+list(i)]
# #     model, acctrain, acctest = LDA(df_)
# #     combinaciones.append([i, acctrain, acctest])
# #     progreso = len(combinaciones)*100/len(list_combinations)
# #     if progreso % 5 == 0:
# #         print(f'Progreso: {progreso}%')

# #%%
# combs = pd.DataFrame(combinaciones, columns=['Variables', 'acctrain', 'acctest'])


# #%%
# combs[combs['acctrain']==combs['acctrain'].max()]
# #%%
# combs[combs['acctest']==combs['acctest'].max()]
# #%%
# combs[['acctrain','acctest']].mean(axis=1).sort_values()

# #%% C9
# model = LDA(# %%
# df[['id']+list(combs.loc[378366,'Variables'])].drop(columns=['ORIENT','INSOL','DISTDREN']))

#%%
# Modelo elegido y construcción de la ecuación

eq = ''
for feat,coef in zip(model.feature_names_in_,model.coef_.flatten()):
    if coef < 0: signo = ''
    else: signo = '+'
    eq += f'{signo}{coef:.2f}*"{feat}_norm@1" '
eq

