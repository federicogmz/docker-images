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

pd.options.display.float_format = '{:,.2f}'.format


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


ruta_insumos = '/Users/federicogomez/Documents/Work/2024_pomca_sanbartolo/Susceptibilidad/Insumos/'

ruta_muestra = '/Users/federicogomez/Documents/Work/2024_pomca_sanbartolo/Inventario/Muestra.shp'
ruta_dem = ruta_insumos + 'dem.tif'
ruta_ugs = '/Users/federicogomez/Documents/Work/2024_pomca_sanbartolo/Susceptibilidad/Ugs_preliminar/Ugs_preliminar.shp'

ugs = gpd.read_file(ruta_ugs)

#%%

pomca = POMCA(ruta_dem, ugs)


muestra_gdf = gpd.read_file(ruta_muestra)

filtered_points = muestra_gdf[muestra_gdf['id'] == 1]

variables_vector = [
    'Ugs', # UGS str
    'gmf', # SUBUNIDAD float->int
    'cobert', # NOMENCLAT float->int
    'usos', # USO_ACT float->int
    'geo', # NOMENCLAT str
]

rasterized_layers = {}
gdf_dict = {}
for var in variables_vector:
    gdf = gpd.read_file(f"{ruta_insumos}{var}.shp")

    joined = gpd.sjoin(filtered_points, gdf, how="inner", predicate="within")
    
    density = joined.groupby(joined.index_right).size() / gdf.area
    
    max_density = density.max()
    gdf['weight'] = density / max_density

    gdf['weight'] = gdf['weight'].fillna(0)
    gdf.to_file(f"{ruta_insumos}{var}_codificada.shp")

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


#%%

import numpy as np
import xarray as xr
import rioxarray
import matplotlib.pyplot as plt

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
for var in list(rasters.keys()):
    print(var)
    cleaned_raster = clean_and_normalize_raster(rasters[var], var)
    cleaned_raster = cleaned_raster.transpose('band', 'y', 'x')
    cleaned_raster.rio.to_raster(f"{ruta_insumos}{var}_norm.tif")

#%%
# Histogramas
folder = '/home/fdgmz/Documents/san_bartolo/Susceptibilidad/Figuras'
for col in df.iloc[:,1:].columns:

    data = df[col]
    sns.displot(data, kde=True, bins=40)
    plt.xlabel(''); plt.ylabel('Conteo'); plt.suptitle(col)
    # plt.savefig(f'{folder}/{col}.png')

#%%Test de normalidad

import pandas as pd
from scipy import stats
from scipy.stats import skew, boxcox
import seaborn as sns
import matplotlib.pyplot as plt

normalidad = []

for col in df.iloc[:, 1:].columns:
    data = df[col].dropna()

    args = stats.norm.fit(data)
    pvalue = stats.kstest(data, 'norm', args=args).pvalue

    # if pvalue > 0.05:
    #     print(f'{col} es normal: p-value={pvalue:.2f}')
    # else:
    #     print(f'{col} no es normal: p-value={pvalue:.2f}')

    if skew(data) > 0.1:
        sesgo = 'positivo'
        transformada = 'lognormal'
        # print(f'{col} tiene sesgo positivo. Se transforma a través de lognormal')
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
        # print(f'{col} tiene sesgo negativo. Se transforma a través del cubo')
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

# plt.savefig(f'{folder}/correlacion.png')
#%%
#T-Test

import pandas as pd
import numpy as np
from scipy.stats import ttest_ind, f_oneway

contraste = []
for col in df.iloc[:, 1:].columns:
    media_in = np.mean(df[col][df['id'] == 1])
    std_in = np.std(df[col][df['id'] == 1])
    media_es = np.mean(df[col][df['id'] == 0])
    std_es = np.std(df[col][df['id'] == 0])
    ttest_pvalue = ttest_ind(df[col][df['id'] == 1], df[col][df['id'] == 0]).pvalue
    
    f_val, anova_pvalue = f_oneway(df[col][df['id'] == 1], df[col][df['id'] == 0])
    
    contraste.append([
        col.upper(), 
        media_es, std_es, 
        media_in, std_in,
        f_val, anova_pvalue,
        ttest_pvalue
    ])

contraste_df = pd.DataFrame(contraste, columns=[
    'VARIABLE', 'MEDIA ES', 'STD',
    'MEDIA IN', 'STD',
    'F', 'Sig.', 'TTest'
])

contraste_df

#%%
import pandas as pd
import numpy as np
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, cohen_kappa_score, f1_score
from sklearn.feature_selection import RFE
from scipy.stats import f
from itertools import product
import warnings
warnings.filterwarnings("ignore")

#Análisis discriminante
from statsmodels.stats.outliers_influence import variance_inflation_factor
from statsmodels.tools.tools import add_constant


def calculate_vif(df):
    df_with_const = add_constant(df)
    vif = pd.DataFrame()
    vif['Variable'] = df_with_const.columns
    vif['VIF'] = [variance_inflation_factor(df_with_const.values, i) for i in range(df_with_const.shape[1])]
    vif = vif[vif['Variable']!='const']
    vif = vif.sort_values('VIF', ascending=False)
    return vif

def iterative_vif(df, threshold=5):
    while True:
        vif = calculate_vif(df)
        if (vif['VIF'] <= threshold).all():
            break
        max_vif_var = vif.sort_values('VIF', ascending=False).iloc[0]['Variable']
        df = df.drop(columns=[max_vif_var])
        print(f'Removiendo {max_vif_var}')
    return df

def wilks_lambda_test(lda, X, y):
    wilks_lambda = np.prod(1 - lda.explained_variance_ratio_)
    n, p = X.shape
    df1 = p * (len(np.unique(y)) - 1)
    df2 = n - 1 - (p - df1) / 2
    f_stat = ((1 - wilks_lambda) / wilks_lambda) * (df2 / df1)
    p_value = f.sf(f_stat, df1, df2)
    return wilks_lambda, p_value

def find_best_lda_model(df, metric):
    X = df.drop('id', axis=1)
    y = df['id']

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)
    results = []

    for num_features in range(1, X_train.shape[1] + 1):
        rfe = RFE(LinearDiscriminantAnalysis(solver='svd'), n_features_to_select=num_features)
        X_train_rfe = rfe.fit_transform(X_train, y_train)
        X_test_rfe = rfe.transform(X_test)
        
        groups = [X_train_rfe[y_train == label] for label in np.unique(y_train)]
        cov_matrices = [np.cov(group, rowvar=False) for group in groups if group.shape[0] > 1]
        
        if len(cov_matrices) < 2:
            continue

        lda = LinearDiscriminantAnalysis(solver='svd')
        lda.fit(X_train_rfe, y_train)
        
        wilks_lambda, p_value = wilks_lambda_test(lda, X_train_rfe, y_train)

        y_train_pred = lda.predict(X_train_rfe)
        y_test_pred = lda.predict(X_test_rfe)

        train_accuracy = accuracy_score(y_train, y_train_pred)
        test_accuracy = accuracy_score(y_test, y_test_pred)
        kappa = cohen_kappa_score(y_test, y_test_pred)
        f1 = f1_score(y_test, y_test_pred)
        y_all = np.concatenate([y_train, y_test])
        y_all_pred = np.concatenate([y_train_pred, y_test_pred])

        overall_accuracy = accuracy_score(y_all, y_all_pred)

        results.append({
            'num_features': num_features,
            'train_accuracy': train_accuracy,
            'test_accuracy': test_accuracy,
            'overall_accuracy': overall_accuracy,
            'kappa': kappa,
            'f1_score': f1,
            'features': X.columns[rfe.support_].tolist(),
            'wilks_lambda': wilks_lambda,
            'p_value': p_value,
            'model': lda
        })

    sorted_results = sorted(results, key=lambda x: -x[metric])
    return sorted_results

results = []
sorted_results = find_best_lda_model(df, metric='f1_score')
results.extend(sorted_results[:1])

df_new = df.drop(columns=['pends', 'densfract'])
df_new = iterative_vif(df_new)
sorted_results = find_best_lda_model(df_new, metric='f1_score')

results.extend(sorted_results[:4])

resultados = pd.DataFrame(results)

#%%

X = df_new.drop('id', axis=1)
y = df_new['id']

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)
rfe = RFE(LinearDiscriminantAnalysis(solver='svd'), n_features_to_select=7)
X_train_rfe = rfe.fit_transform(X_train, y_train)
X_test_rfe = rfe.transform(X_test)

lda = LinearDiscriminantAnalysis(solver='svd')
lda.fit(X_train_rfe, y_train)

import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay, roc_curve, auc

y_pred = lda.predict(X_test_rfe)
y_proba = lda.predict_proba(X_test_rfe)[:, 1]

fpr, tpr, thresholds = roc_curve(y_test, y_proba)
roc_auc = auc(fpr, tpr)

plt.figure(figsize=(8,8))
plt.plot(fpr, tpr, color='dodgerblue', lw=2, label=f'Área bajo la curva: {roc_auc*100:.1f}')
plt.plot([0, 1], [0, 1], color='grey', linestyle='--')
plt.xlabel('Tasa de Falsos Positivos (FPR)')
plt.ylabel('Tasa de Verdaderos Positivos (TPR)')
plt.title('Curva ROC')
plt.legend(loc='lower right')
plt.show()

cm = confusion_matrix(y_test, y_pred)
cm_normalized = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis] * 100

plt.figure(figsize=(8,8))
disp = ConfusionMatrixDisplay(confusion_matrix=cm_normalized, display_labels=lda.classes_)
disp.plot(cmap=plt.cm.Blues)
plt.title('Matriz de Confusión (%)')
plt.xlabel('Etiqueta Predicha')
plt.ylabel('Etiqueta Verdadera')
plt.show()



#%%
# Modelo elegido y construcción de la ecuación

model = resultados.at[1, 'model']
features = resultados.at[1, 'features']

eq = ''
for feat, coef in zip(features, model.coef_.flatten()):
    if coef < 0:
        signo = ''
    else:
        signo = '+' if eq else ''  # Evita un signo '+' al inicio de la ecuación
    eq += f'{signo}{coef:.2f}*"{feat}@1" '

eq



#%%
import matplotlib.pyplot as plt
import numpy as np
import rioxarray
import geopandas as gpd
import seaborn as sns

ruta_muestra = '/Users/federicogomez/Documents/Work/2024_pomca_sanbartolo/Inventario/Muestra.shp'
susceptibilidad = '/Users/federicogomez/Documents/Work/2024_pomca_sanbartolo/Susceptibilidad/Suscept.tif'
raster = rioxarray.open_rasterio(susceptibilidad)
raster = raster.where(raster!=-3.4028234663852886e+38)

col_min = raster.min().values
col_max = raster.max().values
raster = (raster - col_min) / (col_max - col_min)

# raster.rio.to_raster('/Users/federicogomez/Documents/Work/2024_pomca_sanbartolo/Susceptibilidad/Suscept_norm.tif')

data = raster.values.flatten()
data = data[~np.isnan(data)]

#%%
muestra_gdf = gpd.read_file(ruta_muestra)
muestra_gdf = muestra_gdf.explode(index_parts=False).reset_index(drop=True)

muestra_gdf['suscept_value'] = muestra_gdf.geometry.apply(
    lambda geom: raster.sel(x=geom.x, y=geom.y, method='nearest').values[0]
)

muestra_inestable = muestra_gdf[muestra_gdf['id'] == 1]
muestra_estable = muestra_gdf[muestra_gdf['id'] == 0]

k = 15

p1 = np.percentile(muestra_inestable['suscept_value'].dropna(), k)
p2 = np.percentile(muestra_estable['suscept_value'].dropna(), 100-k)

print(f'Valor de p1: {p1:.3f}')
print(f'Valor de p2: {p2:.3f}')
