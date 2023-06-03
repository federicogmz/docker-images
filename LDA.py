#%%
import pandas as pd
from scipy import stats
import numpy as np
from scipy.stats import boxcox, skew, ttest_ind
import matplotlib.pyplot as plt
import ast
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.model_selection import RepeatedStratifiedKFold, GridSearchCV, train_test_split
from sklearn.metrics import confusion_matrix, cohen_kappa_score, accuracy_score

pd.options.display.float_format = '{:,.2f}'.format

#%%
#Leer archivos
df = pd.read_csv('/Volumes/GoogleDrive-111862222221066300789/Shared drives/Proyectos(Fede)/POMCA/Susceptibilidad/Insumos/df.csv')
df = df.fillna(0)

#%%
#Renombrar columnas
df.columns = ['id', 'DEM', 'PEND', 'PENDS', 'RUGOS',
       'CURVAR', 'ORIENT', 'INSOL', 'ACUENCA', 'LONG',
       'DISTVIAS', 'DISTDREN', 'DENSDREN', 'DISTFALLA', 'UGS', 'GMF',
       'COBERT', 'ESPESOR', 'USOS', 'GEO', 'DENSFRACT']
df = df[['id', 'DEM', 'PEND', 'PENDS', 'RUGOS',
       'CURVAR', 'ORIENT', 'INSOL', 'ACUENCA', 'LONG',
       'DISTVIAS', 'DISTDREN', 'DENSDREN', 'DISTFALLA', 'DENSFRACT', 
       'UGS', 'GMF', 'COBERT', 'ESPESOR', 'USOS', 'GEO']]

#%%
#Normalizar variables

df.iloc[:, 1:] = StandardScaler().fit_transform(df.iloc[:, 1:])
df.iloc[:, 1:] = df.iloc[:, 1:] + np.abs(df.iloc[:, 1:].min()) + 1

#%%
#Histogramas
folder = '/Volumes/GoogleDrive-111862222221066300789/Shared drives/Proyectos(Fede)/POMCA/Susceptibilidad/Figuras/Hist'
for col in df.iloc[:,1:].columns:

    data = df[col]
    sns.displot(data, kde=True, bins=40)
    plt.xlabel(''); plt.ylabel('Conteo'); plt.suptitle(col)
    plt.savefig(f'{folder}/{col}.png')

#%%Test de normalidad

normalidad = []

for col in df.iloc[:,1:].columns:
    
    data = df[col]
   
    pvalue = stats.kstest(data, 'norm').pvalue

    # if pvalue > 0.05: print(f'{col} es normal: p-value={pvalue:.2f}')
    # else: print(f'{col} no es normal: p-value={pvalue:.2f}')
    if skew(data) > 0.1:
        sesgo = 'positivo'
        transformada = 'lognormal'
        #print(f'{col} tiene sesgo positivo.Se transforma a través de lognormal')
        data = boxcox(data, 0)
        nuevo_p = stats.kstest(data, 'norm').pvalue
        #if pvalue > 0.01: print(f'{col} transformada es normal: p-value={pvalue:.2f}')
        #else: print(f'{col} transformada no es normal: p-value={pvalue:.2f}')
        #sns.displot(data, kde=True, bins=40)
        #plt.title=f'Datos no transformados {col}'
    elif skew(data) < 0.1:
        sesgo = 'negativo'
        transformada = 'cubo'
        #print(f'{col} tiene sesgo negativo.')
        #print('Se transforma a través del cubo')
        data = data ** 3
        nuevo_p = stats.kstest(data, 'norm').pvalue
        #if pvalue > 0.01: print(f'{col} transformada es normal: p-value={pvalue:.2f}')
        #else: print(f'{col} transformada no es normal: p-value={pvalue:.2f}')
        #sns.displot(data, kde=True, bins=40)
        #plt.title=f'Datos no transformados {col}'
    #plt.show()
    else: sesgo = 'normal'
    
    normalidad.append([col,pvalue,sesgo,transformada,nuevo_p])

normalidad = pd.DataFrame(normalidad, columns=['Variable','p-value','Sesgo','Transformada','p-value'])
normalidad

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
    square = True);

ax.set_xticklabels(
    ax.get_xticklabels(),
    rotation=90);

print('Dan con alta correlación positiva UGS-Espesor y DensFract-Geo\nDan con alta correlación negativa Insol-Pend')
#%%
#T-Test

contraste = []
for col in df.iloc[:,1:].columns:
    media_in = np.mean(df[col][df['id']==1])
    std_in = np.std(df[col][df['id']==1])
    media_es = np.mean(df[col][df['id']==0])
    std_es = np.std(df[col][df['id']==0])
    #print(f'{col}\nInestables: media:{media_in}, std:{std_in}')
    #print(f'{col}\nEstables: media:{media_es}, std:{std_es}')
    pvalue = ttest_ind(df[col][df['id']==1], df[col][df['id']==0]).pvalue
    #if pvalue > 0.05: print(f'Los grupos de {col} son iguales: p-value={pvalue:.2f}')
    #else: print(f'Los grupos de {col} son diferentes: p-value={pvalue:.2f}')
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

#%% C2
#Sacando INSOL por correlación con pendiente
df_ = df.drop(columns=['INSOL'])
model = LDA(df_)

#%% C3
#Sacando INSOL por correlación con pendiente y ESPESOR Y DENSFRACT por alta correlación y ser derivadas de UGS y GEO
df__ = df.drop(columns=['INSOL','ESPESOR','DENSFRACT'])
model = LDA(df__)

#%% C4
#Variables con medias estadísticamente diferentes y y ESPESOR Y DENSFRACT por alta correlación y ser derivadas de UGS y GEO
df___ = df.drop(columns=['DEM','RUGOS','CURVAR','ORIENT','INSOL',
                    'ACUENCA','LONG','DISTDREN','ESPESOR','DENSFRACT'])
model = LDA(df___)

#%%


#%%
#Encuentra todas las combinaciones posibles y analiza LDA
def LDA(df):
    X_train, X_test, y_train, y_test = train_test_split(df.iloc[:,1:], df.iloc[:,0:1], test_size=0.2, random_state=1)
    model = LinearDiscriminantAnalysis(solver='svd').fit(X_train, y_train)
    acctrain = accuracy_score(y_train, model.predict(X_train))
    acctest = accuracy_score(y_test, model.predict(X_test))
    
    return model, acctrain, acctest

from itertools import combinations

sample_list = df.columns[1:]
list_combinations = list()
for n in range(2, len(sample_list) + 1):
    list_combinations += list(combinations(sample_list, n))

combinaciones = []
for i in list_combinations:
    df_ = df[['id']+list(i)]
    model, acctrain, acctest = LDA(df_)
    combinaciones.append([i, acctrain, acctest])
    progreso = len(combinaciones)*100/len(list_combinations)
    if progreso % 5 == 0:
        print(f'Progreso: {progreso}%')

#%%
combs = pd.DataFrame(combinaciones, columns=['Variables', 'acctrain', 'acctest'])
combs.to_pickle('FuncDiscr.pkl')

#%%
#Mejores combinaciones
combs = pd.read_pickle('FuncDiscr.pkl')

#%%
combs[combs['acctrain']==combs['acctrain'].max()]
#%%
combs[combs['acctest']==combs['acctest'].max()]
#%%
combs[['acctrain','acctest']].mean(axis=1).sort_values()

#%% C5
model = LDA(df[['id']+list(combs.loc[1001190,'Variables'])])
#%% C6
model = LDA(df[['id']+list(combs.loc[568622,'Variables'])])
#%% C7
model = LDA(df[['id']+list(combs.loc[378366,'Variables'])])
#%% C8
model = LDA(df[['id']+list(combs.loc[991939,'Variables'])])

#%% C9
model = LDA(# %%
df[['id']+list(combs.loc[378366,'Variables'])].drop(columns=['ORIENT','INSOL','DISTDREN']))

#%%
# Modelo elegido y construcción de la ecuación

eq = ''
for feat,coef in zip(model.feature_names_in_,model.coef_.flatten()):
    if coef < 0: signo = ''
    else: signo = '+'
    eq += f'{signo}{coef:.2f}*{feat} '
eq

#%%
#Construcción de mapa de susceptibilidad
import rasterio as rio
import rasterio.mask
import numpy as np
import fiona
from osgeo import gdal
#%%

ruta = '/Volumes/GoogleDrive-111862222221066300789/Shared drives/Proyectos(Fede)/POMCA/Susceptibilidad/Insumos/'

#%% Cortar todos los raster con el límite de cuenca
for raster in ['Slope.tif','Aspecto.tif','Insolacion.tif','DistDren.tif',
                'DensDren.tif','DistFalla.tif','gmf.tif','coberturas.tif','usos.tif']:
    gdal.Warp(srcDSOrSrcDSTab=ruta + raster,
                destNameOrDestDS='/Volumes/GoogleDrive-111862222221066300789/Shared drives/Proyectos(Fede)/POMCA/Susceptibilidad/Final/'+raster,
                cutlineDSName=ruta +'Limite_cuenca_real.shp',
                cropToCutline=True)
#%%
ruta = '/Volumes/GoogleDrive-111862222221066300789/Shared drives/Proyectos(Fede)/POMCA/Susceptibilidad/Final/'

def read(path):
    raster = rasterio.open(ruta+path).read(1)
    raster = np.where(raster==np.nan,0,raster)
    raster = np.where(raster==raster[0,0],0,raster)
    raster = (raster - np.nanmean(raster)) / np.nanstd(raster)
    raster = raster + np.abs(np.nanmin(raster)) + 1
    return raster

#%%
pend = read('Slope.tif')
orient = read('Aspecto.tif')
insol = read('Insolacion.tif')
distdren = read('DistDren.tif')
densdren = read('DensDren.tif')
distfalla = read('DistFalla.tif')
gmf = read('gmf.tif')
cobert = read('coberturas.tif')
usos = read('usos.tif')

#%%

suscept = 0.29*pend +0.06*orient +0.25*insol -0.19*distdren +0.04*densdren -0.22*distfalla +0.39*gmf -0.16*cobert +0.67*usos
suscept = np.where(suscept==suscept[0,0],np.nan,suscept)

meta = rasterio.open(ruta+'Slope.tif').meta.copy()
meta.update(compress='lzw', nodata=0)

#Exports raster output file and assign metadata
with rasterio.open(ruta+'Suscept.tif', 'w+', **meta) as out:
    out.write_band(1, suscept)
