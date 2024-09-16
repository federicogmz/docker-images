# %%

import geopandas as gpd
import xarray as xr
from rasterstats import zonal_stats
import numpy as np

cuencas_path = "/Users/federicogomez/Documents/Work/2024_pomca_sanbartolo/Susceptibilidad_AvT/Cuencas/Cuencas.shp"
dtm_path = "/Users/federicogomez/Documents/Work/2024_pomca_sanbartolo/Insumos/MDT/MDT_V3.tif"

cuencas = gpd.read_file(cuencas_path)
dtm = xr.open_dataarray(dtm_path, mask_and_scale=True)[0].rio.reproject(cuencas.crs)

cuencas["Area (A)"] = cuencas.geometry.area / 1e6

stats = zonal_stats(
    cuencas, dtm.values, affine=dtm.rio.transform(), stats=["min", "max"]
)

cuencas["Relief (H)"] = [(stat["max"] - stat["min"]) / 1000 for stat in stats]


def process_perimeter_points(cuencas, dtm, num_points=200):
    def generate_perimeter_points(polygon, num_points):
        perimeter = polygon.boundary
        distances = np.linspace(0, perimeter.length, num_points)
        points = [perimeter.interpolate(distance) for distance in distances]
        return gpd.GeoDataFrame(geometry=gpd.GeoSeries(points), crs=cuencas.crs)

    def extract_elevation(points_gdf, dtm):
        dtm_rio = dtm.rio.write_crs(cuencas.crs)
        points_gdf["elevation"] = [
            dtm_rio.sel(x=point.x, y=point.y, method="nearest").item()
            for point in points_gdf.geometry
        ]
        return points_gdf

    def compute_distances_to_min_elevation(points_gdf):
        min_elevation_idx = points_gdf["elevation"].idxmin()
        min_elevation_point = points_gdf.loc[min_elevation_idx, "geometry"]
        points_gdf["distance_to_min_elevation"] = points_gdf.geometry.distance(
            min_elevation_point
        )
        return points_gdf

    for idx, geom in cuencas.iterrows():
        points_gdf = generate_perimeter_points(geom.geometry, num_points)
        points_with_elevation = extract_elevation(points_gdf, dtm)
        points_with_distances = compute_distances_to_min_elevation(
            points_with_elevation
        )
        max_distance_id = points_with_distances["distance_to_min_elevation"].idxmax()
        max_distance = points_with_distances.loc[
            max_distance_id, "distance_to_min_elevation"
        ]
        cuencas.at[idx, "Lb (km)"] = max_distance / 1000

    return cuencas


cuencas = process_perimeter_points(cuencas, dtm)

cuencas["Releif Ratio (Rh)"] = cuencas["Relief (H)"] / cuencas["Lb (km)"]
cuencas["Melton (M)"] = cuencas["Relief (H)"] / np.sqrt(cuencas["Area (A)"])

cuencas.to_file(cuencas_path, engine='fiona')

# %%

import numpy as np
import pandas as pd
import seaborn as sns
import geopandas as gpd
from sklearn.svm import SVC
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler


def train_and_plot_models(df, combinations):
    df['FFR'] = (df['FFR'] == 'T').astype(int)
    X = df[['A (sq. km)', 'Su', 'M', 'Rh']]
    y = df['FFR']
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    models = {}
    fig, axes = plt.subplots(1, len(combinations), figsize=(18, 6))
    fig.suptitle('Clasificación Mejorada para Diferentes Combinaciones de Variables', fontsize=16, fontweight='bold')
    sns.set(style="whitegrid")
    
    for i, ((x_col, y_col), ax) in enumerate(zip(combinations, axes)):
        X_plot = X_scaled[:, [X.columns.get_loc(x_col), X.columns.get_loc(y_col)]]
        model = SVC(kernel='rbf', C=10, gamma=0.1)
        model.fit(X_plot, y)
        models[i] = model
        
        xx, yy = np.meshgrid(np.linspace(X_plot[:, 0].min()-1, X_plot[:, 0].max()+1, 100),
                             np.linspace(X_plot[:, 1].min()-1, X_plot[:, 1].max()+1, 100))
        Z = model.predict(np.c_[xx.ravel(), yy.ravel()])
        Z = Z.reshape(xx.shape)
        
        ax.contourf(xx, yy, Z, alpha=0.3, cmap='coolwarm')
        ax.scatter(X_plot[y==0, 0], X_plot[y==0, 1], c='blue', label='No torrencial', edgecolor='k', s=50)
        ax.scatter(X_plot[y==1, 0], X_plot[y==1, 1], c='red', label='Torrencial', edgecolor='k', s=50)
        ax.set_xlabel(x_col)
        ax.set_ylabel(y_col)
        ax.legend()
    
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.show()
    return models, scaler

# Use the function
df = pd.read_csv("db/torrencialidad.csv")
combinations = [('A (sq. km)', 'Rh'), ('Su', 'Rh'), ('M', 'Rh')]
models, scaler = train_and_plot_models(df, combinations)

def classify_and_plot_new_data(df, combinations, models, scaler):
    X = df[['A (sq. km)', 'Su', 'M', 'Rh']]
    X_scaled = scaler.transform(X)
    
    fig, axes = plt.subplots(1, len(combinations), figsize=(18, 6))
    # fig.suptitle('Clasificación Mejorada para Diferentes Combinaciones de Variables', fontsize=16, fontweight='bold')
    sns.set(style="whitegrid")
    
    for i, ((x_col, y_col), ax) in enumerate(zip(combinations, axes)):
        X_plot = X_scaled[:, [X.columns.get_loc(x_col), X.columns.get_loc(y_col)]]
        
        xx, yy = np.meshgrid(np.linspace(X_plot[:, 0].min()-1, X_plot[:, 0].max()+1, 100),
                             np.linspace(X_plot[:, 1].min()-1, X_plot[:, 1].max()+1, 100))
        Z = models[i].predict(np.c_[xx.ravel(), yy.ravel()])
        Z = Z.reshape(xx.shape)
        
        ax.contourf(xx, yy, Z, alpha=0.3, cmap='coolwarm')
        
        # Predict classes for new data
        y_pred = models[i].predict(X_plot)
        ax.scatter(X_plot[y_pred==0, 0], X_plot[y_pred==0, 1], c='blue', label='No torrencial', edgecolor='k', s=50)
        ax.scatter(X_plot[y_pred==1, 0], X_plot[y_pred==1, 1], c='red', label='Torrencial', edgecolor='k', s=50)
        
        ax.set_xlabel(x_col)
        ax.set_ylabel(y_col)
        ax.legend()
        
        # Add classification to the dataframe
        df[f'Clasificacion_{x_col}_{y_col}'] = y_pred
    
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.show()
    
    return df

cuencas_path = "/Users/federicogomez/Documents/Work/2024_pomca_sanbartolo/Susceptibilidad_AvT/Cuencas/Cuencas.shp"
cuencas = gpd.read_file(cuencas_path)
cuencas = cuencas.rename(
    columns={
        'Area (A)': 'A (sq. km)',
        'Releif Rat': 'Rh',
        'Melton (M)': 'M',
    }
)

cuencas_clasificadas = classify_and_plot_new_data(cuencas, combinations, models, scaler)

# cuencas_clasificadas.to_file('/Users/federicogomez/Documents/Work/2024_pomca_sanbartolo/Susceptibilidad_AvT/Cuencas/Cuencas_clasificadas.shp', engine='fiona')