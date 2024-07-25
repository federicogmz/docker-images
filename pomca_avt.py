# %%

import geopandas as gpd
import xarray as xr
from rasterstats import zonal_stats
import numpy as np

cuencas_path = "/home/fdgmz/Documents/san_bartolo/Avt/Cuencas.shp"
dtm_path = "/home/fdgmz/Documents/san_bartolo/Avt/MDT_V3.tif"

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

# %%

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# Database de cuencas torrenciales y no torrenciales
df = pd.read_csv("Torrencialidad/torrencialidad.csv")
df_cuencas = cuencas

# Set the seaborn style and context
sns.set(context="poster", style="whitegrid")

# Create subplots with 1 row and 3 columns
fig, axes = plt.subplots(1, 3, figsize=(20, 7), sharey=True)

# Plot 1 - A vs Rh
sns.scatterplot(
    x="A (sq. km)",
    y="Rh",
    hue="FFR",
    style="FFR",
    markers={"nT": "o", "T": "X"},
    data=df[df["FFR"].isin(["nT", "T"])],
    s=60,
    palette=["red", "blue"],
    alpha=0.8,
    ax=axes[0],
    legend=False,
)
sns.scatterplot(
    x="Area (A)",
    y="Releif Ratio (Rh)",
    data=df_cuencas,
    s=100,
    color="darkgreen",
    alpha=0.5,
    ax=axes[0],
    markers=".",
)
axes[0].grid(alpha=0.4)

# Plot 2 - Su vs Rh
sns.scatterplot(
    x="Su",
    y="Rh",
    hue="FFR",
    style="FFR",
    markers={"nT": "o", "T": "X"},
    data=df[df["FFR"].isin(["nT", "T"])],
    s=60,
    palette=["red", "blue"],
    alpha=0.8,
    ax=axes[1],
    legend=False,
)
sns.scatterplot(
    x="Su",
    y="Releif Ratio (Rh)",
    data=df_cuencas,
    s=100,
    color="darkgreen",
    alpha=0.5,
    ax=axes[1],
    markers=".",
)
axes[1].grid(alpha=0.4)

# Plot 3 - M vs Rh
sns.scatterplot(
    x="M",
    y="Rh",
    hue="FFR",
    style="FFR",
    markers={"nT": "o", "T": "X"},
    data=df[df["FFR"].isin(["nT", "T"])],
    s=60,
    palette=["red", "blue"],
    alpha=0.8,
    ax=axes[2],
    legend=False,
)
sns.scatterplot(
    x="Melton (M)",
    y="Releif Ratio (Rh)",
    data=df_cuencas,
    s=100,
    color="darkgreen",
    alpha=0.5,
    ax=axes[2],
    markers=".",
)

# Create empty scatter plot artists for legend items
legend_markers = [
    Line2D(
        [0],
        [0],
        marker="o",
        color="w",
        markerfacecolor="blue",
        markersize=15,
        alpha=0.8,
    ),
    Line2D(
        [0], [0], marker="X", color="w", markerfacecolor="red", markersize=15, alpha=0.8
    ),
    Line2D(
        [0],
        [0],
        marker=".",
        color="w",
        markerfacecolor="green",
        markersize=25,
        alpha=0.8,
    ),
]

# Add custom legend to the third subplot
axes[2].legend(
    legend_markers, ["No torrencial", "Torrencial", "Cuencas estudio"], fontsize=15
)

axes[2].grid(alpha=0.4)
plt.tight_layout()
plt.show()

# %%
