#%%

import os
import rasterio
import numpy as np
import matplotlib.pyplot as plt

# Define paths
base_path = '/Users/federicogomez/Downloads/Agriban/Hidrología_hidráulica/Depth_Vel_Tr'
output_path = os.path.join(base_path, 'Amenaza')

# Ensure output directory exists
os.makedirs(output_path, exist_ok=True)

# Define return periods and corresponding files
return_periods = {
    'Tr2V1': 2, 'Tr5V1': 5, 'Tr10V1': 10, 'Tr25V1': 25, 'Tr50V1': 50, 'Tr100V1': 100,
}

# Define classification thresholds
thresholds = {
    'low': {'velocity': 0.5, 'depth': 0.4},
    'medium': {'velocity': 2, 'depth': 0.8, 'product': 0.5}
}

def classify_risk(velocity, depth):
    product = velocity * depth
    if (velocity <= thresholds['low']['velocity'] and depth <= thresholds['low']['depth']):
        return 1  # Low risk
    elif (velocity < thresholds['medium']['velocity'] and depth < thresholds['medium']['depth'] and product <= thresholds['medium']['product']):
        return 2  # Medium risk
    else:
        return 3  # High risk

#%%
# Procesar cada periodo de retorno
for period, years in return_periods.items():
    # Leer los rasters de profundidad y velocidad
    depth_file = os.path.join(base_path, period, 'Depth (Max).Terrain.MDT_2024_30cm_Predio_modelo.tif')
    velocity_file = os.path.join(base_path, period, 'Velocity (Max).Terrain.MDT_2024_30cm_Predio_modelo.tif')
    
    with rasterio.open(depth_file) as depth_src, rasterio.open(velocity_file) as velocity_src:
        depth = depth_src.read(1)
        velocity = velocity_src.read(1)
        
        # Obtener las máscaras de NoData
        depth_mask = depth == depth_src.nodata
        velocity_mask = velocity == velocity_src.nodata
        combined_mask = depth_mask | velocity_mask
        
        # Aplicar la clasificación
        risk_map = np.vectorize(classify_risk)(velocity, depth)
        
        # Aplicar la máscara de NoData
        risk_map[combined_mask] = depth_src.nodata

        # Guardar el raster clasificado
        profile = depth_src.profile
        profile.update(dtype=risk_map.dtype, compress='lzw')
        output_file = os.path.join(output_path, f'{period}_risk_map.tif')
        with rasterio.open(output_file, 'w', **profile) as dst:
            dst.write(risk_map, 1)

#%%

# Initialize the weighted sum array
weighted_sum_depth = np.zeros_like(rasterio.open(os.path.join(base_path, 'Tr2V1', 'Depth (Max).Terrain.MDT_2024_30cm_Predio_modelo.tif')).read(1))
weighted_sum_velocity = np.zeros_like(rasterio.open(os.path.join(base_path, 'Tr2V1', 'Velocity (Max).Terrain.MDT_2024_30cm_Predio_modelo.tif')).read(1))
total_weight = 0

# Define classification thresholds
thresholds = {
    'low': {'velocity': 0.5, 'depth': 0.4},
    'medium': {'velocity': 2, 'depth': 0.8, 'product': 0.5}
}

def classify_risk(velocity, depth):
    product = velocity * depth
    if (velocity <= thresholds['low']['velocity'] and depth <= thresholds['low']['depth']):
        return 1  # Low risk
    elif (velocity < thresholds['medium']['velocity'] and depth < thresholds['medium']['depth'] and product <= thresholds['medium']['product']):
        return 2  # Medium risk
    else:
        return 3  # High risk

weighted_sum_depth = np.zeros_like(rasterio.open(os.path.join(base_path, 'Tr2V1', 'Depth (Max).Terrain.MDT_2024_30cm_Predio_modelo.tif')).read(1))
weighted_sum_velocity = np.zeros_like(rasterio.open(os.path.join(base_path, 'Tr2V1', 'Velocity (Max).Terrain.MDT_2024_30cm_Predio_modelo.tif')).read(1))
total_weight = 0
nodata_mask = np.ones_like(weighted_sum_depth, dtype=bool)

# Process each return period
for period, years in return_periods.items():
    # Read the depth and velocity rasters
    depth_file = os.path.join(base_path, period, 'Depth (Max).Terrain.MDT_2024_30cm_Predio_modelo.tif')
    velocity_file = os.path.join(base_path, period, 'Velocity (Max).Terrain.MDT_2024_30cm_Predio_modelo.tif')
    with rasterio.open(depth_file) as depth_src, rasterio.open(velocity_file) as velocity_src:
        depth = depth_src.read(1)
        profile = depth_src.profile
        velocity = velocity_src.read(1)
        
        # Update the weighted sum arrays, ignoring NoData values
        valid_mask = (depth_src.read_masks(1) > 0) & (velocity_src.read_masks(1) > 0)
        nodata_mask &= ~valid_mask
        weighted_sum_depth[valid_mask] += depth[valid_mask] * years
        weighted_sum_velocity[valid_mask] += velocity[valid_mask] * years
        total_weight += years

# Calculate the final weighted average rasters
final_depth = np.where(nodata_mask, np.nan, weighted_sum_depth / total_weight)
final_velocity = np.where(nodata_mask, np.nan, weighted_sum_velocity / total_weight)

# Classify the final risk map
final_risk = np.zeros_like(final_depth, dtype=np.uint8)
for row in range(final_depth.shape[0]):
    for col in range(final_depth.shape[1]):
        if np.isfinite(final_depth[row, col]) and np.isfinite(final_velocity[row, col]):
            final_risk[row, col] = classify_risk(final_velocity[row, col], final_depth[row, col])

# Set NoData values in the final risk map
final_risk[nodata_mask] = 0

# Save the final risk map
final_risk_file = os.path.join(output_path, 'final_risk_map.tif')
with rasterio.open(final_risk_file, 'w', **profile) as dst:
    dst.write(final_risk, 1)

#%%

from scipy import ndimage

# Apply median filter to smooth small isolated pixels
final_risk = ndimage.median_filter(final_risk, size=3)

# Filter by minimum mapping unit
label_image, num_labels = ndimage.label(final_risk > 0)
pixel_area = 0.3 * 0.3  # Resolution of 0.3 m
min_mapping_area = 4  # Minimum mapping area of 1.5 m2
areas = ndimage.sum(final_risk > 0, label_image, range(num_labels + 1)) * pixel_area
small_object_mask = areas < min_mapping_area
final_risk[small_object_mask[label_image]] = 0

# Save the final risk map
final_risk_file = os.path.join(output_path, 'final_risk_map_amc.tif')
with rasterio.open(final_risk_file, 'w', **depth_src.profile) as dst:
    dst.write(final_risk, 1)

# %%

import numpy as np
import matplotlib.pyplot as plt

# Crear un grid de valores de velocidad y profundidad
velocity_values = np.linspace(0, 3, 300)
depth_values = np.linspace(0, 1.2, 300)
Velocity, Depth = np.meshgrid(velocity_values, depth_values)

# Clasificar el riesgo según los umbrales dados
def classify_risk(velocity, depth):
    product = velocity * depth
    if velocity <= 0.5 and depth <= 0.4:
        return 1  # Bajo riesgo
    elif velocity < 2 and depth < 0.8 and product <= 0.5:
        return 2  # Riesgo medio
    else:
        return 3  # Alto riesgo

# Aplicar la clasificación al grid
Risk = np.vectorize(classify_risk)(Velocity, Depth)

# Crear la gráfica con áreas sombreadas según el riesgo
plt.figure(figsize=(8, 6))
cmap = plt.cm.get_cmap('RdYlGn_r', 3)  # Colores: Verde (bajo), Amarillo (medio), Rojo (alto)
plt.contourf(Velocity, Depth, Risk, levels=[0.5, 1.5, 2.5, 3.5], colors=['green', 'yellow', 'red'], alpha=0.7)

# Añadir etiquetas y formato
cbar = plt.colorbar(ticks=[1,2,3], aspect=50)
cbar.ax.set_yticklabels(['Baja', 'Media', 'Alta'])
plt.title('Clasificación de Amenaza por Inundación', fontsize=14)
plt.xlabel('Velocidad (m/s)', fontsize=12)
plt.ylabel('Profundidad (m)', fontsize=12)
plt.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()
plt.show()
# %%
