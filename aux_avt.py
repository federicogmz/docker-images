import pandas as pd

notebook = r'/mnt/c/Users/Fede/OneDrive - Universidad Nacional de Colombia/Work/2023/Supia/Workingfolder/QmaxSupia.xlsx'
df = pd.read_excel(notebook, sheet_name='Sheet1')
# Create and export tables
folder_path = "/mnt/c/Users/Fede/OneDrive - Universidad Nacional de Colombia/Work/2023/Supia/AvT/Top"

for index, row in df.iterrows():
    q1 = row['Q1']
    v1 = 10
    q3 = row['Q3']
    v3 = 10
    
    # Create the table
    table_data = {'T': range(0, 310, 10),
                  'Q1': [q1] * 31,
                  'V1': [v1] * 31,
                  'Q2': [0] * 31,
                  'V2': [0] * 31,
                  'Q3': [q3] * 31,
                  'V3': [v3] * 31}
    
    table_df = pd.DataFrame(table_data)
    
    # Export the table to CSV
    file_name = f'{folder_path}/{int(row["Nombre"])}.csv'
    table_df.to_csv(file_name, sep='\t', index=False, float_format='%.2f')

#%%
# Script to generate individual shell scripts

from pyproj import Proj, transform

# Define the source and target coordinate reference systems
source_crs = Proj(init='epsg:9377')
target_crs = Proj(init='epsg:32618')

# Input coordinates
hydro_coords = {
    1: '429073.479,606606.742',
    2: '433831.921,603075.079',
    3: '424045.941,600956.618',
    4: '424673.753,607478.887',
    5: '424710.547,600171.898',
    6: '425286.645,599711.642',
    7: '426223.945,599564.848',
    8: '431039.614,601342.015',
    9: '429906.525,600251.532',
    10: '433489.455,600875.457',
    11: '427176.659,606063.293',
    12: '426393.405,604961.478',
    13: '431134.199,602953.697',
    14: '429909.536,603112.014',
    15: '432673.08,598436.945',
    16: '432228.66,596836.593',
    17: '428536.161,606667.402',
    18: '425286.255,611065.805',
    19: '423802.249,606626.601',
    20: '433349.362,601686.933'
}

# # Convert coordinates
# converted_coords = {}
# for key, value in hydro_coords.items():
#     x, y = map(float, value.split(','))
#     converted_x, converted_y = transform(source_crs, target_crs, x, y)
#     converted_coords[key] = f'{converted_x:.8f},{converted_y:.8f}'

# # Print the converted coordinates
# for key, value in converted_coords.items():
#     print(f'{key}: {value}')

template = """\
r.in.gdal input=DTM.tif output=dem

g.region -s raster=dem zoom=dem
g.region -p

r.avaflow prefix=m{prefijo} phases=m elevation=dem controls=0,1,0,1,1,0 friction=30,10,0,0,0,0,0.035 basal=-7,0.01,0,0,0.01 hydrograph={csv_file} hydrocoords={coords},20,-9999 time=300,6000
"""

folder_path = "/mnt/c/Users/Fede/OneDrive - Universidad Nacional de Colombia/Work/2023/Supia/AvT"

for numero, coords in hydro_coords.items():
    script_content = template.format(prefijo=int(numero), csv_file=f'{int(numero)}.csv', coords=coords)
    script_filename = f'{folder_path}/{int(numero)}.sh'
    
    with open(script_filename, 'w') as script_file:
        script_file.write(script_content)

# %%
import os
import re
from osgeo import gdal, osr
import numpy as np

# Define the folder path
folder_path = "/mnt/c/Users/Fede/OneDrive - Universidad Nacional de Colombia/Work/2023/Supia/AvT/Top/resultados"

# Function to find and process the ASCII files
def process_ascii_files(parent_folder):
    # Use regex to extract the m and following numbers
    match = re.search(r'/m(\d+)_results$', parent_folder)

    if match:
        m_number = match.group(1)
        ascii_folder = os.path.join(parent_folder, f"m{m_number}_ascii")
    else:
        print("Pattern not found in the path.")

    ascii_files = [f for f in os.listdir(ascii_folder) if f.endswith('.asc') and f.startswith(f"m{m_number}_hflow_max")]
    
    if not ascii_files:
        print(f"No files found in {ascii_folder}")
        return None
    else:
        print(f'Archivos encontrados: {len(ascii_files)}')

    all_data = []
    for ascii_file in ascii_files:
        file_path = os.path.join(ascii_folder, ascii_file)
        raster = gdal.Open(file_path)
        band = raster.GetRasterBand(1)
        data = band.ReadAsArray()
        all_data.append(data)
    
    # Filter out values lower than 0.1 and set them to NaN
    max_raster = np.maximum.reduce(all_data)
    max_raster[max_raster < 0.1] = np.nan
    
    # Read vfront raster
    vfront_file = os.path.join(ascii_folder, f"m{m_number}_vhmax.asc")
    vfront_raster = gdal.Open(vfront_file)
    vfront_band = vfront_raster.GetRasterBand(1)
    vfront_data = vfront_band.ReadAsArray()
    
    # Compute intensity raster
    intensity_raster = max_raster * vfront_data ** 2
    
    return intensity_raster, raster.GetGeoTransform(), raster.GetProjection()

# Iterate through parent folders (m1_results, m2_results, ..., m20_results)
for i in [9,11]:
    parent_folder = os.path.join(folder_path, f"m{i}_results")
    
    # Process ASCII files and compute intensity raster
    intensity_raster, geotransform, projection = process_ascii_files(parent_folder)
    
    if intensity_raster is not None:
        # Set the CRS to EPSG 32618
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(32618)
        projection_wkt = srs.ExportToWkt()
        
        # Write the intensity raster as a GeoTIFF file
        output_file = os.path.join(folder_path, f"Idf{i}.tif")
        driver = gdal.GetDriverByName("GTiff")
        output_raster = driver.Create(output_file, intensity_raster.shape[1], intensity_raster.shape[0], 1, gdal.GDT_Float32)
        
        # Set geotransform
        output_raster.SetGeoTransform(geotransform)
        
        # Set projection
        output_raster.SetProjection(projection_wkt)  # Use the WKT string for EPSG 32618
        
        # Write array data
        output_raster.GetRasterBand(1).WriteArray(intensity_raster)
        
        # Flush the cache and close the dataset
        output_raster.FlushCache()
        output_raster = None
