import os
import re
from typing import Any
import whitebox
import rasterio
import rioxarray
import subprocess
import numpy as np
import xarray as xr
import richdem as rd
from osgeo import gdal
import geopandas as gpd
from rasterio.plot import show
import matplotlib.pyplot as plt
from geocube.api.core import make_geocube

np.seterr(divide='ignore', invalid='ignore')

class geohazards():

    def __init__(self):
        
        return
    
    def preprocessDEM(self, demPath):

        return
    
    def resample(input, resolution, export=False, output=''):

        inputRaster = rioxarray.open_rasterio(input)

        # Resample the dataset to the new resolution
        rasterResample = inputRaster.rio.reproject(dst_crs=inputRaster.rio.crs, resolution=resolution)[0]
        
        if export:
            rasterResample.rio.to_raster(output)

        return rasterResample

    def exportASCII(self, data, output, fmt='%f'):

        if len(data.shape)==3:
            data = data.squeeze()

        # Convert DataArray to NumPy array
        data_array = data.values

        # Get the shape of the array
        rows, cols = data.sizes['y'], data.sizes['x']

        # Write the ASCII file
        with open(output, 'w') as f:
            # Write the header information
            f.write('ncols {}\n'.format(cols))
            f.write('nrows {}\n'.format(rows))
            f.write('xllcorner {}\n'.format(data.x.values[0]))
            f.write('yllcorner {}\n'.format(data.y.values[0]))
            f.write('cellsize {}\n'.format(float(data.coords['x'].diff('x').values[0])))
            f.write('NODATA_value -9999\n')

            data_array = np.where(np.isnan(data_array), -9999, data_array)

            # Write the data with custom formatting
            np.savetxt(f, data_array, fmt=fmt, delimiter=' ')

        print(f'Raster exported succesfully to {output}.')

    def Catani(self, demPath, geoPath, hmin=0.1, hmax=3.0):

        """
        Function to get soil thickness from model S Catani et. al (2010)
        """

        #Imports dem and opens it with rasterio
        dem = xr.open_dataarray(demPath, mask_and_scale=True)
        dem_ = rd.LoadGDAL(demPath)
        slope = rd.TerrainAttribute(dem_, attrib='slope_radians')
        slope_rad = xr.DataArray(slope, coords=[dem.coords['y'],dem.coords['x']])

        gdf = gpd.read_file(geoPath)

        if isinstance(hmin, str):
            hmin = make_geocube(vector_data=gdf, measurements = [hmin], 
                            like = dem, fill = np.nan)[hmin]
        else: hmin = hmin

        if isinstance(hmax, str):
            hmax = make_geocube(vector_data=gdf, measurements = [hmax], 
                            like = dem, fill = np.nan)[hmax]
        else: 
            hmax = hmax

        #Calculates variables with Tangent
        tan_slope = np.tan(slope_rad)
        tan_slope_max = np.tan(np.nanmax(slope_rad))
        tan_slope_min = np.tan(np.nanmin(slope_rad))

        #Calculates soil thickness
        catani = hmax * (1 - ( (tan_slope - tan_slope_min) / (tan_slope_max - tan_slope_min) ) * (1 - (hmin / hmax)) )

        catani = xr.DataArray(catani, coords=[dem.coords['y'], dem.coords['x']])

        catani = catani.where(~np.isnan(dem), np.nan).squeeze()

        #Returns soil thickness rasterio file
        return catani
    
    def Slope(self, demPath):

        dem = xr.open_dataarray(demPath, mask_and_scale=True)

        dem_ = rd.LoadGDAL(demPath)
        slope = rd.TerrainAttribute(dem_, attrib='slope_degrees')

        slope = xr.DataArray(slope, coords=[dem.coords['y'],dem.coords['x']])

        return slope

    def flowDirection(self, demPath):

        # Initialize the WhiteboxTools
        wbt = whitebox.WhiteboxTools()

        # Run the FlowPointerEpsilon tool to compute flow direction
        wbt.d8_pointer(
            demPath,
            f'{os.path.dirname(demPath)}/flowdir.tif',
            esri_pntr=True,
            callback=self.non_verbose_callback
        )

        flowdir = xr.open_dataarray(f'{os.path.dirname(demPath)}/flowdir.tif')
        flowdir = flowdir.where(flowdir!=0, 1)
        os.remove(f'{os.path.dirname(demPath)}/flowdir.tif')
        
        return flowdir

    def non_verbose_callback(self, progress):
        """ Function to make WhiteBoxTools non verbose"""
        pass


class TRIGRS(geohazards):

    def __init__(self):

        geohazards.__init__(self)

    def __call__(self, demPath, geoPath, outPath):

        self.dem, self.slope, self.zonas, self.zs, self.flowdir = self.Insumos(demPath, geoPath, outPath)
        
        error = self.GridMatch(self.rows, self.cols, outPath)

        if not error:
            n, imax, rows, cols, nwf = self.TopoIndex(self.rows, self.cols, outPath)

        # if n==0:
        #     self.TRIGRS(outPath, imax, rows, cols, nwf)

    def Insumos(self, demPath, geoPath, outPath):

        # DEM
        dem = xr.open_dataarray(demPath, mask_and_scale=True)

        print('Performing hydrological corrections on DEM: fill and breach depressions.')
        
        # Initialize the WhiteboxTools
        wbt = whitebox.WhiteboxTools()
        wbt.fill_depressions(
            demPath,
            f'{os.path.dirname(demPath)}/filldem.tif',
            fix_flats=True,
            flat_increment=0.1,
            callback=self.non_verbose_callback
        )
        wbt.breach_depressions(
            f'{os.path.dirname(demPath)}/filldem.tif',
            f'{os.path.dirname(demPath)}/filldem.tif',
            callback=self.non_verbose_callback
        )
        dem = xr.open_dataarray(f'{os.path.dirname(demPath)}/filldem.tif')
        os.remove(f'{os.path.dirname(demPath)}/filldem.tif')

        # Pendiente
        slope = self.Slope(demPath)

        # Zonas
        gdf = gpd.read_file(geoPath)
        gdf.geometry = gdf.geometry.map(lambda x: x.buffer(1).union(x))
        gdf = gdf.explode().reset_index(drop=True)

        zonas = make_geocube(vector_data=gdf, measurements = ['Zona'], 
                        like=dem, fill = np.nan)['Zona']

        # Set the corresponding pixels in zs to NaN
        zonas = zonas.where(~np.isnan(dem))

        # Espesor
        zs = self.Catani(demPath, geoPath)

        # Dirección de flujo
        flowdir = self.flowDirection(demPath)

        dem = dem.where(~np.isnan(zonas))
        zs = zs.where(~np.isnan(zonas))
        flowdir = flowdir.where(~np.isnan(zonas))
        slope = slope.where(~np.isnan(zonas))

        if dem.isnull().sum().values == zs.isnull().sum().values == flowdir.isnull().sum().values == slope.isnull().sum().values == zonas.isnull().sum().values:
            print('\nNaN match in all rasters.\n')
            print('----- Exporting ASCII. -----')

            self.exportASCII(dem, f'{outPath}/dem.asc')
            self.exportASCII(slope, f'{outPath}/slope.asc')
            self.exportASCII(zonas, f'{outPath}/zonas.asc', fmt='%d')
            self.exportASCII(zs, f'{outPath}/zs.asc')
            self.exportASCII(flowdir, f'{outPath}/flowdir.asc', fmt='%d')

            print('----- ASCII exported succesfully. -----\n')

        self.rows, self.cols = dem.sizes['y'], dem.sizes['x']
        self.gdf = gdf

        return dem, slope, zonas, zs, flowdir

    def GridMatch(self, rows, cols, outPath):

        # Create the content string with variables
        content = 'number of grid files to test\n'
        content += '5\n'
        content += 'rows,columns\n'
        content += f'{rows},{cols}\n'
        content += 'name of input files (text string, up to 255 characters per line; one file per line, list master grid first)\n'
        content += 'dem.asc\n'
        content += 'zonas.asc\n'
        content += 'slope.asc\n'
        content += 'flowdir.asc\n'
        content += 'zs.asc\n'
        content += '*** Note, Flow-direction grids need additional processing beyond the capabilities of GridMatch.'

        # Write the content to the file
        with open(f'{outPath}/gm_in.txt', 'w') as file:
            file.write(content)

        print("GridMatch input created successfully.")
        print(f"----- Executing gridmatch.exe. -----")

        os.chdir(outPath)
        subprocess.run([f'{outPath}/gridmatch.exe'])

        print(f"----- gridmatch.exe finished. -----\n")

        # Results
        mismatches = {}

        with open(f'{outPath}/GridMatchLog.txt', 'r') as file:
            lines = file.readlines()

        # Iterate over the grid results section and extract the grid names and number of mismatches
        for i, line in enumerate(lines):
            if ' Results for grid' in line:
                grid_name = lines[i+1].strip()
            if ' Number of mismatches found:' in line:
                mismatches_line = lines[i]
                mismatches_count = int(mismatches_line.split(':')[1].strip())
                mismatches[grid_name] = mismatches_count

        # Print the mismatches count for each grid
        for grid_name, count in mismatches.items():
            print(f"Mismatches found for {grid_name}: {count}")

        error = 1

        if all(value == 0 for value in mismatches.values()):
            print('GridMatch did not found mismatches. You can procede.\n')
            error = 0
            
        return error

    def TopoIndex(self, rows, cols, outPath):

        iterations = 500

        # Create the content string with variables
        content = 'Name of project (up to 255 characters)\n'
        content += 'project\n'
        content += 'Rows, Columns, flow-direction numbering scheme (ESRI=1, TopoIndex=2)\n'
        content += f'{rows},{cols}, 1\n'
        content += 'Exponent, Number of iterations\n'
        content += f'-1, {iterations}\n'
        content += 'Name of elevation grid file\n'
        content += 'dem.asc\n'
        content += 'Name of direction grid\n'
        content += 'flowdir.asc\n'
        content += 'Save listing of D8 downslope neighbor cells? Enter T (.true.) or F (.false.)\n'
        content += 'T\n'
        content += 'Save grid of D8 downslope neighbor cells? Enter T (.true.) or F (.false.)\n'
        content += 'T\n'
        content += 'Save cell index number grid? Enter T (.true.) or F (.false.)\n'
        content += 'T\n'
        content += 'Save list of cell number and corresponding index number? Enter T (.true.) or F (.false.)\n'
        content += 'T\n'
        content += 'Save flow-direction grid remapped from ESRI to TopoIndex? Enter T (.true.) or F (.false.)\n'
        content += 'T\n'
        content += 'Name of folder to store output?\n'
        content += 'tpx\\\n'
        content += 'ID code for output files? (8 characters or less)\n'
        content += 'project\n'

        # Write the content to the file
        with open(f'{outPath}/tpx_in.txt', 'w') as file:
            file.write(content)

        print("TopoIndex input created successfully.\n")
        print(f"----- Executing TopoIndex.exe with {iterations} iterations. -----\n")

        os.chdir(outPath)
        if not os.path.exists(f'{outPath}/tpx'):
            os.makedirs(f'{outPath}/tpx')
        subprocess.run([f'{outPath}/TopoIndex.exe'])

        print(f"\n----- TopoIndex.exe finished. -----\n")

        # Results

        with open(f'{outPath}/TopoIndexLog.txt', 'r') as file:
            text = file.read()

            success = 'TopoIndex finished normally' in text
            error = 'Corrections did not converge' in text

            if success:
                print("TopoIndex finished normally. You can procede.\n")
                imax, rows, cols, nwf = re.findall(r'\b\d+\b', re.search(r'Data cells, Rows, Columns, Downslope cells\n(.*?)\n', text, re.DOTALL).group(1))
                n = 0
            if error:
                last_line = text.splitlines()[-1].strip()
                matches = re.findall(r'\d+', text.splitlines()[-2:-1][0])
                n = int(matches[-1])
                print(f"TopoIndex did not finished succesfully.\n{last_line}\n")
                print(f"Cells that do not converge: {n}\n")
                imax = rows = cols = nwf = 0
            
        return n, imax, rows, cols, nwf

    def TRIGRS(self, outPath, imax, rows, cols, nwf):

        # Define the variable values
        project_name = 'project'
        tx = 100
        nmax = -10
        nzs = 10
        mmax = 200
        nper = 1
        zmin = 0.001
        uww = 9.8e3
        t = 28800
        zones = self.gdf['Zona'].max()
        zmax = depth = -3
        rizero = 2.59e-07
        min_slope_angle = 0

        cohesion = [13000, 9000, 16000, 13000, 16000, 16000]
        phi = [31, 35, 31, 31, 31, 31]
        uws = [19800, 19800, 19100, 19800, 19100, 19100]
        diffus = [2.68e-06, 5.00e-07, 2.47e-04, 2.68e-06, 2.47e-04, 2.47e-04]
        k_sat = [2.68e-08, 5.00e-09, 2.47e-06, 2.68e-08, 2.47e-06, 2.47e-06]
        theta_sat = [0.50, 0.45, 0.40, 0.51, 0.40, 0.40]
        theta_res = [0.28, 0.22, 0.13, 0.27, 0.13, 0.13]
        alpha = [4.54, 3.84, 3.88, 2.53, 3.88, 3.88]

        cri = 6.28611e-06
        capt = [0, 28800]

        slope_file = f'{outPath}/slope.asc'
        zone_file = f'{outPath}/zonas.asc'
        depth_file = f'{outPath}/zs.asc'
        init_depth_file = f'{outPath}/zs.asc'
        infil_rate_file = "none"
        rainfall_files = "none"

        runoff_receptor_file = "tpx\\TIdscelGrid_Charras.txt"
        runoff_order_file = "tpx\\TIcelindxList_Charras.txt"
        runoff_cell_list_file = "tpx\\TIdscelList_Charras.txt"
        runoff_weighting_file = "tpx\\TIwfactorList_Charras.txt"
        output_folder = "Resultados\\"
        output_suffix = "P"

        save_runoff_grid = False
        save_factor_of_safety_grid = True
        save_depth_of_safety_grid = False
        save_pore_pressure_grid = False
        save_infil_rate_grid = False
        save_unsat_zone_flux_grid = False
        save_pressure_head_flag = 0
        num_output_times = 1
        output_times = 28800
        skip_other_timesteps = False
        use_analytic_solution = True
        estimate_positive_pressure_head = True
        use_psi0 = True
        log_mass_balance = True
        flow_direction = "gener"
        add_steady_background_flux = True

        # Generate the text
        content = "Name of project (up to 255 characters)\n"
        content += f"{project_name}\n"
        content += "imax, row, col, nwf, tx, nmax\n"
        content += f"{imax}, {rows}, {cols}, {nwf}, {tx}, {nmax}\n"
        content += "nzs, mmax, nper, zmin, uww, t, zones\n"
        content += f"{nzs}, {mmax}, {nper}, {zmin}, {uww}, {t}, {zones}\n"
        content += "zmax, depth, rizero, Min_Slope_Angle (degrees)\n"
        content += f"{zmax}, {depth}, {rizero}, {min_slope_angle}\n"

        for i in range(zones):
            content += f"zone,{i+1}\n"
            content += "cohesion,phi,uws,diffus,K-sat,Theta-sat,Theta-res,Alpha\n"
            content += f"{cohesion[i]},{phi[i]},{uws[i]},{diffus[i]},{k_sat[i]},{theta_sat[i]},{theta_res[i]},{alpha[i]}\n"

        content += f"cri(1), {cri}\n"
        content += "capt(1), capt(2), ..., capt(n), capt(n+1)\n"
        content += f"{','.join(map(str, capt))}\n"
        content += "File name of slope angle grid (slofil)\n"
        content += f"{slope_file}\n"
        content += "File name of property zone grid (zonfil)\n"
        content += f"{zone_file}\n"
        content += "File name of depth grid (zfil)\n"
        content += f"{depth_file}\n"
        content += "File name of initial depth of water table grid (depfil)\n"
        content += f"{init_depth_file}\n"
        content += "File name of initial infiltration rate grid (rizerofil)\n"
        content += f"{infil_rate_file}\n"
        content += "List of file name(s) of rainfall intensity for each period, (rifil())\n"
        content += f"{rainfall_files}\n"
        content += "File name of grid of D8 runoff receptor cell numbers (nxtfil)\n"
        content += f"{runoff_receptor_file}\n"
        content += "File name of list of defining runoff computation order (ndxfil)\n"
        content += f"{runoff_order_file}\n"
        content += "File name of list of all runoff receptor cells (dscfil)\n"
        content += f"{runoff_cell_list_file}\n"
        content += "File name of list of runoff weighting factors (wffil)\n"
        content += f"{runoff_weighting_file}\n"
        content += "Folder where output grid files will be stored (folder)\n"
        content += f"{output_folder}\n"
        content += "Identification code to be added to names of output files (suffix)\n"
        content += f"{output_suffix}\n"
        content += "Save grid files of runoff? Enter T (.true.) or F (.false.)\n"
        content += f"{str(save_runoff_grid)}\n"
        content += "Save grid of minimum factor of safety? Enter T (.true.) or F (.false.)\n"
        content += f"{str(save_factor_of_safety_grid)}\n"
        content += "Save grid of depth of minimum factor of safety? Enter T (.true.) or F (.false.)\n"
        content += f"{str(save_depth_of_safety_grid)}\n"
        content += "Save grid of pore pressure at depth of minimum factor of safety? Enter T (.true.) or F (.false.)\n"
        content += f"{str(save_pore_pressure_grid)}\n"
        content += "Save grid files of actual infiltration rate? Enter T (.true.) or F (.false.)\n"
        content += f"{str(save_infil_rate_grid)}\n"
        content += "Save grid files of unsaturated zone basal flux? Enter T (.true.) or F (.false.)\n"
        content += f"{str(save_unsat_zone_flux_grid)}\n"
        content += "Save listing of pressure head and factor of safety (\"flag\")? (Enter -2 detailed, -1 normal, 0 none)\n"
        content += f"{str(save_pressure_head_flag)}\n"
        content += "Number of times to save output grids\n"
        content += f"{num_output_times}\n"
        content += "Times of output grids\n"
        content += f"{','.join(map(str, output_times))}\n"
        content += "Skip other timesteps? Enter T (.true.) or F (.false.)\n"
        content += f"{str(skip_other_timesteps)}\n"
        content += "Use analytic solution for fillable porosity? Enter T (.true.) or F (.false.)\n"
        content += f"{str(use_analytic_solution)}\n"
        content += "Estimate positive pressure head in rising water table zone (i.e. in lower part of unsat zone)? Enter T (.true.) or F (.false.)\n"
        content += f"{str(estimate_positive_pressure_head)}\n"
        content += "Use psi0=-1/alpha? Enter T (.true.) or F (.false.) (False selects the default value, psi0=0)\n"
        content += f"{str(use_psi0)}\n"
        content += "Log mass balance results? Enter T (.true.) or F (.false.)\n"
        content += f"{str(log_mass_balance)}\n"
        content += "Flow direction (enter \"gener\", \"slope\", or \"hydro\")\n"
        content += f"{flow_direction}\n"
        content += "Add steady background flux to transient infiltration rate to prevent drying beyond the initial conditions during periods of zero infiltration?\n"
        content += f"{str(add_steady_background_flux)}\n"

        return

class SHALSTAB(geohazards):

    def __init__(self, demPath, geoPath, q, shalstabPath='', criticalRainPath='', exportShalstab=False, exportCriticalRain=False):

        geohazards.__init__(self)

        #Imports rasters and opens them with rasterio
        dem = rioxarray.open_rasterio(demPath, mask_and_scale=True)

        dem_ = rd.LoadGDAL(demPath)
        accum = rd.FlowAccumulation(dem_, method='D8')
        accum = xr.DataArray(accum, coords=[dem.coords['y'],dem.coords['x']])

        gdf = gpd.read_file(geoPath)

        #Rasterize from geology
        c = make_geocube(vector_data=gdf, measurements = ['C'], 
                        like=dem, fill = np.nan)['C']
        p = make_geocube(vector_data=gdf, measurements = ['Phi'], 
                        like=dem, fill = np.nan)['Phi']
        p = np.radians(p)
        g = make_geocube(vector_data=gdf, measurements = ['Gamma'], 
                        like=dem, fill = np.nan)['Gamma']
        k = make_geocube(vector_data=gdf, measurements = ['K'], 
                        like = dem, fill = np.nan)['K']

        #Calculates slope
        slope = rd.TerrainAttribute(dem_, attrib='slope_radians')
        slope = xr.DataArray(slope, coords=[dem.coords['y'],dem.coords['x']])
        slope_rad = slope.where(slope!=-9999.)

        #Calculates zs from Catani
        zs = self.Catani(demPath, geoPath)

        #Unit weight of water
        gammaw = 9.81

        #Left from unconditional conditions
        izq_un = np.tan(slope_rad)

        #Right from unstable condition
        der_unstable = np.tan(p) + (c / (g * zs * np.cos(slope_rad)**2))

        #Right from stable condition
        der_stable = (1 - (gammaw/g)) * np.tan(p) + (c / (g * zs * np.cos(slope_rad)**2))

        #Left from main equation
        izq = accum / dem.rio.resolution()[0]

        #Right from main equation
        der = ((k * (zs * np.cos(slope_rad)) * np.sin(slope_rad)) / (q/1000)) * ((g / gammaw) * (1 - np.tan(slope_rad) / np.tan(p)) + (c / (gammaw * zs * np.cos(slope_rad)**2 * np.tan(p)))) 

        #Shalstab stability categories
        shalstab = np.where(izq_un >= der_unstable, 2, # Unconditionally unstable
                    np.where(izq_un < der_stable, 1, # Unconditionally stable
                    np.where(izq > der, 3, # Unstable
                    np.where(izq <= der, 4, np.nan)))) # Stable

        shalstab = xr.DataArray(shalstab, coords=[dem.coords['y'],dem.coords['x']])
        shalstab.rio.write_nodata(-9999, inplace=True)

        if exportShalstab:
            shalstab.rio.to_raster(shalstabPath)

        #Calculates critical rainfall
        criticalRain = 1000 * (k * zs * np.cos(slope_rad) * np.sin(slope_rad)) * (dem.rio.resolution()[0] / accum) * ((g / gammaw) * (1 - (np.tan(slope_rad) / np.tan(p))) + c / (gammaw * zs * np.cos(slope_rad)**2 * np.tan(p)))
        criticalRain = criticalRain.where(criticalRain>=0)

        criticalRain = np.where(izq_un >= der_unstable, -2, 
                       np.where(izq_un < der_stable, -1, criticalRain))
        criticalRain = xr.DataArray(criticalRain, coords=[dem.coords['y'],dem.coords['x']])

        criticalRain.rio.write_nodata(-9999, inplace=True)

        if exportCriticalRain:
            criticalRain.rio.to_raster(criticalRainPath)
        
        totalCeldas = shalstab.to_series().sum()

        try: 
            incondEstables = shalstab.to_series().value_counts()[1]
        except:
            incondEstables = 0
        try:
            incondInestables = shalstab.to_series().value_counts()[2]
        except:
            incondInestables = 0
        try:
            inestables = shalstab.to_series().value_counts()[3]
        except:
            inestables = 0
        try:
            estables = shalstab.to_series().value_counts()[4]
        except:
            estables = 0

        stabilityReport = f'Incondicionalmente estables: {incondEstables*100/totalCeldas:.2f}%\n'
        stabilityReport += f'Incondicionalmente inestables: {incondInestables*100/totalCeldas:.2f}%\n'
        stabilityReport += f'Inestables: {inestables*100/totalCeldas:.2f}%\n'
        stabilityReport += f'Estables: {estables*100/totalCeldas:.2f}%'

        shalstab.attrs['Reporte'] = stabilityReport

        print(stabilityReport)

        self.shalstab = shalstab
        self.criticalRain = criticalRain

        return 
    
    def __call__(self):

        return self.shalstab, self.criticalRain

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
