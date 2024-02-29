import os
import re
import whitebox
import rasterio
import rioxarray
import statistics
import subprocess
import numpy as np
import pandas as pd
import xarray as xr
import richdem as rd
import geopandas as gpd
from scipy import ndimage
from matplotlib import colors
from pysheds.grid import Grid
from rasterio.plot import show
import matplotlib.pyplot as plt
from geocube.api.core import make_geocube

np.seterr(divide='ignore', invalid='ignore')

class geohazards():

    def __init__(self):
        
        return
    
    def preprocess_dem(self, dem_path):

        print('--- Preprocessing DTM. ---')

        with rasterio.open(dem_path) as src:
            data = src.read(1, masked=True)
            nodata_mask = data.mask

            dilated_mask = ndimage.binary_dilation(nodata_mask)

            data_filled = data.copy()
            data_filled[nodata_mask] = np.mean(data[dilated_mask])

            profile = src.profile
            with rasterio.open(dem_path, 'w', **profile) as dst:
                dst.write(data_filled, 1)

        print('--- DTM successfully proccessed. ---')

        dem = xr.open_dataarray(dem_path)

        return dem
    
    def resample(self, input, resolution, export=False, output=''):

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

    def Catani(self, dem_path, geo_path, hmin=0.1, hmax=3.0):

        """
        Function to get soil thickness from model S Catani et. al (2010)
        """

        #Imports dem and opens it with rasterio
        dem = xr.open_dataarray(dem_path, mask_and_scale=True)
        dem_ = rd.LoadGDAL(dem_path)
        slope = rd.TerrainAttribute(dem_, attrib='slope_radians')
        slope_rad = xr.DataArray(slope, coords=[dem.coords['y'],dem.coords['x']])

        gdf = gpd.read_file(geo_path)

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
    
    def Slope(self, dem_path):

        dem = xr.open_dataarray(dem_path, mask_and_scale=True)

        dem_ = rd.LoadGDAL(dem_path)
        slope = rd.TerrainAttribute(dem_, attrib='slope_degrees')

        slope = xr.DataArray(slope, coords=[dem.coords['y'],dem.coords['x']])

        return slope

    def flowdir(self, dem_path):

        dem_xr = xr.open_dataarray(dem_path)

        # Create a PySheds grid
        grid = Grid.from_raster(dem_path)
        dem = grid.read_raster(dem_path)

        # Return the preprocessed DEM array
        flow_direction = grid.flowdir(dem)
        
        # Convert to xarray DataArray
        flow_direction = xr.DataArray(flow_direction, coords=[dem_xr.coords['y'],dem_xr.coords['x']])

        flow_direction = flow_direction.where(flow_direction>0, np.nan)

        flow_direction.rio.set_crs(dem_xr.rio.crs, inplace=True)

        return flow_direction

    def non_verbose_callback(self, progress):
        """ Function to make WhiteBoxTools non verbose"""
        pass

class TRIGRS(geohazards):

    def __init__(self):

        geohazards.__init__(self)

    def __call__(self, dem_path, geo_path, out_path, hora, cri, fosm=False):

        self.out_path = out_path

        self.dem, self.slope, self.zonas, self.zs, self.fdir = self.Insumos(dem_path, geo_path)
        self.zones = self.gdf['Zona'].max()

        error = self.GridMatch(self.rows, self.cols)

        if not error:
            n, imax, rows, cols, nwf = self.TopoIndex(self.rows, self.cols)

        if n==0:
            if not os.path.exists(f'{self.out_path}/Resultados'):
                os.makedirs(f'{self.out_path}/Resultados')

            cohesion = self.gdf['Cohesion'].values
            friccion = self.gdf['Friccion'].values
            gamma = self.gdf['Gamma'].values
            ks = self.gdf['Ks'].values

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

            # Corre para valores medios
            print(f'Corriendo para valores medios.')
            self.tr_in_creation(imax, rows, cols, nwf, hora, cri, C_means, phi_means, uws_means, ks_means, output_suffix='M')
            self.TRIGRS_main()

            def read_result(file_path):
                data = np.genfromtxt(file_path, skip_header=6, delimiter=' ')
                data = np.where(data == -9999, np.nan, data)
                data = np.where(data >= 10, 10, data)
                return data
            
            dem = rioxarray.open_rasterio(dem_path, mask_and_scale=True)
            fs_mean = read_result(f'{out_path}/Resultados/TRfs_min_M_1.txt')
            fs_mean = xr.DataArray(fs_mean, coords=[dem.coords['y'],dem.coords['x']])
            fs_mean.rio.write_nodata(-9999, inplace=True)
            fs_mean.rio.write_crs(dem.rio.crs, inplace=True)
            fs_mean.rio.to_raster(f'{out_path}/Resultados/FS.tif')
            
            
            if fosm:

                print(f'Corriendo para valores modificados.')

                parameters = {'C':C_means, 'P':phi_means, 'G':uws_means, 'K':ks_means}

                for key in parameters.keys():
                    params = parameters.copy()
                    params[key] = [ float(i)*1.1 for i in params[key]]

                    self.tr_in_creation(imax, rows, cols, nwf, hora, cri, *list(params.values()), output_suffix=key)
                    self.TRIGRS_main()

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

                # Reading results
                FS_medio = read_result(f'{out_path}/Resultados/TRfs_min_M_1.txt')
                FS_cohe = read_result(f'{out_path}/Resultados/TRfs_min_C_1.txt')
                FS_phi = read_result(f'{out_path}/Resultados/TRfs_min_P_1.txt')
                FS_uws = read_result(f'{out_path}/Resultados/TRfs_min_G_1.txt')
                FS_ks = read_result(f'{out_path}/Resultados/TRfs_min_K_1.txt')

                zonas = np.genfromtxt(f'{out_path}/zonas.asc', skip_header=6, delimiter=' ')
                zonas = np.where(zonas == -9999, np.nan, zonas)

                def dFS(FS, mean_list):
                    for i in range(1, self.zones + 1):
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
                    for i in range(1, self.zones + 1):
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

                Ind_conf = xr.DataArray(Ind_conf, coords=[dem.coords['y'],dem.coords['x']])
                Ind_conf.rio.write_nodata(-9999, inplace=True)
                Ind_conf.rio.write_crs(dem.rio.crs, inplace=True)
                Ind_conf.rio.to_raster(f'{out_path}/Resultados/Ind_Conf.tif')
            
        return

    def Insumos(self, dem_path, geo_path):

        # DEM
        dem = xr.open_dataarray(dem_path)
        fdir = self.flowdir(dem_path)

        # Pendiente
        slope = self.Slope(dem_path)

        # Zonas
        gdf = gpd.read_file(geo_path)
        gdf.geometry = gdf.geometry.map(lambda x: x.buffer(1).union(x))
        gdf = gdf.explode(index_parts=True).reset_index(drop=True)

        zonas = make_geocube(vector_data=gdf, measurements = ['Zona'], like=dem, fill = np.nan)['Zona']

        # Set the corresponding pixels in zs to NaN
        zonas = zonas.where(~np.isnan(fdir))

        # Espesor
        zs = self.Catani(dem_path, geo_path)

        dem = dem.where(~np.isnan(zonas))
        zs = zs.where(~np.isnan(zonas))
        fdir = fdir.where(~np.isnan(zonas))
        slope = slope.where(~np.isnan(zonas))

        if dem.isnull().sum().values == zs.isnull().sum().values == fdir.isnull().sum().values == slope.isnull().sum().values == zonas.isnull().sum().values:
            print('\nNaN match in all rasters.\n')
            print('----- Exporting ASCII. -----')

            self.exportASCII(dem, f'{self.out_path}/dem.asc')
            self.exportASCII(slope, f'{self.out_path}/slope.asc')
            self.exportASCII(zonas, f'{self.out_path}/zonas.asc', fmt='%d')
            self.exportASCII(zs, f'{self.out_path}/zs.asc')
            self.exportASCII(fdir, f'{self.out_path}/flowdir.asc', fmt='%d')

            print('----- ASCII exported succesfully. -----\n')

            self.rows, self.cols = dem.sizes['y'], dem.sizes['x']
            self.gdf = gdf

        else:

            print(f'Dem NaN: {dem.isnull().sum().values}')
            print(f'Slope NaN: {slope.isnull().sum().values}')
            print(f'Zs NaN: {zs.isnull().sum().values}')
            print(f'FlowDir NaN: {fdir.isnull().sum().values}')
            print(f'Zonas NaN: {zonas.isnull().sum().values}')

            raise Exception('ASCII do not exported. NaN do not match.')
            
        return dem, slope, zonas, zs, fdir

    def GridMatch(self, rows, cols):

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
        with open(f'{self.out_path}/gm_in.txt', 'w') as file:
            file.write(content)

        print("GridMatch input created successfully.")
        print(f"----- Executing gridmatch.exe. -----")

        os.chdir(self.out_path)
        subprocess.run([
            f'{self.out_path}/gridmatch.exe'],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL
        )

        print(f"----- gridmatch.exe finished. -----")

        # Results
        mismatches = {}

        with open(f'{self.out_path}/GridMatchLog.txt', 'r') as file:
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

    def TopoIndex(self, rows, cols):

        iterations = 100

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
        with open(f'{self.out_path}/tpx_in.txt', 'w') as file:
            file.write(content)

        print("TopoIndex input created successfully.\n")
        print(f"----- Executing TopoIndex.exe with {iterations} iterations. -----")

        os.chdir(self.out_path)
        if not os.path.exists(f'{self.out_path}/tpx'):
            os.makedirs(f'{self.out_path}/tpx')
        subprocess.run(
            [f'{self.out_path}/TopoIndex.exe'],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL
        )

        print(f"----- TopoIndex.exe finished. -----")

        # Results

        with open(f'{self.out_path}/TopoIndexLog.txt', 'r') as file:
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

    def tr_in_creation(self, imax, rows, cols, nwf, hora, i, cohesion, friccion, gamma, ks, output_suffix):

        # tr_in file creation
        # Define the variable values
        tx = 1
        nmax1 = 30
        nzs = 10
        mmax2 = 0
        nper = 1
        zmin = 0.01
        uww = 9.8e3
        t = hora * 3600
        zones = self.zones
        zmax = depth = -3
        rizero = 7.7e-08
        min_slope_angle = 0
        cri = i / (1000 * 3600)  # Intensidad en m/s 

        cohesion = cohesion * 1000
        phi = friccion
        gamma = gamma * 1000
        k_sat = ks
        diffus = ks * 100
        theta_sat = [-0.23] * zones
        theta_res = [-0.48] * zones
        alpha = [-0.06] * zones

        capt = [0, t]

        slope_file = f'{self.out_path}/slope.asc'
        zone_file = f'{self.out_path}/zonas.asc'
        depth_file = f'{self.out_path}/zs.asc'
        init_depth_file = f'{self.out_path}/zs.asc'
        infil_rate_file = "none"
        rainfall_files = "none"

        runoff_receptor_file = "tpx\\TIdscelGrid_project.txt"
        runoff_order_file = "tpx\\TIcelindxList_project.txt"
        runoff_cell_list_file = "tpx\\TIdscelList_project.txt"
        runoff_weighting_file = "tpx\\TIwfactorList_project.txt"
        output_folder = "Resultados\\"

        save_runoff_grid = False
        save_factor_of_safety_grid = True
        save_depth_of_safety_grid = False
        save_pore_pressure_grid = False
        save_infil_rate_grid = False
        save_unsat_zone_flux_grid = False
        save_pressure_head_flag = 0
        num_output_times = 1
        output_times = t
        skip_other_timesteps = False
        use_analytic_solution = True
        estimate_positive_pressure_head = True
        use_psi0 = True
        log_mass_balance = True
        flow_direction = "gener"
        add_steady_background_flux = True

        # Generate the text
        content = "Name of project (up to 255 characters)\n"
        content += "Project\n"
        content += "imax, row, col, nwf, tx, nmax\n"
        content += f"{imax}, {rows}, {cols}, {nwf}, {tx}, {nmax1}\n"
        content += "nzs, mmax, nper, zmin, uww, t, zones\n"
        content += f"{nzs}, {mmax2}, {nper}, {zmin}, {uww}, {t}, {zones}\n"
        content += "zmax, depth, rizero, Min_Slope_Angle (degrees)\n"
        content += f"{zmax}, {depth}, {rizero}, {min_slope_angle}\n"

        for i in range(zones):
            content += f"zone,{i+1}\n"
            content += "cohesion,phi,uws,diffus,K-sat,Theta-sat,Theta-res,Alpha\n"
            content += f"{cohesion[i]},{phi[i]},{gamma[i]},{diffus[i]},{k_sat[i]},{theta_sat[i]},{theta_res[i]},{alpha[i]}\n"

        content += "cri(1), cri(2), ..., cri(nper)\n"
        content += f"{cri}\n"
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
        content += f"{output_times}\n"
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

        # Write the content to the file
        with open(f'{self.out_path}/tr_in.txt', 'w') as file:
            file.write(content)

        print("--- TRIGRS input created successfully. ---")

        return

    def TRIGRS_main(self):

        print(f"----- Executing TRIGRS.exe. -----")

        subprocess.run(
            [f'{self.out_path}/TRIGRS.exe'],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL
        )

        print(f"----- TRIGRS.exe finished. -----")

class SHALSTAB(geohazards):

    def __init__(self, dem_path, geo_path, geoColumns, q, zsPath='', shalstabPath='', criticalRainPath='', exportZs=False, exportShalstab=False, exportCriticalRain=False):

        geohazards.__init__(self)

        #Imports rasters and opens them with rasterio
        dem = rioxarray.open_rasterio(dem_path, mask_and_scale=True)

        dem_ = rd.LoadGDAL(dem_path)
        accum = rd.FlowAccumulation(dem_, method='D8')
        accum = xr.DataArray(accum, coords=[dem.coords['y'],dem.coords['x']])

        gdf = gpd.read_file(geo_path)

        cohesion, friccion, gamma, permeabilidad = geoColumns

        #Rasterize from geology
        c = make_geocube(vector_data=gdf, measurements = [cohesion], 
                        like=dem, fill = np.nan)[cohesion]
        p = make_geocube(vector_data=gdf, measurements = [friccion], 
                        like=dem, fill = np.nan)[friccion]
        p = np.radians(p)
        g = make_geocube(vector_data=gdf, measurements = [gamma], 
                        like=dem, fill = np.nan)[gamma]
        k = make_geocube(vector_data=gdf, measurements = [permeabilidad], 
                        like = dem, fill = np.nan)[permeabilidad]

        #Calculates slope
        slope = rd.TerrainAttribute(dem_, attrib='slope_radians')
        slope = xr.DataArray(slope, coords=[dem.coords['y'],dem.coords['x']])
        slope_rad = slope.where(slope!=-9999.)

        #Calculates zs from Catani
        zs = self.Catani(dem_path, geo_path)

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
        shalstab.rio.write_crs(dem.rio.crs, inplace=True)

        if exportShalstab:
            shalstab.rio.to_raster(shalstabPath)
        if exportZs:
            zs.rio.to_raster(zsPath)

        #Calculates critical rainfall
        criticalRain = 1000 * (k * zs * np.cos(slope_rad) * np.sin(slope_rad)) * (dem.rio.resolution()[0] / accum) * ((g / gammaw) * (1 - (np.tan(slope_rad) / np.tan(p))) + c / (gammaw * zs * np.cos(slope_rad)**2 * np.tan(p)))
        criticalRain = criticalRain.where(criticalRain>=0)

        criticalRain = np.where(izq_un >= der_unstable, -2, # Unconditionally unstable
                       np.where(izq_un < der_stable, -1, # Unconditionally stable
                       criticalRain))
        criticalRain = xr.DataArray(criticalRain, coords=[dem.coords['y'],dem.coords['x']])

        criticalRain.rio.write_nodata(-9999, inplace=True)

        if exportCriticalRain:
            criticalRain.rio.to_raster(criticalRainPath)

        gdf = gpd.read_file(geo_path)
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

        self.zs = zs
        self.shalstab = shalstab
        self.criticalRain = criticalRain

        return 
    
    def __call__(self):

        return self.shalstab, self.criticalRain

class FS(geohazards):

    def __init__(self):

        geohazards.__init__(self)

    def __call__(self, dem_path, geo_path, zw_path, c, phi, gammas, output='./FS.tif'):

        dem = rasterio.open(dem_path)
        zw = rasterio.open(zw_path)

        #Calculates slope
        slope = self.Slope(dem_path)
        slope_array = slope.read(1)

        #Sets no data to np.nan
        slope_array[slope_array == slope.nodata] = np.nan

        #Converts slope from degrees to radians
        slope_rad = np.radians(slope_array)

        gammaw = 9.81

        C = make_geocube(vector_data=geo_path, measurements = [c], like=dem, fill = np.nan)[c]
        phi = make_geocube(vector_data=geo_path, measurements = [phi], like=dem, fill = np.nan)[phi]
        gammas = make_geocube(vector_data=geo_path, measurements = [gammas], like=dem, fill = np.nan)[gammas]

        FS = C + (gammas - gammaw) * zw * (np.cos(slope_rad)**2) * np.tan(phi) / gammas * zw * np.sin(slope_rad) * np.cos(slope_rad)

        #Copies metadata from dem
        meta = dem.meta.copy()
        meta.update(compress='lzw', nodata=-9999)

        #Exports raster output file and assign metadata
        with rasterio.open(output, 'w+', **meta) as out:
            out.write_band(1, FS)
        
        #Returns soil thickness rasterio file
        return rasterio.open(output)