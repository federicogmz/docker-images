#%%
import numpy as np
import statistics, os
import matplotlib.pyplot as plt 
from matplotlib import colors
import rasterio as rio
import pandas as pd

def TRIGRS(tr, i, hora, zonas, prefix, zmax, depth, export, FOSM):

    prefijo = f'{prefix}{tr}'
    cri=(i)/(1000*3600) #Intensidad en m/s
    capt_int=hora*3600 #Tiempo en s

    TopoIndexLog=open('TopoIndexLog.txt', "r")
    list_of_lines_tx = TopoIndexLog.readlines()
    list_of_lines_tx[-5].split('  ')
    imax=list_of_lines_tx[-5].split('  ')[1]
    row=list_of_lines_tx[-5].split('  ')[2]
    col=list_of_lines_tx[-5].split('  ')[3]
    nwf=list_of_lines_tx[-5].split('  ')[4][:-1]

    #Para separar cada fila en una lista
    tr_in_file = open('tr_in.txt', "r")
    list_of_lines = tr_in_file.readlines()

    n=zonas
    n_zones=list_of_lines[5].split(',')
    n_zones[-1]=str(n)+'\n'
    list_of_lines[5]=','.join(n_zones)

    topo=list_of_lines[3].split(',')
    topo[0]=imax
    topo[1]=row
    topo[2]=col
    topo[3]=nwf
    list_of_lines[3]=','.join(topo)

    #Use a positive value of alpha to use the unsaturated
    sat='Sat'
    Alpha_1=-0.06
    Alpha_2=-0.06
    Alpha_3=-0.06
    Alpha_4=-0.06
    Alpha_5=-0.06
    Alpha_6=-0.06
    Alpha_7=-0.06
    Alpha_8=-0.06
    Alpha_9=-0.06

    list_of_lines[12+(3*(n-1))]=str(cri)+'\n'
    capt_intt=list_of_lines[14+(3*(n-1))].split(',')
    capt_intt[-1]=str(capt_int)+'\n'
    list_of_lines[14+(3*(n-1))]=','.join(capt_intt)
    capt_intt_2=list_of_lines[5].split(',')
    capt_intt_2[-2]=str(capt_int)
    list_of_lines[5]=','.join(capt_intt_2)
    list_of_lines[56+(3*(n-1))]=str(capt_int)+'\n'

    zmax_depth=list_of_lines[7].split(',')
    zmax_depth[:2]=[str(zmax),str(depth)]
    list_of_lines[7]=','.join(zmax_depth)

    for i in np.arange(1,n+1):
        alpha_int=list_of_lines[10+(3*(i-1))].split(',')
        alpha_int[-1]=str(locals()["Alpha_"+str(i)])+'\n'
        list_of_lines[10+(3*(i-1))]=','.join(alpha_int)
        
    "Theta_res"
    Theta_res_1=-0.48
    Theta_res_2=-0.48
    Theta_res_3=-0.48
    Theta_res_4=-0.48
    Theta_res_5=-0.48
    Theta_res_6=-0.48
    Theta_res_7=-0.48
    Theta_res_8=-0.48
    Theta_res_9=-0.48

    for i in np.arange(1,n+1):
        Theta_res_int=list_of_lines[10+(3*(i-1))].split(',')
        Theta_res_int[-2]=str(locals()["Theta_res_"+str(i)])
        list_of_lines[10+(3*(i-1))]=','.join(Theta_res_int)

    "Theta-sat"
    Theta_sat_1=-0.23
    Theta_sat_2=-0.23
    Theta_sat_3=-0.23
    Theta_sat_4=-0.23
    Theta_sat_5=-0.23
    Theta_sat_6=-0.23
    Theta_sat_7=-0.23
    Theta_sat_8=-0.23
    Theta_sat_9=-0.23

    for i in np.arange(1,n+1):
        Theta_sat_int=list_of_lines[10+(3*(i-1))].split(',')
        Theta_sat_int[-3]=str(locals()["Theta_sat_"+str(i)])
        list_of_lines[10+(3*(i-1))]=','.join(Theta_sat_int)

    "C"
    variation_C=(40/100)

    C_mean_1=20000
    C_mean_2=30000
    C_mean_3=20000
    # C_mean_4=16000
    # C_mean_5=7000
    # C_mean_6=7000
    # C_mean_7=7000
    # C_mean_8=20000
    # C_mean_9=10000

    for i in np.arange(1,n+1):
        C_mean=locals()["C_mean_"+str(i)]
        locals()["C_vector_"+str(i)]=[C_mean,C_mean+(C_mean*variation_C),C_mean-(C_mean*variation_C)]
        locals()["C_var_"+str(i)]=statistics.variance(locals()["C_vector_"+str(i)])
        locals()["C_std_"+str(i)]=statistics.stdev(locals()["C_vector_"+str(i)])

    "Phi"
    variation_phi=(13/100)

    phi_mean_1=36
    phi_mean_2=30
    phi_mean_3=33
    # phi_mean_4=29
    # phi_mean_5=23
    # phi_mean_6=23
    # phi_mean_7=23
    # phi_mean_8=32
    # phi_mean_9=27

    for i in np.arange(1,n+1):
        phi_mean=locals()["phi_mean_"+str(i)]
        locals()["phi_vector_"+str(i)]=[phi_mean,phi_mean+(phi_mean*variation_phi),phi_mean-(phi_mean*variation_phi)]
        locals()["phi_var_"+str(i)]=statistics.variance(locals()["phi_vector_"+str(i)])
        locals()["phi_std_"+str(i)]=statistics.stdev(locals()["phi_vector_"+str(i)])

    "Uws"
    variation_uws=(7/100)

    uws_mean_1=21*1000
    uws_mean_2=18.5*1000
    uws_mean_3=20.5*1000
    # uws_mean_4=16.7*1000
    # uws_mean_5=17.35*1000
    # uws_mean_6=17.35*1000
    # uws_mean_7=17.35*1000
    # uws_mean_8=19.25*1000
    # uws_mean_9=19.25*1000

    for i in np.arange(1,n+1):
        uws_mean=locals()["uws_mean_"+str(i)]
        locals()["uws_vector_"+str(i)]=[uws_mean,uws_mean+(uws_mean*variation_uws),uws_mean-(uws_mean*variation_uws)]
        locals()["uws_var_"+str(i)]=statistics.variance(locals()["uws_vector_"+str(i)])
        locals()["uws_std_"+str(i)]=statistics.stdev(locals()["uws_vector_"+str(i)])

    "KS"
    variation_ks=(90/100)

    ks_mean_1=5.00E-3
    ks_mean_2=5.00E-7
    ks_mean_3=5.00E-6
    # ks_mean_4=1.26E-6
    # ks_mean_5=6.22E-7
    # ks_mean_6=6.22E-7
    # ks_mean_7=6.22E-7
    # ks_mean_8=1.00E-8
    # ks_mean_9=1.47E-7

    for i in np.arange(1,n+1):
        ks_mean=locals()["ks_mean_"+str(i)]
        locals()["ks_vector_"+str(i)]=[ks_mean,ks_mean+(ks_mean*variation_ks),ks_mean-(ks_mean*variation_ks)]
        locals()["ks_var_"+str(i)]=statistics.variance(locals()["ks_vector_"+str(i)])
        locals()["ks_std_"+str(i)]=statistics.stdev(locals()["ks_vector_"+str(i)])

    '''Parámetros a variar'''
    for j in np.arange(1,n+1):
        locals()["C_"+str(j)]=[str(i) for i in [locals()["C_mean_"+str(j)],locals()["C_mean_"+str(j)]+(0.1*locals()["C_mean_"+str(j)])]]
        locals()["phi_"+str(j)]=[str(i) for i in [locals()["phi_mean_"+str(j)],locals()["phi_mean_"+str(j)]+(0.1*locals()["phi_mean_"+str(j)])]]
        locals()["uws_"+str(j)]=[str(i) for i in [locals()["uws_mean_"+str(j)],locals()["uws_mean_"+str(j)]+(0.1*locals()["uws_mean_"+str(j)])]]
        locals()["ks_"+str(j)]=[str(i) for i in [locals()["ks_mean_"+str(j)],locals()["ks_mean_"+str(j)]+(0.1*locals()["ks_mean_"+str(j)])]]

    '''Ciclo for para correr trigrs'''

    print(f'Corre para valores medios.')

    for w in np.arange(1,n+1):
        locals()["Zone"+str(w)]=list_of_lines[10+(3*(w-1))].split(',')
        locals()["Zone"+str(w)]=list(locals()["Zone"+str(w)])
        locals()["Zone"+str(w)][0]=locals()["C_"+str(w)][0] #Cohesión
        locals()["Zone"+str(w)][1]=locals()["phi_"+str(w)][0] #Ángulo de fricción
        locals()["Zone"+str(w)][2]=locals()["uws_"+str(w)][0] #Peso específico 
        locals()["Zone"+str(w)][4]=locals()["ks_"+str(w)][0] #Peso específico 
        locals()["Zone"+str(w)][3]=str(float(locals()["ks_"+str(w)][0])*100) #Diffusion
        list_of_lines[10+(3*(w-1))]=','.join(locals()["Zone"+str(w)])        

    #Para cambiar el nombre de los archivos que se van a generar
    list_of_lines[38+(3*(n-1))]=f'{prefijo}\n'

    # Sobre escribir el archivo txt
    a_file = open('tr_in.txt', "w")
    a_file.writelines(list_of_lines)
    a_file.close()

    # ejecutar trigrs
    os.system('"'+"TRIGRS.exe"+'"')

    if FOSM:
            
        print('Corre para valores variados.')

        for w in np.arange(1,n+1):
            locals()["Zone"+str(w)]=list_of_lines[10+(3*(w-1))].split(',')
            locals()["Zone"+str(w)]=list(locals()["Zone"+str(w)])
                
        for a,b,c,d,e in zip([1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1],[f'{prefijo}_C\n',f'{prefijo}_Phi\n',f'{prefijo}_Uws\n',f'{prefijo}_Ks\n']):
            for w in np.arange(1,n+1):
                locals()["Zone"+str(w)][0]=locals()["C_"+str(w)][a]
                locals()["Zone"+str(w)][1]=locals()["phi_"+str(w)][b]
                locals()["Zone"+str(w)][2]=locals()["uws_"+str(w)][c]
                locals()["Zone"+str(w)][4]=locals()["ks_"+str(w)][d]
                locals()["Zone"+str(w)][3]=str(float(locals()["ks_"+str(w)][d])*100)
                list_of_lines[10+(3*(w-1))]=','.join(locals()["Zone"+str(w)])
                        
            list_of_lines[38+(3*(n-1))]=e
            
            # # Sobre escribir el archivo txt
            a_file = open('tr_in.txt', "w")
            a_file.writelines(list_of_lines)
            a_file.close()
                
            #ejecutar trigrs
            os.system('"'+"TRIGRS.exe"+'"')

        '''Val medio'''
        FS_1_medio = np.genfromtxt(f'./Resultados/TRfs_min_{prefijo}_1.txt', skip_header=6, delimiter=' ')
        FS_1_medio =np.where(FS_1_medio==-9999,np.nan,FS_1_medio)
        FS_1_medio =np.where(FS_1_medio==10,10,FS_1_medio)

        '''Cambio solo de la cohesion'''
        FS_1_cohe = np.genfromtxt(f'./Resultados/TRfs_min_{prefijo}_C_1.txt', skip_header=6, delimiter=' ')
        FS_1_cohe=np.where(FS_1_cohe==-9999,np.nan,FS_1_cohe)
        FS_1_cohe=np.where(FS_1_cohe>=9.8,10,FS_1_cohe)

        '''Cambio solo del ángulo de fricción'''
        FS_1_phi = np.genfromtxt(f'./Resultados/TRfs_min_{prefijo}_Phi_1.txt', skip_header=6, delimiter=' ')
        FS_1_phi=np.where(FS_1_phi==-9999,np.nan,FS_1_phi)
        FS_1_phi=np.where(FS_1_phi>=9.8,10,FS_1_phi)


        '''Cambio solo del peso específico'''
        FS_1_uws = np.genfromtxt(f'./Resultados/TRfs_min_{prefijo}_Uws_1.txt', skip_header=6, delimiter=' ')
        FS_1_uws=np.where(FS_1_uws==-9999,np.nan,FS_1_uws)
        FS_1_uws=np.where(FS_1_uws>=9.8,10,FS_1_uws)

        '''Cambio solo del KS'''
        FS_1_ks = np.genfromtxt(f'./Resultados/TRfs_min_{prefijo}_Ks_1.txt', skip_header=6, delimiter=' ')
        FS_1_ks=np.where(FS_1_ks==-9999,np.nan,FS_1_ks)
        FS_1_ks=np.where(FS_1_ks>=9.8,10,FS_1_ks)

        '''Mapa de zonas'''
        zonas = np.genfromtxt('zonas.asc', skip_header=6, delimiter=' ')
        zonas =np.where(zonas ==-9999,np.nan,zonas )

        '''FOSM'''

        dFS_1_cohe=FS_1_cohe-FS_1_medio
        dFS_1_phi=FS_1_phi-FS_1_medio
        dFS_1_uws=FS_1_uws-FS_1_medio
        dFS_1_ks=FS_1_ks-FS_1_medio

        for j in np.arange(1,n+1):
            if j==1:
                locals()["dFS_"+str(hora)+"_dx_cohe"]=np.where(zonas==j,dFS_1_cohe/(0.1*locals()["C_mean_"+str(j)]),dFS_1_cohe)
            else:
                locals()["dFS_"+str(hora)+"_dx_cohe"]=np.where(zonas==j,dFS_1_cohe/(0.1*locals()["C_mean_"+str(j)]),locals()["dFS_"+str(hora)+"_dx_cohe"])

        for j in np.arange(1,n+1):
            if j==1:
                locals()["dFS_"+str(hora)+"_dx_phi"]=np.where(zonas==j,dFS_1_phi/(0.1*locals()["phi_mean_"+str(j)]),dFS_1_phi)
            else:
                locals()["dFS_"+str(hora)+"_dx_phi"]=np.where(zonas==j,dFS_1_phi/(0.1*locals()["phi_mean_"+str(j)]),locals()["dFS_"+str(hora)+"_dx_phi"])

        for j in np.arange(1,n+1):
            if j==1:
                locals()["dFS_"+str(hora)+"_dx_uws"]=np.where(zonas==j,dFS_1_uws/(0.1*locals()["uws_mean_"+str(j)]),dFS_1_uws)
            else:
                locals()["dFS_"+str(hora)+"_dx_uws"]=np.where(zonas==j,dFS_1_uws/(0.1*locals()["uws_mean_"+str(j)]),locals()["dFS_"+str(hora)+"_dx_uws"])

        for j in np.arange(1,n+1):
            if j==1:
                locals()["dFS_"+str(hora)+"_dx_ks"]=np.where(zonas==j,dFS_1_ks/(0.1*locals()["ks_mean_"+str(j)]),dFS_1_ks)
            else:
                locals()["dFS_"+str(hora)+"_dx_ks"]=np.where(zonas==j,dFS_1_ks/(0.1*locals()["ks_mean_"+str(j)]),locals()["dFS_"+str(hora)+"_dx_ks"])

        for j in np.arange(1,n+1):
            if j==1:
                locals()["VF_"+str(hora)+"_C"]=np.where(zonas==j,(locals()["dFS_"+str(hora)+"_dx_cohe"]**(2))*locals()["C_var_"+str(j)],locals()["dFS_"+str(hora)+"_dx_cohe"]) 
            else:
                locals()["VF_"+str(hora)+"_C"]=np.where(zonas==j,(locals()["dFS_"+str(hora)+"_dx_cohe"]**(2))*locals()["C_var_"+str(j)],locals()["VF_"+str(hora)+"_C"])

        for j in np.arange(1,n+1):
            if j==1:
                locals()["VF_"+str(hora)+"_phi"]=np.where(zonas==j,(locals()["dFS_"+str(hora)+"_dx_phi"]**(2))*locals()["phi_var_"+str(j)],locals()["dFS_"+str(hora)+"_dx_phi"])
            else:
                locals()["VF_"+str(hora)+"_phi"]=np.where(zonas==j,(locals()["dFS_"+str(hora)+"_dx_phi"]**(2))*locals()["phi_var_"+str(j)],locals()["VF_"+str(hora)+"_phi"])

        for j in np.arange(1,n+1):
            if j==1:
                locals()["VF_"+str(hora)+"_uws"]=np.where(zonas==j,(locals()["dFS_"+str(hora)+"_dx_uws"]**(2))*locals()["uws_var_"+str(j)],locals()["dFS_"+str(hora)+"_dx_uws"])
            else:
                locals()["VF_"+str(hora)+"_uws"]=np.where(zonas==j,(locals()["dFS_"+str(hora)+"_dx_uws"]**(2))*locals()["uws_var_"+str(j)],locals()["VF_"+str(hora)+"_uws"])

        for j in np.arange(1,n+1):
            if j==1:
                locals()["VF_"+str(hora)+"_ks"]=np.where(zonas==j,(locals()["dFS_"+str(hora)+"_dx_ks"]**(2))*locals()["ks_var_"+str(j)],locals()["dFS_"+str(hora)+"_dx_ks"])
            else:
                locals()["VF_"+str(hora)+"_ks"]=np.where(zonas==j,(locals()["dFS_"+str(hora)+"_dx_ks"]**(2))*locals()["ks_var_"+str(j)],locals()["VF_"+str(hora)+"_ks"])

        std_1_Fs=np.sqrt(locals()["VF_"+str(hora)+"_C"]+locals()["VF_"+str(hora)+"_phi"]+locals()["VF_"+str(hora)+"_uws"]+locals()["VF_"+str(hora)+"_ks"]) #OJO si hay datos faltantes, revisar esto

        VF_1=locals()["VF_"+str(hora)+"_C"]+locals()["VF_"+str(hora)+"_phi"]+locals()["VF_"+str(hora)+"_uws"]+locals()["VF_"+str(hora)+"_ks"]

        Ind_conf_1= (FS_1_medio-1)/std_1_Fs 
        Ind_conf_1=np.where(Ind_conf_1<0,-9999,Ind_conf_1)
        Ind_conf_1=np.where(Ind_conf_1==np.math.inf,9999,Ind_conf_1)

        ind = np.arange(4)  
        width = 0.35    
        porc=pd.DataFrame([np.nanmean(locals()["VF_"+str(hora)+"_C"]/VF_1),np.nanmean(locals()["VF_"+str(hora)+"_phi"]/VF_1),np.nanmean(locals()["VF_"+str(hora)+"_uws"]/VF_1), np.nanmean(locals()["VF_"+str(hora)+"_ks"]/VF_1)])
        porc.plot.bar(legend=False)
        plt.ylabel('Porcentaje de la varianza')
        plt.xticks(ind, ('Cohesión', 'Ángulo de\nfricción', 'Peso\nEspecífico', 'Ks'), rotation=0)
        plt.tight_layout()
        plt.show()

        c_map = colors.ListedColormap(['red', 'yellow', 'green'])
        bounds = [0, 1.0, 2.5, 10]
        norm = colors.BoundaryNorm(bounds, c_map.N)
        plt.imshow(Ind_conf_1, cmap=c_map, norm=norm)

        if export:

            file=rio.open('fill.asc')
            raster=file.read(1)

            meta=file.profile
            file_transform = meta['transform']
            file_crs = meta['crs']

            with rio.open(f'./Resultados/Ind_conf_{prefijo}.tif', 'w', 
                            driver='Gtiff',height=raster.shape[0],width=raster.shape[1],count=1,
                            dtype='float64',nodata=-999,crs=file_crs ,transform=file_transform ) as dst:
                dst.write(Ind_conf_1,1)
#%%

os.chdir(r"G:\Unidades compartidas\Proyectos(Fede)\Tarso\Urbano\TRIGRS")
trs = [2,233,5,10,25,50,100,300,500]
ints = [13.080,14.015,17.107,19.766,23.153,25.561,28.058,32.409,34.344]
zonas = 9
hora = 4
zmax = -3
depth = -3

#Valores medios
for i,tr in zip(ints,trs):
    TRIGRS(tr, i, hora, zonas, prefix='M', zmax=zmax, depth=depth, FOSM=False, export=False)

#%%
os.chdir("./Mohan/Codigo/")
tr=2
i=13.080
zonas = 9
hora = 4
zmax = -3
depth = -3
TRIGRS(tr, i, hora, zonas, prefix='M', zmax=zmax, depth=depth, FOSM=False, export=False)
#%% Amenaza
tr = 100
i = 114.6750106
zonas = 3
hora = 30
zmax = -3
depth = -3
os.chdir(r"G:\Unidades compartidas\Proyectos(Fede)\Tarso\Urbano\TRIGRS")
TRIGRS(tr, i, hora, zonas, prefix='U', zmax=zmax, depth=depth, FOSM=True, export=True)

#%%
import rasterio

dem = rasterio.open('./Mohan/Codigo/dem.txt')
meta = dem.meta.copy()
meta.update(compress='lzw', nodata=-9999)

txt = rasterio.open(f'./Mohan/Codigo/out/TRz_at_fs_min_M2_1.txt').read(1)
with rasterio.open(f'./Mohan/Codigo/out/FS/TRz_at_fs_min_M2_1.tif', 'w+', **meta) as out:
    out.write_band(1,txt)


# for file in os.listdir('./Boton/Codigo/out/FS/'):
#     if file.startswith('TRfs_min'):
#         txt = rasterio.open(f'./Boton/Codigo/out/FS/{file}').read(1)
#         with rasterio.open(f'./Boton/Codigo/out/FS/{file[:-4]}.tif', 'w+', **meta) as out:
#             out.write_band(1,txt)