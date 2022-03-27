#%%
import numpy as np
import statistics 
import os
import matplotlib.pyplot as plt
import rasterio as rio
import math
import pandas as pd

ruta_in=r'G:/Unidades compartidas/Proyectos(Fede)/Tarso/Amenaza/TRIGRS/'

#%%

#Para leer 
with open(ruta_in+'tr_in.txt', "r+") as f:
     old = f.read() # read everything in the file

#Para separar cada fila en una lista
tr_in_file = open(ruta_in+'tr_in.txt', "r")
list_of_lines = tr_in_file.readlines()

n=3 #NUMERO DE ZONAS

# Zona 1: GW, Zona 2: SW, Zona 3: SC

"C"
variation_C=(40/100)
C_mean_1=20000
C_vector_1=[C_mean_1,C_mean_1+(C_mean_1*variation_C),C_mean_1-(C_mean_1*variation_C)]
C_var_1=statistics.variance(C_vector_1)
C_std_1=statistics.stdev(C_vector_1)

C_mean_2=20000
C_vector_2=[C_mean_2,C_mean_2+(C_mean_2*variation_C),C_mean_2-(C_mean_2*variation_C)]
C_var_2=statistics.variance(C_vector_2)
C_std_2=statistics.stdev(C_vector_2)

C_mean_3=30000
C_vector_3=[C_mean_3,C_mean_3+(C_mean_3*variation_C),C_mean_3-(C_mean_3*variation_C)]
C_var_3=statistics.variance(C_vector_3)
C_std_3=statistics.stdev(C_vector_3)

"Phi"
variation_phi=(13/100)
phi_mean_1=36
phi_vector_1=[phi_mean_1,phi_mean_1+(phi_mean_1*variation_phi),phi_mean_1-(phi_mean_1*variation_phi)]
phi_var_1=statistics.variance(phi_vector_1)
phi_std_1=statistics.stdev(phi_vector_1)

phi_mean_2=33
phi_vector_2=[phi_mean_2,phi_mean_2+(phi_mean_2*variation_phi),phi_mean_2-(phi_mean_2*variation_phi)]
phi_var_2=statistics.variance(phi_vector_2)
phi_std_2=statistics.stdev(phi_vector_2)

phi_mean_3=30
phi_vector_3=[phi_mean_3,phi_mean_3+(phi_mean_3*variation_phi),phi_mean_3-(phi_mean_3*variation_phi)]
phi_var_3=statistics.variance(phi_vector_3)
phi_std_3=statistics.stdev(phi_vector_3)

"Uws"
variation_uws=(7/100)
uws_mean_1=21000
uws_vector_1=[uws_mean_1,uws_mean_1+(uws_mean_1*variation_uws),uws_mean_1-(uws_mean_1*variation_uws)]
uws_var_1=statistics.variance(uws_vector_1)
uws_std_1=statistics.stdev(uws_vector_1)

uws_mean_2=20500
uws_vector_2=[uws_mean_2,uws_mean_2+(uws_mean_2*variation_uws),uws_mean_2-(uws_mean_2*variation_uws)]
uws_var_2=statistics.variance(uws_vector_2)
uws_std_2=statistics.stdev(uws_vector_2)

uws_mean_3=18500
uws_vector_3=[uws_mean_3,uws_mean_3+(uws_mean_3*variation_uws),uws_mean_3-(uws_mean_3*variation_uws)]
uws_var_3=statistics.variance(uws_vector_3)
uws_std_3=statistics.stdev(uws_vector_3)


"KS"
variation_ks=(90/100)
ks_mean_1=5.0E-03
ks_vector_1=[ks_mean_1,ks_mean_1+(ks_mean_1*variation_ks),ks_mean_1-(ks_mean_1*variation_ks)]
ks_var_1=statistics.variance(ks_vector_1)
ks_std_1=statistics.stdev(ks_vector_1)

ks_mean_2=5.0E-06
ks_vector_2=[ks_mean_2,ks_mean_2+(ks_mean_2*variation_ks),ks_mean_2-(ks_mean_2*variation_ks)]
ks_var_2=statistics.variance(ks_vector_2)
ks_std_2=statistics.stdev(ks_vector_2)

ks_mean_3=5.0E-07
ks_vector_3=[ks_mean_3,ks_mean_3+(ks_mean_3*variation_ks),ks_mean_3-(ks_mean_3*variation_ks)]
ks_var_3=statistics.variance(ks_vector_3)
ks_std_3=statistics.stdev(ks_vector_3)


'''Parámetros a variar'''
C_1=[str(i) for i in [C_mean_1,C_mean_1+(0.1*C_mean_1)]]
phi_1=[str(i) for i in [phi_mean_1,phi_mean_1+(0.1*phi_mean_1)]]
uws_1=[str(i) for i in [uws_mean_1,uws_mean_1+(0.1*uws_mean_1)]]
ks_1=[str(i) for i in [ks_mean_1,ks_mean_1+(0.1*ks_mean_1)]]

C_2=[str(i) for i in [C_mean_2,C_mean_2+(0.1*C_mean_2)]]
phi_2=[str(i) for i in [phi_mean_2,phi_mean_2+(0.1*phi_mean_2)]]
uws_2=[str(i) for i in [uws_mean_2,uws_mean_2+(0.1*uws_mean_2)]]
ks_2=[str(i) for i in [ks_mean_2,ks_mean_2+(0.1*ks_mean_2)]]

C_3=[str(i) for i in [C_mean_3,C_mean_3+(0.1*C_mean_3)]]
phi_3=[str(i) for i in [phi_mean_3,phi_mean_3+(0.1*phi_mean_3)]]
uws_3=[str(i) for i in [uws_mean_3,uws_mean_3+(0.1*uws_mean_3)]]
ks_3=[str(i) for i in [ks_mean_3,ks_mean_3+(0.1*ks_mean_3)]]

'''Ciclo for para correr trigrs'''
contador=0
for j in np.arange(2): #una iteración para el promedio y cambiando los otros
    contador=contador+1
    if contador==1:    
        print(contador)
        Zone1=list_of_lines[10].split(',')
        Zone1=list(Zone1)
        Zone1[0]=C_1[0] #Cohesión
        Zone1[1]=phi_1[0] #Ángulo de fricción
        Zone1[2]=uws_1[0] #Peso específico 
        Zone1[4]=ks_1[0] #Peso específico 
        list_of_lines[10]=','.join(Zone1)
        
        Zone2=list_of_lines[13].split(',') #CAMBIA!!
        Zone2=list(Zone2)
        Zone2[0]=C_2[0] #Cohesión
        Zone2[1]=phi_2[0] #Ángulo de fricción
        Zone2[2]=uws_2[0] #Peso específico 
        Zone2[4]=ks_2[0] #Peso específico 
        list_of_lines[13]=','.join(Zone2) #CAMBIA!!
        
        Zone3=list_of_lines[16].split(',') #CAMBIA!!
        Zone3=list(Zone3)
        Zone3[0]=C_3[0] #Cohesión
        Zone3[1]=phi_3[0] #Ángulo de fricción
        Zone3[2]=uws_3[0] #Peso específico 
        Zone3[4]=ks_3[0] #Peso específico 
        list_of_lines[16]=','.join(Zone3)  #CAMBIA!!
        
        #Para cambiar el nombre de los archivos que se van a generar
        list_of_lines[38+(3*(n-1))]='S_5_30\n'  #MODIFICAR MANUALMENTE
        
        # Sobre escribir el archivo txt
        a_file = open(ruta_in+'tr_in.txt', "w")
        a_file.writelines(list_of_lines)
        a_file.close()
          
        # ejecutar trigrs
        os.chdir(ruta_in)
        os.system(r'"G:/Unidades compartidas/Proyectos(Fede)/Tarso/Amenaza/TRIGRS/TRIGRS.exe"')    

    elif contador>1:
        Zone1=list_of_lines[10].split(',')
        Zone1=list(Zone1)
        Zone2=list_of_lines[13].split(',')
        Zone2=list(Zone2)
        Zone3=list_of_lines[16].split(',')
        Zone3=list(Zone3)
        
        for a,b,c,d,e in zip([1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1],['Cohe_4h\n','Phi_4h\n','Uws_4h\n','Ks_4h\n']):
           Zone1[0]=C_1[a]
           Zone2[0]=C_2[a] 
           Zone3[0]=C_3[a]  
            
           Zone1[1]=phi_1[b]
           Zone2[1]=phi_2[b] 
           Zone3[1]=phi_3[b] 
            
           Zone1[2]=uws_1[c]
           Zone2[2]=uws_2[c] 
           Zone3[2]=uws_3[c]

           Zone1[4]=ks_1[d]
           Zone2[4]=ks_2[d] 
           Zone3[4]=ks_3[d]
            
           list_of_lines[10]=','.join(Zone1)
           list_of_lines[13]=','.join(Zone2)
           list_of_lines[16]=','.join(Zone3)
            
           list_of_lines[38+(3*(n-1))]=e
            
            # # Sobre escribir el archivo txt
           a_file = open(ruta_in+'tr_in.txt', "w")
           a_file.writelines(list_of_lines)
           a_file.close()
              
            #ejecutar trigrs
           os.chdir(ruta_in)
           os.system(r'"G:/Unidades compartidas/Proyectos(Fede)/Tarso/Amenaza/TRIGRS/TRIGRS.exe"')

'''Val medio'''
FS_1_medio = np.genfromtxt(ruta_in + '/Resultados/TRfs_min_S_5_30_1.txt', skip_header=6, delimiter=' ')
FS_1_medio =np.where(FS_1_medio==-9999,np.nan,FS_1_medio)
FS_1_medio =np.where(FS_1_medio==10,10,FS_1_medio)
# plt.imshow(FS_1_medio)


'''Cambio solo de la cohesion'''
FS_1_cohe = np.genfromtxt(ruta_in + '/Resultados/TRfs_min_Cohe_4h_1.txt', skip_header=6, delimiter=' ')
FS_1_cohe=np.where(FS_1_cohe==-9999,np.nan,FS_1_cohe)
FS_1_cohe=np.where(FS_1_cohe>=9.8,10,FS_1_cohe)

'''Cambio solo del ángulo de fricción'''
FS_1_phi = np.genfromtxt(ruta_in + '/Resultados/TRfs_min_phi_4h_1.txt', skip_header=6, delimiter=' ')
FS_1_phi=np.where(FS_1_phi==-9999,np.nan,FS_1_phi)
FS_1_phi=np.where(FS_1_phi>=9.8,10,FS_1_phi)


'''Cambio solo del peso específico'''
FS_1_uws = np.genfromtxt(ruta_in + '/Resultados/TRfs_min_uws_4h_1.txt', skip_header=6, delimiter=' ')
FS_1_uws=np.where(FS_1_uws==-9999,np.nan,FS_1_uws)
FS_1_uws=np.where(FS_1_uws>=9.8,10,FS_1_uws)


# '''Cambio solo del KS'''
FS_1_ks = np.genfromtxt(ruta_in + '/Resultados/TRfs_min_ks_4h_1.txt', skip_header=6, delimiter=' ')
FS_1_ks=np.where(FS_1_ks==-9999,np.nan,FS_1_ks)
FS_1_ks=np.where(FS_1_ks>=9.8,10,FS_1_ks)



'''Mapa de zonas''' #para dividir sobre cada valor de cada zona
# "G:/Unidades compartidas/ETCR(9)/3.Mov. Masa/11.TRIGRS/Plancha/Prueba_Daissy/zone.txt"
zonas = np.genfromtxt(f'{ruta_in}/zonas.asc', skip_header=6, delimiter=' ')
zonas =np.where(zonas ==-9999,np.nan,zonas )


'''FOSM'''

dFS_1_cohe=FS_1_cohe-FS_1_medio
dFS_1_phi=FS_1_phi-FS_1_medio
dFS_1_uws=FS_1_uws-FS_1_medio
dFS_1_ks=FS_1_ks-FS_1_medio

dFS_1_dx_cohe=np.where(zonas==1.,dFS_1_cohe/(0.1*C_mean_1),dFS_1_cohe)
dFS_1_dx_cohe=np.where(zonas==2.,dFS_1_cohe/(0.1*C_mean_2),dFS_1_dx_cohe)
dFS_1_dx_cohe=np.where(zonas==3.,dFS_1_cohe/(0.1*C_mean_3),dFS_1_dx_cohe)

dFS_1_dx_phi=np.where(zonas==1.,dFS_1_phi/(0.1*phi_mean_1),dFS_1_phi)
dFS_1_dx_phi=np.where(zonas==2.,dFS_1_phi/(0.1*phi_mean_2),dFS_1_dx_phi)
dFS_1_dx_phi=np.where(zonas==3.,dFS_1_phi/(0.1*phi_mean_3),dFS_1_dx_phi)

dFS_1_dx_uws=np.where(zonas==1.,dFS_1_uws/(0.1*uws_mean_1),dFS_1_uws)
dFS_1_dx_uws=np.where(zonas==2.,dFS_1_uws/(0.1*uws_mean_2),dFS_1_dx_uws)
dFS_1_dx_uws=np.where(zonas==3.,dFS_1_uws/(0.1*uws_mean_3),dFS_1_dx_uws)

dFS_1_dx_ks=np.where(zonas==1.,dFS_1_ks/(0.1*ks_mean_1),dFS_1_ks)
dFS_1_dx_ks=np.where(zonas==2.,dFS_1_ks/(0.1*ks_mean_2),dFS_1_dx_ks)
dFS_1_dx_ks=np.where(zonas==3.,dFS_1_ks/(0.1*ks_mean_3),dFS_1_dx_ks)

VF_1_C=np.where(zonas==1.,(dFS_1_dx_cohe**(2))*C_var_1,dFS_1_dx_cohe) 
VF_1_C=np.where(zonas==2.,(dFS_1_dx_cohe**(2))*C_var_2,VF_1_C)
VF_1_C=np.where(zonas==3.,(dFS_1_dx_cohe**(2))*C_var_3,VF_1_C)

VF_1_phi=np.where(zonas==1.,(dFS_1_dx_phi**(2))*phi_var_1,dFS_1_dx_phi)
VF_1_phi=np.where(zonas==2.,(dFS_1_dx_phi**(2))*phi_var_2,VF_1_phi)
VF_1_phi=np.where(zonas==3.,(dFS_1_dx_phi**(2))*phi_var_3,VF_1_phi)

VF_1_uws=np.where(zonas==1.,(dFS_1_dx_uws**(2))*uws_var_1,dFS_1_dx_uws)
VF_1_uws=np.where(zonas==2.,(dFS_1_dx_uws**(2))*uws_var_2,VF_1_uws)
VF_1_uws=np.where(zonas==3.,(dFS_1_dx_uws**(2))*uws_var_3,VF_1_uws)


VF_1_ks=np.where(zonas==1.,(dFS_1_dx_ks**(2))*ks_var_1,dFS_1_dx_ks)
VF_1_ks=np.where(zonas==2.,(dFS_1_dx_ks**(2))*ks_var_2,VF_1_ks)
VF_1_ks=np.where(zonas==3.,(dFS_1_dx_ks**(2))*ks_var_3,VF_1_ks)

std_1_Fs=np.sqrt(VF_1_C+VF_1_phi+VF_1_uws+VF_1_ks) #OJO si hay datos faltantes, revisar esto

VF_1=VF_1_C+VF_1_phi+VF_1_uws+VF_1_ks




Ind_conf_1= (FS_1_medio-1)/std_1_Fs 
Ind_conf_1=np.where(Ind_conf_1<0,-9999,Ind_conf_1)
Ind_conf_1=np.where(Ind_conf_1==np.math.inf,9999,Ind_conf_1)
np.unique(Ind_conf_1[~np.isnan(Ind_conf_1)])
len(np.where(Ind_conf_1==math.inf)[0])
np.nanmax(Ind_conf_1)

file=rio.open(f'{ruta_in}/dem.asc')
raster=file.read(1)

meta=file.profile
file_transform = meta['transform']
file_crs = meta['crs']

with rio.open(ruta_in+'/Resultados/Ind_conf_Sat_30min.tif', 'w', 
                driver='Gtiff',height=raster.shape[0],width=raster.shape[1],count=1,
                dtype='float64',nodata=-999,crs=file_crs ,transform=file_transform ) as dst:
      dst.write(Ind_conf_1,1);


FS_1_medio_sin=np.where(FS_1_medio==10,np.nan,FS_1_medio)
plt.hist(FS_1_medio_sin[~np.isnan(FS_1_medio_sin)], bins=100)
plt.axvline(x=1, color='red', ls='--')
plt.xlabel('Factor de seguridad')
plt.ylabel('Frecuencia')
plt.tight_layout()
plt.savefig(ruta_in+'/Resultados/FS_medio_Sat_30min.jpg')



ind = np.arange(4)  
width = 0.35    
porc=pd.DataFrame([np.nanmean(VF_1_C/VF_1),np.nanmean(VF_1_phi/VF_1),np.nanmean(VF_1_uws/VF_1), np.nanmean(VF_1_ks/VF_1)])
porc.plot.bar(legend=False)
plt.ylabel('Porcentaje de la varianza')
plt.xticks(ind, ('Cohesión', 'Ángulo de\nfricción', 'Peso\nEspecífico', 'Ks'), rotation=0)
plt.tight_layout()
plt.savefig(ruta_in+'/Resultados/Var_FS_Sat_30min.jpg')


"Número de celdas que fallan"
Falla=[(FS_1_medio<1).sum(),
        (FS_1_cohe<1).sum(),
        (FS_1_phi<1).sum(),
        (FS_1_uws<1).sum(),
        (FS_1_ks<1).sum()]
Falla=pd.DataFrame(Falla)
Falla.index=['FS_Medio','FS_C','FS_Phi','FS_Uws','FS_Ks']
print(Falla)

