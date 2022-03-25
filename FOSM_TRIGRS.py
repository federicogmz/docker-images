#%%
import numpy as np
import statistics 
import os

ruta_in='G:/Unidades compartidas/Proyectos(Fede)/Tarso/Amenaza/TRIGRS/'

#Para leer 
with open(ruta_in+'tr_in.txt', "r+") as f:
     old = f.read() # read everything in the file

#Para separar cada fila en una lista
tr_in_file = open(ruta_in+'tr_in.txt', "r")
list_of_lines = tr_in_file.readlines()

n=3 #NUMERO DE ZONAS

# Zona 1: GW, Zona 2: SW, Zona 3: SC

#%%
"C"
variation_C=(100/100)
C_mean_1=5000
C_vector_1=[C_mean_1,C_mean_1+(C_mean_1*variation_C),C_mean_1-(C_mean_1*variation_C)]
C_var_1=statistics.variance(C_vector_1)
C_std_1=statistics.stdev(C_vector_1)

variation_C=(100/100)
C_mean_2=5000
C_vector_2=[C_mean_2,C_mean_2+(C_mean_2*variation_C),C_mean_2-(C_mean_2*variation_C)]
C_var_2=statistics.variance(C_vector_2)
C_std_2=statistics.stdev(C_vector_2)

variation_C=(50/100)
C_mean_3=20000
C_vector_3=[C_mean_3,C_mean_3+(C_mean_3*variation_C),C_mean_3-(C_mean_3*variation_C)]
C_var_3=statistics.variance(C_vector_3)
C_std_3=statistics.stdev(C_vector_3)

"Phi"
variation_phi=(11/100)
phi_mean_1=36
phi_vector_1=[phi_mean_1,phi_mean_1+(phi_mean_1*variation_phi),phi_mean_1-(phi_mean_1*variation_phi)]
phi_var_1=statistics.variance(phi_vector_1)
phi_std_1=statistics.stdev(phi_vector_1)

variation_phi=(12/100)
phi_mean_2=33
phi_vector_2=[phi_mean_2,phi_mean_2+(phi_mean_2*variation_phi),phi_mean_2-(phi_mean_2*variation_phi)]
phi_var_2=statistics.variance(phi_vector_2)
phi_std_2=statistics.stdev(phi_vector_2)

variation_phi=(17/100)
phi_mean_3=30
phi_vector_3=[phi_mean_3,phi_mean_3+(phi_mean_3*variation_phi),phi_mean_3-(phi_mean_3*variation_phi)]
phi_var_3=statistics.variance(phi_vector_3)
phi_std_3=statistics.stdev(phi_vector_3)

"Uws"
variation_uws=(5/100)
uws_mean_1=21000
uws_vector_1=[uws_mean_1,uws_mean_1+(uws_mean_1*variation_uws),uws_mean_1-(uws_mean_1*variation_uws)]
uws_var_1=statistics.variance(uws_vector_1)
uws_std_1=statistics.stdev(uws_vector_1)

variation_uws=(10/100)
uws_mean_2=20500
uws_vector_2=[uws_mean_2,uws_mean_2+(uws_mean_2*variation_uws),uws_mean_2-(uws_mean_2*variation_uws)]
uws_var_2=statistics.variance(uws_vector_2)
uws_std_2=statistics.stdev(uws_vector_2)

variation_uws=(8/100)
uws_mean_3=18500
uws_vector_3=[uws_mean_3,uws_mean_3+(uws_mean_3*variation_uws),uws_mean_3-(uws_mean_3*variation_uws)]
uws_var_3=statistics.variance(uws_vector_3)
uws_std_3=statistics.stdev(uws_vector_3)


"KS"
variation_ks=(100/100)
ks_mean_1=5.0E-03
ks_vector_1=[ks_mean_1,ks_mean_1+(ks_mean_1*variation_ks),ks_mean_1-(ks_mean_1*variation_ks)]
ks_var_1=statistics.variance(ks_vector_1)
ks_std_1=statistics.stdev(ks_vector_1)

ks_mean_2=5.0E-04
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
#%%
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
        list_of_lines[38+(3*(n-1))]='Val_4h\n'  #MODIFICAR MANUALMENTE
        
        # Sobre escribir el archivo txt
        a_file = open(ruta_in+'tr_in.txt', "w")
        a_file.writelines(list_of_lines)
        a_file.close()
          
        # ejecutar trigrs
        os.chdir(ruta_in)
        os.system(f'{ruta_in}TRIGRS.exe')    

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
           os.system(f'{ruta_in}TRIGRS.exe')       