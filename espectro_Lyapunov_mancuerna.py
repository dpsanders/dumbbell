# -*- coding: utf-8 -*-
"""
Created on Sat May 24 16:09:58 2014

@author: rafael
"""

f1 = open('exponents.dat', 'w')


from Gram_Schmidt import GS 
import Clase_Mancuerna as MC

from time import time as Time

import numpy as np
import matplotlib.pyplot as plt
import lyapunov_mancuerna as lym


from numpy import pi, log

la = np.linalg
ar = np.array
rd = np.random.rand


rcm = [1.,1.] # posición inicial del centro de masa (CM)

l0 = 0.2 # la mancuerna mide 2l

m0 = 1. # y masa TOTAL m0

vcm = [2.1, 0.6] # velocidad inicial del CM

theta = pi/6.  #inclinación inicial de la mancuerna 30°

omega = 2. # velocidad angular inicial

#posiciones de las paredes
w1 = [-3.,-2.]
w2 = [3.,2.]

manc = MC.Mancuerna()


manc.reset()

# vectores iniciales de desplazamiento ALEATORIOS; la primera entrada son los vectores ya ortogonalizados
# y la segunda ortoNORMalizados

dispvecs, normdispvecs = GS(rd(6,6))  

dispvecs, normdispvecs

N=5000  ### número de colisiones que se quieren

col_count = 0.

exponents_col = np.zeros((N, 6))
exponents_time = np.zeros((N, 6))
lengths = np.zeros((N, 6))

col_times = []

times = []

t_T = 0.


start_time = Time()

for k in xrange(N): # se calculan los 6 exponentes de la k-ésima colisión, de las N en total

    col_count +=1 # para contar el número de colisiones
    
    # la primera entrada es el vector en el espacio fase JUSTO ANTES de la colisión
    # la segunda es el vector transformado por el mapeo de la colisión (DESPUÉS de la colisión)
    # la tercera y cuarta especifican la orientación de la pared de colisión y qué partícula chocó
    # la última es el tiempo de colisión

    Gbef, Gaft, wall, part, tc = manc.collision()
    
    t_T += tc #tiempo total que lleva el sistema en evolución
    
    times.append(t_T) # lista de tiempos de colisión
    
    col_times.append(tc) # lista de los tiempos (intervalos) ENTRE colisiones

    # calcular la evolución de los vectores de desplazamiento hasta el momento de la colisión
    dispvecsbef = []

    for v in normdispvecs:
        v = lym.evol_Gamma(v, tc, l0, m0)
        dispvecsbef.append(v)

    dispvecsbef = ar(dispvecsbef)


    # calcular cómo se transforman los vectores de desplazamiento con la colsión

    dispvecsaft = []
    for dG in dispvecsbef:
        dispvecsaft.append( lym.collision_displacement_vector_map(Gbef, dG, part, wall, l0, m0) )
    
    dispvecsaft = ar(dispvecsaft)


    # ortogonalizar los vectores de desplazamiento, ya transformados;
    # calcular sus normas y LUEGO ortonormalizarlos

    orthovecs, normdispvecs = GS(dispvecsaft)
    
    
    for j, vv in enumerate(orthovecs):        
             
        # suma de los logaritmos hasta la k-esima colisión
        if k>=1:
            lengths[k,j] = lengths[k-1,j] + log(la.norm(vv))    
        
        elif k<1:
           lengths[0,j] = log(la.norm(vv))
        
        exponents_col[k,j] = lengths[k,j]/(k+1.)
        exponents_time[k,j] = lengths[k,j]/(1.*t_T)
                                      
    
    f1.write(str(t_T)+' '+str(tc)+' ')
    for j in range(6):
        f1.write(str(exponents_col[k,j])+' '+str(exponents_time[k,j])+' ')
    f1.write(str(manc.cc)+'\n')
        
    x1, y1 = manc.part(1)[0:2]
    x2, y2 = manc.part(2)[0:2]    
    
    # para monitorear qué tan rápido va el código
    if ((k+1)%500) ==0. :
        
        # analizar la evolución de la suma de los exponentes    
        suma1 = 0.
        suma2 = 0.
        for i in range(6):
            suma1 += exponents_time[k, i]
            suma2 += exponents_col[k,i]
        print 'suma de exponentes; '+str(manc.cc)+' colisiones'
        print suma1, suma2, 1.*suma1/suma2
        print '***'
    
        
        current_time = Time()
        minutes = str(int((current_time-start_time)/60.))
        secs = int(round((current_time-start_time)%60))
        
        if secs<10:
            secs = '0'+str(secs)
        else:
            secs=str(secs)
        
        
        print k+1, manc.x, manc.y, manc.cc
        print x1, y1, x2, y2
        print 'elapsed time = '+minutes+':'+secs
        print '\n'
    
    # monitorear si la partícula permanece DENTRO del billar
    if abs(x1)>w2[0] or abs(x2)>w2[0] or abs(y1)>w2[1] or abs(y2)>w2[1] :
        current_time = Time()
        minutes = int((current_time-start_time)/60.)
        secs = round((current_time-start_time)%60)
        
        print k+1, x1, y1, x2, y2, manc.cc
        print 'elapsed time = '+str(minutes)+':'+str(secs)
        print 'ya valió'
        print '\n'
        break
    
print '\n'

print manc

print '\n'

print 'tiempo de movimiento de la mancuerna = ', t_T

current_time = Time()
minutes = str(int((current_time-start_time)/60.))
secs = int(round((current_time-start_time)%60))
        
if secs<10:
     secs = '0'+str(secs)
else:
     secs=str(secs)
print 'elapsed time = '+minutes+':'+secs

suma1 = 0.
suma2 = 0.
n = manc.cc-1
for i in range(6):
    suma1 += exponents_time[n, i]
    suma2 += exponents_col[n,i]


print'\n*****\n'
print 'Espectro de Lyapunov'
for i in range(6):
    print exponents_time[-1,i], exponents_col[-1,i]


print'\n*****\n'
print 'Suma de los exponentes'
print suma1, suma2, suma1/suma2


f1.close()

print'\n*****\n'
print 'Done!'
print'\n*****\n'
tcut=0
tf = k

colors = ['cyan', 'blue', 'red', 'green', 'purple', 'yellow']

plt.figure(figsize=(15,15))

for i in range(6):
    plt.plot(exponents_col[tcut:tf, i], color=colors[i])
    plt.plot(exponents_time[tcut:tf, i], color=colors[i], ls='--')

plt.ylim(-2., 2.);
plt.savefig('Espectro-Lyapunov-Mancuerna.pdf')



