# -*- coding: utf-8 -*-
"""
Created on Sat May 24 16:08:12 2014

@author: rafael
"""

#from Gram_Schmidt import GS 
#import Clase_Mancuerna as MC

import numpy as np

from numpy import sin, cos, pi, tan, absolute, sqrt, dot, log

la = np.linalg
ar = np.array
rd = np.random.rand

# la mancuerna mide 2*l0
l0=0.2
# y tiene una masa m0
m0=1.
def coefficients(theta, part, wall, l=l0):
    #choque con paredes horizontales
    if wall=='h':
        if part==1:
            csca = 1./(cos(theta)) # choca la partícula 1
        else:
            csca = -1./(cos(theta))  # choca la partícula 2
    # choque con paredes verticales
    elif wall=='v':
        if part==1:
            csca = -1./(sin(theta))  # choca la partícula 1
        else:
            csca = 1./(sin(theta))  # choca la partícula 2
            
    a = 2./(1. + csca**2) - 1. #adimensional
    b = 2.*l*csca/(1. + csca**2)  # unidades de [l]
    c = b/(l**2)  # unidades de [1/l]
    
    return a,b,c


def collision_map(gamma, part, wall, l=l0):
    # con esta función se define el mapeo que describe una colisión de la partícula 'part' 
    # con una pared, 'wall'
    # gamma es un vector en el espacio fase, compuesto por entradas gamma = (qq, pp)
    # donde qq y pp son los vectores de posición y momentos dados a su vez por
    # qq = (x_CM, y_CM, theta); y pp = (p_x, p_y, L)
    qq = gamma[0:3]
    x, y, theta = qq
    pp = gamma[3:]
    px, py, L = pp
    
    # La acción del mapeo se expresa en términos de la matriz B, definida más abajo.
    #Los coeficientes a,b y c se obtienen con la función coefficients
    
    a, b, c = coefficients(theta, part, wall, l)
    
    # ahora hay que definir el vector normal dependiendo de si la colisión es
    #con las paredes verticales u horizontales
    if wall=='h':
        n1 = 0.; n2 = 1.
    elif wall=='v':
        n1=1.; n2=0.
    
    # Definición de la matriz B
    B = ar([[n2**2 + a*n1**2, n1*n2*(a-1.), -c*n1], 
            [n1*n2*(a-1), n1**2 + a*n2**2, -c*n2],
            [-b*n1, -b*n2, -a]])
    
    pp = dot(B, pp)
    
    return np.hstack([qq, pp])





def collision_map_derivative(gamma, vec, part, wall, l=l0):
    # Esta función calcula la derivada del mapeo (evaluado en Gamma)
    # La función regresa el producto punto de la derivada del mapeo con el vector vec.
    # Para usarla es necesario definir la matriz A descrita abajo 
    # y esto involucra encontrar el ángulo alpha en términos de theta,
    # su derivada y las derivadas de los coeficientes
    
    qq = gamma[0:3]; x, y, theta = qq
    pp = gamma[3:]; px, py ,L = pp
    
    vecqq = vec[0:3]
    vecpp = vec[3:]
    
    if part==1:
        deralpha = -1.
        if wall=='h':
            alpha = pi/2. - theta
        elif wall=='v':
            alpha = - theta
    
    if part==2:
        deralpha = 1.
        if wall=='h':
            alpha = theta - pi/2.
        elif wall=='v':
            alpha = theta
    
    a, b, c = coefficients(theta, part, wall, l)
    
    # definición de las derivadas de los coefficientes. Se calcularon usando Mathematica
    dera = deralpha*(8.*sin(2*alpha))/((cos(2*alpha)-3.)**2)
    derb = deralpha*(8.*l*(cos(alpha)**3))/((cos(2*alpha)-3.)**2)
    derc = derb/(l**2)
    
    
    if wall=='h':
        n1 = 0.; n2 = 1.
    elif wall=='v':
        n1=1.; n2=0.
    
    # Definición de la matriz B
    B = ar([[n2**2 + a*(n1**2), n1*n2*(a-1.), -c*n1], 
            [n1*n2*(a-1), n1**2 + a*(n2**2), -c*n2],
            [-b*n1, -b*n2, -a]])
    
    # Definición de la matriz A    
    A = ar([[0., 0., n1*(dera*(px*n1*+py*n2) - L*derc)],
            [0., 0., n2*(dera*(n1*px+ n2*py) - L*derc)],
            [0., 0., -L*dera - derb*(n1*px + n2*py)]])
            
    # transformación del vector dpp después de la colisión
    vecpp = dot(A, vecqq) + dot(B, vecpp)
    
    return np.hstack([vecqq, vecpp])

# esta función es la evolución del vector Gamma en el tiempo;
# i.e. es la parte continua de la dinámica, mientras la mancuerna se mueve dentro de la mesa
def evol_function(Gamma, l=l0, m =m0):
    I = m*(l**2)
    px, py, L = Gamma[3:]; vx, vy = 1.*px/m, 1.*py/m
    omega = 1.*L/I
    
    return ar([vx, vy, omega, 0, 0, 0])


# con esta función se calcula el tiempo de desfase entre la colisión de una
# trayectoria de referencia y una satélite.
def delta_tau(deltaGamma, Gamma, part, wall, l=l0, m=m0):
    I = m*(l**2)
    dx, dy, dtheta  = deltaGamma[0:3]
    px, py, L = Gamma[3:]; vx, vy = 1.*px/m, 1.*py/m
    theta = Gamma[2]; omega = 1.*L/I
    
    if part==1:
        dxi = dx - l*sin(theta)*dtheta
        dyi = dy + l*cos(theta)*dtheta
        vxi = vx - l*omega*sin(theta)
        vyi = vy + l*omega*cos(theta)
        
    if part ==2:
        dxi = dx + l*sin(theta)*dtheta
        dyi = dy - l*cos(theta)*dtheta
        vxi = vx + l*omega*sin(theta)
        vyi = vy - l*omega*cos(theta)
    
    dqqi = ar([dxi, dyi])
    vvi = ar([vxi, vyi])
    
    if wall =='v': # colisión con paredes verticales
        n = ar([1, 0])
    if wall =='h': #colisión con paredes horizontales
        n = ar([0, 1])
        
    return  -1.*dot(dqqi, n)/dot(vvi, n)


def collision_displacement_vector_map(Gamma, deltaGamma, part, wall, l=l0, m=m0):
    
    Gp = collision_map(Gamma, part, wall) # vector después de la colisión
    
    DMdG = collision_map_derivative(Gamma, deltaGamma, part, wall, l) # primer término de dG'
    FG = evol_function(Gamma, l, m)  # F(Gamma)
    DMF = collision_map_derivative(Gamma, FG, part, wall, l) # primer término del segundo sumando
    FGp = evol_function(Gp, l, m) # F(Gamma')
    
    dT = delta_tau(deltaGamma, Gamma, part, wall, l) #delta tau
    
    return DMdG + (DMF - FGp)*dT


#==============================================================================
# Las siguientes dos funciones describen cómo evoluciona Gamma y deltaGamma, respectivamente.
# Por lo tanto, no son sino la integral de la función 'evol_function' definida más arriba.
# Dado que en este caso el movimiento es bastante sencillo --MRU superpuesto a un movimiento circular
# uniforme-- las funciones son muy fáciles de definir e implementar.
#==============================================================================

def evol_Gamma(Gamma0, t, l=l0, m=m0):
    I = m*(l**2)
    qq = Gamma0[0:3]
    pp = Gamma0[3:]
    vx, vy = (1./m)*pp[0:2]
    omega = 1.*pp[-1]/I
    qq += t* ar([vx, vy, omega])
    
    return np.hstack([qq, pp])


def evol_deltaGamma(deltaGamma0, t, l=l0, m=m0):
    I = m*(l**2)
    dqq = deltaGamma0[0:3]
    dpp = deltaGamma0[3:]
    dvx, dvy = (1./m)*dpp[0:2]
    omega = 1.*dpp[-1]/I
    dqq += t*ar([dvx, dvy, omega])
    
    return np.hstack([dqq, dpp])
