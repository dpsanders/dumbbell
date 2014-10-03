# -*- coding: utf-8 -*-
"""
Created on Thu Jun  5 20:44:08 2014

@author: rafael
"""



import numpy as np

from numpy import cos, sin, pi, tan, arctan, sqrt, absolute

ar = np.array
rd = np.random.rand

rcm = [1.,1.] # posición inicial del centro de masa (CM)

l0 = 0.2 # la mancuerna mide 2l0

vcm = [2.1, 0.6] # velocidad inicial del CM

theta = pi/6.  #inclinación inicial de la mancuerna 30°

omega = 2. # velocidad angular inicial

#posiciones de las paredes
wall1 = [-3.,-2.]
wall2 = [3.,2.]


class Mancuerna(object):
    
    def __init__(self, x=1.0, y=1., vx=2.1, vy=0.6, l=l0, theta=pi/6, omega=2., cc=0):
        self.x, self.y  = x, y
        self.vx, self.vy = vx, vy
        self.l, self.theta, self.omega = l, theta, omega
        self.cc = cc
        #valores iniciales
        self.x0, self.y0 = self.x, self.y
        self.vx0, self.vy0 = self.vx, self.vy
        self.theta0, self.omega0 = self.theta, self.omega
        
    def __repr__(self):
        return "Rcm(%s, %s), Vcm(%s, %s), theta=%s, omega=%s" % \
    (self.x, self.y, self.vx, self.vy, self.theta, self.omega)
    
    
    def reset(self):
        self.x, self.y = self.x0, self.y0
        self.vx, self.vy = self.vx0, self.vy0
        self.theta, self.omega = self.theta0, self.omega0
        self.cc = 0
    
    def mover(self, dt):
        self.x += self.vx*dt
        self.y += self.vy*dt
        self.theta += self.omega*dt
        self.theta = np.mod(self.theta, 2*pi)
    
        # define la posición y velocidad las partículas
    def part(self, p, t=0):
        x, y = self.x, self.y
        vx, vy = self.vx, self.vy
        theta, omega, l = self.theta, self.omega, self.l
        
        thtm = np.mod(theta + omega*t, 2*pi)        
        
        if p==1:
            x1 = x + vx*t + l*cos(thtm)
            y1 = y + vy*t + l*sin(thtm) 
        
            vx1 = vx - l*omega*sin(thtm)
            vy1 = vy + l*omega*cos(thtm)
            return ar([x1, y1, vx1, vy1])
        
        elif p==2:
            x2 = x + vx*t - l*cos(thtm)
            y2 = y + vy*t - l*sin(thtm) 
        
            vx2 = vx + l*omega*sin(thtm)
            vy2 = vy - l*omega*cos(thtm)
        
            return ar([x2, y2, vx2, vy2])
        else:
            print "No hay tantas particulas wey"
    
#==============================================================================
#     Obtener una aproximación al tiempo de colisión
#==============================================================================
    def t_approx(self, w1, w2, direc):
        # 'direc' especifica se se quiere estimar el tiempo de colisión
        # en la dirección 'x': 'direc=0'; o en la dirección 'y'
        #  'direc = 1'.
        
       # posición de la partícula 1 y 2 en la dirección especificada
       r1 = self.part(1)[direc]
       r2 = self.part(2)[direc]
                
       t_elapsed = 0.
       dt = 0.005
        
       while w1 < r1 < w2 and w1 < r2 < w2:
            t_elapsed += dt
            
            r1 = self.part(1, t_elapsed)[direc]
            r2 = self.part(2, t_elapsed)[direc]      
        
       return t_elapsed - dt, dt
    
    
    ####
    #  colisión  con CUALQUIERA de las paredes VERTICALES
    ####
    
    def vcol(self, w1=wall1[0], w2=wall2[0]):
        
        vx, vy = self.vx, self.vy
        theta, omega, l = self.theta, self.omega, self.l
        
        if vx ==0.:
            return ar([-1, "No hay velocidad en x"])
        
        else:
        
            err = 1E-10
            tbig = 50000.
            t11 = [tbig, tbig]; t12 = [tbig, tbig]
            t21 = [tbig, tbig]; t22 = [tbig, tbig]
            
            # tiempo aproximado de colisión en la dirección x ('0')
            tapprox, dt = self.t_approx(w1, w2, 0)
            
            
            n=0
            tx1w1, tx1w2 = tapprox, tapprox
            x1w1, vx1w1 = self.part(1, tx1w1)[[0, 2]]
            x1w2, vx1w2 = self.part(1, tx1w2)[[0, 2]]
        
            tx2w1, tx2w2 = tapprox, tapprox
            x2w1, vx2w1 = self.part(2, tx2w1)[[0, 2]]
            x2w2, vx2w2 = self.part(2, tx2w2)[[0, 2]]
                
            # condición de colisión: una de las partículas está a menos de 'err' de alguna pared 
            # y el tiempo para dicha colisión está en el intervalo (tapprox, tapprox + dt)
            
            cond11 = ( abs(x1w1-w1)<err and tapprox < tx1w1 < tapprox + dt )
            cond12 = ( abs(x1w2-w2)<err and tapprox < tx1w2 < tapprox + dt )
            cond21 = ( abs(x2w1-w1)<err and tapprox < tx2w1 < tapprox + dt )
            cond22 = ( abs(x2w2-w2)<err and tapprox < tx2w2 < tapprox + dt )
            
            cond= ( cond11 or cond12 or cond21 or cond22 )
            
            # aplicar Newton mientras NO se satisfaga 'cond'
            while not(cond):
                tx1w1 = tx1w1 - (x1w1-w1)/vx1w1
                tx1w2 = tx1w2 - (x1w2-w2)/vx1w2
                x1w1, vx1w1 = self.part(1, tx1w1)[[0, 2]]
                x1w2, vx1w2 = self.part(1, tx1w2)[[0, 2]]
            
                tx2w1 = tx2w1 - (x2w1-w1)/vx2w1
                tx2w2 = tx2w2 - (x2w2-w2)/vx2w2
                x2w1, vx2w1 = self.part(2, tx2w1)[[0, 2]]
                x2w2, vx2w2 = self.part(2, tx2w2)[[0, 2]]
                
                #actualizar condición
                d11 = abs(x1w1-w1); d12 = abs(x1w2-w2);
                d21 = abs(x2w1-w1); d22 = abs(x2w2-w2);
                cond11 = ( d11<err and tapprox < tx1w1 < tapprox + dt )
                cond12 = ( d12<err and tapprox < tx1w2 < tapprox + dt )
                cond21 = ( d21<err and tapprox < tx2w1 < tapprox + dt )
                cond22 = ( d22<err and tapprox < tx2w2 < tapprox + dt )
            
                cond= ( cond11 or cond12 or cond21 or cond22 )
                        
                n+=1
                if n>300:
                    tx1w1 =-1.; tx1w2 =-1.; tx2w1 =-1.; tx2w2 =-1.;
                    break
                    
            # elegir cuál es el tiempo correcto
            if cond11:
                t11.append(tx1w1)
            
            if cond12:
                t12.append(tx1w2)
            
            if cond21:
                t21.append(tx2w1)
            
            if cond22:
                t22.append(tx2w2) 
                
            # elegir el tiempo mínimo
            tx1w1 = np.min([t for t in t11 if t >0]); tx1w2 = np.min([t for t in t12 if t >0])
            tx2w1 = np.min([t for t in t21 if t >0]); tx2w2 = np.min([t for t in t22 if t >0])
                    
            times = ar([t for t in [tx1w1, tx1w2, tx2w1, tx2w2] if t>0])
            tx = np.min(times)
            
            tht = np.mod(theta + omega*tx, 2*pi)
            u = vx
            oi = omega
            
            if tx == tx1w1 or tx == tx1w2: #colisionó la partícula 1
                csca = -1./(sin(tht))
                
                # transformaciones para las velocidades después de la colisión
                vx = -u + 2.*( u- l*oi*csca )/(1. + csca**2)
                omega = oi -2.*( l*oi + u*csca )/(l*( 1. + csca**2))
                    
                return ar([tx, vx, vy, omega, 1, tx1w1, tx1w2, tx2w1, tx2w2])
                            
            if tx== tx2w1 or tx == tx2w2: # colisionó la partícula 2
                csca = 1./(sin(tht))
                
                # transformaciones para las velocidades después de la colisión
                vx = -u + 2.*( u- l*oi*csca )/(1. + csca**2)
                omega = oi -2.*( l*oi + u*csca )/(l*( 1. + csca**2))
                    
                return ar([tx, vx, vy, omega, 2, tx1w1, tx1w2, tx2w1, tx2w2])
        
        
        
        
        
    ####
    #  colisión  con CUALQUIERA de las paredes HORIZONTALES
    ####
    
    def hcol(self, w1=wall1[1], w2=wall2[1]):
        
        vx, vy = self.vx, self.vy
        theta, omega, l = self.theta, self.omega, self.l
        
        if vy ==0.:
            return ar([-1, "No hay velocidad en y"])   
        
        err = 1E-10
        tbig = 50000.
        t11 = [tbig, tbig]; t12 = [tbig, tbig]
        t21 = [tbig, tbig]; t22 = [tbig, tbig]
        
        # tiempo aproximado de colisión en la dirección y ('1')
        tapprox, dt = self.t_approx(w1, w2, 1) 
        
            
        
        n=0
        ty1w1, ty1w2 = tapprox, tapprox
        y1w1, vy1w1 = self.part(1, ty1w1)[[1, 3]]
        y1w2, vy1w2 = self.part(1, ty1w2)[[1, 3]]
    
        ty2w1, ty2w2 = tapprox, tapprox
        y2w1, vy2w1 = self.part(2, ty2w1)[[1, 3]]
        y2w2, vy2w2 = self.part(2, ty2w2)[[1, 3]]
        
        # condición de colisión: una de las partículas está a menos de 'err' de alguna pared 
        # y el tiempo para dicha colisión está en el intervalo (tapprox, tapprox + dt)
        
        cond11 = ( abs(y1w1-w1)<err and tapprox < ty1w1 < tapprox + dt )
        cond12 = ( abs(y1w2-w2)<err and tapprox < ty1w2 < tapprox + dt )
        cond21 = ( abs(y2w1-w1)<err and tapprox < ty2w1 < tapprox + dt )
        cond22 = ( abs(y2w2-w2)<err and tapprox < ty2w2 < tapprox + dt )
        
        cond= ( cond11 or cond12 or cond21 or cond22 )    
        
        # aplicar Newton mientras NO se satisfaga la 'cond'
            
        while not(cond):
            ty1w1 = ty1w1 - (y1w1-w1)/vy1w1
            ty1w2 = ty1w2 - (y1w2-w2)/vy1w2
            y1w1, vy1w1 = self.part(1, ty1w1)[[1, 3]]
            y1w2, vy1w2 = self.part(1, ty1w2)[[1, 3]]
        
            ty2w1 = ty2w1 - (y2w1-w1)/vy2w1
            ty2w2 = ty2w2 - (y2w2-w2)/vy2w2
            y2w1, vy2w1 = self.part(2, ty2w1)[[1, 3]]
            y2w2, vy2w2 = self.part(2,ty2w2)[[1, 3]]
            
            #actualizar condición
            d11 = abs(y1w1-w1); d12 = abs(y1w2-w2);
            d21 = abs(y2w1-w1); d22 = abs(y2w2-w2);
            cond11 = ( d11<err and tapprox < ty1w1 < tapprox + dt )
            cond12 = ( d12<err and tapprox < ty1w2 < tapprox + dt )
            cond21 = ( d21<err and tapprox < ty2w1 < tapprox + dt )
            cond22 = ( d22<err and tapprox < ty2w2 < tapprox + dt )
        
            cond= ( cond11 or cond12 or cond21 or cond22 )
                       
                    
            n+=1
            if n>300:
                ty1w1 =-1.; ty1w2 =-1.; ty2w1 =-1.; ty2w2 =-1.;
                break
                
        # elegir cuál es el tiempo correcto
        if cond11:
            t11.append(ty1w1); t12.append(tbig); t21.append(tbig); t22.append(tbig)
        
        if cond12:
            t11.append(tbig); t12.append(ty1w2); t21.append(tbig); t22.append(tbig)
        
        if cond21:
            t11.append(tbig); t12.append(tbig); t21.append(ty2w1); t22.append(tbig)
        
        if cond22:
            t11.append(tbig); t12.append(tbig); t21.append(tbig); t22.append(ty2w2)
        
        # elegir el tiempo mínimo                           
        ty1w1 = np.min([t for t in t11 if t >0]); ty1w2 = np.min([t for t in t12 if t >0])
        ty2w1 = np.min([t for t in t21 if t >0]); ty2w2 = np.min([t for t in t22 if t >0])        
                
        times = ar([t for t in [ty1w1, ty1w2, ty2w1, ty2w2] if t>0])
        ty = np.min(times)
                
        
        tht = np.mod(theta + omega*ty, 2*pi)
        u = vy
        oi = omega
        
                
        if ty == ty1w1 or ty == ty1w2: #colisionó la partícula 1
            csca = 1./(cos(tht))
            
            # transformaciones para las velocidades después de la colisión
            vy = -u + 2.*( u- l*oi*csca )/(1. + csca**2)
            omega = oi -2.*( l*oi + u*csca )/(l*( 1. + csca**2))
                
            return ar([ty, vx, vy, omega, 1, ty1w1, ty1w2, ty2w1, ty2w2])
                                    
        if ty== ty2w1 or ty == ty2w2: # colisionó la partícula 2
            csca = -1./(cos(tht))
            
            # transformaciones para las velocidades después de la colisión
            vy = -u + 2.*( u- l*oi*csca )/(1. + csca**2)
            omega = (oi -2.*( l*oi + u*csca )/(l*( 1. + csca**2)))
                                
            return ar([ty, vx, vy, omega, 2, ty1w1, ty1w2, ty2w1, ty2w2])
        
         
        
     # se actualiza el valor de las PARTÍCULAS, velocidad, theta y omega después de la colisión   
    def collision_parts(self, w1x=wall1[0], w2x=wall2[0], w1y=wall1[1], w2y=wall2[1]): 
        self.cc += 1
        tx = self.vcol(w1x, w2x)[0]
        ty = self.hcol(w1y, w2y)[0]
               
        if tx <= ty and tx>0.: # colisión con paredes verticales
            tc, vxn, vyn, omen, p = self.vcol(w1x, w2x)[0:5]
            wall = 'x'
        
        if ty<tx and ty>0.:
            tc, vxn, vyn, omen, p = self.hcol(w1y, w2y)[0:5]
            wall = 'y'
        
        self.mover(tc)
        self.vx, self.vy, self.omega = vxn, vyn, omen   
        x1, y1 = self.part(1)[0:2]
        x2, y2 = self.part(2)[0:2]
        
         
        vM = sqrt(self.vx**2 + self.vy**2 + (self.l*self.omega)**2)
       
        
        d1 = 0.001*abs(self.x -w1x)/vM
        d2 = 0.001*abs(self.x -w2x)/vM
        d3 = 0.001*abs(self.y -w1y)/vM
        d4 = 0.001*abs(self.y -w2y)/vM
                
        dt = min(d1, d2, d3, d4, 0.0001)                   
        
        self.mover(dt)
        
        
                        
        return x1, y1, x2, y2, tc, "p=%i" %p, wall, self.cc
        
            
        
                    
        
                       
    def collision(self,  w1x=wall1[0], w2x=wall2[0], w1y=wall1[1], w2y=wall2[1]):
        self.cc += 1   
            #########
            ## Con esta función se obtiene el punto en el espacio fase justo antes de la colisión,
            ## justo después de la colisión, la partícula que choca y la orientación de la pared con la que choca.
        
        tx = self.vcol(w1x, w2x)[0]
        ty = self.hcol(w1y, w2y)[0]
        xi, yi, thetai = self.x, self.y, self.theta
        vxi, vyi, omegai = self.vx, self.vy, self.omega
        li = self.l
        
        if tx <= ty and tx>0.: # colisión con paredes verticales
            tc, vxn, vyn, omen, p = self.vcol(w1x, w2x)[0:5]
            wall = 'v'
            
        if ty <= tx and ty>0.: # colisión con paredes horizontales
            tc, vxn, vyn, omen, p = self.hcol(w1y, w2y)[0:5]
            wall = 'h'   
        
        ### valores de Gamma ANTES de la colisión
        self.mover(tc)
        Gammabefore = ar([xi, yi, thetai, vxi, vyi, omegai*(li**2) ])
        
        ### valores de Gamma DESPUÉS de la colisión
        self.vx, self.vy, self.omega = vxn, vyn, omen           
        Gammaafter = ar([self.x, self.y, self.theta, vxn, vyn, omen*(li**2)])
           
        
        vM = sqrt(self.vx**2 + self.vy**2 + (self.l*self.omega)**2)
        
        
        d1 = 0.001*abs(self.x -w1x)/vM
        d2 = 0.001*abs(self.x -w2x)/vM
        d3 = 0.001*abs(self.y -w1y)/vM
        d4 = 0.001*abs(self.y -w2y)/vM
                
        dt = min(d1, d2, d3, d4, 0.0001)           
        
        
        self.mover(dt)

      
        return Gammabefore, Gammaafter, wall, p, tc
    
        
        
    
    #calcular tiempo de colisión
    def ct(self, w1x=wall1[0], w2x=wall2[0], w1y=wall1[1], w2y=wall2[1]): 
        tcx = self.vcol(w1x, w2x)[0]
        tcy = self.hcol(w1y, w2y)[0]
        return min(tcx, tcy)
    
    # obtener las posiciones de las 2 masas
    def positions(self, T, dt=0.05, w1x=wall1[0], w2x=wall2[0], w1y=wall1[1], w2y=wall2[1]): 
        tc = self.ct(w1x, w2x, w1y, w2y)
        t=0
        x1, y1 = self.part(1, t)[0:2]
        x2, y2 = self.part(2, t)[0:2]
        posm1 = ar([x1, y1])
        posm2 = ar([x2, y2])
        times = ar([t])
        
        # no hay colisiones en tiempos menores a T
        if T < tc:
            
            n = int(T/dt)
            for t in np.linspace(dt, T, n):
                x1, y1 = self.part(1, t)[0:2]
                x2, y2 = self.part(2, t)[0:2]
                posm1 = np.vstack((posm1, ar([x1,y1])))
                posm2 = np.vstack((posm2, ar([x2,y2])))
                times = np.vstack((times, ar([t])))
                
            return times, posm1, posm2
        
         # T es mayor que el tiempo de colisión
        else:
            while t < T:
                nt = 0.
                while nt < tc and t<T:
                    x1, y1 = self.part(1, nt)[0:2]
                    x2, y2 = self.part(2, nt)[0:2]
                    posm1 = np.vstack((posm1, ar([x1,y1])))
                    posm2 = np.vstack((posm2, ar([x2,y2])))
                    times = np.vstack((times, ar([t])))
                    t+= dt
                    nt+= dt
                    #print t
                self.collision_parts(w1x, w2x, w1y, w2y)
                tc = self.ct(w1x, w2x, w1y, w2y)
            
            return times, posm1, posm2
            
        