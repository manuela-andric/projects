# -*- coding: utf-8 -*-
"""
Created on Sun May  7 17:24:25 2017

@author: luka
"""
import numpy as np
from sympy import Symbol, Matrix, sin, cos, nsolve, plot

#generiranje matrice B
def generirajMatricuB(W,kc,Ic,km,Im):
    B=1
    for i in reversed(range(len(W))):
        B*=Matrix([
        [1, kc/km[i]*sin(W[i]), -Ic/Im[i]*(kc/km[i])**2*(1-cos(W[i])), -Ic/Im[i]*(kc/km[i])**3*(W[i]-sin(W[i]))],
        [0, cos(W[i]), -(Ic/Im[i])*(kc/km[i])*sin(W[i]), -Ic/Im[i]*(kc/km[i])**2*(1-cos(W[i]))],
        [0, (Im[i]/Ic)*(km[i]/kc)*sin(W[i]), cos(W[i]), kc/km[i]*sin(W[i])],
        [0, 0, 0, 1]
        ])
    
    return B        

#deklaracija
P = Symbol('P')
w = Symbol('w')
Io = Symbol('Io')
w1 = 0.2721*Symbol('w')
w2 = 0.2108 * Symbol('w')
w3 = 0.1782*Symbol('w')
    
#vrijednosti po segmentima
Im = [1.5*Io, 2.5*Io, 3.5*Io]
km = [np.sqrt(1/1.5), np.sqrt(1/2.5), np.sqrt(1/3.5)]
kc = km[0] #definirarti krutost prvog segmenta
Ic = Im[0] #definirati moment tromosti prvog segmenta
W=[w1,w2,w3] #vektor ovisan o broju segmenata

#stvori matricu B    
B = generirajMatricuB(W,kc,Ic,km,Im)
    
#jednadzba u matrici B na B22
B22 = B.row(1)[1]
  
#rijesi numerickom metodom
result_numeric = nsolve(B22, w, 1)
print("Rezultat numerickom metodom: w = {}".format(result_numeric))
Sila=result_numeric**2
print('Kritična sila: Pcr= {}*EI/L**2'.format(Sila))
   
#graf jednadzbe B22
plot(B22).save('./plots/B22.png')



