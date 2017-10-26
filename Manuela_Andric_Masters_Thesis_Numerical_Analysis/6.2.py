# -*- coding: utf-8 -*-
"""
Created on Sun May  7 17:24:25 2017

@author: Manuela Andric
"""

from sympy import Symbol, Matrix, sin, cos, linsolve
from pprint import pprint
import numpy as np

#generiraj matrice A
def stvoriMatriceA(W,kc,Ic,km,Im):
    A = []
    for i in range(len(W)):
        temp=Matrix([
        [1, kc/km[i]*sin(W[i]), -Ic/Im[i]*(kc/km[i])**2*(1-cos(W[i])), -Ic/Im[i]*(kc/km[i])**3*(W[i]-sin(W[i]))],
        [0, cos(W[i]), -(Ic/Im[i])*(kc/km[i])*sin(W[i]), -Ic/Im[i]*(kc/km[i])**2*(1-cos(W[i]))],
        [0, (Im[i]/Ic)*(km[i]/kc)*sin(W[i]), cos(W[i]), kc/km[i]*sin(W[i])],
        [0, 0, 0, 1]
        ])
        A.append(temp)
    
    return A  
    
#deklaracija
E = 2.1*10**8
F = 100
L = 15
Io = 2.2930*10**(-4)
H = 10

#broj presjeka
N = 24

#deklaracija simbola jednadzbe
v = 0
ro = 0
M0 = Symbol("M0")
V0 = Symbol("V0")

#vrijednosti po segmentima
Im = []
km = []
W = []
delta = 1 + (3/N/2)
for i in range(N):
    Im.append(delta*Io)
    km.append(np.sqrt(1/delta)*np.sqrt(F/(E*Io)))
    W.append(km[i]*L/N)
    delta = delta + (3/N)    
kc = km[0]
Ic = Im[0]  

#stvori matrice A
A = []
A = stvoriMatriceA(W,kc,Ic,km,Im)

#jednadzbe po segmentima
S = []
for i in reversed(range(N)):
    temp = i
    mat = A[temp]
    for j in range(N-i):
        if temp != i:
            mat = mat * A[temp]
        temp += 1
    if i != 0:
        S.append(mat * Matrix([v,ro,M0,V0]))
    else:
        S.append(mat * Matrix([v,ro,M0,V0])-Matrix([0,0,0,H]))

#racunanje M0 i V0 te substitucija rezultata
res = linsolve((S[N-1][2],S[N-1][3]),(M0,V0))
(M0, V0) = next(iter(res))
print(M0/kc)

for i in range(N):
    S[i] = S[i].subs({"M0": M0, "V0": V0})
    S[i][0] = S[i][0]/(E*Ic*kc**3)
    S[i][1] = S[i][1]/(E*Ic*kc**2)
    S[i][2] = S[i][2]/kc
    pprint(S[i])

