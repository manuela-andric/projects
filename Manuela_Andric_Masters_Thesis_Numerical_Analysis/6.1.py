# -*- coding: utf-8 -*-
"""
Created on Sat May 13 17:38:47 2017

@author: Manuela Andric
"""

from sympy import Symbol, Matrix, solve_linear_system  
    
#generiraj matricu (dimenzije min>6)
def stvoriMatricu(dimenzije, lamb, F, W, c):
    pocetna = lamb
    M=Matrix.zeros(dimenzije, dimenzije+1)
    delta = pocetna-1  
    pozicija = dimenzije-4
    
    if(dimenzije<6):
        print("Dimenzije min 6*6!")        
        return 0
    
    for i in range(dimenzije):
        #iznimke pravila za prvi red
        if(i==0):
            M[i,dimenzije-3] = pocetna
            M[i,dimenzije-2] = F/c-2*pocetna
            M[i,dimenzije-1] = pocetna-F/c
            M[i,dimenzije] = (i+1)*W*lamb/c
        #iznimke pravila za predzadnji red
        elif(i==dimenzije-2):
            M[i,0] = F/c-(delta*(dimenzije-1)+1)*2
            M[i,1] = (delta*(dimenzije-1)+1)
            M[i,dimenzije-1] = -F/c
            M[i,dimenzije] = (i+1)*W*lamb/c
        #iznimke pravila za zadnji red
        elif(i==dimenzije-1):
            M[i,0] = 8
            M[i,dimenzije-1] = -F/c
            M[i,dimenzije] = (i+1)*W*lamb/c
        #ostali redovi
        else:
            M[i, pozicija-i+1] = delta*(i+1)+1
            M[i, pozicija-i+2] = F/c - (delta*(i+1)+1)*2
            M[i, pozicija-i+3] = delta*(i+1)+1
            M[i, dimenzije-1] = -F/c   
            M[i,dimenzije] = (i+1)*W*lamb/c
        
    return M    


#deklaracija pocetnih uvjeta
F = 100
W = 10
E = 2.1*10**8
L = 15
Io = 2.2930*10**(-4)
dimenzije = 24
lamb = L/dimenzije
c = (E*Io)/(lamb**2)
system=stvoriMatricu(dimenzije, lamb, F, W, c)

#deklaracija simbola u listu
v = []
for i in range(dimenzije):
    v.insert(i, Symbol('v{}'.format(i+1))) 

#izracunaj sustav
A=solve_linear_system(system,v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11],v[12],v[13],v[14],v[15],v[16],v[17],v[18],v[19],v[20],v[21],v[22],v[23])

#izracunaj M vrijedonosti
M = []
for i in range(dimenzije):
    if i == 0:
        M.insert(i, F*(-A[v[dimenzije-1]])-W*lamb*dimenzije)
    else:
        M.insert(i, F*A[v[i-1]]-A[v[dimenzije-1]]-(dimenzije-i)*W*lamb)
        

print(M)
print(A)