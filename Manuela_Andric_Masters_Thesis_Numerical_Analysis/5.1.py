# -*- coding: utf-8 -*-
"""
Created on Sat May 13 17:38:47 2017

@author: Manuela Andric
"""

from sympy import Symbol, Matrix, nsolve, plot, solveset, pprint, S   
    
#generiraj matricu (dimenzije min>6)
def stvoriMatricu(dimenzije, pocetna):
    x = Symbol("x")
    M=Matrix.zeros(dimenzije)
    delta = pocetna-1  
    pozicija = dimenzije-4   
    
    if(dimenzije<6):
        print("Dimenzije min 6*6!")        
        return 0
    
    for i in range(dimenzije):
        #iznimke pravila za prvi red
        if(i==0):
            M[i,dimenzije-3] = pocetna
            M[i,dimenzije-2] = x-2*pocetna
            M[i,dimenzije-1] = pocetna-x
        #iznimke pravila za predzadnji red
        elif(i==dimenzije-2):
            M[i,0] = x-(delta*(dimenzije-1)+1)*2
            M[i,1] = (delta*(dimenzije-1)+1)
            M[i,dimenzije-1] = -x
        #iznimke pravila za zadnji red
        elif(i==dimenzije-1):
            M[i,0] = 8
            M[i,dimenzije-1] = -x
        #ostali redovi
        else:
            M[i, pozicija-i+1] = delta*(i+1)+1
            M[i, pozicija-i+2] = x - (delta*(i+1)+1)*2
            M[i, pozicija-i+3] = delta*(i+1)+1
            M[i, dimenzije-1] = -x    
        
    return M    


#deklaracija simbola
x = Symbol("x")  
l = Symbol("l")
E = Symbol("E")
I = Symbol("I")

#deklaracija pocetnih uvjeta
lmbd = l/24
M=stvoriMatricu(24, 1.125)

#print(M)

#izracunaj
eq=M.det()
rez=nsolve(eq,x,-1) 
#plot(eq)   
#pprint(solveset(M.det(), x, domain=S.Reals))
print(rez*(E*I/lmbd**2))    