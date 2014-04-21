# -*- coding: utf-8 -*-
import scipy.constants as consts
from ncorr_stdev import *
# cw107

btz = consts.physical_constants["Boltzmann constant"]
c   = consts.physical_constants["speed of light in vacuum"]
h   = consts.physical_constants["Planck constant"]

filters = {'green' : './cw107/green_filter',
           'orange': './cw107/orange_filter',
           'blue'  : './cw107/blue_filter'}
lamb = {'green' : 0.500*10**-6,          # nm
        'orange': 0.589*10**-6,          # nm
        'blue'  : 0.458*10**-6}          # nm
dlamb = 20*(10**-9)                      # nm
color = 'orange'

#~ reading data
with open("./cw107/r_od_t", "r") as input_file:
    #~ zaleznosc temperatury zarowki od jej oporu
    opor = zip(*(double(line.strip().replace(",", ".").split()) for line in input_file))
    R_eq = regresja(opor[0], opor[1], 0)

with open(filters[color], "r") as input_file:
    data = zip(*(line.strip().split() for line in input_file))

#~  M-3850 for U&I, U722A for i
U   = [double(eval(i)) for i in data[0]]        # V
U_n = [0.003*j +0.01 for j in U]                # V
I   = [double(eval(i)) for i in data[1]]        # A
I_n = [0.015*j +5*0.01 for j in I]              # A
sc  = [double(eval(i)) for i in data[2]]
pre_i = [double(eval(j)) for j in data[3]]
i   = [j[0]*j[1] for j in zip(sc, pre_i)]       # A
i_n = [0.02*j for j in sc]                      # A

A, Us, Is, iss, B = sym.symbols('A Us Is iss B')
funx = 1.0/(A*Us/Is + B)
funy = sym.ln(iss)

def funkcja(pom, j):
    slown = dict(A=R_eq[0], 
                 B=R_eq[2],
                 Us=U[j],
                 Is=I[j],
                 iss=i[j])
    
    aaax = [ [A, R_eq[1], 0],
             [B, R_eq[3], 0],
             [Us, U_n[j]*p3, 0],
             [Is, I_n[j]*p3, 0],
             [iss, i_n[j]*p3, 0] ]
    aaay = [ [iss, i_n[j]*p3, 0] ]
    
    ax = nsk(aaax, funx, slown, 0.0)
    bx = ax.get_val()
    ay = nsk(aaay, funy, slown, 0.0)
    by = ay.get_val()
    return [bx[0], bx[1], by[0], by[1]]

x, dx, y, dy = ciag(range(len(U)), funkcja)
wspy = reg_tls(x, dx, y, dy, 0)

#                               kalkulacja h
slown = dict(K=btz[0], 
             lam=lamb[color],
             C=c[0],
             B=-wspy[0])

B, K, lam, C = sym.symbols('B K lam C')
fun = B*K*lam/C

aaa = [ [B, wspy[1], 0],
        [K, btz[2], 0],
        [lam, dlamb*p3, 0] ]

a = nsk(aaa, fun, slown, 0.0)
print "Filtr: "+str(color)+" "+str(lamb[color])+" nm"
print "Obliczone h:     ",
a.frmat("J s")
print "Wartosc tablicowa:   "+str(h[0])+" ± "+str(h[2])+h[1]
