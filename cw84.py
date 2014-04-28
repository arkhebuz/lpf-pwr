# -*- coding: utf-8 -*-
from ncorr_stdev import *

#       dokladnosci
dok_x = 2.0/1000        # mm → m
dok_lam = 10.0**(-9.0)      # nm → m

# Odleglosc l
L1 = 20.0/100           # cm → m
L2 = 60.0/100           # cm → m

# Polozenia prazkow lewy-prawy (1 rzad)
d_n1 = [ [0.204, 0.296],
         [0.196, 0.304],
         [0.187, 0.313],
         [0.181, 0.319],
         [0.174, 0.325] ]       # m
# (2 rzad)
d_n2 = [ [0.153, 0.347],
         [0.137, 0.363],
         [0.120, 0.380],
         [0.105, 0.395],
         [0.088, 0.413] ]       # m

# Lambda, m
lam = [420*10**(-9), 480*10**(-9), 550*10**(-9), 600*10**(-9), 660*10**(-9)]

# Obliczenie stalej siatki
xl = [m[0] for m in d_n1] + [n[0] for n in d_n2]
xp = [m[1] for m in d_n1] + [n[1] for n in d_n2]

lam2 = [m for m in lam]*2
n = [1,1,1,1,1, 2,2,2,2,2]

ds, ns, lams, L1s, L2s, xls, xps = sym.symbols('ds ns lams L1s L2s xls xps')
fun =  (2*ns*lams*(sym.sqrt((L2s-L1s)**2 + ((xps-xls)/2)**2)) )/(xps-xls)
aaa = [ [lams, dok_lam*p3, 0],
        [L1s, dok_x*p3, 0],
        [L2s, dok_x*p3, 0],
        [xls, dok_x*p3, 0],
        [xps, dok_x*p3, 0] ]

def fun1(__, k):
    slown = dict(ns=1.0*n[k],
                 lams=lam2[k],
                 L1s = L1,
                 L2s = L2,
                 xls = xl[k],
                 xps = xp[k])

    a = nsk(aaa, fun, slown, 0)
    return [a.get_val()[0], a.get_val()[1], 0, 0]

d, d_niep, __, __ = ciag(range(0, 10), fun1)

stala = mean_wgh(d, d_niep, 0)
print "Stala siatki d:     ",
frmat(stala, "m")


####################### obliczenia filtrow
l1p = 0.100
l2p = 0.500

# Filtr "IF 625Ca", polozenia prazkow
ca  = [ [0.178, 0.322],     # rzad 1
        [0.096, 0.406] ]    # rzad 2
# Filtr ?
uknw = [ [0.195, 0.305],
         [0.139, 0.361] ]
# Filtr "700"
nir = [ [0.166, 0.334] ]
# ↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓
set_filtr = ca

fun = ds*(xps-xls)/(2*ns*sym.sqrt( (L2s - L1s)**2 + ((xps-xls)/2)**2 ))
aaa.append([ds, stala[1], 0])

def fun2(flt, n):
    slown = dict(ns=n+1.0,
                 ds = stala[0],
                 L1s = l1p,
                 L2s = l2p,
                 xls = flt[0],
                 xps = flt[1])

    a = nsk(aaa, fun, slown, 0)
    return a.get_val()[0], a.get_val()[1], 0, 0

flt_lam, flt_n, __, __ = ciag(set_filtr,fun2)
print "Dlugosc fali filtra:",

if len(flt_lam) > 1:
    mn = mean_wgh(flt_lam, flt_n)
    frmat(mn,"m")
else:
    frmat([flt_lam[0], flt_n[0]], "m")
