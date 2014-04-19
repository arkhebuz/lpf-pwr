# -*- coding: utf-8 -*-
from ncorr_stdev import *
 
# Calculates the density of a pipe fragment
 
h = 38.1/1000
d1 = 16.1/1000
d2 = 11.95/1000
m = 9.12/1000
 
dok_suw = 0.05/1000
dok_wag = 0.01/1000
 
slown = dict(hs=h, d1s=d1, d2s=d2, ms=m)
 
ms, hs, d1s, d2s = sym.symbols('ms hs d1s d2s')
fun = 4*ms/(hs*sym.pi*(d1s**2 - d2s**2))
 
aaa = [ [hs, dok_suw*p3, 0],
[d1s, dok_suw*p3, 0],
[d2s, dok_suw*p3, 0],
[ms, dok_wag*p3, 0] ]
 
a = nsk(aaa, fun, slown, 0)
print a.get_val()
a.frmat('kg/m^3')
a.wspolczynniki("p")
