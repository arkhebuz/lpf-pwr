# -*- coding: utf-8 -*-
from ncorr_stdev import *

#          k    rl     rp
pomsy = [ [5, 46.35, 41.49],
          [6, 46.52, 41.34],
          [7, 46.67, 41.17],
          [8, 46.81, 41.04],
          [9, 46.90, 40.91],
          [10, 47.10, 40.79],
          [11, 47.18, 40.70],
          [12, 47.31, 40.60],
          [13, 47.42, 40.48],
          [14, 47.50, 40.33],
          [15, 47.59, 40.25],
          [16, 47.74, 40.15],
          [17, 47.86, 40.07],
          [18, 47.95, 39.97],
          [19, 48.03, 39.89],
          [20, 48.12, 39.79] ]

#grone = [1, 45.56, 45.43]
centers = [(pomsy[i][1]+pomsy[i][2])/2 for i in range(len(pomsy))]
dz = 0.01/1000                                     # mm → m
dok_sr = dz + (max(centers)-min(centers))/2000     # mm → m
lam = [589*10**(-9), 20*10**-9]                    # [λ, ∆λ] nm → m
dK = 3
alfa = 0

#               init
rl, rp, k, l = sym.symbols('rl rp k l')
funy = ((rl-rp)**2)/4.0
funx = k*l

aaay = [ [rl, dok_sr*p3, 0],
         [rp, dok_sr*p3, 0]]
aaax = [ [l, lam[1]*p3, 0] ]

def funkcjaxy(pom, iteracja):
    slown = dict(rl=pom[1]*10**-3,
                 rp=pom[2]*10**-3,
                 k=pom[0]+dK,
                 l=lam[0])

    a = nsk(aaay, funy, slown, alfa)
    b1 = a.get_val()
    bbb = nsk(aaax, funx, slown, alfa)
    b2 = bbb.get_val()
    #      [  x,    dx,     y,    dy   ]
    return [b2[0], b2[1], b1[0], b1[1] ]

x, dx, y, dy = ciag(pomsy, funkcjaxy)
reg_zw = regresja(x, y, 0)
reg_tls = reg_tls(x, dx, y, dy, 0)

#~ srednia ciagu estymat czastkowych
fun = ((rl-rp)**2)/(4*k*l)
aaa = [aaay[0], aaay[1], aaax[0]]

def funkcja(pom, i):
    slown = dict(rl=pom[1]*0.001,
                 rp=pom[2]*0.001,
                 k=pom[0]+dK,       # k+2.815
                 l=lam[0])

    a = nsk(aaa, fun, slown, alfa)
    b = a.get_val()
    return [b[0], b[1], 0, 0]

x, dx, __, __ = ciag(pomsy, funkcja)
sr_waz = mean_wgh(x, dx)

print("                     A (R)                 dA (dR)               "
      "B                         dB                       R^2")
print "Reg. zwykla:     ", reg_zw
print "Reg. TLS:        ", reg_tls
print "Sr. wazona:      ", sr_waz
