# -*- coding: utf-8 -*-
from numpy import pi, mean, std, double, longdouble as ldbl
# log10 from numpy breaks on double/longdouble
from math import log10, floor, ceil
from scipy.stats import linregress, t, uniform
import scipy.odr.odrpack as odrpack
import sympy as sym

# Misc.
p3 = ldbl(3.0)**(-0.5)
fun  = sym.symbols('fun', cls=sym.Function)
funx = sym.symbols('funx', cls=sym.Function)
funy = sym.symbols('funy', cls=sym.Function)


def frmat(liczba, jednostka):
    """Rounds the uncertainlity properly, to 1 or 2 significant digits."""
    wrt, niep = liczba

    p = int(floor(log10(niep)))
    pot = 10**p

    n_rnd = ceil(niep/pot)*pot
    dgt = 0
    if (n_rnd/niep) >= 1.1:
        n_rnd = ceil(10*niep/pot)*(pot/10)
        dgt = 1

    p1 = int(floor(log10(n_rnd)))
    p2 = int(floor(log10(abs(wrt))))
    w_rnd = round(wrt, dgt-p1)

    r = p2-p1+dgt
    if r < 0:   r = 0

    print ("%."+str(r)+"E") % (w_rnd), "±", ("%."+str(dgt)+"E") % (n_rnd), jednostka


def mean_wgh(wlk, wlk_niep, eng_out=0):
    """Calculates weighted mean of values with different
    variances and its standard deviation.
    """
    # https://en.wikipedia.org/wiki/Weighted_mean#Dealing_with_variance
    # http://www.physicsforums.com/showthread.php?t=612633

    Wi = [(1.0/(i**2)) for i in wlk_niep]
    xWi = [i[0]*i[1] for i in zip(wlk, Wi)]
    V1 = sum(Wi)
    X_sr = sum(xWi)/V1

    N_mns_1 = len(wlk)-1
    if N_mns_1 == 0:    N_mns_1 = 1

    varint = 1.0/V1
    wi_t_resid_sq = [i[1]*(i[0] - X_sr)**2 for i in zip(wlk, Wi)]
    varext = (varint / N_mns_1)*sum(wi_t_resid_sq)

    if eng_out:
        V2 = sum([i**2 for i in Wi])
        varone = ((V1/(V1**2 - V2))*sum(wi_t_resid_sq))**0.5
        print("inside( int: "+str((varint)**0.5)+", ext: "+str((varext)**0.5)+
              ", ext/int: "+str((varext/varint)**0.5)+", "+str(varone)+" )")

    return [X_sr, max((varint)**0.5, (varext)**0.5)]


def regresja(x, y, eng_out=0):
    """Normal linear regression"""
    slope, intercept, r_value, p_value, std_err_A = linregress(x, y)
    std_err_B = std_err_A*(sum([i**2 for i in x])/len(x))**0.5
    if eng_out is True: print slope, std_err_A, intercept, std_err_B, r_value**2
    return [slope, std_err_A, intercept, std_err_B, r_value**2]


def reg_tls(x, x_n, y, y_n, eng_out=0):
    """Total least squares linear regression"""
    x = double(x)
    y = double(y)
    x_n = double(x_n)
    y_n = double(y_n)
    reg = regresja(x, y, 0)     # initial guess

    def f(B, z):
        return B[0]*z + B[1]
    linear = odrpack.Model(f)
    mydata = odrpack.RealData(x, y, sx=x_n, sy=y_n)
    myodr = odrpack.ODR(mydata, linear, beta0=[reg[0], reg[2]], maxit=100)
    myoutput = myodr.run()
    A, B = myoutput.beta
    dA, dB = myoutput.sd_beta

    if eng_out is True: myoutput.pprint()
    return [A, dA, B, dB]


def ciag(data_lst, funkcja):
    """
    funkcja - function returning [x, dx, y, dy]
    """
    out = [double(funkcja(data_lst[eN], eN)) for eN in range(len(data_lst))]
    x  = [i[0] for i in out]
    dx = [i[1] for i in out]
    y  = [i[2] for i in out]
    dy = [i[3] for i in out]

    return x, dx, y, dy


class nsk(object):
    def __init__(self, dane, fun, slown, alfa=0):
        """
        Calculates standard uncertainty of indirect measurements
        (non-correlated input quantities assumed)

        dane = [ [C1, B1, a1], [C2, B2, a2], ... ]
        C - SymPy symbol
        B - B-type standard uncertainty
        a - Values list for calculating A-type standard uncertainty

        fun - SymPy function
        slown - ditctionary with values for SymPy symbols

        For SymPy objects, evaluation with 33 digits of precision is used.
        For SciPy/NumPy functions, Numpy's longdouble is used where possible,
        with fallback to standard double when it breaks.
        """
        self.precision = 33
        self.dane = dane
        self.fun = fun
        self.slown = {k : sym.Float(str(v), self.precision) for k,v in slown.iteritems()}
        alfa = double(alfa)
        suma = []
        self.pochodne = []

        if alfa == 0:
            kp_u = ldbl(1)
        else:
            kp_u = ldbl(uniform.ppf(1-alfa, 0, p3**-1))

        for wielkosc in range(len(dane)):
            poch = fun.diff(dane[wielkosc][0])
            Cx = poch.evalf(subs=slown, n=self.precision)       # 33 digits of precision...

            if dane[wielkosc][2] == 0:
                niep = kp_u*ldbl(dane[wielkosc][1])
            else:
                eNka = len(dane[wielkosc][2])
                if alfa == 0:
                    kp_t = ldbl(1)
                else:
                    kp_t = ldbl(t.ppf(1-alfa/2, eNka))
                niep_A = std(ldbl(dane[wielkosc][2]), ddof=1) / (eNka)**0.5
                niep_B = ldbl(dane[wielkosc][1])
                niep = ((kp_t*niep_A)**2 + (kp_u*niep_B)**2)**0.5

            self.pochodne.append([dane[wielkosc][0], poch, Cx])
            suma.append((Cx*sym.S(niep))**2)

        self.u = (sum(suma))**0.5

    def get_val(self):
        """Returns values as sympy.core.numbers.Float"""
        wartosc = self.fun.evalf(subs=self.slown, n=self.precision)
        niepewnosc = self.u
        return [wartosc, niepewnosc]

    def frmat(self, jednostka):
        frmat(self.get_val(), jednostka)

    def wspolczynniki(self, nazwa_funkcji):
        """Prints the sensitivity coefficients"""
        for poch in range(len(self.pochodne)):
            rown1 = sym.Eq(sym.Derivative(nazwa_funkcji, self.pochodne[poch][0]),
                           self.pochodne[poch][1])

            rown2 = sym.Eq(sym.symbols("C"+str(self.pochodne[poch][0])), rown1)

            sym.pprint(sym.Eq(rown2, self.pochodne[poch][2]), use_unicode=True)
            print ""

