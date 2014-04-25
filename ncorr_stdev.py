# -*- coding: utf-8 -*-
from numpy import dtype, finfo, mean, std, double
from numpy import errstate as numpyerrstate, longdouble as ldbl
# log10 from numpy breaks on double/longdouble
from math import log10, floor, ceil
from scipy.stats import linregress, t, uniform
import scipy.odr.odrpack as odrpack
import sympy as sym

# Set SymPy's precision little higher than longdouble's
SymPy_dgt_increment = 3
SymPy_precision = finfo(ldbl).precision + SymPy_dgt_increment

# Misc.
p3 = ldbl(3.0)**(-0.5)
fun  = sym.symbols('fun', cls=sym.Function)
funx = sym.symbols('funx', cls=sym.Function)
funy = sym.symbols('funy', cls=sym.Function)


def frmat(number, unit):
    """Rounds the uncertainlity properly, to 1 or 2 significant digits.

    number - [value, uncertainlity] tuple
    unit - unit, a string
    """
    wrt, niep = number

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

    print ("%."+str(r)+"E") % (w_rnd), "±", ("%."+str(dgt)+"E") % (n_rnd), unit


def mean_wgh(val_lst, sderr_lst, eng_out=0):
    """Calculates weighted mean of values with different
    variances and its standard deviation.

    val_lst - list of values
    sderr_lst - list of corresponding standard deviations
    """
    # https://en.wikipedia.org/wiki/Weighted_mean#Dealing_with_variance
    # http://www.physicsforums.com/showthread.php?t=612633

    wlk = ldbl(val_lst)
    wlk_niep = ldbl(sderr_lst)
    with numpyerrstate(divide='raise'):
        try:
            Wi = [(1.0/(i**2)) for i in wlk_niep]
            xWi = [i[0]*i[1] for i in zip(wlk, Wi)]
            V1 = sum(Wi)
            N_mns_1 = abs(len(wlk))-1

            X_sr = sum(xWi)/V1
            varint = 1.0/V1
            wi_t_resid_sq = [i[1]*(i[0] - X_sr)**2 for i in zip(wlk, Wi)]
            varext = (varint/N_mns_1)*sum(wi_t_resid_sq)
            if eng_out:
                V2 = sum([i**2 for i in Wi])
                varone = ((V1/(V1**2 - V2))*sum(wi_t_resid_sq))**0.5
                print("inside mean_wgh: int: "+str((varint)**0.5)+
                      ", ext: "+str((varext)**0.5)+
                      ", ext/int: "+str((varext/varint)**0.5)+
                      ", "+str(varone)+" ;")
            return [X_sr, max((varint)**0.5, (varext)**0.5)]

        except (ZeroDivisionError, FloatingPointError) as mean_wgh_error:
            print ("Whoops, check your input data. You might have too short l"
                "ists (len<2) or standard deviation equal to zero somewhere.")
            raise mean_wgh_error


def regresja(x, y, eng_out=0):
    """Normal linear regression"""
    x = ldbl(x)
    y = ldbl(y)
    slope, intercept, r_value, p_value, std_err_A = ldbl(linregress(x, y))
    std_err_B = std_err_A*(sum([i**2 for i in x])/len(x))**0.5
    if eng_out is True:
        print slope, std_err_A, intercept, std_err_B, r_value**2, p_value
    return [slope, std_err_A, intercept, std_err_B, r_value**2]


def reg_tls(x, dx, y, dy, eng_out=0):
    """Total least squares linear regression"""
    # odrpack breaks on longdbl/sympy.float
    x = double(x)
    y = double(y)
    x_n = double(dx)
    y_n = double(dy)
    reg = double(regresja(x, y, 0))     # initial guess

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


def ciag(data_lst, function):
    """data_lst - list of values
    function - function taking data_lst[i] and returning [x, dx, y, dy]
    """
    out = [ldbl(function(data_lst[eN], eN)) for eN in range(len(data_lst))]
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

        For SciPy/NumPy functions, Numpy's longdouble is used where possible,
        with fallback to standard double when it breaks.
        For SymPy objects, evaluation with (Numpy's longdouble
        precision +3) digits is used.
        """
        self.precision = SymPy_precision
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
            Cx = poch.evalf(subs=slown, n=self.precision)

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
        """See frmat function"""
        frmat(self.get_val(), jednostka)

    def wspolczynniki(self, function_name):
        """Prints the sensitivity coefficients"""
        for poch in range(len(self.pochodne)):
            rown1 = sym.Eq(sym.Derivative(function_name, self.pochodne[poch][0]),
                           self.pochodne[poch][1])

            rown2 = sym.Eq(sym.symbols("C"+str(self.pochodne[poch][0])), rown1)

            sym.pprint(sym.Eq(rown2, self.pochodne[poch][2]), use_unicode=True)
            print ""


if __name__ == '__main__':
    typ_e = dtype(ldbl)
    size = dtype(ldbl).itemsize
    precision = finfo(ldbl).precision
    print ("The NumPy's long double on you machine seems to "
           "be called {0}, \nwith {1} bytes of size and "
           "{2} digits of precision.").format(typ_e, size, precision)
    print ("SymPy precision is thus set a little higher, at {0} digits."
            ).format(SymPy_precision)
