# -*- coding: utf-8 -*-
from numpy import pi, mean, std, longdouble, double
from math import log10, floor, ceil
from scipy.stats import linregress, t, uniform
import scipy.odr.odrpack as odrpack
import sympy as sym

# Misc.
p3 = 3.0**(-0.5)
fun = sym.symbols('fun', cls=sym.Function)
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
	if r < 0:	r = 0
	
	print ("%."+str(r)+"E") % (w_rnd),"±", ("%."+str(dgt)+"E") % (n_rnd), jednostka

def mean_wgh(wlk, wlk_niep, eng_out):
	"""Calculates weighted mean of values with different variances and its standard deviation."""
	# https://en.wikipedia.org/wiki/Weighted_mean#Dealing_with_variance
	# http://www.physicsforums.com/showthread.php?t=612633
	
	Wi = [ (1.0/(i**2)) for i in wlk_niep]
	xWi = [ i[0]*i[1] for i in zip(wlk, Wi)]
	V1 = sum(Wi)
	X_sr = sum(xWi)/V1
	
	N_mns_1 = len(wlk)-1
	if N_mns_1 == 0:	N_mns_1 = 1
	
	varint = 1.0/V1
	varext = (varint / N_mns_1)*sum([i[1]*(i[0] - X_sr)**2 for i in zip(wlk, Wi)])
	
	if eng_out: 
		V2 = sum([i**2 for i in Wi])
		varone = ( (V1/(V1**2 - V2))*sum([i[1]*(i[0] - X_sr)**2 for i in zip(wlk, Wi)]) )**0.5
		print "inside( int: "+str((varint)**0.5)+", ext: "+str((varext)**0.5)+", ext/int: "+str((varext/varint)**0.5)+", "+str(varone)+" )"
	return [X_sr, max((varint)**0.5, (varext)**0.5)]

def regresja(x, y, eng_out):
	"""Normal linear regression"""
	slope, intercept, r_value, p_value, std_err_A = linregress(x, y)
	std_err_B = std_err_A*(sum([i**2 for i in x])/len(x))**0.5
	if eng_out == True: print slope, std_err_A, intercept, std_err_B, r_value**2
	return [slope, std_err_A, intercept, std_err_B, r_value**2]

def reg_tls(x, x_n, y, y_n, eng_out):
	"""Total least squares linear regression"""
	x = double(x)
	y = double(y)
	x_n = double(x_n)
	y_n = double(y_n)
	reg = regresja(x, y, 0) 	# initial guess
	
	def f(B, z):
		return B[0]*z + B[1]
	linear = odrpack.Model(f)
	mydata = odrpack.RealData(x, y, sx=x_n, sy=y_n)
	
	myodr = odrpack.ODR(mydata, linear, beta0=[reg[0], reg[2]], maxit=100)
	myoutput = myodr.run()
	if eng_out == True: myoutput.pprint()
	A, B = myoutput.beta
	dA, dB =  myoutput.sd_beta
	return [ A, dA, B, dB]

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
	def __init__(self, dane, fun, slown, alfa):
		"""
		Calculates standard uncertainty of indirect measurements (non-correlated input quantities assumed)
		
		dane = [ [C1, B1, a1], [C2, B2, a2], ... ]
		C - SymPy symbol 
		B - B-type standard uncertainty
		a - Values list for calculating A-type standard uncertainty
		"""
		self.dane = dane
		self.fun = fun
		self.slown = slown
		
		suma = []
		self.pochodne = []
		
		if alfa != 0:	kp_u = uniform.ppf( 1-alfa, 0, p3**-1) 
		else:			kp_u = 1
		
		for wielkosc in range(len(dane)):
			poch = fun.diff(dane[wielkosc][0])
			#print poch, poch.evalf(subs=slown)
			Cx = longdouble(poch.evalf(subs=slown))
			self.pochodne.append([dane[wielkosc][0], poch, Cx])
			
			if dane[wielkosc][2] != 0:
				eNka = len(dane[wielkosc][2])
				
				if alfa != 0:		kp_t = t.ppf( 1-alfa/2, eNka)
				else:				kp_t = 1
				
				niep_A = std(dane[wielkosc][2], ddof=1) / (eNka)**0.5
				niep = ( (kp_t*niep_A)**2 + (kp_u*dane[wielkosc][1])**2)**0.5
			else:
				niep = kp_u*dane[wielkosc][1]
			suma.append((Cx*niep)**2)
		
		self.u = (sum(suma))**0.5
	
	def get_val(self):
		wartosc = self.fun.evalf(subs=self.slown)
		niepewnosc = self.u
		return [wartosc, niepewnosc]
	
	def frmat(self, jednostka):
		frmat(self.get_val(), jednostka)
	
	def wspolczynniki(self, nazwa_funkcji):
		"""Prints the sensitivity coefficients"""
		for poch in range(len(self.pochodne)):
			rown1 = sym.Eq(sym.Derivative(nazwa_funkcji, self.pochodne[poch][0]), 
						   self.pochodne[poch][1])
			sym.symbols("C"+str(self.pochodne[poch][0]))
			rown2 = sym.Eq(sym.symbols("C"+str(self.pochodne[poch][0])), rown1)
			
			sym.pprint(sym.Eq(rown2, self.pochodne[poch][2]), use_unicode=True)
			print ""

