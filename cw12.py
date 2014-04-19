# -*- coding: utf-8 -*-
from ncorr_stdev import *

#										dokladnosci
dok_sow = 0.05/1000 	# mm → m
dok_sr = 0.01/1000	# mm → m
dok_przym = 1.0/1000 	# mm → m
dok_wag = 1.0/1000	# g → kg
dok_stop = 0.2 		# s

#										pomiary
n = 50.0
dl_drut = 614.0/1000	# mm → m
d_tarczy = 140.1/1000	# mm → m
sr_drutu = 0.59/1000	# mm → m

tbt = [5*60+41.35, 5*60+41.86]				# s
tzt = [6*60+46.79, 6*60+47.24, 6*60+47.08]		# s
tbt_sr = mean(tbt)					# s
tzt_sr = mean(tzt)					# s

masa = [0.2473, 0.2474, 0.2475, 0.2474, 0.2474]		# g → kg
m_sr = mean(masa)					# kg

#										definiowanie funkcji
slown = dict(ms=m_sr, 
	     ls=dl_drut, 
	     ss=d_tarczy, 
	     ns=n, 
	     ds=sr_drutu, 
	     t2s=tzt_sr, 
	     t1s=tbt_sr)

ms, ls, ss, ns, ds, t2s, t1s = sym.symbols('ms ls ss ns ds t2s t1s')
fun = 16*sym.pi*ms*ls*(ss**2)*(ns**2)/((ds**4)*(t2s**2 - t1s**2))

#										pakowanie danych
aaa = [ [ms, dok_wag*p3, masa], 
		[ls, dok_przym*p3, 0], 
		[ss, dok_sow*p3, 0], 
		[ds, dok_sr*p3, 0], 
		[t1s, dok_stop*p3, tbt], 
		[t2s, dok_stop*p3, tzt] ]

a = nsk(aaa, fun, slown, 0.0)
print a.get_val()
a.frmat('GPa')
#a.wspolczynniki("G")

