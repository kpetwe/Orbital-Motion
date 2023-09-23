import numpy as np
from math import *
import matplotlib.pyplot as plt
from scipy.fft import fft, ifft
from scipy.signal import find_peaks
from scipy.stats import entropy as ent

from planetary_motion import OrbitalMotion as om
from planetary_motion import Coordinates


d = om.date(19, 4.0, 1990)
# Example from Tutorial
# https://www.stjarnhimlen.se/comp/tutorial.html#5
d = om.date(19, 4.0, 1990)
sun_arr = om.getSunData(d = 19, m = 4, y = 1990, ct = "rect", csys = "eclip")

print(sun_arr)
merc_arr = om.getPlanData(d = 19, m = 4, y = 1990, p = "saturn", ct="sphere", geoc=False)
# Heliocentric
print(merc_arr)
#merc_geo = om.togeo(sun_arr, merc_arr)
#print(merc_geo)
#mer_geo_eq = Coordinates.ecEq(merc_geo[0], merc_geo[1], merc_geo[2], sun_arr[11])
#print(mer_geo_eq[0], mer_geo_eq[1], mer_geo_eq[2])
#sphere = Coordinates.sphere(mer_geo_eq[0], mer_geo_eq[1], mer_geo_eq[2])
#print(sphere)

"""
sun_arr:
x, y, z
(0.8810475240419581, 0.4820994287140291, 0.0, 
xeq, yeq, zeq
0.8810475240419581, 0.44231332298204745, 0.19177795357906033, 
RA, Decl, r
26.658077679334305, 11.008374735025573, 1.0043229554216402, 
L, M, oblecl
26.838831863999985, 104.06528413449996, 23.4405623709)


merc_arr: xeclip  yeclip  zeclip long lat r 
(-0.3678208693993782, 0.06108452909848423, 0.03869908805497192, 170.57086510954457, 5.925527266740486, 0.37486148252012297)
rect equat
0.5132266546425799 0.48296234260038473 0.25158260486587974
ra decl r
(43.25988423931299, 19.64595158381238, 0.7482967529508375)

"""