from math import *

"""
    A series of miscellaneous functions for 
    calculating trigonmetric values using degrees instead of radians. 
    Simplifies planetary calculations.
"""

# Methods/Consts covered by Python:
    # pi = 3.141592653589793
    # radians(deg) -> degrees to radians
    # degrees(rad) -> radians to degrees
    # np.cbrt(x) -> cube root function

# Conversions between radians and degrees for various functions
sind = lambda x: sin(radians(x))
cosd = lambda x: cos(radians(x))
tand = lambda x: tan(radians(x))
asind = lambda x: degrees(asin(x))
acosd = lambda x: degrees(acos(x))
atand = lambda x: degrees(atan(x))
atan2d = lambda y,x: degrees(atan2(y,x))

# Reduces angle to be between 0 and 360 degrees
def rev(x):
    return x - floor(x/360.0) * 360.0