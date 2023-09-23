# Planetary-Motion
Contains the planetary-motion package used for calculating the planet locations and generating data for the Analyzing Divination Systems as Cryptographic Random Number Generators WXML project. 
The equations and processes used in this package are courtesy of the "Computing planetary positions" tutorial by Paul Schlyter and can be found [here](https://www.stjarnhimlen.se/comp/tutorial.html). 

### planetary-motion package
The package contains the following files:
- **Coordinates.py** (Used for coordinate system conversions and for deriving Ecliptic coordinates from the primary orbital elements.)
- **OrbitalElements.py** (Provides the primary orbital elements of a planet for a given date)
- **OrbitalMotion.py** (Returns the location of a planet on a given date using a specified coordinate system)
- **Perturbations.py** (Revises the coordinates of the Moon, Jupiter, and Saturn due to perturbations caused by gravity)
- **Trig.py** (Used for radian/degree conversions and for keeping outputted angles between 0 and 360 degrees)

### Calculating planetary positions from the package
To calculate the values, we first need to import the OrbitalMotion file into our project:
```python
from planetary_motion import OrbitalMotion as om
```

To calculate the position of the Sun relative to the Earth on the date April 19th, 1990 at 12:00 am Universal Time,
we would use the method getSunData() to do so. <br>
getSunData() has the following parameters: 
- d - requested date - default of 31
- m - requested month - default of 12
- y - requested year - default of 1999
- ut - requested time in Universal Time (as a float value) - default of 0.0
- ct - coordinate type (input "rect" for rectangular coordinates and "sphere" for spherical coordinates) - default "rect"
- csys - coordinate system (input "eclip" for ecliptic coordinates and "equat" for equatorial coordinates) - default "eclip"

```python
sun_arr = om.getSunData(d = 19, m = 4, y = 1990, ut = 0.0, ct = "sphere", csys = "eclip")
```

Similarly, the method getPlanData() will return the requested data for a planet (or the moon). As of this point in time, planetary-motion returns the heliocentric ecliptic
coordinates for Mercury, Venus, Mars, Jupiter, and Saturn, and returns the geocentric ecliptic and equatorial coordinates for these planets and the Moon. <br>
getPlanData() has the following parameters:
- d - requested date - default of 31
- m - requested month - default of 12
- y - requested year - default of 1999
- ut - requested time in Universal Time (as a float value) - default of 0.0
- p - requested planet (acceptable inputs are: "mercury", "venus", "mars", "jupiter", "saturn", and "moon"  and are case-insensitive) - default of "mercury"
- ct - coordinate type (input "rect" for rectangular coordinates and "sphere" for spherical coordinates) - default "rect"
- csys - coordinate system (input "eclip" for ecliptic coordinates and "equat" for equatorial coordinates) - default "equat"
- geoc - boolean = True for a geocentric orbit, False for heliocentric* - default "True"
    - *The Moon's coordinates will always be geocentric, and as of this moment heliocentric will only return ecliptic coordinates

```python
merc_arr = om.getPlanData(d = 19, m = 4, y = 1990, p = "mercury", ct="sphere", geoc=False)
```

Side note: For latitude/longitudinal coordinates, use ecliptic spherical, and for Right-Ascension/Declination coordinates, use equatorial spherical. Also, note that all spherical coordinates are given in degrees and not in hours.

### WXML-Testing notebook
This file contains the graphs and entropy data used in the WXML presentation, along with sample code for using the original OrbitalMotion file. 
