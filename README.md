# Planetary-Motion
Contains the planetary-motion package used for calculating the planet locations and generating data for the Analyzing Divination Systems as Cryptographic Random Number Generators WXML project. 
The equations and processes used in this package are courtesy of the "Computing planetary positions" tutorial by Paul Schlyter and can be found [here](https://www.stjarnhimlen.se/comp/tutorial.html). 

### planetary-motion package
The package contains the following files:
Coordinates.py (Used for coordinate system conversions and for deriving Ecliptic coordinates from the primary orbital elements.)
OrbitalElements.py (Provides the primary orbital elements of a planet for a given date)
OrbitalMotion.py (Returns the location of a planet on a given date using a specified coordinate system)
Perturbations.py (Revises the coordinates of the Moon, Jupiter, and Saturn due to perturbations caused by gravity)
Trig.py (Used for radian/degree conversions and for keeping outputted angles between 0 and 360 degrees)

### WXML-Testing notebook
This file contains the graphs and entropy data used in the WXML presentation, along with sample code for using the original OrbitalMotion file. It is currently being revised for the package.
