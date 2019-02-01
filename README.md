Luce - 2D surface flux patterns on bodies inside an N Body integrator
=

This C++ code produces maps of surface flux (longitude/latitude) received by bodies under the effects of an arbitrary number and arrangement of other bodies, such as stars, planets and moons.

This has been developed from progenitor codes whose work was published in:

Forgan, Mead, Cockell & Raven (2015), "Surface flux patterns on planets in circumbinary systems and potential for photosynthesis", International Journal of Astrobiology, Volume 14, Issue 3, pp. 465-478 (DOI:10.1017/S147355041400041X)

Brown, Mead, Forgan Raven & Cockell (2014), "Photosynthetic potential of planets in 3 : 2 spin-orbit resonances", International Journal of Astrobiology, Volume 13, Issue 4, pp. 279-289 (DOI: 10.1017/S1473550414000068)

If you plan to use this code for your own research, or if you would like to contribute to this repository then please get in touch with a brief description of what you would like to do.  I will seek co-authorship for any subsequent publications.

Features:
--------
* Simple orbital setup routines or direct cartesian vector input positions for bodies
* 2D flux and darkness time recorded on any number of bodies
* 4th Order Hermite N Body integration (shared variable timestep)
* Library of Python 2.7 plotting scripts 
* Library of example parameter setups to run

Future Features/Wishlist:
-
* Planetary Illumination

Requirements:
-------------
* C++ compiler (g++ recommended) and Makefile software (e.g. gmake)
* Python for plotting scripts (scripts developed in Python 2.7) - dependencies include numpy, matplotlib, scipy

The code reads in a single input parameter file, which contains a set of global parameters for all bodies in the simulation, along with specific parameters for each
body included in the simulation. Parameter files can either specify the initial positions of all bodies, or the initial Keplerian orbits of all bodies.

Further details of the parameter file structure can be found in the userguide in `\docs`, and example parameter files are given in `\paramfiles`

Once compiled, the code is executed with the command

`> ./luce input.params`

The code was originally developed using the eclipse CDT.  We now recommend using the Makefile in `\src` to compile with g++.
