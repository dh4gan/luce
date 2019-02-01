# Useful functions for reading info file, etc

import numpy as np
from string import split

h = 6.63e-27
c = 2.99e10
k = 1.38e-16

hck = h*c/k

def wavelengthToRGB(wavelength):
    '''Function converts a wavelength (nm)in the visible spectrum to an RGB colour'''
    red =0.0
    green = 0.0
    blue = 0.0
    
    # If wavelength in UV or beyond, make it blue
    if(wavelength<380.0):
        red = 0.0
        green = 0.0
        blue = 1.0
        print "Wavelength ", wavelength, " too blue for visible RGB!"

    if(wavelength>781.0):
        red = 1.0
        green = 0.0
        blue = 0.0
        print "Wavelength ", wavelength, " too red for visible RGB!"


    if(wavelength>=380.0 and wavelength < 440.0):
        red = (wavelength -440)/(440-380)
        green = 0.0
        blue =1.0
    elif(wavelength >=440 and wavelength < 490):
        red = 0.0
        green = ((wavelength)-440)/(490-440)        
        blue = 1.0            
    elif(wavelength>=490 and wavelength < 510):
        red = 0.0
        green = 1.0
        blue = -(wavelength-510)/(510-490)
    elif(wavelength>=510 and wavelength < 580):
        red = (wavelength-510)/(580-510)
        green=1.0
        blue= 0.0
    elif(wavelength>=580 and wavelength < 645):
        red = 1.0
        green = -(wavelength-645)/(645-580)
        blue = 0.0
    elif(wavelength>=645 and wavelength< 781):
        red = 1.0
        green = 0.0
        blue = 0.0

    rgb = (red,green,blue)

    return rgb

def calc_planck_function(Teff, lambda_min, lambda_max, npoints, norm=False):
    '''Calculates the planck function in linear space from lambda_min to lambda max
    (in cm), normalised to unity if norm=True'''
    
    wavelength = np.linspace(lambda_min, lambda_max, npoints, endpoint=True)
    function = np.zeros(npoints)
    
    function[:] = np.exp(hck/(wavelength[:]*Teff))-1
    for i in range(npoints):
        if function[i]!=0.0:
            function[i] = 1.0/function[i]                        
        
    function[:] = function[:]*(2.0*h*c*c)/pow(wavelength[:],5)
    
    if norm==True:
        dl = (lambda_max-lambda_min)/float(npoints)
        integral= 0.0
        for i in range(npoints):
            integral = integral + function[i]*dl
        
        function[:] = function[:]/integral
    
  
    return wavelength,function


def read_infofile(prefix):
    infofile = prefix +'.info'

    print "Reading information file ",infofile
    f =  open(infofile, 'r')

    line = f.readline()
    nfiles=int(line)

    print "There are ", nfiles, " files"

    line = f.readline()
    nstars = int(line)

    starname = []
    starradius = []
    startemp = []

    starcolor = np.zeros((nstars,3))

    for s in range(nstars):
        line = f.readline()
        starname.append(line.strip())
        
        line = f.readline()
        numbers = split(line)
    
        starradius.append(float(numbers[0]))
        startemp.append(float(numbers[1]))
        wavelength =float(numbers[2])*1e7
        rgb = wavelengthToRGB(wavelength)        
        starcolor[s,:] = rgb
    
    line = f.readline()
    fluxmax = float(line)
        
    f.close()
    
    print "There are ",nstars, " stars"
    print "Names: ",starname
    print "Radii: ", starradius
    print "Temperatures: ", startemp
    print "Colours ", starcolor
    print "Maximum Flux: ", fluxmax
    
    return nfiles, nstars, starname, starradius, startemp, starcolor, fluxmax