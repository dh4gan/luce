# Written by D Forgan, 17/12/2014
# Reads in output dumps from nbody_2DFlux (C++ code)
# Produces 2D maps of each timestep (PNG)

import numpy as np
import matplotlib.pyplot as plt
import infofile
from string import split
from os import system
#from sys import exit

pi = 3.1415926585

longcol = 0
latcol = 1
fluxcol = 2
darkcol = 3


# Read in input parameters

prefix = raw_input("What is the file prefix? ")
planetname = raw_input("What is the planet name? ")

nfiles, nstars, starname, starradius, startemp, starcolor,fluxmax = infofile.read_infofile(prefix)

moviechoice = raw_input("Make an animated gif at end? (y/n) ")
deletechoice = 'n'
if(moviechoice=='y'):
    deletechoice = raw_input("Delete .png files? (y/n) ")

nzeros = int(np.log10(nfiles))+1

# Loop over files

for i in range(nfiles):

    # Create filename - how many zeros needed?
    num = str(i+1)
    k=np.log10(i+1)
    while (k<nzeros): 
        num = "0"+num
        k+=1     
        
    
    inputfile = prefix +'_'+planetname+'_'+num+'.flux'    
    fluxfile = 'flux_'+prefix+'_'+planetname+'.'+num+'.png'    
       
    # Read in header - time, position data etc

    f = open(inputfile, 'r')

    line = f.readline()

    numbers = split(line)

    time=float(numbers[0])
    nlat = int(numbers[1])
    nlong = int(numbers[2])
        
    f.close()
    
    print 'File ', str(i+1),' Time  ', time, " yr"
    
    # Read in rest of file
    
    data = np.genfromtxt(inputfile, skiprows=1)
        
    # Reshape to fit 2D array
    
    latitude = data[:,latcol].reshape(nlat,nlong)*180.0/pi
    longitude= data[:,longcol].reshape(nlat,nlong)*180.0/pi
    flux = data[:,fluxcol].reshape(nlat,nlong)
    darkness = data[:,darkcol].reshape(nlat,nlong)

    # Plot 2D maps of this timestep    
    
    # Flux
    
    fig1 = plt.figure(1)
    ax = fig1.add_subplot(111)
    ax.set_xlabel('Longitude (degrees)')
    ax.set_ylabel('Latitude (degrees)')
    plt.pcolor(longitude,latitude,flux, cmap='spectral',vmin = 0.0, vmax = fluxmax)
    plt.colorbar()

    plt.savefig(fluxfile, format= 'png')
    plt.clf()

# end of loop


# Command for converting images into gifs - machine dependent

convertcommand = '/opt/ImageMagick/bin/convert '
#convertcommand = '/usr/bin/convert '

# Create movie if requested
if(moviechoice=='y'):
    print 'Creating animated gif of flux pattern, filename fluxmovie.gif'
    system(convertcommand +'-delay 10 flux_'+prefix+'*.png fluxmovie.gif')

    if(deletechoice=='y'):
        print 'Deleting png files'
        system('rm '+prefix+'_'+planetname+'*.png')                        
        

print 'Complete'
