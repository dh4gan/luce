# Written by D Forgan, 18/12/2014
# Reads in .location files from nbody_2Dflux (C++ code)
# Produces curves for data at specific longitude and latitude

import numpy as np
import matplotlib.pyplot as plt
import infofile

from string import split
from os import system
import sys

pi = 3.1415926585
solrad = 6.9e10
AU = 1.496e13

# Columns of .location file

timecol = 0
xcol =1
ycol = 2
zcol = 3
longcol = 4
latcol = 5
starfluxcol = 6
altcol = 7
azcol = 8
hourcol = 9

# Columns of .position file

timecol = 0
distcol = 1
anomcol = 2
xcol = 3
ycol = 4
zcol = 5

starname = []
starradius= []
startemp = []
starlambda = []
starcolor = []

# Read in info file to obtain input parameters

prefix = raw_input("What is the file prefix? ")
planetname = raw_input("What is the planet name? ")

nfiles, nstars, starname, starradius, startemp, starcolor,fluxmax = infofile.read_infofile(prefix)

# Read in movie input parameters
   
moviechoice = raw_input("Make an animated gif at end? (y/n) ")
if(moviechoice=='y'):
    deletechoice = raw_input("Delete .png files? (y/n) ")

#moviechoice = 'y'
#deletechoice = 'y'

nzeros = int(np.log10(nfiles))

# Create planck functions for both stars

lambda_min = 1.0e2 # minimum/maximum in nm
lambda_max = 1.0e5

lambda_min = lambda_min/1.0e7 # convert these to cm
lambda_max = lambda_max/1.0e7

nspectrum = 1000
dl = (lambda_max-lambda_min)/float(nspectrum)

planck= np.zeros((nstars,nspectrum))

for s in range(nstars):
    wavelengths,holder = infofile.calc_planck_function(startemp[s], lambda_min, lambda_max, nspectrum)
    planck[s,:] = holder[:]

wavelengths[:] = wavelengths[:]*1.0e7
lambda_min = lambda_min*1.0e7
lambda_max = lambda_max*1.0e7
dl = dl*1.0e7

# Arrays to hold flux, altitude and azimuth data

time = np.zeros(nfiles)
flux = np.zeros((nstars,nfiles))
fluxtot = np.zeros(nfiles)
darkness = np.zeros(nfiles)
SFD = np.zeros((nstars,nspectrum))
SFDtot = np.zeros(nspectrum)

altitude = np.zeros((nstars,nfiles))
azimuth = np.zeros((nstars,nfiles))
hourangle = np.zeros((nstars,nfiles))

height= np.zeros((nstars,nfiles))
horizontal = np.zeros((nstars,nfiles))
distance = np.zeros((nstars,nfiles))
starsize = np.zeros((nstars,nfiles))

# Holds positions of all objects
separations = np.zeros((nstars,3,nfiles))

# Read position data

for s in range(nstars):
    filename = prefix+'_'+planetname+'_'+starname[s]+'.location'
    data = np.genfromtxt(filename, usecols=(xcol,ycol,zcol))
    separations[s,:,:] = data.transpose()

# Calculate separation between planet and each star

for s in range(nstars):
    for i in range(nfiles):
        sep = 0.0
        for k in range(3):
            sep =sep + (separations[s,k,i])**2
            
        distance[s,i] = sep
 
# Loop over location files

for s in range(nstars):
    
    locationfile = prefix+'_'+planetname+'_'+starname[s]+'.location'
                    
    # Read in location file
    # First header - time, position data etc

    data = np.genfromtxt(locationfile)

    time = data[:,timecol]
    mylong = data[0,longcol]
    mylat = data[0,latcol]
    
    flux[s,:] = data[:,starfluxcol]
    
    altitude[s,:] = data[:,altcol]
    azimuth[s,:] = data[:,azcol]
    hourangle[s,:] = data[:,hourcol]*180.0/pi                    
        
    horizontal[s,:] = np.cos(altitude[s,:])*np.sin(azimuth[s,:])
    height[s,:] = np.sin(altitude[s,:])
    starsize[s,:] = 20.0*starradius[s]/(distance[s,:]*distance[s,:])
    
        
fluxtot = np.zeros(len(time))

for s in range(nstars):
    fluxtot[:] = fluxtot[:]+flux[s,:]
# end of loop

# Plot flux

fluxfile = 'latlong_'+str(mylat)+'_'+str(mylong)+'_flux_'+prefix+'_'+planetname+'.png'    

fig1 = plt.figure()
ax = fig1.add_subplot(111)
ax.set_xlabel('Time (yr)')
ax.set_ylabel('Flux ($W \, m^{-2}$)')
#ax.set_yscale('log')
#ax.set_ylim(0.0, fluxmax)
ax.plot(time,fluxtot, label = 'Total')
for s in range(nstars):
    ax.plot(time, flux[s,:],label = starname[s])
ax.legend(loc='upper right')

plt.savefig(fluxfile, format= 'png')

# Plot altitude

altfile = 'latlong_'+str(mylat)+'_'+str(mylong)+'_altitude_'+prefix+'_'+planetname+'.png'  

fig1 = plt.figure()
ax = fig1.add_subplot(111)
ax.set_xlabel('Time (yr)')
ax.set_ylabel('Altitude (degrees)')
for s in range(nstars):
    plt.plot(time,altitude[s,:], label=starname[s])
ax.legend()

plt.savefig(altfile, format= 'png')

# Plot azimuth

azfile = 'latlong_'+str(mylat)+'_'+str(mylong)+'_azimuth_'+prefix+'_'+planetname+'.png'  

fig1 = plt.figure()
ax = fig1.add_subplot(111)
ax.set_xlabel('Time (yr)')
ax.set_ylabel('azimuth (degrees)')
for s in range(nstars):
    plt.plot(time,azimuth[s,:], label=starname[s])
ax.legend()

plt.savefig(azfile, format= 'png')

# Hour Angle

hourfile = 'latlong_'+str(mylat)+'_'+str(mylong)+'_hourangle_'+prefix+'_'+planetname+'.png'  

fig1 = plt.figure()
ax = fig1.add_subplot(111)
ax.set_xlabel('Time (yr)')
ax.set_ylabel('Hour Angle (degrees)')
for s in range(nstars):
    plt.plot(time,hourangle[s,:], label=starname[s])
ax.legend()

plt.savefig(hourfile, format= 'png')

# If user wants to see movies, then do this stuff here

if(moviechoice == 'y'):
    
    for itime in range(len(time)):
        
        # Create filename - how many zeros needed?
        num = str(itime+1)
        k=np.log10(itime+1)
        while (k<nzeros): 
            num = "0"+num
            k+=1        
    
        
        skyoutput = 'skypos_'+prefix+num+'_latlong_'+str(mylat)+'_'+str(mylong)+'.png'
        
        print 'Plotting Files for snapshot ',num
        
        # Plot sky position for this timestep
    
        fig2 = plt.figure(1)
        ax = fig2.add_subplot(111)
        ax.set_xlabel('Horizontal Position')
        ax.set_ylabel('Height')
        ax.set_ylim(0,1)
        ax.set_xlim(-1,1)
    
        for s in range(nstars):    
            plt.scatter(horizontal[s,itime],height[s,itime], marker='o',s=starsize[s,itime],c=starcolor)
        
        plt.savefig(skyoutput, format= 'png')
        plt.clf()

        # Create SFD for this position from planck functions, separations and total flux
    
        SFDoutput = 'sfd_'+prefix+num+'_latlong_'+str(mylat)+'_'+str(mylong)+'.png'        
        SFDtot[:] = 0.0
    
        for s in range(nstars):
            SFD[s,:] = pi*starradius[s]*starradius[s]*planck[s,:]/(distance[s,itime]*distance[s,itime])                


            if(np.sum(SFD[s,:])!=0.0):
                SFD[s,:] = SFD[s,:]*flux[s,itime]/(np.sum(SFD[s,:])*dl)
            else:
                SFD[s,:] = 0.0
            
            SFDtot[:] = SFDtot[:] +SFD[s,:]
        
    
        # Plot the SFD
    
        fig3 = plt.figure(2)
        ax = fig3.add_subplot(111)
        ax.set_xlabel('Wavelength (nm)')
        ax.set_ylabel('Spectral Flux Density ($W\,m^{-2} \,nm^{-1}$)')            
        ax.set_yscale('log')
        ax.set_xscale('log')
        ax.set_ylim(1e-10,1e1)
        ax.set_xlim(lambda_min, lambda_max)
        
        for s in range(nstars):        
            ax.plot(wavelengths,SFD[s,:], label=starname[s])
        
        plt.plot(wavelengths,SFDtot,label='Total')
        ax.legend()
        plt.savefig(SFDoutput, format= 'png')

        plt.clf()



    print 'All timesteps plotted'
    # Command for converting images into gifs - machine dependent

    convertcommand = '/opt/ImageMagick/bin/convert '
    #convertcommand = '/usr/bin/convert '

    # Create movie if requested
    
    print 'Creating animated gif of sky pattern, filename skymovie.gif'
    system(convertcommand +'skypos_'+prefix+'*_latlong_'+str(mylat)+'_'+str(mylong)+'.png skymovie.gif')
    
    print 'Creating animated gif of SFD, filename sfdmovie.gif'
    system(convertcommand +'sfd_'+prefix+'*_latlong_'+str(mylat)+'_'+str(mylong)+'.png sfdmovie.gif')
    
    
    if(deletechoice=='y'):
        print 'Deleting png files'
        system('rm skypos_'+prefix+'*.png')
        system('rm sfd_'+prefix+'*.png')                    


print 'Complete'

