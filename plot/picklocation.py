# Written by D Forgan, 12/8/2013
# Reads in output dumps from multistar_2D_planet_flux (C++ code)
# Produces curves for data at specific longitude and latitude

import numpy as np
import matplotlib.pyplot as plt
import infofile
from string import split
from os import system

pi = 3.1415926585
solrad = 6.9e10
AU = 1.496e13

# Columns of .flux file

longcol = 0
latcol = 1
fluxcol = 2
darkcol = 3

# columns of .sky file

starfluxcol=2
altcol = 3
azcol = 4
hourcol = 5

# columns of .position files

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
positions = np.zeros((nstars+1,3,nfiles))

# Read .position files

for s in range(nstars):
    filename = prefix+'_'+starname[s]+'.position'
    data = np.genfromtxt(filename, usecols=(xcol,ycol,zcol))
    positions[s,:,:] = data.transpose()

filename = prefix+'_Planet.position'
data = np.genfromtxt(filename, usecols=(xcol,ycol,zcol))
positions[nstars,:,:] = data.transpose()

# Calculate separation between planet and each star

for s in range(nstars):
    for i in range(nfiles):
        sep = 0.0
        for k in range(3):
            sep =sep + (positions[nstars,k,i] - positions[s,k,i])**2
            
        distance[s,i] = sep
 
# Loop over files

for i in range(nfiles):

    # Create filename - how many zeros needed?
    num = str(i+1)
    k=np.log10(i+1)
    while (k<nzeros): 
        num = "0"+num
        k+=1        
    
    fluxfile = prefix +num+'.flux'
    skyfile = []
    for star in starname:        
        skyfile.append(prefix+'_'+star+'_'+num+'.sky')
        
            
    # Read in flux file
    # First header - time, position data etc

    f = open(fluxfile, 'r')

    line = f.readline()

    numbers = split(line)

    time[i]=float(numbers[0])
    nlat = int(numbers[1])
    nlong = int(numbers[2])
    
    
    f.close()
    
    print 'File ', str(i+1),' Time  ', time[i]

    # Read in rest of file
    
    data = np.genfromtxt(fluxfile, skiprows=1)
    
    if(i==0):
        # Pick from possible latitudes
        
        lat_possibles = np.unique(data[:,latcol])
                
        print "Fix value of latitude: here are the choices"
        for j in range (len(lat_possibles)):
            print '(',j,')', lat_possibles[j]
        latselect = input("Enter integer corresponding to desired value: ")        
        mylat= lat_possibles[latselect]
        
        long_possibles = np.unique(data[:,longcol])
        print "Fix value of longitude: here are the choices"
        for j in range (len(long_possibles)):
            print '(',j,')', long_possibles[j]
    
        longselect = input("Enter integer corresponding to desired value: ")
        mylong= long_possibles[longselect]
                                
        
    # Select flux only from the correct latitude and longitude    
    
    myentry = data[data[:,latcol]==mylat]
    myentry = myentry[myentry[:,longcol]==mylong]        
        
    fluxtot[i] = myentry[0,fluxcol]
    darkness[i] = myentry[0,darkcol]
    
           
    # Now do sky positions
    
    for s in range(nstars):
        filename = skyfile[s]
        
        f = open(filename, 'r')

        line = f.readline()

        numbers = split(line)

        time[i]=float(numbers[0])
        nlat = int(numbers[1])
        nlong = int(numbers[2])
    
    
        f.close()
    
        #print 'File '+filename+ ' read' 

        # Read in rest of file
    
        data = np.genfromtxt(filename, skiprows=1)
        
        # Select altitude, azimuth, hour angle only from the correct latitude and longitude    
    
        myentry = data[data[:,latcol]==mylat]
        myentry = myentry[myentry[:,longcol]==mylong]        
        
        flux[s,i] = myentry[0,starfluxcol]
        altitude[s,i] = myentry[0,altcol]
        azimuth[s,i] = myentry[0,azcol]
        hourangle[s,i] = myentry[0,hourcol]*180.0/pi
    
        skyoutput = 'skypos_'+prefix+num+'_latlong_'+str(mylat)+'_'+str(mylong)+'.png'        
        
        horizontal[s,i] = np.cos(altitude[s,i])*np.sin(azimuth[s,i])
        height[s,i] = np.sin(altitude[s,i])
        starsize[s,i] = 20.0*starradius[s]/(distance[s,i]*distance[s,i])

        # print horizontal[s,i], height[s,i], starsize[s,i]
        
    if(moviechoice == 'y'):
        # Plot sky position for this timestep
    
        fig2 = plt.figure(2)
        ax = fig2.add_subplot(111)
        ax.set_xlabel('Horizontal Position')
        ax.set_ylabel('Height')
        ax.set_ylim(0,1)
        ax.set_xlim(-1,1)
    
        for s in range(nstars):    
            plt.scatter(horizontal[s,i],height[s,i], marker='o',s=starsize[s,i],c=starcolor)
        
        plt.savefig(skyoutput, format= 'png')
        plt.clf()

    # Create SFD for this position from planck functions, separations and total flux
    
    SFDoutput = 'sfd_'+prefix+num+'_latlong_'+str(mylat)+'_'+str(mylong)+'.ng'        
    SFDtot[:] = 0.0
    
    for s in range(nstars):
        SFD[s,:] = pi*starradius[s]*starradius[s]*planck[s,:]/(distance[s,i]*distance[s,i])                


        if(np.sum(SFD[s,:])!=0.0):
            SFD[s,:] = SFD[s,:]*flux[s,i]/(np.sum(SFD[s,:])*dl)
        else:
            SFD[s,:] = 0.0
            
        SFDtot[:] = SFDtot+SFD[s,:]
        
    
    #Plot the SFD
    
    fig1 = plt.figure(1)
    ax = fig1.add_subplot(111)
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


# end of loop

# Plot flux

fluxfile = 'latlong_'+str(mylat)+'_'+str(mylong)+'_flux_'+prefix+'.png'    

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



# Plot darkness

darkfile = 'latlong_'+str(mylat)+'_'+str(mylong)+'_darkness_'+prefix+'.png'   

fig1 = plt.figure()
ax = fig1.add_subplot(111)
ax.set_xlabel('Time (yr)')
ax.set_ylabel('Period of Darkness (yr)')
plt.plot(time,darkness)

plt.savefig(darkfile, format= 'png')


# Plot altitude

altfile = 'latlong_'+str(mylat)+'_'+str(mylong)+'_altitude_'+prefix+'.png'  

fig1 = plt.figure()
ax = fig1.add_subplot(111)
ax.set_xlabel('Time (yr)')
ax.set_ylabel('Altitude (degrees)')
for s in range(nstars):
    plt.plot(time,altitude[s,:], label=starname[s])
ax.legend()

plt.savefig(altfile, format= 'png')

# Plot azimuth

azfile = 'latlong_'+str(mylat)+'_'+str(mylong)+'_azimuth_'+prefix+'.png'  

fig1 = plt.figure()
ax = fig1.add_subplot(111)
ax.set_xlabel('Time (yr)')
ax.set_ylabel('azimuth (degrees)')
for s in range(nstars):
    plt.plot(time,azimuth[s,:], label=starname[s])
ax.legend()

plt.savefig(azfile, format= 'png')

# Hour Angle

hourfile = 'latlong_'+str(mylat)+'_'+str(mylong)+'_hourangle_'+prefix+'_.png'  

fig1 = plt.figure()
ax = fig1.add_subplot(111)
ax.set_xlabel('Time (yr)')
ax.set_ylabel('Hour Angle (degrees)')
for s in range(nstars):
    plt.plot(time,hourangle[s,:], label=starname[s])
ax.legend()

plt.savefig(hourfile, format= 'png')

# Plot height and horizontal on same graph

# azfile = 'hh_'+prefix+'_latlong_'+str(mylat)+'_'+str(mylong)+'.png'    
# 
# fig1 = plt.figure()
# ax = fig1.add_subplot(111)
# ax.set_xlabel('Time (yr)')
# ax.set_ylabel('x(y) position')
# 
#     plt.plot(time,height, color='red')
#     plt.plot(time,horizontal, color='blue')
# 
# plt.savefig(azfile, format= 'png')

# Command for converting images into gifs - machine dependent

#convertcommand = '/opt/ImageMagick/bin/convert '
convertcommand = '/usr/bin/convert '

# Create movie if requested
if(moviechoice=='y'):
    print 'Creating animated gif of sky pattern, filename skymovie.gif'
    system(convertcommand +'skypos_'+prefix+'*_latlong_'+str(mylat)+'_'+str(mylong)+'.png skymovie.gif')
    
    print 'Creating animated gif of SFD, filename sfdmovie.gif'
    system(convertcommand +'sfd_'+prefix+'*_latlong_'+str(mylat)+'_'+str(mylong)+'.png sfdmovie.gif')
    
    
    if(deletechoice=='y'):
        print 'Deleting png files'
        system('rm skypos_'+prefix+'*.png')
        system('rm sfd_'+prefix+'*.png')                    


print 'Complete'


