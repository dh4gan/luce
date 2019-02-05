/*
 * main.cpp
 *
 *  Created on: Jan 9, 2014
 *      Author: dh4gan
 *
 *	Reads in parameter files and runs N Body code
 *	where some objects have surface flux patterns computed in tandem
 */

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "System.h"
#include "Star.h"
#include "Planet.h"
#include "PlanetSurface.h"
#include "parFile.h"

#include <chrono>
#include <time.h>
#include <fstream>
#include <sstream>

using namespace std;

int main(int argc, char* argv[])
    {

    double G = Gmau;
    double unit2sec = year/twopi;

    int i, fileType, nTime;
    int snapshotNumber=0;

    double tStop;
    double timeunit,timeyr;
    double dtyr,dtunit, dtflux;
    double tSnap, tMax, dtmax;


    Vector3D body_i_position;
    Vector3D body_i_velocity;

    System nBodySystem;
    vector<Body*> BodyArray;
    vector<Body*> Bodies;

    parFile input;
    FILE * outputfile;
        
        
        printf("  \n");
        printf("%s",screenBar.c_str());
        printf("\tLuce - 2D surface flux patterns on bodies inside an N Body integrator\n");
        printf("\t\t\tVersion: %s\n", VERSION);
        printf("\t\t\tCompiled: %s\n", __DATE__);
        printf("\t\t\tgit commit: %s \n", GIT_HASH);
        printf("%s",screenBar.c_str());
        printf("  \n");
        
        // Record start time
        std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
        

        // Read in parameters file
        
        if (argc == 2)
        {
            string fileString = string(argv[1]);
            
            printf("\tReading file %s\n", fileString.c_str());
            input.readFile(fileString);
        }
        else
        {
            
            input.readFile();
            
        }
        
        // Check and display input parameters
        input.checkParameters();  //TODO - rewrite parFile::checkParameters and displayParameters
        input.displayParameters();

        // Retrieve parameter data from parFile object
        
        tMax = input.getDoubleVariable("MaximumTime");
        tSnap = input.getDoubleVariable("SnapshotTime");
        
        bool restart = input.getBoolVariable("Restart");
        bool fullOutput = input.getBoolVariable("fullOutput");
        string systemName = input.getStringVariable("SystemName");
        

    nTime = int(tMax/tSnap)+1;

    if(restart)
	{
		cout << "Restart - Using vector data from nbody output" << endl;
	}
        
        
        printf("Creating bodies \n");
        
        //First loop through each of the bodies in the system and set them up
        //adding them to the BodyArray
        for (i = 0; i < input.getIntVariable("Number_Bodies"); i++)
        {
            
            printf("Creating Body %s \n",input.getStringVariable("BodyName",i).c_str());
            if (input.getStringVariable("BodyType",i).compare("Star")==0)
            {
                BodyArray.push_back(new Star(input, i, G));
            }
            else if (input.getStringVariable("BodyType",i) == "Planet")
            {
                BodyArray.push_back(new Planet(input, i, G));
            }
            else if(input.getStringVariable("BodyType",i) == "PlanetSurface")
            {
                BodyArray.push_back(new PlanetSurface(input, i, G));
                
              /*     if(input.getBoolVariable("Restart"))
               *
               * {
               *     cout << "Reading Temperature data for World " << BodyArray.back()->getName() << endl;
               *     snapshotNumber = BodyArray.back()->getRestartParameters();
                }*/
                
            }
            
        }
        
        
        printf("Setting up system %s \n", systemName.c_str());
        
        nBodySystem = System(systemName, BodyArray);
        
        // If the System is created from orbital parameters, set up vectors here
        vector<int> orbitCentre;
        
        for (int i=0; i<input.getIntVariable("Number_Bodies"); i++)
        
        {
            orbitCentre.push_back(input.getIntVariable("OrbitCentre",i));
            
        }
        nBodySystem.setHostBodies(orbitCentre);
        
        if(input.getStringVariable("ParType").compare("Orbital")==0 and input.getBoolVariable("Restart")==false)
        {
            nBodySystem.setupOrbits(orbitCentre);
        }
        
        // Calculate the system's initial properties
        nBodySystem.calcInitialProperties();
        
        // Switch Planetary Illumination on/off
        nBodySystem.setIllumination(input.getBoolVariable("PlanetaryIllumination"));
        
        // Set up the outputs
        
        if (restart and snapshotNumber !=0)
        {
            outputfile = fopen(input.getStringVariable("NBodyOutput").c_str(), "a");
        }
        else
        {
            outputfile = fopen(input.getStringVariable("NBodyOutput").c_str(), "w");
            fprintf(outputfile, "Number of Bodies, %i \n", input.getIntVariable("Number_Bodies"));
        }


    nBodySystem.initialise2DFluxOutput(systemName);
    nBodySystem.setFluxOutput(fullOutput);

    if(fullOutput)
    {
    	printf("Run will produce full output \n");
    }

    // Now loop over snap shots, outputting the system data each time

    tStop = 0.0;
    tMax = tMax * twopi; // Convert maximum time to code units
    tSnap = tSnap*twopi; // Convert snapshot time to code units

    dtmax = 0.1*tSnap;

    // Timesteps will be calculated in NBody units, and converted back to years for Flux calculation

    // Calculate the N Body timestep, which has a maximum minimum LEBM timestep for all worlds and NBody timestep

    nBodySystem.calcNBodyTimestep(dtmax);
    dtunit = nBodySystem.getTimestep();

    dtyr = dtunit/twopi;
    timeunit = 0.0;

    printf("System set up: Running \n");
    while (timeunit < tMax)
	{
	tStop = timeunit + tSnap;

	dtflux = 0.0;

	while (timeunit < tStop)
	    {

	    // Evolve the NBody particles for the minimum timestep in code units
	    nBodySystem.evolveSystem(dtunit);

	    timeunit = timeunit + dtunit;
	    dtflux = dtflux + dtyr;

	    // Recalculate the minimum timestep
	    nBodySystem.calcNBodyTimestep(dtmax);
        
	    dtunit = nBodySystem.getTimestep();
	    timeyr = timeunit/twopi;
	    dtyr = dtunit/twopi;

	    }

	printf("Time: %+.4E yr, N-Body Timestep: %+.4E years, %+.4E units\n",timeyr, dtyr, dtunit);

	// Calculate 2D Fluxes
	nBodySystem.calc2DFlux(timeyr, dtflux);

	// Output data to files
	snapshotNumber++;
	timeyr = timeunit/twopi;

	// N Body data goes to a single file
	nBodySystem.outputNBodyData(outputfile, timeyr, orbitCentre);

	// 2D Flux data goes to separate files for each World in the System
	nBodySystem.output2DFluxData(snapshotNumber, timeyr, systemName);

	}

    // Close the N Body file

    fclose(outputfile);

    // Write integrated Data to files

    nBodySystem.outputIntegratedFluxData();
    nBodySystem.outputInfoFile(snapshotNumber);
        
    // Simulation has now ended
    // Write elapsed runtime to the screen
        
    std::chrono::high_resolution_clock::time_point finish = std::chrono::high_resolution_clock::now();
        
    std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double> >(finish - start);
        
    printf("%s",screenBar.c_str());
    printf("Run %s complete \n", nBodySystem.getName().c_str());
    printf("Wall Clock Runtime: %f s \n", time_span.count());
    printf("%s",screenBar.c_str());
        

    return 0;
    }

