/*
 * main.cpp
 *
 *  Created on: Jan 9, 2014
 *      Author: dh4gan
 *
 *	Reads in parameter files and runs N Body code
 *	where some objects have climate modelling done in tandem
 *	via LEBM modelling
 * // TODO - change description text
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

#include <fstream>
#include <sstream>
using namespace std;

int main(int argc, char* argv[])
    {

    double G = 1;
    double pi = 3.141592654;
    double twopi = 2.0*pi;
    double year =  3.1556926e7;

    int i, fileType;
    int snapshotNumber=0;

    double tStop;
    double timeunit,timeyr;
    double dtyr,dtunit;
    double tSnap, tMax, dtmax;


    Vector3D body_i_position;
    Vector3D body_i_velocity;

    System nBodySystem;
    vector<Body*> BodyArray;
    vector<Body*> Bodies;

    parFile input;
    FILE * outputfile;

    printf("  \n");
    printf("*********************************************** \n");
    printf("    NBODY 2D Flux Code \n ");
    printf("    Date Created : 4th December 2014 \n");
    printf("*********************************************** \n");
    printf("  \n");

    // Read in parameters file

    if (argc == 2)
	{
	string fileString = string(argv[1]);
	fileType = input.readParFile(fileString);
	}
    else
	{
	fileType = input.readParFile();
	if (fileType > 1)
	    {
	    return -1;
	    }
	}

    // Record parameter data

    tMax = input.maximumTime;
    tSnap = input.snapshotTime;

    if(input.restart)
	{
		cout << "Restart - Using vector data from nbody output" << endl;
	}

    //First loop through each of the bodies in the system and set them up
    //adding them to the BodyArray
    for (i = 0; i < input.number_bodies; i++)
	{
	if (fileType == 0 or input.restart)
	    {

	    cout << "setting up body from vectors" << endl;

	    body_i_position = input.getBodyPosition(i);
	    body_i_velocity = input.getBodyVelocity(i);


	    // If the Body is a Star, add a Star Object

	    if (input.BodyTypes[i] == "Star")
		{
		BodyArray.push_back(
			new Star(input.BodyNames[i], input.BodyTypes[i],
				input.Mass[i], input.Radius[i], body_i_position,
				body_i_velocity, input.luminosity[i],input.effectiveTemperature[i], input.nLambda));
		}

	    // If the Body is a Planet, add a Planet Object
	    if (input.BodyTypes[i] == "Planet")
		{
		BodyArray.push_back(
			new Planet(input.BodyNames[i], input.BodyTypes[i],
				input.Mass[i], input.Radius[i], body_i_position,
				body_i_velocity, input.albedo[i]));
		}

	    // If the Body is a PlanetSurface, add a PlanetSurface Object and set up surface grids

	    if (input.BodyTypes[i] == "PlanetSurface")
		{
		// Code will halt if initial temperature zero
		// this stops incomplete params files running successfully


		BodyArray.push_back(
			new PlanetSurface(input.BodyNames[i],
				input.BodyTypes[i], input.Mass[i],
				input.Radius[i], body_i_position,
				body_i_velocity, input.number_bodies,
				input.nLatitude, input.nLongitude,
				input.rotationPeriod[i], input.obliquity[i]));

		}


	    }
	else if (fileType == 1 and input.restart==false)
	    {
	   	printf("setting up body with orbital parameters \n");

	    // If the Body is a Star, add a Star Object
	    if (input.BodyTypes[i] == "Star")
		{

		BodyArray.push_back(
			new Star(input.BodyNames[i], input.BodyTypes[i],
				input.Mass[i], input.Radius[i],
				input.semiMajorAxis[i], input.eccentricity[i],
				input.inclination[i], input.longAscend[i],
				input.Periapsis[i], input.meanAnomaly[i], G,
				input.totalMass, input.luminosity[i], input.effectiveTemperature[i], input.nLambda));

		}

	    // If the Body is a Planet, add a Planet Object
	    if (input.BodyTypes[i] == "Planet")
		{
		BodyArray.push_back(
			new Planet(input.BodyNames[i], input.BodyTypes[i],
				input.Mass[i], input.Radius[i],
				input.semiMajorAxis[i], input.eccentricity[i],
				input.inclination[i], input.longAscend[i],
				input.Periapsis[i], input.meanAnomaly[i], G,
				input.totalMass, input.albedo[i]));

		}

	    // If the Body is a World, add a World Object and set up LEBM
	    if (input.BodyTypes[i] == "PlanetSurface")
		{


		BodyArray.push_back(
			new PlanetSurface(input.BodyNames[i],
				input.BodyTypes[i], input.Mass[i],
				input.Radius[i], input.semiMajorAxis[i],
				input.eccentricity[i], input.inclination[i],
				input.longAscend[i], input.Periapsis[i],
				input.meanAnomaly[i], G, input.totalMass,
				input.number_bodies, input.nLatitude,
				input.nLongitude, input.rotationPeriod[i],
				input.obliquity[i]));
		}

	    }

	printf("body %s set up \n", BodyArray.back()->getName().c_str());

	}

    // Set up System object using BodyArray

    printf("Setting up system %s \n", input.SystemName.c_str());

    nBodySystem = System(input.SystemName, BodyArray);

    // If the System is created from orbital parameters, set up vectors here

    if(fileType ==1 and input.restart==false)
	{
	nBodySystem.setupOrbits(input.orbitCentre);
	}

    // Calculate its initial properties
    nBodySystem.calcInitialProperties();
    nBodySystem.setHostBodies(input.orbitCentre);

    // Switch Planetary Illumination on/off
    nBodySystem.setIllumination(input.illumination);

    // Set up the outputs

    if (input.restart and snapshotNumber !=0)
	{
	outputfile = fopen(input.NBodyFile.c_str(), "a");
	}
    else
	{
	outputfile = fopen(input.NBodyFile.c_str(), "w");
	fprintf(outputfile, "Number of Bodies, %i \n", input.number_bodies);
	}


    // Now loop over snap shots, outputting the system data each time

    tStop = 0.0;
    tMax = tMax * twopi; // Convert maximum time to code units
    tSnap = tSnap*twopi; // Convert snapshot time to code units

    dtmax = 0.1*input.snapshotTime;

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
	while (timeunit < tStop)
	    {

	    timeyr = timeunit/twopi;
	    // Evolve the LEBMs in the system for the minimum timestep in seconds
	    nBodySystem.calc2DFlux(timeyr, dtyr);

	    // Evolve the NBody particles for the minimum timestep in code units
	    nBodySystem.evolveSystem(dtunit);

	    timeunit = timeunit + dtunit;

	    // Recalculate the minimum timestep
	    nBodySystem.calcNBodyTimestep(dtmax);
	    dtunit = nBodySystem.getTimestep();

	    dtyr = dtunit/twopi;


	    }

	printf("Time: %+.4E yr, Combined Timestep: %+.4E s, %+.4E units\n",timeyr, dtyr, dtunit);
	// Output data to files
	snapshotNumber++;
	timeyr = timeunit/twopi;

	// N Body data goes to a single file
	nBodySystem.outputNBodyData(outputfile, timeyr, input.orbitCentre);

	printf("N Body data written \n");
	// 2D Flux data goes to separate files for each World in the System

	nBodySystem.output2DFluxData(snapshotNumber, timeyr);

	printf("Flux data written \n");

	}

    //Close the file before returning

    fclose(outputfile);


    // Write integrated Data to files



    // Write .info file

    string fileString = input.SystemName+".info";
    FILE *infoFile = fopen(fileString.c_str(), "w");

    double globalFluxMax= 0.0;
    int nStars = nBodySystem.countStars();

    	fprintf(infoFile,"%i \n", snapshotNumber);
    	fprintf(infoFile,"%i \n", nStars);

    for (int s = 0; s < nBodySystem.getBodyCount(); s++)
	{
	if (BodyArray[s]->getType() == "Star")
	    {
	    fprintf(infoFile, "%s \n", BodyArray[s]->getName().c_str());
	    fprintf(infoFile, "%+.4E %+.4E %+.4E \n", BodyArray[s]->getRadius(),
		    BodyArray[s]->getTeff(),
		    BodyArray[s]->calculatePeakWavelength());
	    }

	if (BodyArray[s]->getType() == "PlanetSurface")
	    {
	    if (BodyArray[s]->getFluxMax() > globalFluxMax)
		{
		globalFluxMax = BodyArray[s]->getFluxMax();
		}

	    }
	}

    	fprintf(infoFile, "%+.4E \n", globalFluxMax);
    	fclose(infoFile);



    return 0;
    }

