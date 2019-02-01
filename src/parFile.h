/*
 * parFile.h
 *
 *  Created on: Sep 23, 2013
 *      Author: davidharvey
 */

#ifndef PARFILE_H_
#define PARFILE_H_

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "System.h"

#include <fstream>
#include <sstream>
using namespace std;

class parFile {
public:
	parFile( );
	parFile(string name);

	string NBodyFile;
	string parFileName;
	string SystemName;
	string fileType;

	// TODO - need to put prepicked location data in as parameters too
	bool restart;
	bool illumination;
	bool tidal;
	bool fullOutput;

	int number_bodies;
	int snapshotNumber;
	int nLongitude, nLatitude, nLambda;
	double snapshotTime;
	double maximumTime;
	double systemTime;
	double totalMass;

	vector<string> parameters;
	vector<string> BodyNames;
	vector<string> BodyTypes;

	vector<double> Mass;
	vector<double> Radius;

	vector<double> x_position;
	vector<double> y_position;
	vector<double> z_position;

	vector<double> x_velocity;
	vector<double> y_velocity;
	vector<double> z_velocity;

	vector<double> semiMajorAxis;
	vector<double> eccentricity;
	vector<double> inclination;
	vector<double> longAscend;
	vector<double> Periapsis;
	vector<double> meanAnomaly;

	vector<double> luminosity;
	vector<double> effectiveTemperature;
	vector<double> albedo;

	vector<double> rotationPeriod;
	vector<double> obliquity;

	double longTrack;
	double latTrack;

	vector<int> orbitCentre;

	Vector3D getBodyPosition(int index);
	Vector3D getBodyVelocity(int index);

	int readParFile();
	int readParFile(string fileName);

	void readPosFile();
	void readOrbFile();
	int parType();
	int parType(string fileName);

	void setupRestartPositions();


};


#endif /* PARFILE_H_ */