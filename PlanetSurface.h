/*
 * PlanetSurface.h
 *
 *  Created on: Nov 28, 2014
 *      Author: dh4gan
 */

#include <iostream>
#include "Body.h"

using namespace std;

#ifndef PLANETSURFACE_H_
#define PLANETSURFACE_H_

class PlanetSurface: public Body {
public:

	PlanetSurface();
	PlanetSurface(string &namestring, string &typestring, double &m, double &rad, Vector3D &pos,
			Vector3D &vel, int &nstar, int &nlat, int &nlong, double &spin,
			double &obliq);
	PlanetSurface(string &namestring, string &typestring, double &m, double &rad, double semimaj,
			double ecc, double inc, double trueAnom, double longascend,
			double argper, double G, double totalMass, int &nstar, int &nlat, int &nlong, double &spin,
			double &obliq);
	virtual ~PlanetSurface();

	/* Accessors */


	int getNLongitude(){return nLongitude;}
	int getNLatitude(){return nLatitude;}
	int getNStars(){return nStars;}
	int getLongPick(){return iLongPick;}
	int getLatPick(){return iLatPick;}

	double getPSpin(){return Pspin;}
	double getObliquity(){return obliquity;}


	void setNLongitude(int nlong){nLongitude=nlong;}
	void setNLatitude(int nlat){nLatitude=nlat;}
	void setNStars(int s){nStars=s;}

	void setPSpin(double spin){Pspin=spin;}
	void setObliquity(double obliq){obliquity = obliq;}

	// Standard cloning method
	virtual PlanetSurface* Clone() { return new PlanetSurface(*this); }

	// Other Methods

	void initialiseArrays();
	void initialiseOutputVariables(string prefixString, vector<Body*> stars);

	void resetFluxTotals();

	void findSurfaceLocation(double &longitude, double &latitude);
	void calcLongitudeOfNoon(Body* &star, int &istar);

	void calcFlux(int &istar, Body* &star, double &eclipseFraction, double &time, double &dt,double &fluxmax, double &fluxunit);
	void calcIntegratedQuantities(double &dt);

	void writeFluxFile(int &snapshotNumber, double &time);
	void writeSkyFile(FILE* outputSky, int &istar, double &time);
	void writeToLocationFiles(double &time);

	void writeIntegratedFile();

	static const int nStarMax = 10;
	static const int nLatMax = 500;
	static const int nLongMax = 500;

protected:

	int nStars;
	int nLatitude;
	int nLongitude;

	int iLatPick;
	int iLongPick;

	double Pspin;
	double obliquity;

	double noon[nStarMax];

	double longitude[nLongMax];
	double latitude[nLatMax];

	double flux[nStarMax][nLongMax][nLatMax];
	double altitude[nStarMax][nLongMax][nLatMax];
	double azimuth[nStarMax][nLongMax][nLatMax];

	double hourAngle[nStarMax][nLongMax];

	double fluxtot[nLongMax][nLatMax];
	double integratedflux[nLongMax][nLatMax];
	double darkness[nLongMax][nLatMax];

	FILE *integratedFile, *fluxFile, *locationFile[nStarMax];

};

#endif /* PLANETSURFACE_H */

