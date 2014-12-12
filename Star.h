/*
 * Star.h
 *
 *  Created on: Nov 8, 2012
 *      Author: dh4gan
 */

#include "Body.h"

using namespace std;

#ifndef STAR_H_
#define STAR_H_

class Star: public Body {
public:

    Star();

    // Position Constructors



    Star(string &namestring, string &typestring, double &m, double &rad,
	    Vector3D &pos, Vector3D &vel);



    //Without Luminosity
    Star(string &namestring, string &typestring, double &m, double &rad,
	    Vector3D &pos, Vector3D &vel, double &T, int &n);

    // With Luminosity
    Star(string &namestring, string &typestring, double &m, double &rad,
    	    Vector3D &pos, Vector3D &vel, double &lum, double &T, int &n);

    // Orbital Constructors

    // Without Luminosity
    Star(string &namestring, string &typestring, double &m, double &rad,
	    double semimaj, double ecc, double inc, double longascend,
	    double argper, double meananom, double G, double totalMass,
	    double &T, int &n);

    // With Luminosity
    Star(string &namestring, string &typestring, double &m, double &rad,
	    double semimaj, double ecc, double inc, double longascend,
	    double argper, double meananom, double G, double totalMass,
	    double &lum, double &T, int &n);

	virtual ~Star();

	void setLuminosity(double lum){luminosity = lum;}
	double getLuminosity() {return luminosity;}

	void setTeff(double T){Teff=T;}
	double getTeff() {return Teff;}

	void setLambdaMin(double l){lambda_min=l;}
	double getLambdaMin() {return lambda_min;}

	void setLambdaMax(double l){lambda_max=l;}
	double getLambdaMax() {return lambda_max;}

	void setNLambda(int n){nlambda = n;}
	int getNLambda(){return nlambda;}

	void setInnerHZ(double r){innerHZ = r;}
	double getInnerHZ(){return innerHZ;}

	void setOuterHZ(double r){outerHZ = r;}
	double getOuterHZ(){return innerHZ;}

	vector<double> getILambda(){return I_lambda;}

	// Standard cloning method
	virtual Star* Clone() { return new Star(*this); }

	// Calculation Methods

	void calcLuminosityStefanBoltzmann();
	void calculateBlackbodySpectrum();
	double calculatePeakWavelength();
	void calculateSingleHZ();

protected:

	double luminosity; // Luminosity of Star IN SOLAR UNITS
	double Teff;
	double lambda_min;
	double lambda_max;
	double innerHZ;
	double outerHZ;

	int nlambda;
	vector<double> I_lambda;

};

#endif /* STAR_H */

