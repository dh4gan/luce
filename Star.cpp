/*
 * Star.cpp
 *
 *  Created on: 13 Sep 2012
 *      Author: dhf
 */

#include "Star.h"
#include "RadiationConstants.h"
#include <iostream>
/* Constructors and Destructor */

Star::Star() :
	Body() {
}
Star::Star(string &namestring, double &m, double &rad, Vector3D  &pos, Vector3D  &vel) :
	Body(namestring, m, rad, pos, vel) {
	nlambda = 1000;
	I_lambda = vector<double> (nlambda, 0.0);
	Teff = 5000.0;
	   cout << "A star is born  " << Teff << endl;

}

Star::Star(string &namestring, double &m, double &rad,Vector3D  &pos, Vector3D &vel, double &T, int &n) :
	Body(namestring, m, rad, pos, vel) {
	Teff = T;
	nlambda = n;
	I_lambda = vector<double> (nlambda, 0.0);

}

Star::~Star() {
}

/* Calculation Methods */

void Star::calculateBlackbodySpectrum() {
	/* Written by D Forgan
	 *
	 * This method calculates a simple Planck Function using the Star's effective temperature
	 *
	 */

	int i;
	double dlambda, lambda;
	/* Loop over all wavelengths */

	dlambda = (lambda_max - lambda_min) / float(nlambda);

	for (i = 0; i < nlambda; i++) {
		// Calculate lambda

		lambda = lambda_min + i * dlambda;

		// Calculate Planck Function at this lambda

		I_lambda[i] = exp(hplanck * c / (lambda * k_B * Teff)) - 1.0;
		I_lambda[i] = 2.0 * hplanck * c * c / (pow(lambda, 5) * I_lambda[i]);

	}

}

double Star::calculatePeakWavelength()
{
	/* Written by D Forgan
	 *
	 * This method returns the peak wavelength of emission according to Wien's Law
	 * (in centimetres)
	 *
	 */

	if(Teff > 0.0)
	{
		return 2.8977685e-1/Teff;
	}
	else
	{
		return 0.0;
	}
}

void Star::calculateSingleHZ() {
	/* Written by D Forgan
	 *
	 * This method calculates the distances from the star to the Inner and Outer Habitable Zone radii
	 * based on Effective Temperature (see Underwood et al (2003))
	 *
	 */

	double Sinner, Souter;

	// Calculate critical inner and outer HZ flux (in Solar Units)

	Sinner = 4.19e-8 * Teff * Teff - 2.139e-4 * Teff + 1.268;
	Souter = 6.19e-9 * Teff * Teff - 1.319e-5 * Teff + 0.2341;

	// Convert this flux into a distance

	innerHZ = sqrt(luminosity/Sinner);
	outerHZ = sqrt(luminosity/Souter);

}
