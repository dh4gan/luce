/*
 * PlanetSurface.cpp
 *
 *  Created on: Nov 28 2014
 *      Author: dh4gan
 */

#include "PlanetSurface.h"
#include <stdio.h>

PlanetSurface::PlanetSurface() :
		Body() {
	nStars = 1;
	nLatitude = 100;
	nLongitude = 100;
	Pspin = 2.7e-3;
	obliquity = 0.5123;

	initialiseArrays();

}
PlanetSurface::PlanetSurface(string &namestring, string &typestring, double &m,
		double &rad, Vector3D &pos, Vector3D &vel, int &nstar, int &nlat,
		int &nlong, double &spin, double &obliq) :
		Body(namestring, typestring, m, rad, pos, vel) {
	nStars = nstar;
	nLatitude = nlat;
	nLongitude = nlong;
	Pspin = spin;
	obliquity = obliq;

	initialiseArrays();

}

PlanetSurface::PlanetSurface(string &namestring, string &typestring, double &m,
		double &rad, double semimaj, double ecc, double inc, double trueAnom,
		double longascend, double argper, double G, double totalMass,
		int &nstar, int &nlat, int &nlong, double &spin, double &obliq) :
		Body(namestring, typestring, m, rad, semimaj, ecc, inc, trueAnom,
				longascend, argper, G, totalMass) {

	nStars = nstar;
	nLatitude = nlat;
	nLongitude = nlong;
	Pspin = spin;
	obliquity = obliq;

	initialiseArrays();

}

PlanetSurface::~PlanetSurface() {

}

void PlanetSurface::initialiseArrays() {
	/**
	 * Written 01/12/14 by dh4gan
	 * Initialises 1D, 2D and 3D arrays to store information
	 *
	 */

	// Initialise array values
	double pi = 3.14159265285;
	double dlat = pi / float(nLatitude);
	double dlong = 2.0 * pi / float(nLongitude);

	// 1D arrays
	for (int j = 0; j < nLongitude; j++) {
		longitude[j] = j * dlong;
	}

	for (int k = 0; k < nLatitude; k++) {
		latitude[k] = k * dlat;
	}

	// 2D arrays

	for (int j = 0; j < nLongMax; j++) {

		for (int istar = 0; istar < nStarMax; istar++) {
			hourAngle[istar][j] = 0.0;
		}

		for (int k = 0; k < nLatMax; k++) {
			fluxtot[j][k] = 0.0;
			integratedflux[j][k] = 0.0;
			darkness[j][k] = 0.0;

		}
	}

	//3D arrays

	for (int istar = 0; istar < nStarMax; istar++) {
		for (int j = 0; j < nLongMax; j++) {

			for (int k = 0; k < nLatMax; k++) {
				flux[istar][j][k] = 0.0;
				altitude[istar][j][k] = 0.0;
				azimuth[istar][j][k] = 0.0;

			}

		}

	}

	cout << "Arrays initialised " << endl;

}

void PlanetSurface::resetFluxTotals() {
	/*
	 * Written 1/12/14 by dh4gan
	 * This resets the flux arrays so they can be used in the next timestep
	 *
	 */

	for (int j = 0; j < nLongitude; j++) {
		for (int k = 0; k < nLatitude; k++) {
			fluxtot[j][k] = 0.0;
			for (int istar = 0; istar < nStars; istar++) {
				flux[istar][j][k] = 0.0;

			}
		}
	}

}

// Simple function to stop acos becoming infinite if abs(x) > 1.0
static double safeAcos(double x) {
	if (x < -1.0)
		x = -1.0;
	else if (x > 1.0)
		x = 1.0;
	return acos(x);
}

// Simple function to find array indices corresponding with input longitude and latitude
void PlanetSurface::findSurfaceLocation(double &longitude, double &latitude) {

	int j;
	double longtry, lattry, delta_long, delta_lat;
	double pi = 3.14159265285;
	double dlat = pi / float(nLatitude);
	double dlong = 2.0 * pi / float(nLongitude);

	delta_long = 1.0e30;
	delta_lat = 1.0e30;

	for (j = 0; j < nLongitude; j++) {
		longtry = j * dlong;

		if (fabs(longitude - longtry) < delta_long) {
			delta_long = fabs(longitude - longtry);
			iLongPick = j;
		}

	}

	for (j = 0; j < nLatitude; j++) {
		lattry = j * dlat;

		if (fabs(latitude - lattry) < delta_lat) {
			delta_lat = fabs(latitude - lattry);
			iLatPick = j;
		}

	}

}

void PlanetSurface::calcLongitudeOfNoon(Body* &star, int &istar) {
	/*
	 * Written 01/12/14 by dh4gan
	 * Calculates the initial longitude of noon of a star on the planet's surface
	 *
	 */

	// Calculate current noon longitudes for each star on planet's surface
	Vector3D planetpos = getPosition();
	Vector3D starpos = star->getPosition();

	// Find relative vector from planet to star
	planetpos = planetpos.relativeVector(starpos);

	// Decompose this to find longitude of noon

	noon[istar] = atan2(planetpos.elements[1], planetpos.elements[0]);

	cout << "Star " << star->getName()
			<< " is initially overhead for longitude " << noon[istar] << endl;

}

void PlanetSurface::writeFluxFile(FILE* outputFlux, double &time)

{
	/**
	 * Written 01/12/14 by dh4gan
	 * Writes a snapshot of the total flux to file
	 *
	 */

	fprintf(outputFlux, "%+.4E %i %i \n", time, nLatitude, nLongitude);
	for (int j = 0; j < nLongitude; j++) {
		for (int k = 0; k < nLatitude; k++) {
			fprintf(outputFlux, "%+.4E %+.4E %+.4E %+.4E \n", longitude[j],
					latitude[k], fluxtot[j][k], darkness[j][k]);
		}
	}
	fflush(outputFlux);
	fclose(outputFlux);

}

void PlanetSurface::writeSkyFile(FILE* outputSky, int &istar, double &time) {

	fprintf(outputSky, "%+.4E %i %i \n", time, nLatitude, nLongitude);
	for (int j = 0; j < nLongitude; j++) {
		for (int k = 0; k < nLatitude; k++) {
			fprintf(outputSky, "%+.4E  %+.4E  %+.4E %+.4E %+.4E %+.4E \n",
					longitude[j], latitude[k], flux[istar][j][k],
					altitude[istar][j][k], azimuth[istar][j][k],
					hourAngle[istar][j]);
		}
	}
	fflush(outputSky);
	fclose(outputSky);

}

void PlanetSurface::writeToLocationFile(FILE* outputLocation, int &istar,
		int &iLongPick, int &iLatPick, double &time) {
	/*
	 * Written 01/12/14 by dh4gan
	 * writes prepicked location data for a timestep
	 */

	fprintf(outputLocation,
			"%+.4E  %+.4E  %+.4E  %+.4E  %+.4E  %+.4E  %+.4E \n", time,
			longitude[iLongPick], latitude[iLatPick],
			flux[istar][iLongPick][iLatPick],
			altitude[istar][iLongPick][iLatPick],
			azimuth[istar][iLongPick][iLatPick], hourAngle[istar][iLongPick]);
	fflush(outputLocation);

}

void PlanetSurface::writeIntegratedFile(FILE *outputFlux) {

	fprintf(outputFlux, "%i %i \n", nLatitude, nLongitude);

	for (int j = 0; j < nLongitude; j++) {
		for (int k = 0; k < nLatitude; k++) {
			fprintf(outputFlux, "%+.4E  %+.4E  %+.4E  %+.4E \n", longitude[j],
					latitude[k], integratedflux[j][k], darkness[j][k]);
		}
	}
	fflush(outputFlux);
	fclose(outputFlux);
}

void PlanetSurface::calcFlux(int &istar, Body* &star, double &eclipseFraction,
		double &time, double &dt, double &fluxmax, double &fluxunit) {
	/*
	 * Written 1/12/14 by dh4gan
	 * Calculates the flux received from one star over the entire surface
	 */

	double rdotn, declination, long_apparent;
	double magpos, lstar, fluxtemp;
	Vector3D pos, unitpos, planetpos, starpos;
	Vector3D zvector(0.0, 0.0, 1.0);

	double pi = 3.14159265285;

	// Get position of planet relative to star
	planetpos = getPosition();
	starpos = star->getPosition();
	pos = (starpos).relativeVector(planetpos);
	magpos = pos.magVector();
	unitpos = pos.unitVector();
	lstar = star->getLuminosity();

	// Declination of the Sun - angle between planet's position vector and equator (at noon)

	Vector3D decVector(unitpos.elements[0], unitpos.elements[1],
			unitpos.elements[2]);

	// Rotate this vector if planet has non-zero obliquity
	if (obliquity != 0.0) {
		decVector.rotateX(obliquity);
	}

	rdotn = unitpos.dotProduct(decVector);
	declination = safeAcos(rdotn);

	// Now begin calculation over latitude and longitude points

	// Loop over longitude and latitude

	for (int j = 0; j < nLongitude; j++) {
#pragma omp parallel default(none) \
shared(j,longitude,latitude,hourAngle,flux,nLatitude)\
shared(noon,altitude,azimuth,time,obliquity,nStars) \
shared(fluxsol,eclipseFraction,darkness,integratedflux,dt) \
	private(k,s,fluxtemp) \
	reduction(max: fluxmax)
		{
			// Rotate planet according to its spin period

			long_apparent = fmod(
					longitude[j] - noon[istar] + 2.0 * pi * time / Pspin,
					2.0 * pi);

			// Calculate hour angle - distance between current longitude and noon

			// Distance between longitude and noon = hourAngle
			// hour angle between -180 and +180

			Vector3D longSurface(cos(long_apparent), sin(long_apparent), 0.0);

			rdotn = unitpos.dotProduct(longSurface);
			hourAngle[istar][j] = safeAcos(rdotn);

			if ((unitpos.crossProduct(longSurface)).dotProduct(zvector) > 0.0) {
				hourAngle[istar][j] = -hourAngle[istar][j];
			}

#pragma omp for schedule(runtime) ordered*/
			for (int k = 0; k < nLatitude; k++) {

				// construct surface normal vector at this longitude and latitude

				Vector3D surface(sin(latitude[k]) * cos(long_apparent),
						sin(latitude[k]) * sin(long_apparent),
						cos(latitude[k]));

				surface = surface.unitVector();

				// If necessary, rotate surface vector by obliquity

				if (obliquity != 0.0) {
					surface.rotateX(obliquity);
				}

				// take the dot product with the unit position vector

				rdotn = unitpos.dotProduct(surface);

				// Calculate fluxes
				// if position.surface is less than zero, long/lat location is not illuminated

				if (rdotn > 0.0) {

					fluxtemp = lstar * rdotn / (4.0 * pi * magpos * magpos);
				}

				flux[istar][j][k] = flux[istar][j][k]
						+ fluxtemp * (1.0 - eclipseFraction) * fluxunit;

				fluxtot[j][k] = fluxtot[j][k] + flux[istar][j][k];
				//cout << s << "  "<<flux[j][k]<< "  " << eclipseFraction[s] << endl;

				if (fluxtot[j][k] > fluxmax) {
					fluxmax = fluxtot[j][k];
				}

				// Calculate altitude and azimuthal position on sky, and angular size
				// Formulae not exactly as used normally
				// As we measure latitudes from 0 to 180, not -90 to 90, some sines have become cosines, and vice versa
				// (some extra plus and minus signs as a result)

				//altitude = cos(declination) * cos(latitude) * cos(hourAngle) + sin(declination) * sin(latitude);  // These calculations assume lat->-90,90 deg

				altitude[istar][j][k] = -cos(declination)
						* cos(hourAngle[istar][j]) * sin(latitude[k])
						+ sin(declination) * cos(latitude[k]); // These calculations assume lat->0,180 deg
				altitude[istar][j][k] = asin(altitude[istar][j][k]);

				// Azimuth angle measured from north, clockwise

				if (cos(altitude[istar][j][k]) * sin(latitude[k]) != 0.0) {

					//azimuth = (sin(declination) - sin(altitude)*cos(latitude))/(cos(altitude)*cos(latitude)); // These calculations assume lat->-90, 90 deg
					azimuth[istar][j][k] = (sin(altitude[istar][j][k])
							* sin(latitude[k]) - sin(declination))
							/ (cos(altitude[istar][j][k]) * sin(latitude[k])); // These calculations assume lat->0, 180 deg
					azimuth[istar][j][k] = safeAcos(azimuth[istar][j][k]);
				} else {
					azimuth[istar][j][k] = 0.0;
				}

				// If hour angle positive (afternoon) azimuth is > 180

				if (hourAngle[istar][j] > 0.0) {
					azimuth[istar][j][k] = 2.0 * pi - azimuth[istar][j][k];
				}

			}
		}

	}
}

void PlanetSurface::calcIntegratedQuantities(double &dt) {
	/*
	 * Written 03/12/14 by dh4gan
	 * Updates the darkness array based on the total flux received at a given timestep
	 *
	 */

	for (int j = 0; j < nLongitude; j++) {
#pragma omp parallel default(none) \
	shared(j,longitude,latitude,hourAngle,flux,nLatitude)\
	shared(noon,altitude,azimuth,time,obliquity,nStars) \
	shared(fluxsol,eclipseFraction,darkness,integratedflux,dt) \
		private(k,s,fluxtemp) \
		reduction(max: fluxmax)
		{
#pragma omp for schedule(runtime) ordered*/
			for (int k = 0; k < nLatitude; k++) {

				integratedflux[j][k] = integratedflux[j][k]
						+ fluxtot[j][k] * dt;

				// If flux zero, add to darkness counter
				if (fluxtot[j][k] < 1.0e-6) {

					darkness[j][k] = darkness[j][k] + dt;
				}
			}
		}
	}

}

