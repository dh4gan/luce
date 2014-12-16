/*
 * Constants.h
 *
 *  Created on: Sep 16, 2012
 *      Author: dhf
 *      NOTE ALL CONSTANTS ARE IN SI UNITS
 */

#ifndef CONSTANTS_H_
#define CONSTANTS_H_

const double pi = 3.141592;
const double twopi = 2.0*pi;
const double piby2 = pi/2.0;
const double Gsi = 6.67e-11;
const double Gmau = 1.0;  // Value of G for solar mass-AU units (time units...2pi units= 1 year)
const double AU = 1.496e11;
const double msol = 1.99e30;
const double lsol = 3.826e26;
const double rsol = 6.955e8;

const double fluxsol = lsol/(AU*AU);
const double fluxsolcgs = fluxsol*1000;
//double fluxsol = 1;



#endif /* CONSTANTS_H_ */
