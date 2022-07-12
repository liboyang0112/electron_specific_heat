#include "constants.h"
#include <iostream>
#include "Experiment.h"
#include "formulas.h"
#include "Plot.h"
#include <fstream>
#include <sstream>
#include "CalcMaterial.h"
#include <gsl/gsl_poly.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_filter.h>
#include <gsl/gsl_sf_fermi_dirac.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_math.h>
using namespace std;


int main(int argc, char const *argv[])
{
	Plot plot;
	plot.init(".");
	double risingTime = 0.5; // in us
	double droppingTime = 0.01; // in us
	double duration = 2; // in us
	double lifeLife= 800;  // excited state life time, in us
	double maxLight = 1; // brightest
	double pulseTime = 10; // time of pulse arrival
	double saturation = 1; // pump max
	double conversionRate = 5; // light -> pump
	std::vector<double> time;
	std::vector<double> excitingNumber;
	std::vector<double> light;
	std::vector<double> pulse;
	double startTime = 0;
	double endTime = 4;
	double interval = (endTime-startTime)/5000;
	bool pulseCame = 0;
	light.push_back(0);
	excitingNumber.push_back(0);
	time.push_back(startTime);
	pulse.push_back(0);
	for (double t = startTime+interval; t <= endTime; t+=interval)
	{
		time.push_back(t);
		double prevlight = light.back();
		if(t<=risingTime) light.push_back(prevlight+maxLight/risingTime*interval);
		else if(t>risingTime && t<=risingTime+duration) light.push_back(maxLight);
		else if(t>risingTime+duration && t<=risingTime+duration+droppingTime) light.push_back(prevlight-maxLight/droppingTime*interval);
		else if(t>risingTime+duration+droppingTime) light.push_back(0);
		double prevNumber = excitingNumber.back();
		prevNumber += prevlight*conversionRate*interval*(1-prevNumber/saturation);
		prevNumber -= prevNumber*interval/lifeLife;
		if(prevNumber > saturation) prevNumber = saturation;
		if(pulseCame==0 && t>pulseTime) {
			prevNumber=0;
			pulseCame = 1;
			pulse.push_back(1);
		}else{
			pulse.push_back(0);
		}
		excitingNumber.push_back(prevNumber);
	}
	plot.plotTripple("PumpLight","excited atom number","pulse",time,light,excitingNumber,pulse,"time/us","fibre.png");
	return 0;
}
