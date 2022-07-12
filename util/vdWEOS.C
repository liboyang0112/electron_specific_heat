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
	double Tc = 647.5; //K
	double rho0= 1000;  //kg/m3   //Affect the value of a and b
	double Tref = 273.15; //K   //Affect the value of a and b
	double molV = 22.4e-3; // m3/mol
	double atm = 1e5; //Pa
	double molmass = 18e-3; // kg/mol
	double moldensity = rho0/molmass; //mol/m3
	double r = 0;//atm/R/Tref/moldensity approx 0;
	double A = 27./8*Tc/Tref;
	double y = -sqrt(pow(r+A,2)-4*A)+r+A;
	double x = 1-y/2/A;
	double b = x/moldensity;//0.03049e-3; m3/mol
	double a = 27.*R*b*Tc/8;
	double Ttest = Tref;
	double ggamma = 4./3;

struct param{
};

double func(double x, void* p){
	double ret = x*R*Ttest/(1-b*x)-a*x*x-1e5;
	//printf("ret=%e\n", ret);
	return x*R*Ttest/(1-b*x)-a*x*x;
}

int main(int argc, char const *argv[])
{
	Plot plot;
	plot.init(".");



	printf("r= %e, a=%e, b=%e, maxDensity=%e, moldensity=%e\n", r,a,b,molmass/b, moldensity);
	struct param par;
	const gsl_root_fsolver_type * Tg = gsl_root_fsolver_brent;
	gsl_root_fsolver * s = gsl_root_fsolver_alloc (Tg);
	gsl_function F;
	F.function = &func;
	F.params = &par;

	double rhosolved = 0;
	double rhotmp = 0;
	std::vector<double> density;
	std::vector<double> pressure;
	std::vector<double> vol;
	std::vector<double> soundSpeed;
	double max = 1.001*rho0/molmass;
	for (double i = 1e-4; i <= max+1e-4; i+=max/100)
	{
		vol.push_back(1./i);
		density.push_back(i*b);
		pressure.push_back(i*R*Ttest/(1-b*i)-a*i*i);
		soundSpeed.push_back(sqrt(fabs(ggamma*Ttest*R/pow(1-b*i,2)-2*a*i)/molmass));
	}
	plot.plotSingle("Pressure",density,pressure,"density","dEOS.pdf");
	plot.plotSingle("volumn",vol,pressure,"density","vEOS.pdf");
	plot.plotSingle("C",density,soundSpeed,"density","soundSpeed.pdf");
	gsl_root_fsolver_set (s, &F, 40000,1./b-1e-5);
	while(1){
		gsl_root_fsolver_iterate(s);
		rhosolved=gsl_root_fsolver_root(s);
		if(fabs(rhosolved-rhotmp)<1e-100) break;
		rhotmp=rhosolved;
	}

	printf("solved rho = %f\n", rhosolved*molmass);
	printf("solved c1 = %f\n", sqrt(fabs(ggamma*Ttest*R/pow(1-b*rhosolved,2)-2*a*rhosolved)/molmass));
	printf("gas c1 = %f\n", sqrt(fabs(ggamma*Ttest*R/pow(1-b*0,2)-2*a*0)/molmass));
	printf("c2 = %f\n", sqrt(fabs(ggamma-2+rhosolved*2*b)/pow(1-rhosolved*b,2)*R*Ttest/molmass));
	return 0;
}