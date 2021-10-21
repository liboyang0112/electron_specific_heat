#include "Experiment.h"
#include "Plot.h"
#include <iostream>
const double eV = 11600;
int main(int argc, char const *argv[])
{
	Experiment exp;

	Plot plot;
	plot.init(".");
	exp.setEnvironment(1e-1*eV, 4000*eV, 10*eV, 1e5);
	//exp.setEnvironment(3e-2*eV, 10*eV, 0.0010*eV, 1e5);
	exp.setLaser(1e-6, 20e-6, 250e-15, 1e5, 20);
	exp.initMaterial(Cu,"DOSCAR");
	std::vector<double> rates;
	std::vector<double> escapeds;
	std::vector<double> Ts;
	for (double T_e = 1000; T_e < 1000000; T_e+=1000)
	{
		Ts.push_back(T_e);
		escapeds.push_back(exp.getNElectronEscaped(T_e));
		rates.push_back(exp.getElectronEscapeRate(T_e));
	}
	plot.plotSingle("Electron emission rate / Electron$\\cdot$s$^{-1}$",Ts, rates,"T/K","emissionRate.pdf");
	plot.plotSingle("Electron escaped",Ts,escapeds,"T/K","escaped.pdf");

	return 0;
}