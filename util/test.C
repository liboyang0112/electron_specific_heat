#include "constants.h"
#include <iostream>
#include "Experiment.h"
#include "formulas.h"
#include "Plot.h"
#include <fstream>
#include <sstream>
#include "CalcMaterial.h"
using namespace std;
int main(int argc, char const *argv[])
{
	Plot plot;
	plot.init(".");
/*	
	Experiment exp;
	exp.workPiece = Cu;

	exp.workPiece.electronTemprature=1200;
	exp.workPiece.latticeTemprature=600;
	cout<<getFermiV(Cu.rho,Cu.n)<<endl;
	//cout<<getFermiEnergy_eV(Ag.rho,Ag.n)<<endl;
	//cout<<getFermiEnergy_eV(Au.rho,Au.n)<<endl;
	double R, Alpha, n, k, ep1, ep2;
	exp.laserWaveLength = 1e-6;
	std::vector<double> ns, ks, Rs, Alphas, wls, ep1s, ep2s;
	//for (double i = 0.02e-6; i < 2e-6; i+=0.01e-6)
	//{
	//	exp.laserWaveLength = i;
		exp.getRandAlpha(R, Alpha, n, k, ep1, ep2);
	//	wls.push_back(i);
	//	ns.push_back(n);
	//	ks.push_back(k);
	//	Alphas.push_back(Alpha);
	//	Rs.push_back(R);
	//	ep1s.push_back(ep1);
	//	ep2s.push_back(ep2);
	//}
	//plot.plotDouble("n","k",wls,ns,ks,"lambda","nk.pdf");
	//plot.plotDouble("ep1","ep2",wls,ep1s,ep2s,"lambda","epsilon.pdf");
	//plot.plotSingle("R",wls,Rs,"lambda","R.pdf");
	//plot.plotSingle("Alpha",wls,Alphas,"lambda","Alpha.pdf");
//
	//getRandAlpha(-7.3421, 3.5438, lambda, R, Alpha);
	cout<<"R: "<<R<<", alpha: "<<Alpha<<"test: "<<8.02946e+13*2*PI*h/q_e/2/PI/10.83<<endl;
*/
/*
	CalcMaterial cv;
	cv.readDOS("DOSCAR");
	string filename = "Cu.yml";
	std::ifstream file(filename);
	if (!file)
	{
			std::cerr << "File " << filename << " could not be opened" << std::endl;
			exit(0);
	}
	std::string inputline;
	while(getline(file,inputline))
	{
		if(inputline.find("data:")!=std::string::npos) break;
	}
	double wl,n,k;
	std::vector<double> ns, ks, Rs, Alphas, wls, ep1s, ep2s,incs;
	while(getline(file,inputline))
	{
		std::stringstream stream(inputline);
		stream>>wl>>n>>k;
		Alphas.push_back(4*PI*k/wl);
		wls.push_back(-h*c/wl/10/q_e*1e6);
		incs.push_back(h*c/wl/10/q_e*1e6/cv.DOS_incidence(1,1e-20,h*c/wl/10/q_e*1e6));
		Rs.push_back((pow(n-1,2)+k*k) / (pow(n+1,2)+k*k));
		ns.push_back(n);
		ks.push_back(k);
	}
	plot.plotSingle("Alpha",wls,Alphas,"reached energy","Alpha.pdf");
	plot.plotSingle("R",wls,Rs,"reached energy","R.pdf");
	plot.plotSingle("incs",wls,incs,"reached energy","incs.pdf");

	std::vector<double> vZvsV;
	std::vector<double> vUbarvsV;
	std::vector<double> vZvsT[190];
	std::vector<double> vUbarvsT[190];
	std::vector<double> vT;
	for (double vol = 1; vol < 20; vol+=1)
	{
		if(vol==1) readLine("V_"+std::to_string(vol),"T",vT);
		readLine("V_"+std::to_string(vol),"Z",vZvsT[int(vol-1)]);
		readLine("V_"+std::to_string(vol),"Ubar",vUbarvsT[int(vol-1)]);
	}

	plot.plotTripple("Z, rV=1","Z, rV=2","Z, rV=3",vT,vZvsT[0],vZvsT[1],vZvsT[2],"T/eV","ZdiffV.pdf");
	plot.plotTripple("Ubar, rV=1","Ubar, rV=2","Ubar, rV=3",vT,vUbarvsT[0],vUbarvsT[1],vUbarvsT[2],"T/eV","UbardiffV.pdf");
	*/

double k_B = 3.285e4;
double rho0_water=1./0.5824;
double rho0_f = 1.3897e-3;
double phi_fluid_initial = 0.6;
double alpha_water=7.438e4;
double ratio = 1./(1-rho0_f/rho0_water);
double ret = rho0_f*k_B*phi_fluid_initial*ratio-alpha_water*rho0_f*rho0_f;
printf("%e\n", ret);
return 0;
}