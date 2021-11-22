#include "CalcMaterial.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <sstream>
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
#include <vector>
#include <map>
#include "gsl/gsl_errno.h"
#include "constants.h"
#include "formulas.h"
#include <complex.h>
#include "Faddeeva.h"
struct multiGaussianFit_params { std::vector<double> &x; std::vector<double> &y ;};

double gaussianFunc(const gsl_vector *p, void *param){
	struct multiGaussianFit_params *params = (struct multiGaussianFit_params*)param;
	auto &x = params->x;
	auto &y = params->y;
	double func = 0;
	double error2 = 0;
	for (u_int j = 0; j < x.size(); ++j)
	{
		func = 0;
		for (u_int i = 0; i < p->size/2; ++i)
		{
			func+=gsl_vector_get(p,i*2)*gsl_ran_gaussian_pdf(x.at(j),gsl_vector_get(p,i*2+1));
		}
		error2+=pow((func-y.at(j)),2);
	}
	return error2;
}

double erfFunc(const gsl_vector *p, void *param){
	struct multiGaussianFit_params *params = (struct multiGaussianFit_params*)param;
	auto &x = params->x;
	auto &y = params->y;
	double func = 0;
	double error2 = 0;
	double norm = 0;
	for (u_int i = 0; i < p->size/2; ++i)
	{
		norm+=gsl_vector_get(p,i*2);
	}
	norm = 1.5-norm;
	for (u_int j = 0; j < x.size(); ++j)
	{
		func = 0;
		for (u_int i = 0; i < p->size/2; ++i)
		{
			func+=gsl_vector_get(p,i*2)*erf(x.at(j)/gsl_vector_get(p,i*2+1));
		}
		func+=norm*erf(x.at(j)/gsl_vector_get(p,p->size-1));
		error2+=pow((func-y.at(j)),2);
	}
	return error2;
}

double gaussianPlusCosExpFunc(const gsl_vector *p, void *param){
	struct multiGaussianFit_params *params = (struct multiGaussianFit_params*)param;
	auto &x = params->x;
	auto &y = params->y;
	double func = 0;
	double error2 = 0;
	for (u_int j = 0; j < x.size(); ++j)
	{
		func = 0;
		for (u_int i = 0; i < p->size/2-1; ++i)
		{
			func+=gsl_vector_get(p,i*2)*gsl_ran_gaussian_pdf(x.at(j),gsl_vector_get(p,i*2+1));
		}
		double a = gsl_vector_get(p,p->size-3);
		double b = gsl_vector_get(p,p->size-2);
		double c = gsl_vector_get(p,p->size-1);
		func+=c*( a*sin(a*x.at(j)) + b*cos(a*x.at(j)) ) * exp(-b*x.at(j));
		error2+=pow((func-y.at(j)),2);
	}
	return error2;
}

void multiGaussianFit(std::vector<double> &v, std::vector<double> &x, std::vector<double> &y){
	u_int size = v.size();
	gsl_vector *gv = gsl_vector_alloc(size);
	for (u_int i = 0; i < size; ++i)
	{
		gsl_vector_set (gv, i, v.at(i));
	}
	auto ss = gsl_vector_alloc (size);
	gsl_vector_set_all (ss, 0.2);
	gsl_multimin_fminimizer *mn = gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex2,size);
	gsl_multimin_function minex_func;
	struct multiGaussianFit_params params = { x, y };
	minex_func.n = size;
	minex_func.f = erfFunc;
	minex_func.params = &params;
	gsl_multimin_fminimizer_set (mn, &minex_func, gv, ss);
	int status;
	double mnsize;
	int iter = 0;
	do
	{
		iter++;
		status = gsl_multimin_fminimizer_iterate(mn);
		if (status)
		break;
		mnsize = gsl_multimin_fminimizer_size (mn);
		status = gsl_multimin_test_size (mnsize, 1e-2);
		
		if (status == GSL_SUCCESS)
		{
		printf ("converged to minimum at\n");
		}
		printf ("%10.3e %10.3e f() = %7.3f size = %.3f \n",
		gsl_vector_get (mn->x, 0),
		gsl_vector_get (mn->x, 1),
		mn->fval, mnsize);
	} while (status == GSL_CONTINUE);
	for (u_int i = 0; i < size; ++i)
	{
		v[i] = gsl_vector_get(mn->x,i);
	}
	gsl_vector_free(gv);
	gsl_vector_free(ss);
	gsl_multimin_fminimizer_free (mn);
}


std::vector<double> smooth(std::vector<double> v){
	u_int size = v.size();
	gsl_vector *gv = gsl_vector_alloc(size);
	gsl_vector *gvs = gsl_vector_alloc(size);
	for (u_int i = 0; i < size; ++i)
	{
		gsl_vector_set(gv, i, v[i]);
	}

	gsl_filter_gaussian_workspace *fws = gsl_filter_gaussian_alloc(10);
	gsl_filter_gaussian(GSL_FILTER_END_PADVALUE, 0.1, 0, gv, gvs, fws);
	gsl_filter_gaussian_free(fws);
	//gsl_filter_median_workspace *fws = gsl_filter_median_alloc(10);
	//gsl_filter_median(GSL_FILTER_END_PADVALUE, gv, gvs, fws);
	//gsl_filter_median_free(fws);
	for (u_int i = 0; i < size; ++i)
	{
		v[i] = gsl_vector_get(gvs, i);
	}
	return v;
}



struct dirac_func_params { double x; double j; };

double dirac_func(double k,	void * p){
	struct dirac_func_params * params = (struct dirac_func_params *)p;
	return pow(k,params->j)/(exp(k-params->x)+1) + pow(k,-params->j-2)/(exp(1./k-params->x)+1);
}

double dirac_int(double x, double j){
	if(x>0) return (j==1./2?gsl_sf_fermi_dirac_half(x)*gsl_sf_gamma(1.5):gsl_sf_fermi_dirac_3half(x)*gsl_sf_gamma(2.5));
	gsl_integration_workspace *iws = gsl_integration_workspace_alloc(500);
	double result, error;

	gsl_function F;
	struct dirac_func_params params = { x, j };

	F.function = &dirac_func;
	F.params = &params;

	gsl_integration_qag(&F,1e-6,1,0,1e-7,100,3,iws,&result,&error);

	gsl_integration_workspace_free(iws);
	return result;
}

void CalcMaterial::initMaterial(material m, std::string filename){
	mat = m;
	readLine("/home/boyang/softwares/electron_specific_heat/config/ionize.txt",mat.name,ionize);

	std::ifstream fileenergyLevel("/home/boyang/softwares/electron_specific_heat/config/energyLevel.txt");
	if (!fileenergyLevel)
	{
			std::cerr << "File " << "/home/boyang/softwares/electron_specific_heat/config/energyLevel.txt" << " could not be opened" << std::endl;
			exit(0);
	}
	std::string inputline;
	std::string matname;
	double altmp,eltmp;
	double neleremain = mat.ele_number;
	if(bareAtomOrbit){
		for (int i = 1; i < bareAtomOrbit+1; ++i)
		{
			if(neleremain<2*i*i) {
				a_l.push_back(neleremain);
				epsilon_l.push_back(-13.6*mat.ele_number*mat.ele_number/i/i);
				a_l.push_back(2*i*i-neleremain);
				epsilon_l.push_back(-13.6*mat.ele_number*mat.ele_number/i/i+1);
			}

			else {
				neleremain-=2*i*i;
				a_l.push_back(2*i*i);
				epsilon_l.push_back(-13.6*mat.ele_number*mat.ele_number/i/i);
			}
		}
	}else if(useIon){
		for (uint i = 0; i < ionize.size(); ++i)
		{
			a_l.push_back(1);
			epsilon_l.push_back(-ionize[i]*1e6/q_e/N_A);
		}
	}else while(getline(fileenergyLevel,inputline)){
		if(inputline.find(mat.name)!=std::string::npos){
			std::stringstream ss(inputline);
			ss>>matname;
			while(ss){
				ss>>altmp>>eltmp; //kJ/mol
				a_l.push_back(altmp);
				epsilon_l.push_back(eltmp);
			}
			break;
		}
	}

	readDOS(filename);
}

double CalcMaterial::computeG(double x, double y, double lambdaOmega2_meV){
	double result;
	for (double i = 0; i < energyEnd; i+=interval)
	{
		result+=pow(DOS.at(int(i/interval)),2)/(exp((i-x)/y)+2+exp((x-i)/y));
	}
	return 2*PI*PI*lambdaOmega2_meV*1e-6*q_e*q_e/h/DOS.at(int(1./interval))/y*result*interval;
}

double CalcMaterial::DOS_incidence(double x, double y, double power, double start, double end){
	if(end == 0) end = energyEnd;

	double result;
	for (double i = start*interval; i+power < end; i+=interval)
	{
		result+=DOS.at(int(i/interval))/(exp((i-x)/y)+1)*DOS.at(int((i+power)/interval))*(1-(1./(exp((i+power-x)/y)+1)));
	}
	return result*interval;
}

double freeDOSInt(double x, double y, double power){
	return 1.5*pow(y,1.5+power)*dirac_int(x/y,0.5+power);
}

double CalcMaterial::DOSInt(double x, double y, double power, double start, double end){
	if(end == 0) end = energyEnd;

	double result = 0;
	for (double i = start*interval; i < end; i+=interval)
	{
		result+=DOS.at(int(i/interval))*pow(i,power)/(exp((i-x)/y)+1);
	}
	return result*interval;
}

double CalcMaterial::latticePlasmaDOSInt(double x, double y, double power){
	double result = 0;
	for (double l = 0; l < 1 + mat.workFunc/fermiEnergy; l+=interval)
	{
		result+=DOS.at(int(l/interval)) * pow(l,power) / (exp((l-1-mat.workFunc/fermiEnergy-x)/y)+1);
	}
	result = result*interval+1.5*pow(y,1.5+power)*dirac_int(x/y,0.5+power);
	int nEle_inner = mat.ele_number-nEle;
	for (uint i = 0; i < a_l.size(); ++i)
	{
		if(nEle_inner<=0) break;
		nEle_inner-=a_l[i];
		result += a_l[i]*pow(epsilon_l[i]/fermiEnergy,power) / (exp((epsilon_l[i]/fermiEnergy-x)/y)+1);
	}
	return result;
}

double CalcMaterial::latticePlasmaDOSIntInner(double x, double y, double power){
	double result = 0;
	for (double l = 0; l < 1 + mat.workFunc/fermiEnergy; l+=interval)
	{
		result+=DOS.at(int(l/interval)) * pow(l,power) / (exp((l-1-mat.workFunc/fermiEnergy-x)/y)+1);
	}
	result *= interval;
	int nEle_inner = mat.ele_number-nEle;
	for (uint i = 0; i < a_l.size(); ++i)
	{
		if(nEle_inner<=0) break;
		nEle_inner-=a_l[i];
		result += a_l[i]*pow(epsilon_l[i]/fermiEnergy,power) / (exp((epsilon_l[i]/fermiEnergy-x)/y)+1);
	}
	return result;
}

double CalcMaterial::DLatticePlasmaDOSIntInner(double x, double y, double xp, double yp, double power){
	double result = 0;
	double e1,e2,de;
	for (double l = 0; l < 1 + mat.workFunc/fermiEnergy; l+=interval)
	{
		e1=exp((l-1-mat.workFunc/fermiEnergy-x)/y);
		e2=exp((l-1-mat.workFunc/fermiEnergy-xp)/yp);
		de = (exp((l-1-mat.workFunc/fermiEnergy-xp)/yp-(l-1-mat.workFunc/fermiEnergy-x)/y)-1)*e1;
		result+=DOS.at(int(l/interval)) * pow(l,power) *(de /(e1+1)/(e2+1));
	}
	result *= interval;
	int nEle_inner = mat.ele_number-nEle;
	for (uint i = 0; i < a_l.size(); ++i)
	{
		if(nEle_inner<=0) break;
		nEle_inner-=a_l[i];
		e1=exp((epsilon_l[i]/fermiEnergy-x)/y);
		e2=exp((epsilon_l[i]/fermiEnergy-xp)/yp);
		de = (exp((epsilon_l[i]/fermiEnergy-xp)/yp-(epsilon_l[i]/fermiEnergy-x)/y)-1)*e1;
		result += a_l[i]*pow(epsilon_l[i]/fermiEnergy,power)*(de /(e1+1)/(e2+1));
	}
	return result;
}

double CalcMaterial::gasPlasmaDOSInt(double x, double y, double power, double Vratio){
	double result = 1.5*pow(y,3./2+power)*dirac_int(x/y,0.5+power)*Vratio;
	for (uint i = 0; i < a_l.size(); ++i)
	{
		result += a_l[i]*pow(epsilon_l[i]/fermiEnergy,power) / (exp((epsilon_l[i]/fermiEnergy-x)/y)+1);
	}
	return result;
}

double CalcMaterial::gasPlasmaDOSIntInner(double x, double y, double power){
	double result = 0;
	for (uint i = 0; i < a_l.size(); ++i)
	{
		result += a_l[i]*pow(epsilon_l[i]/fermiEnergy,power) / (exp((epsilon_l[i]/fermiEnergy-x)/y)+1);
	}
	return result;
}

double CalcMaterial::DGasPlasmaDOSIntInner(double x, double y, double xp, double yp, double power){
	double result = 0;
	double e1,e2,de;
	for (uint i = 0; i < a_l.size(); ++i)
	{
		e1=exp((epsilon_l[i]/fermiEnergy-x)/y);
		e2=exp((epsilon_l[i]/fermiEnergy-xp)/yp);
		de = (exp((epsilon_l[i]/fermiEnergy-xp)/yp-(epsilon_l[i]/fermiEnergy-x)/y)-1)*e1;
		de = a_l[i]*( de/(e1+1)/(e2 +1));
		result +=pow(epsilon_l[i]/fermiEnergy,power)*de;
	}
	return result;
}

double CalcMaterial::DOSIntI(double x, double y, double power, double Vratio){
	if(Vratio==0) Vratio = VratioDefault;
	switch(DOSmode){
		case lattice_plasma:
		return latticePlasmaDOSInt(x, y, power);
		case gas_plasma:
		return gasPlasmaDOSInt(x, y, power, Vratio);
		case free:
		return freeDOSInt(x,y,power);
		default:
		return DOSInt(x, y, power);
	}
}

double CalcMaterial::DOSIntInnerI(double x, double y, double power){
	switch(DOSmode){
		case lattice_plasma:
		return latticePlasmaDOSIntInner(x, y, power);
		case gas_plasma:
		return gasPlasmaDOSIntInner(x, y, power);
		default:
		return 0;
	}
}

double CalcMaterial::DDOSIntInnerI(double x, double y, double xp, double yp, double power){
	switch(DOSmode){
		case lattice_plasma:
		return DLatticePlasmaDOSIntInner(x, y, xp, yp, power);
		case gas_plasma:
		return DGasPlasmaDOSIntInner(x, y, xp, yp, power);
		default:
		return 0;
	}
}

double CalcMaterial::latticePlasmaZ(double x, double y){
	double result = 0;
	//for (double l = 0; l < 1 + mat.workFunc/fermiEnergy; l+=interval)
	//{
	//	result+=DOS.at(int(l/interval)) / (exp((l-1-mat.workFunc/fermiEnergy-x)/y)+1);
	//}
	result = result*interval+1.5*pow(y,1.5)*dirac_int(x/y,0.5);
	return result;
}


double CalcMaterial::gasPlasmaZ(double x, double y, double Vratio){
	if(Vratio==0) Vratio = VratioDefault;
	return 1.5*pow(y,3./2)*dirac_int(x/y,0.5)*Vratio;
}

void CalcMaterial::readDOS(std::string filename){
	printf("reading DOS file: %s\n", filename.c_str());
	DOS.clear();
	std::ifstream file(filename);
	if (!file)
	{
			std::cerr << "File " << filename << " could not be opened" << std::endl;
			exit(0);
	}
	std::string inputline;
	for (int i = 0; i < 6; ++i)
	{
		getline(file,inputline);
	}
	std::stringstream stream(inputline);
	double energyStart, norm, energy, DOS_up, DOS_down, intDOS_up, intDOS_down, totPoints;
	int spinFormat;
	stream>>energyEnd>>energyStart>>totPoints>>fermiEnergy>>norm>>spinFormat;
	bool started = 0;
	interval = (energyEnd-energyStart)/(totPoints-1);
	std::vector<double> plotx;
	plotx.push_back(0);
	DOS_up=DOS_down=0;
	intDOS_up=intDOS_down=0;
	DOS.push_back(DOS_up+DOS_down);
	intDOS.push_back(intDOS_up+intDOS_down);
	for (int i = 0; i < totPoints; ++i){
		getline(file,inputline);
		std::stringstream stream(inputline);
		if(spinFormat) stream>>energy>>DOS_up>>DOS_down>>intDOS_up>>intDOS_down;
		else{
			stream>>energy>>DOS_up>>intDOS_up;
		}if(DOS_up+DOS_down+intDOS_up+intDOS_down==0) {
			if(!started){
				energyStart=energy;
				continue;
			}else{
				energyEnd = energy;
				DOS.push_back(DOS_up+DOS_down);
				intDOS.push_back(intDOS_up+intDOS_down);
				break;
			}
		}
		started = 1;
		DOS.push_back(DOS_up+DOS_down);
		intDOS.push_back(intDOS_up+intDOS_down);
	}
	fermiEnergy-=energyStart;
	energyEnd=(energyEnd-energyStart)/fermiEnergy;
	interval=energyEnd/DOS.size();
	for (u_int i = 0; i < DOS.size()-1; ++i)
	{
		plotx.push_back(plotx.back()+interval*fermiEnergy);
	}
	norm = DOSInt(1, 1e-20, 0, 0, 1);

	for(auto &x : DOS){
		x*=nEle/norm;
	}
	for(auto &x : intDOS){
		x*=nEle/norm;
	}

	std::cout<<"integrated: "<<DOSInt(1, 1e-20, 0, 0, 1)<<std::endl;
	std::cout<<"Fermi Energy: "<<fermiEnergy<<std::endl;
	plot.plotSingle("DOS, ${\\epsilon_F}$="+std::to_string(fermiEnergy)+"eV", plotx, DOS, "${\\epsilon/eV}$", "DOS.pdf");

	return;
}


void CalcMaterial::freeDOS(){
	nEle = 1;
	DOSmode=free;
	dimentionless = 1;
	TUnit="y";
	fermiEnergy = getFermiEnergy_eV(mat);
	printf("create free electron DOS\n");
	DOS.clear();
	for (double i = 0; i < 5; i+=1e-3){
		DOS.push_back(sqrt(i));
	}
	energyEnd=5;
	interval=energyEnd/DOS.size();
	std::vector<double> plotx;
	plotx.push_back(0);
	for (u_int i = 0; i < DOS.size()-1; ++i)
	{
		plotx.push_back(plotx.back()+interval);
	}
	double norm = DOSInt(1, 1e-20, 0, 0, 1);

	for(auto &x : DOS){
		x*=nEle/norm;
	}
	for(auto &x : intDOS){
		x*=nEle/norm;
	}

	std::cout<<"integrated: "<<DOSInt(1, 1e-20, 0, 0, 1)<<std::endl;
	plot.plotSingle("${DOS}$", plotx, DOS, "${l}$", "DOS_free.pdf");

	return;
}

struct my_f_params { double y; CalcMaterial *calc; };

double my_f(double x, void * p)
{
	
	struct my_f_params * params = (struct my_f_params *)p;
	//double val = pow(y,3./2)*dirac_int(x/y,0.5)-2./3;
	//{ // Wrong results without these rubbish, WTF?
	//	std::stringstream ss;
	//	ss<<1;
	//}
	double val = params->calc->DOSIntI(x,params->y) - (params->calc->DOSmode==params->calc->free?1:params->calc->mat.ele_number);
	//double val = params->calc->DOSInt(x,params->y) - params->calc->nEle;
	return val;
}

void CalcMaterial::plotDOSNandF(double x, double y, std::string savename){\

	std::cout<<"plot N and f: "<<savename<<std::endl;
	std::vector<double> epsilon;
	std::vector<double> Dl;
	std::vector<double> f;
	std::vector<double> N;
	double Nint = 0;
	double Uint = 0;
	int nEle_inner = mat.ele_number-nEle;
	if(DOSmode>1){
		for (uint i = 0; i < a_l.size(); ++i)
		{
			if(nEle_inner<=0 && DOSmode==lattice_plasma) break;
			nEle_inner-=a_l[i];
			epsilon.push_back(epsilon_l[i]/fermiEnergy-interval);
			f.push_back(1./(exp((epsilon.back()-x)/y)+1)*100);
			Dl.push_back(0);
			N.push_back(f.back()*Dl.back()/100);
			epsilon.push_back(epsilon_l[i]/fermiEnergy);
			f.push_back(1./(exp((epsilon.back()-x)/y)+1)*100);
			Dl.push_back(a_l[i]/interval);
			N.push_back(f.back()*Dl.back()/100);
			Nint+=N.back()*interval;
			Uint+=i*N.back()*interval;
			epsilon.push_back(epsilon_l[i]/fermiEnergy+interval);
			f.push_back(1./(exp((epsilon.back()-x)/y)+1)*100);
			Dl.push_back(0);
			N.push_back(f.back()*Dl.back()/100);
		}
		if(doPlots) plot.plotTripple("DOS", "electron density", "Occupation (%)", epsilon, Dl, N, f, "l", "inner_"+savename);
	}
	epsilon.clear();
	Dl.clear();
	f.clear();
	N.clear();
	if(DOSmode!=gas_plasma)
	for (double i = (DOSmode==lattice_plasma?-1-mat.workFunc/fermiEnergy:0); i < (DOSmode==lattice_plasma?0:energyEnd); i+=interval)
	{
		epsilon.push_back(i);
		f.push_back(1./(exp((i-x)/y)+1)*(DOSmode==free?1:100));
		Dl.push_back(DOS.at(int((i+(DOSmode==lattice_plasma?1+mat.workFunc/fermiEnergy:0))/interval)));
		N.push_back(f.back()*Dl.back()/(DOSmode==free?1:100));
		Nint+=N.back()*interval;
		Uint+=i*N.back()*interval;
	}

	if(DOSmode>1){
		for (double i = 0; i < energyEnd; i+=interval)
		{
			epsilon.push_back(i);
			f.push_back(1./(exp((i-x)/y)+1)*100);
			Dl.push_back(3./2*sqrt(i));
			N.push_back(f.back()*Dl.back()/100);
			Nint+=N.back()*interval;
			Uint+=i*N.back()*interval;
		}
	}

	std::cout<<"N integral = "<<Nint<<std::endl;
	//std::cout<<"DOS integral = "<<DOSInt(x,y)<<std::endl;
	std::cout<<"DOS integral = "<<DOSIntI(x,y)<<std::endl;
	std::cout<<"U integral = "<<Uint<<std::endl;
	//std::cout<<"UDOS integral = "<<DOSInt(x,y,1)<<std::endl;
	std::cout<<"UDOS integral = "<<DOSIntI(x,y,1)<<std::endl;

	plot.plotTripple("DOS", "electron density", DOSmode==free?"Occupation":"Occupation (%)", epsilon, Dl, N, f, "l", savename);
}

struct rad_int_params { double T_e; double h; CalcMaterial *calc; };

double rad_int_Core(double lambda, void *p){
	double omega = 2*PI*c/lambda;
	struct rad_int_params * params = (struct rad_int_params *)p;
	double alpha = params->calc->getAlphaInt(params->T_e, lambda);
	double R = params->calc->getRInt(params->T_e, lambda);
	//std::cout<<"alpha="<<alpha<<", R="<<R<<", T="<<params->T_e<<", h="<<params->h<<", lambda="<<lambda<<std::endl;
	double result = 2*alpha/PI*exp(-2*alpha/PI*params->h)*(1-R)*2*PI*c/lambda/lambda*
	hbar*pow(omega,3)/2/PI/PI/c/c/(exp(hbar*omega/k_B/params->T_e)-1);
	//std::cout<<"int="<<result<<std::endl;
	return result;
}

int CalcMaterial::calculate(double sT, double eT, double iT)
{
	conductivity_a = mat.rho*N_A/mat.n*q_e*q_e/mat.conductivity/250/k_B;
	conductivity_b = mat.rho*N_A/mat.n*q_e*q_e*mat.rhoee*2/k_B;
	conductivity_c = fermiEnergy*q_e*2/3/k_B;
	std::cout<<"a="<<conductivity_a<<",b="<<conductivity_b<<",c="<<conductivity_c<<std::endl;
	std::vector<double> C_e_lowT;
	std::vector<double> C_e_highT;
	std::vector<double> dC_e;
	std::vector<double> K_es;
	std::vector<double> Ubar;
	startT = sT;
	endT = eT;
	intervalT = iT;
	logstepT = iT;
	C_e.clear();
	C_e_lowT.clear();
	C_e_highT.clear();
	dC_e.clear();
	U_e.clear();
	Gs.clear();
	ys.clear();
	xs.clear();
	Ts.clear();
	mus.clear();
	Z_l.clear();
	TsineV.clear();
	TsinEpsilon.clear();
	Uin.clear();
	C_l.clear();
	double Nk=mat.rho*N_A/mat.n*k_B;
	const gsl_root_fsolver_type * T = gsl_root_fsolver_brent;
	gsl_root_fsolver * s = gsl_root_fsolver_alloc (T);
	gsl_function F;
	double y = 1e-4;
	struct my_f_params params = { y, this };

	F.function = &my_f;
	F.params = &params;

	std::cout<<"endT="<<endT<<std::endl<<"N="<<mat.rho*N_A/mat.n<<std::endl<<"Nk="<<Nk<<std::endl;

	double lastU=0;
	double solution;
	bool plotyet = 0;
	int i = 0;

	double prevX=0;
	double prevY=0;
	//std::vector<double> dUin;
	std::vector<double> dZ_l;
	//dUin.push_back(0);
	dZ_l.push_back(0);
	double Uin0=0;
	for (double T = startT; T < endT;)
	{
		y=T/fermiEnergy/q_e*k_B;
		params.y=y;
		gsl_root_fsolver_set (s, &F, -100000,100000);
		double root;
		solution = 0;
		while(1){
			gsl_root_fsolver_iterate(s);
			root=gsl_root_fsolver_root(s);
			if(fabs(solution-root)<1e-300) break;
			solution=root;
		}
		double U_e_tmp=DOSIntI(solution,y,1)*(dimentionless?1:Nk*fermiEnergy*q_e/k_B);
		if(DOSmode == lattice_plasma) {
			Z_l.push_back(latticePlasmaZ(solution, y));
			if(i>0) {
				Uin.push_back(DOSIntInnerI(solution,y,1)*(dimentionless?1:Nk*fermiEnergy*q_e/k_B)-Uin0);
//				dUin.push_back(DDOSIntInnerI(solution,y,prevX,prevY,1)*(dimentionless?1:Nk*fermiEnergy*q_e/k_B));
//				dZ_l.push_back(DDOSIntInnerI(solution,y,prevX,prevY,0));
			}else{
				Uin0=DOSIntInnerI(solution,y,1)*(dimentionless?1:Nk*fermiEnergy*q_e/k_B);
				Uin.push_back(0);
			}
		}
		else if(DOSmode == gas_plasma) {
			Z_l.push_back(gasPlasmaZ(solution, y));
			if(i>0) {
				Uin.push_back(DOSIntInnerI(solution,y,1)*(dimentionless?1:Nk*fermiEnergy*q_e/k_B)-Uin0);
//				dUin.push_back(DDOSIntInnerI(solution,y,prevX,prevY,1)*(dimentionless?1:Nk*fermiEnergy*q_e/k_B));
//				dZ_l.push_back(DDOSIntInnerI(solution,y,prevX,prevY,0));
			}else{
				Uin0=DOSIntInnerI(solution,y,1)*(dimentionless?1:Nk*fermiEnergy*q_e/k_B);
				Uin.push_back(0);
			}
			
		}
 		U_e.push_back(U_e_tmp);
		if(i>0) C_e.push_back((U_e_tmp-lastU)/intervalT*(dimentionless?fermiEnergy*q_e/k_B:1));
		else C_e.push_back(0);
		C_e_highT.push_back(3./2*(dimentionless?1:Nk));
		//if(i>0) C_e_lowT.push_back(std::min(C_e_lowT.back()+PI*PI/2*intervalT*Nk/fermiEnergy/q_e*k_B,C_e_highT.back()));
		C_e_lowT.push_back(T*PI*PI/2/fermiEnergy/q_e*k_B*(dimentionless?1:Nk));
		if(i>0) dC_e.push_back(i>1?(U_e[i]+U_e[i-2]-2*U_e[i-1])/intervalT/intervalT:0);
		xs.push_back(solution);
		ys.push_back(y);
		Ts.push_back(T);
		TsineV.push_back(T/11600);
		TsinEpsilon.push_back(T/11600/fermiEnergy);
		mus.push_back(solution*fermiEnergy);
		C_l.push_back(3*Nk*debye_int(mat.debyeTemprature/T));
		lastU=U_e_tmp;
		Gs.push_back(computeG(solution,y,mat.lambdaOmega2_meV)*N_A*mat.rho/mat.n/fermiEnergy*q_e/k_B);
		if(doPlots)
			if((y>0.5||T >= endT-intervalT)&&!plotyet){
				plotDOSNandF(solution, y, "NandFAt0.5.pdf");
				plotyet=1;
			}
		i++;
		prevX=solution;
		prevY=y;
		if(logstepT) {
			intervalT = (logstepT-1)*T;
			T+=intervalT;
		}
		else T+=intervalT;
	}
	//Ubar.push_back(0);
	//for (uint i = 0; i < Z_l.size()-1; ++i)
	//{
	//	Ubar.push_back(-dUin[i]/dZ_l[i]);
	//}
	//Ubar[0]=Ubar[1];
	C_e[0]=2*C_e[1]-C_e[2];
	if(C_e[0]<0) C_e[0]=0;
	dC_e[0]=dC_e[1]=2*dC_e[2]-dC_e[3];
	dC_e.push_back(dC_e.back());
	auto smootheddC = smooth(smooth(smooth(dC_e)));
	for (double T = startT; T < endT; T*=logstepT)
	{
		K_es.push_back(getK_e(300,T));
	}

	if(doPlots){
	//plot.plotSingle("${C_v^e/Nk}$", ys, smooth(C_e), "${y}$", "C_e.pdf");
//	plot.plotSingle("$G$", ys, Gs, "${y}$", "G.pdf");
//	plot.plotDouble("${C_v^e/Nk}$", "${C_v^e/Nk,low T approx}$", ys, smooth(C_e), C_e_lowT, "${y}$", "C_e.pdf");
//	plot.plotSingle("${dC_v^e/dy/Nk}$", ys, smootheddC, "${y}$", "dC_V.pdf");
//	plot.plotSingle("${U/T_{F}}Nk$", ys, U_e, "${y}$", "U.pdf");

		if(DOSmode >= 2) plot.plotSingle("Degree of ionization", TsInUnit(), Z_l, TUnit, "Z_l.pdf");
//		if(DOSmode >= 2) plot.plotSingle("d Degree of ionization", TsInUnit(), dZ_l, TUnit, "dZ_l.pdf");
		if(DOSmode >= 2) 
			plot.plotSingle("${C_e}$/J${\\cdot}$K${^{-1}\\cdot}$m${^{-3}}$", TsInUnit(), smooth(C_e), TUnit, "C_e.pdf");

		if(dimentionless) plot.plotSingle("x", ys, xs, "${y}$", "x.pdf");
		else plot.plotSingle("${\\mu}$/eV", TsInUnit(), mus, TUnit, "mu.pdf");
		std::string cunit = dimentionless?"${k_{B}}$":"J${\\cdot}$K${^{-1}\\cdot}$m${^{-3}}$";
		plot.plotTripple("${C_e}$/"+cunit, "${C_{e,FEA}}$/"+cunit, "${C_{e,PG}}$/"+cunit, TsInUnit(), smooth(C_e), C_e_lowT, C_e_highT, TUnit, "C_e_compare.pdf");
		plot.plotSingle("$G$/J${\\cdot}$s$^{-1}{\\cdot}$m${^{-3}\\cdot}$K${^{-1}}$", TsInUnit(), Gs, TUnit, "G.pdf");
		plot.plotSingle("${dC_e/dT}$/J${\\cdot}$K${^{-2}\\cdot}$m${^{-3}}$", TsInUnit(), smootheddC, TUnit, "dC_e.pdf");
		std::string uunit = dimentionless?"${\\epsilon_F}$":"J${\\cdot}$m${^{-3}}$";
		plot.plotSingle("${U}$/"+uunit, TsInUnit(), U_e, TUnit, "U.pdf");
//		plot.plotSingleLogXY("${\\bar{U}}$/"+uunit, TsInUnit(), Ubar, TUnit, "Ubar.pdf");
		plot.plotSingle("${\\bar{U}}$/"+uunit, TsInUnit(), Uin, TUnit, "Uin.pdf");
		plot.plotSingle("${C_l}$/J${\\cdot}$K${^{-1}\\cdot}$m${^{-3}}$",TsInUnit(),C_l,TUnit,"C_l.pdf");
		plot.plotDouble("${C_l}$/J${\\cdot}$K${^{-1}\\cdot}$m${^{-3}}$","${C_e}$/J${\\cdot}$K${^{-1}\\cdot}$m${^{-3}}$", TsInUnit(), C_l, C_e,TUnit,"C_lC_e.pdf");
		plot.plotSingle("${K_e}$/J${\\cdot}$K${^{-1}\\cdot}$m${^{-1}}$s${^{-1}}$",TsInUnit(),K_es,TUnit,"K_e.pdf");
	}
	return 0;
}
int CalcMaterial::calculateOptical(){
	gsl_function F;
	std::string inputline;
	std::string matname;
	std::ifstream file("/home/boyang/softwares/electron_specific_heat/config/BBparam.txt");
	if (!file)
	{
			std::cerr << "File " << "/home/boyang/softwares/electron_specific_heat/config/BBparam.txt" << " could not be opened" << std::endl;
			exit(0);
	}
	double f0, gamma0, f,gamma,omega,sigma;
	std::vector<double> fs, gammas, omegas, sigmas;
	while(getline(file,inputline)){
		if(inputline.find(mat.name)!=std::string::npos){
			std::stringstream ss(inputline);
			ss>>matname>>f0>>gamma0;
			while(ss){
				ss>>f>>gamma>>omega>>sigma;
				if(f==0) break;
				fs.push_back(f);
				gammas.push_back(gamma);
				omegas.push_back(omega);
				sigmas.push_back(sigma);
			}
			break;
		}
	}
	double N_n = mat.rho*N_A/mat.n;
	std::cout<<"Working Piece "<<matname<<":"<<std::endl;
	std::cout<<"Nuclei density: "<< N_n<<std::endl;
	std::cout<<"BB model:"<<std::endl;
	std::cout<<"f0="<<f0<<", gamma0="<<gamma0<<std::endl;
	for (uint i = 0; i < fs.size(); ++i)
	{
		std::cout<<"f="<<fs[i]<<", gamma="<<gammas[i]<<", omega="<<omegas[i]<<", sigma="<<sigmas[i]<<std::endl;
	}

	std::complex<double> chi=0;
	double omegaP = getOmegaP(mat);
	std::complex<double> myI = I;
	omega = 2*PI*c/lambda;
	double E_ph = omega*hbar/q_e;
	for (uint i = 0; i < fs.size(); ++i)
	{
		std::complex<double> aj = pow(E_ph,2)-myI*E_ph*gammas[i];
		aj=sqrt(aj);
		chi+=myI*sqrt(PI/8)*fs[i]*pow(omegaP*hbar/q_e,2)/aj/sigmas[i]*
		(Faddeeva::w((aj-omegas[i])/sqrt(2)/sigmas[i])+Faddeeva::w((aj+omegas[i])/sqrt(2)/sigmas[i]));
		//std::cout<<"chi: "<<chi<<std::endl;
	}
	double param1 = gamma0/(getOhmCollisionRate(mat.conductivity, N_n)+getOhmCollisionRate(1./mat.rhoee/90000, N_n));
	//std::cout<<"param1: "<<param1<<std::endl;
	for (int T_l = startT; T_l < endT; T_l*=logstepT)
	{
		Rs.emplace_back();
		alphas.emplace_back();
		double v_ep=getOhmCollisionRate(mat.conductivity*300/T_l, N_n);
		for (int T_e = startT; T_e < endT; T_e*=logstepT)
		{
			double v_ee=getOhmCollisionRate(1./mat.rhoee/T_e/T_e, N_n);
			double v_pl= 2*PI*getPlasmaCollisionRate(N_n, T_e)*hbar/q_e;
//			std::cout<<"plasma collision rate: "<<v_pl<<", Ohm collision rate: "<<v_ep<<", ee collision rate: "<<v_ee<<std::endl;
			double v_e = param1*(v_ep+v_ee);
			double gamma = (v_e*v_pl)/(v_e+v_pl);
//			std::cout<<"gamma: "<<gamma<<std::endl;

			if(T_l == 300){
				RsInt.emplace_back();
				alphasInt.emplace_back();
				for (double lambdaLoop = startInt; lambdaLoop < endInt; lambdaLoop+=intervalInt)
				{
					double E_phLoop = 2*PI*c/lambdaLoop*hbar/q_e;
					std::complex<double> chiLoop=0;
					for (uint i = 0; i < fs.size(); ++i)
					{
						std::complex<double> aj = pow(E_phLoop,2)-myI*E_phLoop*gammas[i];
						aj=sqrt(aj);
						chiLoop+=myI*sqrt(PI/8)*fs[i]*pow(omegaP*hbar/q_e,2)/aj/sigmas[i]*
						(Faddeeva::w((aj-omegas[i])/sqrt(2)/sigmas[i])+Faddeeva::w((aj+omegas[i])/sqrt(2)/sigmas[i]));
						//std::cout<<"chiLoop: "<<chiLoop<<std::endl;
					}
					std::complex<double> epsilonLoop = 1. - f0*pow(omegaP*hbar/q_e,2)/(E_phLoop*(E_phLoop+gamma*myI));
					std::complex<double> nLoop = sqrt(epsilonLoop+chiLoop);
					double f1Loop=nLoop.real();
					double f2Loop=nLoop.imag();
					RsInt.back().push_back((pow(f1Loop-1,2)+f2Loop*f2Loop) / (pow(f1Loop+1,2)+f2Loop*f2Loop));
					alphasInt.back().push_back(4*PI*f2Loop/lambda);
					//std::cout<<"epsilonD: "<<epsilonLoop<<std::endl;
					//std::cout<<"epsilon: "<<epsilonLoop+chiLoop<<std::endl;
					//std::cout<<"R: "<<RsInt.back().back()<<", Alpha: "<<alphasInt.back().back()<<std::endl;
				}
			}
			std::complex<double> epsilon = 1. - f0*pow(omegaP*hbar/q_e,2)/(E_ph*(E_ph+gamma*myI));
			std::complex<double> n = sqrt(epsilon+chi);
			double f1=n.real();
			double f2=n.imag();
			Rs.back().push_back((pow(f1-1,2)+f2*f2) / (pow(f1+1,2)+f2*f2));
			alphas.back().push_back(4*PI*f2/lambda);
			//std::cout<<"epsilonD: "<<epsilon<<std::endl;
			//std::cout<<"epsilon: "<<epsilon+chi<<std::endl;
			//std::cout<<"n: "<<n<<std::endl;
			//std::cout<<"R: "<<Rs.back().back()<<", Alpha: "<<alphas.back().back()<<std::endl;
			//return 0;
		}
	}
	if(doPlots){
		plot.plotSingle("Reflectivity",TsInUnit(),Rs[0],TUnit,"ReflvsT.pdf");
		plot.plotSingle("Skin Depth / m${^{-1}}$",TsInUnit(),alphas[0],TUnit,"alphavsT.pdf");
	}

	gsl_integration_workspace *iws = gsl_integration_workspace_alloc(5000);
	double result, error;

	struct rad_int_params radparams = {0, 0, this};

	F.function = &rad_int_Core;
	F.params = &radparams;
	std::vector<double> hs;
	std::vector<double> Ps;
	for (double T_e = startT; T_e < endT; T_e*=logstepT)
	{
		radparams.T_e = T_e;
		radP.emplace_back();
		for (double h = 0; h < 5e-8; h+=1e-9)
		{
			if(T_e == startT) hs.push_back(h);
			radparams.h = h;
			gsl_integration_qag(&F,1e-9,endInt,0,1,100,3,iws,&result,&error);
			radP.back().push_back(result);
			if(h == 0) Ps.push_back(result);
		}
	}
	if(doPlots){
		plot.plotSingle("Radiation power at 1000K /J${\\cdot}$(s${^{-1}\\cdot}$m${^{-3}}$)",hs,radP[convert(1000)],"Depth/m","radvsh.pdf");
		plot.plotSingle("Radiation power at surface /J${\\cdot}$(s${^{-1}\\cdot}$m${^{-3}}$)",Ts,Ps,TUnit,"radvsT.pdf");
	}
	gsl_integration_workspace_free(iws);

	return 0;
}

double CalcMaterial::averageZ(CalcMaterial *source){

	std::vector<double> C_e_orig; // electron specific heat ok
	std::vector<double> U_e_orig; // electron internal energy
	std::vector<double> mus_orig; // chemical potential
	std::vector<double> Z_l_orig;
	std::vector<double> Uin_orig;

	for (uint i = 0; i < Z_l.size(); ++i)
	{
		C_e_orig.push_back(C_e[i]);
		U_e_orig.push_back(U_e[i]);
		mus_orig.push_back(mus[i]);
		Z_l_orig.push_back(Z_l[i]);
		Uin_orig.push_back(Uin[i]);
		//double r2 = 1./2;
		//double r2 = sqrt((Z_l[i]*source->Z_l[i]))/mat.ele_number;
		double r2 = Z_l[i]/mat.ele_number;
		double r1 = 1-r2;
		Z_l[i] = r1*Z_l[i] + r2*source->Z_l[i];
		C_e[i] = r1*C_e[i] + r2*source->C_e[i];
		mus[i] = r1*mus[i] + r2*source->mus[i];
		if(i==0) {
			U_e[0] = 0;
			for (int j = 0; j < mat.ele_number; ++j)
			{
				//printf("ionize[%d]=%e\n",j,ionize[j]);
				U_e[0]-=ionize[j]/mat.n*mat.rho*1e6;
			}
		}else{
			U_e[i] = C_e[i]*(Ts[i]-Ts[i-1])+U_e[i-1];
		}
		Uin[i] = r1*Uin[i] + r2*source->Uin[i];
	}
	if(doPlots){
		if(DOSmode >= 2)
		plot.plotTripple("Degree of ionization, avg",
			"Degree of ionization, gas",
			"Degree of ionization, nuclei",
			TsInUnit(), Z_l, Z_l_orig, source->Z_l, TUnit, "Z_l_averaged.pdf");
		plot.plotTripple("${\\mu}$/eV avg", 
			"${\\mu}$/eV gas",
			"${\\mu}$/eV nuclei",
			TsInUnit(), mus, mus_orig, source->mus, TUnit, "mu_averaged.pdf");
		plot.plotTripple("${C_e}$/J${\\cdot}$K${^{-1}\\cdot}$m${^{-3}}$ avg",
			"${C_e}$/J${\\cdot}$K${^{-1}\\cdot}$m${^{-3}}$ gas",
			"${C_e}$/J${\\cdot}$K${^{-1}\\cdot}$m${^{-3}}$ nuclei",
			TsInUnit(), smooth(C_e), smooth(C_e_orig), smooth(source->C_e), TUnit, "C_e_averaged.pdf");
		plot.plotTripple("${U}$/J${\\cdot}$m${^{-3}}$ avg", 
			"${U}$/J${\\cdot}$m${^{-3}}$ gas",
			"${U}$/J${\\cdot}$m${^{-3}}$ nuclei",
			TsInUnit(), U_e, U_e_orig, source->U_e, TUnit, "U_e_averaged.pdf");
		plot.plotTripple("${\\bar{U}}$/J${\\cdot}$m${^{-3}}$ avg", 
			"${\\bar{U}}$/J${\\cdot}$m${^{-3}}$ gas",
			"${\\bar{U}}$/J${\\cdot}$m${^{-3}}$ nuclei",
			TsInUnit(), Uin, Uin_orig, source->Uin, TUnit, "Uin_averaged.pdf");
	}
	return 0;
}

void CalcMaterial::write(std::string filename){
	std::ofstream outfile;
	outfile.open(filename, std::ios::out);
	outfile<<"Te";
	for (uint i = 0; i < TsineV.size(); ++i)
	{
		outfile<<" "<<TsineV[i];
	}
	outfile<<std::endl<<"Z";
	for (uint i = 0; i < Z_l.size(); ++i)
	{
		outfile<<" "<<Z_l[i];
	}
	outfile<<std::endl<<"Uin";
	for (uint i = 0; i < Uin.size(); ++i)
	{
		outfile<<" "<<Uin[i];
	}
	printf("Written to file %s\n", filename.c_str());
	outfile.close();
}