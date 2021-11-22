#include "formulas.h"
#include "constants.h"
#include <algorithm> 
#include <iostream>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_complex.h>
#include <complex.h>
double getFermiEnergy_eV(material mat){
	return pow(3*mat.rho/mat.n/PI,2./3)*pow(h*pow(N_A,1./3)/sqrt(m_e)/sqrt(q_e),2)/8;
}

double getFermiT(material mat){
	return pow(3*mat.rho/mat.n/PI,2./3)*pow(h/sqrt(k_B)*pow(N_A,1./3)/sqrt(m_e),2)/8;
}

double getFermiV(material mat){
	return sqrt(2*getFermiEnergy_eV(mat)*q_e/m_e);
}

double getOmegaP(material mat){
	return sqrt(mat.rho*N_A/mat.n*q_e*q_e/m_e/epsilon_0);
}

double getOmegaP(double n){
	return sqrt(n*q_e*q_e/m_e/epsilon_0);
}

void getEpsilon(double omegaP, double omega, double gamma, double &ep1, double &ep2){
	std::cout<<"omegaP: "<<omegaP<<", omega: "<<omega<<std::endl;
	ep1 = (1-omegaP*omegaP/(pow(omega,2)+gamma*gamma));
	ep2 = (omegaP*omegaP*gamma/(pow(omega,3)+gamma*gamma*omega));
	std::cout<<"epsilon_r: "<<ep1<<"+i"<<ep2<<std::endl;
}

double getPlasmaCollisionRate(double n, double T){
	return 4./3*sqrt(2*PI)*n*pow(q_e,4)/sqrt(m_e)/pow(4*PI*epsilon_0,2)/pow(k_B*T,3./2)*10;
}

double getOhmCollisionRate(double sigma, double n){
	if(n<0) return 0;
	return n*q_e*q_e/m_e/sigma;
}

void getRandAlpha(double ep1, double ep2, double lambda, double &R, double &Alpha, double &f1, double &f2){
	f1 = sqrt( (ep1+sqrt(ep1*ep1+ep2*ep2))/2 );
	f2 = sqrt( (-ep1+sqrt(ep1*ep1+ep2*ep2))/2 );
	std::cout<<"n="<<f1<<", k="<<f2<<std::endl;
	R = (pow(f1-1,2)+f2*f2) / (pow(f1+1,2)+f2*f2);
	Alpha = 4*PI*f2/lambda;
}
void getRandAlphaMat(material mat, double N_e, double rV, double lambda, double T_l, double T_e, 
	double &R, double &Alpha, double &f1, double &f2, double &ep1, double &ep2){
	double param1 = 1;
	double Nn = mat.rho*N_A/mat.n/rV;
	double Z_av = N_e/Nn;
	double A = 0.5e15;
	//double omega = 2*PI*c/lambda;
	//double omegaP = getOmegaP(rho,n);
	//double bmax = sqrt(k_B*T_e/m_e)*std::min(1./omegaP,1./omega);
	//double bmin = -std::min(-Z_av*q_e*q_e/k_B/T_e,-h/2/PI/sqrt(m_e*k_B*T_e));
	//double v_F=getFermiV(rho,n);
	double v_ep=getOhmCollisionRate(mat.conductivity*300/T_l, Nn-N_e);
	double v_ee=getOhmCollisionRate(A/T_e/T_e, Nn-N_e);
	//double v_ep= 26*q_e*q_e*k_B*T_l/h/h/getFermiV(rho,n_n)/epsilon_0;
	double v_pl= getPlasmaCollisionRate(Z_av*N_e, T_e);
	//std::cout<<"plasma collision rate: "<<v_pl<<", Ohm collision rate: "<<v_ep<<", ee collision rate: "<<v_ee<<std::endl;
	double ratio = 1-exp(-Z_av);
	//double v_e = (v_ep+v_ee)*(1-ratio)+v_pl*ratio;
	double v_e = v_ep+v_ee+v_pl;
	double gamma = 2*PI*v_e;
	//std::cout<<"gamma: "<<gamma<<std::endl;
	double omega = 2*PI*c/lambda;
	if(gamma>omega) gamma=omega;
	std::complex<double> epsilon =  1 - pow(getOmegaP(N_e+Nn),2)/(omega*omega+omega*gamma*I);
	ep1=epsilon.real();
	ep2=epsilon.imag();
	std::complex<double> n = sqrt(epsilon);
	f1=n.real();
	f2=n.imag();
	R = (pow(f1-1,2)+f2*f2) / (pow(f1+1,2)+f2*f2);
	Alpha = 4*PI*f2/lambda;
	//printf("R=%f, alpha=%f\n", R, Alpha);
	//getEpsilon(getOmegaP(N_e)/param2, 2*PI*c/lambda, gamma, ep1, ep2);
	//getRandAlpha(ep1, ep2, lambda, R, Alpha, f1, f2);
}

void getRandAlphaMat(material mat, double N_e, double rV, double lambda, double T_l, double T_e, 
	double &R, double &Alpha){
	double f1,f2,ep1,ep2;
	getRandAlphaMat(mat, N_e, rV, lambda, T_l, T_e, R, Alpha, f1, f2, ep1, ep2);
}

double debye_func(double y, void *p){
	return pow(y,3)/(exp(y)-1);
}

double debye_int(double x){
	gsl_integration_workspace *iws = gsl_integration_workspace_alloc(500);
	double result, error;

	gsl_function F;

	F.function = &debye_func;

	gsl_integration_qag(&F,1e-10,x,0,1e-7,100,3,iws,&result,&error);

	gsl_integration_workspace_free(iws);
	return 3*result/pow(x,3);
}