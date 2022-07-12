#include "algorithm.h"
#include "materials.h"
#include "constants.h"
#include "Plot.h"
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
#include "Experiment.h"

struct param{
	double param1;
	//double param2;
	//double param3;
	std::vector<double> *vcZvsT;
	std::vector<double> *vcUinvsT;
	//double T;
	double Nn;
	std::vector<double> *vTs;
};

double func(double x, void* p){
	struct param * par = (struct param*) p;
	//printf("pars %e, %e, %e,%e, %e\n", par->param1,par->param2,par->param3,par->T,interpolate(order_in(*(par->vTs),x+par->T),*(par->vcZvsT)));
	//printf("pars %e\n", par->param1*x+par->param2*interpolate(order_in(*(par->vTs),x+par->T),*(par->vcZvsT))-par->param3);
	//return par->param1*x+par->param2*interpolate(order_in(*(par->vTs),x+par->T),*(par->vcZvsT))-par->param3;
	//return par->param1*x+par->param2*interpolate(order_in(*(par->vTs),x+par->T),*(par->vcZvsT))-par->param3;
	//printf("pars %f, %f, %e, %e\n", x, order_in(*(par->vTs),x),interpolate(order_in(*(par->vTs),x),*(par->vcUinvsT))/par->Nn,par->param1);
	return 1.5*x*q_e*interpolate(order_in(*(par->vTs),x),*(par->vcZvsT))
	+interpolate(order_in(*(par->vTs),x),*(par->vcUinvsT))/par->Nn-par->param1;
}
int main(int argc, char const *argv[])
{
	Plot plot;
	plot.init(".");
	Experiment exp;
	exp.setLaser(1e-6, 10e-6, 250e-15, 1e5, 10);
	//exp.setLaser(1e-6, 20e-6, 10e-12, 1e5, 20);
	//exp.setLaser(1e-6, 20e-6, 10e-9, 1, 2e-3);
	double reflectance = 0.8;
	std::vector<double> vtime;
	std::vector<double> vZ;
	std::vector<double> vrhoe;
	std::vector<double> vp;
	std::vector<double> vv;
	std::vector<double> vvol;
	std::vector<double> vT;
	std::vector<double> vTs;
	std::vector<double> vZvsT[2];
	std::vector<double> vUinvsT[2];
	material mat = Cu;
	double Nn=mat.rho*N_A/mat.n;
	double vol = 1;
	double volstep = 0.1;
	double thickness = 1e-8;
	double uppervol = vol*pow(2,0.1);
	readLine("V_"+std::to_string(vol),"T",vTs);
	readLine("V_"+std::to_string(vol),"Z",vZvsT[0]);
	readLine("V_"+std::to_string(vol),"Uin",vUinvsT[0]);
	readLine("V_"+std::to_string(uppervol),"Z",vZvsT[1]);
	readLine("V_"+std::to_string(uppervol),"Uin",vUinvsT[1]);
	bool currentV = 0;
	double T=0.03; //t=0, T=50eV
	double v=0;
	double A = 1;
	double Z;
	double rhoe;
	double p;
	double a;
	double endt = 1e-9;
	double dt = 1e-13;
	double logstept = pow(10,0.001);
	double absorbed = 0;
	double dl,dU,tmpdU,Uin,Tstep;
	std::vector<double> Us;
	std::vector<double> vcZvsT;
	std::vector<double> vcUinvsT;
	struct param par;
	const gsl_root_fsolver_type * Tg = gsl_root_fsolver_brent;
	gsl_root_fsolver * s = gsl_root_fsolver_alloc (Tg);
	gsl_function F;
	F.function = &func;
	F.params = &par;
	double dp1=0;
	double dp2=0;
	double dTtmp;
	double U=0;
	double alpha;
	std::vector<double> vRef;
	for (double time = 1e-15; time < endt; time+=dt)
	{
		dt = (logstept-1)*time;
		vtime.push_back(time);
		Tstep = order_in(vTs,T);
		vcZvsT.clear();
		interpolate(vZvsT[currentV],vZvsT[!currentV],vcZvsT,1-(uppervol-vol)/volstep);
		Z=interpolate(Tstep,vcZvsT);
		vcUinvsT.clear();
		interpolate(vUinvsT[currentV],vUinvsT[!currentV],vcUinvsT,1-(uppervol-vol)/volstep);
		Uin = interpolate(Tstep,vcUinvsT)/Nn;
		rhoe=Z/vol*Nn;
		p=rhoe*T*q_e;
		a = 3*p/mat.rho/vol/thickness/A;
		dl = v*dt+1./2*a*dt*dt;
		vol+=dl*A/thickness;
		dp2+=1./6*mat.rho*vol*thickness*A*(2*v*a*dt);
		v+=a*dt;
		dp1+=-p*dl*A;
		printf("rhoe=%e\n", rhoe);
		getRandAlphaMat(mat, rhoe, vol, 1e-6, 300, T*eVtoK, reflectance, alpha);
		vRef.push_back(reflectance);
		dU = -p*dl*A+exp.laserStrength(time-2*exp.drilling.pulseWidth,0)*dt*(1-reflectance)*A;
		absorbed+=exp.laserStrength(time-2*exp.drilling.pulseWidth,0)*dt*(1-reflectance);
		U+=dU;
		Us.push_back(U);
		printf("t= %e, rhoe=%e, p=%e, a=%e, dl=%e, dU=%e\n", time, rhoe,p,a,dl,dU/Nn/thickness);
		//par.param1 = 3./2*q_e*Z;
		//par.param2 = Uin+1.5*T*q_e;
		//par.param3 = dU/Nn/thickness+Z*par.param2;
		//par.T = T;
		par.param1 = 3./2*q_e*Z*T+dU/Nn/thickness+Uin;
		par.Nn = Nn;
		par.vcZvsT = &vcZvsT;
		par.vcUinvsT = &vcUinvsT;
		par.vTs = &vTs;
		gsl_root_fsolver_set (s, &F, T*0.5,T*2);
		while(1){
			gsl_root_fsolver_iterate(s);
			T=gsl_root_fsolver_root(s);
			if(fabs(T-dTtmp)<1e-100) break;
			dTtmp=T;
		}
		//T += dT;
		//printf("dT = %e, T= %e\n", dT,T);
		//exit(0);
		vT.push_back(T); 
		vZ.push_back(Z);
		vrhoe.push_back(rhoe);
		vvol.push_back(vol);
		vp.push_back(p);
		vv.push_back(v);
		if(vol>uppervol){
			while(vol>uppervol) {
				volstep = uppervol*(pow(2,0.1)-1);
				uppervol *= pow(2,0.1);
				//printf("uppervol=%e\n", uppervol);
			}
			if(uppervol>=20000) break;
			vZvsT[currentV].clear();
			vUinvsT[currentV].clear();
			readLine("V_"+std::to_string(uppervol),"Z",vZvsT[currentV]);
			readLine("V_"+std::to_string(uppervol),"Uin",vUinvsT[currentV]);
			currentV = !currentV;
		}
	}
	vp[0]=vp[1];
	printf("absorbed energy: %e\n", absorbed);
	printf("dp1: %e, dp2: %e, U: %e\n", dp1,dp2,U);

	plot.plotSingle("Reflectance",vtime,vRef,"t/s","R.pdf");
	plot.plotSingle("T/eV",vtime,vT,"t/s","Tvst.pdf");
	plot.plotSingleLogXY("vol",vtime,vvol,"t/s","vol.pdf");
	plot.plotSingleLogXY("${\\rho_e}$",vtime,vrhoe,"t/s","rhoe.pdf");
	plot.plotSingleLogXY("p/Pa",vtime,vp,"t/s","p.pdf");
	plot.plotSingle("v/m/s",vtime,vv,"t/s","v.pdf");
	plot.plotSingle("Z",vtime,vZ,"t/s","Z.pdf");
	plot.plotSingle("U",vtime,Us,"t/s","U.pdf");
	return 0;
}