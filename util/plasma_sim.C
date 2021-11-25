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
	//printf("pars %e, %e, %e, %e\n", x, order_in(*(par->vTs),x),interpolate(order_in(*(par->vTs),x),*(par->vcUinvsT))/par->Nn,par->param1);
	return 1.5*x*q_e*interpolate(order_in(*(par->vTs),x),*(par->vcZvsT))
	+interpolate(order_in(*(par->vTs),x),*(par->vcUinvsT))/par->Nn-par->param1;
}
int main(int argc, char const *argv[])
{
	Plot plot;
	plot.init("./FV");
		mkdir("FV/T",S_IRWXU);
		mkdir("FV/v",S_IRWXU);
		mkdir("FV/p",S_IRWXU);
		mkdir("FV/rho",S_IRWXU);
		mkdir("FV/mv",S_IRWXU);
		mkdir("FV/e",S_IRWXU);
		mkdir("FV/E",S_IRWXU);
	Experiment experiment;
	experiment.setLaser(1e-6, 10e-6, 250e-15, 1e5, 10);
	//experiment.setLaser(1e-6, 20e-6, 10e-12, 1e5, 20);
	//experiment.setLaser(1e-6, 20e-6, 10e-9, 1, 2e-3);
	double reflectance = 0.8;
	std::vector<double> vtime;
	std::vector<double> vx;
	std::vector<double> vTs;
	std::vector<std::vector<double>*> vvee;  //electron internal energy
	std::vector<std::vector<double>*> vveE;  //electron total energy
	std::vector<std::vector<double>*> vvev;  //electron velosity
	std::vector<std::vector<double>*> vverho;//electron density
	std::vector<std::vector<double>*> vvep;  //electron pressure
	std::vector<std::vector<double>*> vveT;  //electron temperature
	
	std::vector<std::vector<double>*> vvarho;//atom density
	std::vector<std::vector<double>*> vvav;  //atom velosity
	std::vector<std::vector<double>*> vvap;  //atom pressure
	std::vector<std::vector<double>*> vvaT;  //atom temperature
	std::vector<std::vector<double>*> vvZ;   //ionization depth

	std::vector<std::vector<double>*> vvEF;  //electric field

	std::vector<std::vector<double>*> vvZvsVT; //Z(V,Te) function
	std::vector<std::vector<double>*> vvUinvsVT; //U(V,Te) function

	material mat = Cu;
	double Nn=mat.rho*N_A/mat.n;
	double vol = 1;
	double airthickness = 1e-8;
	double matthickness = 1e-7;
	double dx = 1e-10;
	int nmeshpoints = round((airthickness+matthickness)/dx);

	readLine("V_"+std::to_string(vol),"T",vTs);
	std::vector<double> vVs;
	while(vol<20000) {
		vVs.push_back(vol);
		vvZvsVT.push_back(readLine("V_"+std::to_string(vol),"Z"));
		vvUinvsVT.push_back(readLine("V_"+std::to_string(vol),"Uin"));
		vol *= pow(2,0.1);
	}
	double T=0.03; //t=0, T=50eV
	double Z;
	double rhoe;
	double p;
	double dv;
	double logstept = pow(10,0.01);
	double absorbed = 0;
	double Uin,Tstep;
	std::vector<double> Us;
	std::vector<double> vcZvsT;
	std::vector<double> vcUinvsT;
	struct param par;
	const gsl_root_fsolver_type * Tg = gsl_root_fsolver_brent;
	gsl_root_fsolver * s = gsl_root_fsolver_alloc (Tg);
	gsl_function F;
	F.function = &func;
	F.params = &par;
	double dTtmp;
	double alpha;
	std::vector<double> vRef;
	double dp,de;
	double maxreflect=0;

	//Initial condition
	std::vector<double>* tmpvrho = new std::vector<double>();
	std::vector<double>* tmpve = new std::vector<double>();
	std::vector<double>* tmpvE = new std::vector<double>();
	std::vector<double>* tmpvv = new std::vector<double>();
	std::vector<double>* tmpvp = new std::vector<double>();
	std::vector<double>* tmpvT = new std::vector<double>();
	std::vector<double>* tmpvZ = new std::vector<double>();
	for (int i = -airthickness/dx; i < 0; ++i)
	{
		tmpvrho->push_back(0);
		tmpve->push_back(0);
		tmpvE->push_back(0);
		tmpvv->push_back(0);
		tmpvp->push_back(0);
		tmpvZ->push_back(0);
		vx.push_back(i*dx);
		tmpvT->push_back(0.03);
	}
	for (int i = 0; i < matthickness/dx ; ++i)
	{
		tmpvrho->push_back(mat.rho);
		//tmpve->push_back(1./(2./3/0.03/q_e/(mat.rho*dx/mat.n*N_A)));
		//tmpvE->push_back(1./(2./3/0.03/q_e/(mat.rho*dx/mat.n*N_A)));
		tmpve->push_back(0);
		tmpvE->push_back(-mat.bindingE*Nn*q_e*dx);
		tmpvv->push_back(0);
		tmpvp->push_back(0);
		tmpvZ->push_back(0);
		vx.push_back(i*dx);
		tmpvT->push_back(0.03);
	}
	vvrho.push_back(tmpvrho);
	vve.push_back(tmpve);
	vvE.push_back(tmpvE);
	vvv.push_back(tmpvv);
	vvp.push_back(tmpvp);
	vvT.push_back(tmpvT);
	vvZ.push_back(tmpvZ);
	//Start evolution
	std::vector<double> phipl;
	std::vector<double> phipr;
	std::vector<double> phie;
	std::vector<double> phi;
	std::vector<double> tmpvve;
	double CFLmax = 0;
	double CFLmaxpoint;
	double CFLmaxt;
	double endt = 1e-9;
	double dt = 1e-16;
	double recorddt = 1e-14;
	double recorded = 0;
	double ndelete = 0;
	for (double time = 0; time < endt; time+=dt)
	{
		//dt = (logstept-1)*time;
		double Ileft=1;
		double desum = 0;
		tmpvv = new std::vector<double> ();
		tmpve = new std::vector<double> ();
		tmpvE = new std::vector<double> ();
		tmpvT = new std::vector<double> ();
		tmpvp = new std::vector<double> ();
		tmpvZ = new std::vector<double> ();
		tmpvrho = new std::vector<double> ();
		phipl.clear();
		phipr.clear();
		phie.clear();
		phi.clear();
		tmpvve.clear();
		vtime.push_back(time);
		//Left bondary condition
		tmpvv->push_back(0);
		tmpvve.push_back(0);
		tmpvrho->push_back(0);
		tmpvp->push_back(0);
		tmpvT->push_back(0);
		tmpve->push_back(0);
		tmpvE->push_back(0);
		tmpvZ->push_back(0);
		//Find phase boundaries here.
		
		//Calculate fluency
		for (int i = 1; i < nmeshpoints+1 ; ++i)
		{
			//double rhoh = (vvrho.back()->at(i-1)+vvrho.back()->at(i))/2;
			//double vh = rhoh?(vvv.back()->at(i-1)+vvv.back()->at(i))/2/rhoh:0;
			//double ph = (vvp.back()->at(i-1)+vvp.back()->at(i))/2;
			//double eh = (vve.back()->at(i-1)+vve.back()->at(i))/2;
			//phip.push_back(ph*dt+rhoh*vh*vh*dt);
			//phie.push_back(vh*dt*(eh+ph)+0.5*rhoh*vh*vh*vh*dt);
			//phi.push_back(vh*dt*rhoh);
			phipl.push_back(0);
			phipr.push_back(0);
			phie.push_back(0);
			phi.push_back(0);
			for (int j = 1; j >=0 ; --j)
			{
				if(i-j==nmeshpoints) continue;
				double rhol = vvrho.back()->at(i-j);
				bool isSolid = vvZ.back()->at(i-j)<1 && rhol>0.9*mat.rho;
				if(rhol==0) {
					continue;
				}
				double vl = vvv.back()->at(i-j)/rhol/dx;
				if(j==0) vl = -vl;
				double Al = sqrt(q_e*vvT.back()->at(i-j)/2/PI/mat.n*N_A);
				double Bl = Al==0?0:exp(-vl*vl/4/PI/Al/Al);
				double Cl = Al==0?0:(1-erfc(vl/2/sqrt(PI)/Al)/2);
				double gammal = isSolid?0.5*vl*rhol:((Al*Bl+vl*Cl)*rhol);

				double rrho = rhol/mat.rho;
				double pmol = mat.bindingE*Nn*q_e*2*(pow(rrho,4)-pow(rrho,2));
				double emol = mat.bindingE*Nn*q_e*(pow(rrho,4)-2*rrho*rrho);
				double dpl = pmol+isSolid?0.5*vl*fabs(vl)*rhol:(2*PI*Al*Al*Cl*rhol+vl*gammal);
				//if(fabs(i-airthickness/dx)<10 && time >= 319e-15) printf("rrho=%e, pmol=%e, dpl=%e, i=%d, j=%d\n", rrho, pmol,(2*PI*Al*Al*Cl*rhol+vl*gammal),i,j);
				double del = isSolid?gammal*(0.5*vvE.back()->at(i-j)*fabs(vl)):((4*PI*pow(Al,3)*Bl+5*PI*Al*Al*vl*Cl)*rhol+0.5*vl*vl*gammal+emol*gammal/rhol);
				double bouncebackrate = isSolid?0:(i-!j==nmeshpoints?1:std::min(pow(vvrho.back()->at(i-!j)/mat.rho,2./3),1.));
				double r0 = pow(mat.n/N_A/mat.rho,1./3);
				double CFL = gammal/rhol*dt/dx;
				if(CFLmax<CFL){
					CFLmax = CFL;
					CFLmaxt = time;
					CFLmaxpoint = i;
				}
				if(i==nmeshpoints){
					//Right bondary condition
					phipl[phipl.size()-1]=dpl*dt;
					phipr[phipr.size()-1]=-dpl*dt;
				}else {

					if(j==1){
							phipl[phipl.size()-1]+=dpl*dt*(1+bouncebackrate);//*(1-bouncebackrate)+1*dpl*dt*bouncebackrate;
							//phipr[phipl.size()-1]-=dpl*dt*(1-bouncebackrate)+2*dpl*dt*bouncebackrate;
							phie[phie.size()-1]+=del*dt*(1-bouncebackrate);
							phi[phi.size()-1]+=gammal*dt*(1-bouncebackrate);
					}else{
						phipr[phipr.size()-1]-=dpl*dt;//*(1-bouncebackrate)+1*dpl*dt*bouncebackrate;
						//phipl[phipl.size()-1]+=dpl*dt*(1-bouncebackrate)+2*dpl*dt*bouncebackrate;
						phie[phie.size()-1]-=del*dt*(1-bouncebackrate);
						phi[phi.size()-1]-=gammal*dt*(1-bouncebackrate);
					}
				}
				//if(fabs(i-airthickness/dx)<10 && time >= 319e-15)// && (i==100||i==101))
				//	printf("isSolid %d, Al=%e, Bl=%e, Cl=%e, v=%e, gammal=%e, del=%e, dpl=%e, pmol=%e, i=%d, j=%d\n",isSolid, Al, Bl, Cl, vvv.back()->at(i-j)/rhol/dx, gammal*dt,del*dt,dpl*dt,pmol, i,j);
			}//
				//if(abs(i-airthickness/dx)<10 && time >= 319e-15)
				//	printf("gammal=%e, del=%e, i=%d\n", phi.back(),phie.back(),i);
		}
		//Calcualte changes
		for (int i = 0; i < nmeshpoints-1; ++i)
		{
			double tmprho = vvrho.back()->at(i+1);
			Tstep = order_in(vTs,vvT.back()->at(i+1));
			tmpvrho->push_back(tmprho+(phi[i]-phi[i+1])/dx);
			de = phie[i]-phie[i+1];
			double dE = de;
			dv = phipl[i]-phipr[i]-phipl[i+1]+phipr[i+1];
			//if(fabs(i-airthickness/dx)<10 && time>319e-15) printf("dp1=%e, dp2=%e, dp=%e, drho=%e\n", phipl[i]-phipr[i], 
			//	-phipl[i+1]+phipr[i+1], dv, 
			//	tmprho+(phi[i]-phi[i+1])/
			//	dx-(i>airthickness/dx?mat.rho:0));

//			if(dv>0) exit(0);
			double dep = 0;
			double rrho = tmpvrho->back()/mat.rho;
			double rrhob = tmprho/mat.rho;
			double demol = 0;//mat.bindingE*Nn*q_e*(pow(rrho,4)-2*rrho*rrho-pow(rrhob,4)+2*rrhob*rrhob)*dx;
			if(tmprho!=0) dep=((vvv.back()->at(i+1)+dv)/tmpvrho->back()*(vvv.back()->at(i+1)+dv)-(vvv.back()->at(i+1))/tmprho*(vvv.back()->at(i+1)))/2/dx;
			else if(tmpvrho->back()!=0) dep=(vvv.back()->at(i+1)+dv)/2/tmpvrho->back()/dx*(vvv.back()->at(i+1)+dv);
			//if(tmpvrho->back()!=0 && abs(i-airthickness/dx)<10 && time >= 319e-15) 
			//printf("%d, rho=%e, demol=%e, de=%e, dv=%e, v=%e\n", i, tmprho, demol, de, dv, (vvv.back()->at(i+1)+dv)/tmpvrho->back());
			de-=dep;
			de-=demol;
			//if(tmpvrho->back()<1e-10 && de < 0) de=0;

			if(tmpvrho->back()>1e-200){
				vcZvsT.clear();
				interpolate(vvZvsVT,vcZvsT,order_in(vVs,1./tmprho));
				//printf("1./tmprho=%e, vcZvsT.size()=%d, i=%d\n", 1./tmprho, vcZvsT.size(), i);
				Z=interpolate(Tstep,vcZvsT);
				//printf("Z=%e\n", Z);
				vcUinvsT.clear();
				interpolate(vvUinvsVT,vcUinvsT,order_in(vVs,1./tmprho));
				Uin = interpolate(Tstep,vcUinvsT)/Nn;
				//printf("Uin=%e\n", Uin);
				rhoe=Z*tmprho/mat.n*N_A;
				T = vvT.back()->at(i+1);
				if(tmprho==0) T = vvT.back()->at(i+2);
				if(tmprho>1e-20){
					//alpha = 1e8*tmprho;
					getRandAlphaMat(mat, rhoe, mat.rho/tmprho, 1e-6, 300, vvT.back()->at(i+1)*eVtoK, reflectance, alpha);
					if(reflectance>maxreflect) maxreflect = reflectance;
					alpha = tmprho>1?alpha*tmprho/mat.rho:0;
					//if(tmpvrho->back()>8000)
					double abs =experiment.laserStrength(time-2*experiment.drilling.pulseWidth,0)*dt*(1-maxreflect)*dx*alpha*Ileft;
					de += abs;
					dE += abs;
					//printf("Uin=%e, Z=%e, T=%f, de=%e laser=%e\n", Uin,Z,vvT.back()->at(i+1),de);
					//T=2./3*(de+vve.back()->at(i+1))/q_e/(tmpvrho->back()*dx/mat.n*N_A);
					//printf("Z=%e, rho=%e, T=%f, ref=%f, alpha=%f, Ileft=%e\n", Z, tmpvrho->back(), T*eVtoK, reflectance, alpha, Ileft);
					//printf("rho=%e, drho = %e, de=%e, T=%e, dv=%e, dep=%e, i=%d\n", tmprho, (phi[i]-phi[i+1])/dx,de,T,dv,dep, i);
					//if(tmpvrho->back()>2000)
					Ileft *= -(dx*alpha-1);
				}
				if(!isfinite(T)||T<0) exit(0);
				par.param1 = 3./2*q_e*Z*vvT.back()->at(i+1)+de/(tmpvrho->back()*dx/mat.n*N_A)+Uin;
				//if(time>319e-15 && abs(i-airthickness/dx)<10) printf("param1=%e, %d\n", par.param1, i);
				par.Nn = Nn;
				par.vcZvsT = &vcZvsT;
				par.vcUinvsT = &vcUinvsT;
				par.vTs = &vTs;
				gsl_root_fsolver_set (s, &F, T*0.01,T*100);
				while(1){
					gsl_root_fsolver_iterate(s);
					T=gsl_root_fsolver_root(s);
					if(fabs(T-dTtmp)<1e-100) break;
					dTtmp=T;
				}
				//if(time>319e-15 && abs(i-airthickness/dx)<10) printf("T=%e, %d\n", T, i);
				p=rhoe*T*q_e;
			}else{
				T=0;
				Z=0;
				p=0;
			}

			tmpvv->push_back(vvv.back()->at(i+1)+dv);
			tmpve->push_back(vve.back()->at(i+1)+de);
			tmpvve.push_back(tmpvv->back()/tmpvrho->back()/dx);
			tmpvT->push_back(T);
			tmpvp->push_back(p);
			tmpvZ->push_back(Z);
			tmpvE->push_back(vvE.back()->at(i+1)+dE);
			(*tmpvv)[0] = tmpvv->at(1);
			(*tmpvrho)[0] = tmpvrho->at(1);
			(*tmpvp)[0] = tmpvp->at(1);
			(*tmpvT)[0] = tmpvT->at(1);
			(*tmpve)[0] = tmpve->at(1);
			(*tmpvE)[0] = tmpvE->at(1);
			(*tmpvZ)[0] = tmpvZ->at(1);
			tmpvve[0] = tmpvve[1];
		}
		//if(time >= 181e-15) exit(0);
//		plot.plotSingle("rho"+std::to_string(time*1e15),vx,*tmpvrho,"x/m","rho"+std::to_string(time*1e15)+".png");
		//if(time>5000e-15) {
		//printf("CFLmax: %f, at x=%f, t=%e\n", CFLmax, CFLmaxpoint, CFLmaxt);
		if(round((time-recorded)*1e15)>=0 &&time>320e-15) {
			
			if(time>1000e15) recorded=10*recorddt+time;
			if(time>10000e15) recorded=100*recorddt+time;
			else recorded=recorddt+time;
			plot.plotSingle("T"+std::to_string(time*1e15),vx,*tmpvT,"x/m","T/"+std::to_string(time*1e15)+".png");
			plot.plotSingle("v"+std::to_string(time*1e15),vx,tmpvve,"x/m","v/"+std::to_string(time*1e15)+".png");
			plot.plotSingle("p"+std::to_string(time*1e15),vx,*tmpvp,"x/m","p/"+std::to_string(time*1e15)+".png");
			plot.plotSingle("rho"+std::to_string(time*1e15),vx,*tmpvrho,"x/m","rho/"+std::to_string(time*1e15)+".png");
			plot.plotSingle("mv"+std::to_string(time*1e15),vx,*tmpvv,"x/m","mv/"+std::to_string(time*1e15)+".png");
			plot.plotSingle("e"+std::to_string(time*1e15),vx,*tmpve,"x/m","e/"+std::to_string(time*1e15)+".png");
			plot.plotSingle("E"+std::to_string(time*1e15),vx,*tmpvE,"x/m","E/"+std::to_string(time*1e15)+".png");
		}
		vvrho.push_back(tmpvrho);
		vve.push_back(tmpve);
		vvE.push_back(tmpvE);
		vvv.push_back(tmpvv);
		vvp.push_back(tmpvp);
		vvT.push_back(tmpvT);
		vvZ.push_back(tmpvZ);
		if(vvrho.size()>2){
		delete vvrho[ndelete];
		delete vve[ndelete];
		delete vvE[ndelete];
		delete vvv[ndelete];
		delete vvp[ndelete];
		delete vvT[ndelete];
		delete vvZ[ndelete];
		ndelete+=1;
		}
		absorbed+=experiment.laserStrength(time-2*experiment.drilling.pulseWidth,0)*dt*(1-maxreflect);
	}/*
	printf("absorbed energy: %e\n", absorbed);
	printf("dp1: %e, dp2: %e, U: %e\n", dp1,dp2,U);

	plot.plotSingle("Reflectance",vtime,vRef,"t/s","R.png");
	plot.plotSingle("T/eV",vtime,vT,"t/s","Tvst.png");
	plot.plotSingleLogXY("vol",vtime,vvol,"t/s","vol.png");
	plot.plotSingleLogXY("${\\rho_e}$",vtime,vrhoe,"t/s","rhoe.png");
	plot.plotSingleLogXY("p/Pa",vtime,vp,"t/s","p.png");
	plot.plotSingle("v/m/s",vtime,vv,"t/s","v.png");
	plot.plotSingle("Z",vtime,vZ,"t/s","Z.png");
	plot.plotSingle("U",vtime,Us,"t/s","U.png");
	*/
	return 0;
}