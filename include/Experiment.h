#include "materials.h"
#include "formulas.h"
#include "CalcMaterial.h"
#include <string>
#include "constants.h"
#include "Plot.h"
class Experiment
{
public:
	Experiment(){}
	~Experiment(){}
	double temprature;
	double maxTemprature;
	double tempraturePrecision;
	double pressure;
	laser drilling;
	material workPiece;
	CalcMaterial *matData;
	CalcMaterial *gasData;
	double getElectronEscapeRate(double T_e){
		return PI*drilling.FWHM*drilling.FWHM*4*PI*m_e/pow(h,3)*pow(k_B*T_e,2)*
		exp(q_e*(-matData->fermiEnergy-workPiece.workFunc+matData->getMu(T_e))/k_B/T_e);
	}
	double getNElectronEscaped(double T_emax){
		double a = 1.5;
		return k_B*T_emax*epsilon_0/a/q_e/q_e*drilling.FWHM*
		log(1+4*PI*m_e/pow(h,3)*drilling.pulseWidth*PI*drilling.FWHM*a*k_B*T_emax*q_e*q_e/epsilon_0*
			exp((-matData->fermiEnergy-workPiece.workFunc+matData->getMu(T_emax))*q_e/k_B/T_emax));
	}
	double getLaserSource(double x, double y, double depth, double t, double T_e, double T_l){
		return drilling.pulseEnergy*matData->getAlpha(T_e, T_l)*(1-matData->getR(T_e, T_l))*
		exp(-matData->getAlpha(T_e, T_l)*depth-(x*x+y*y)/drilling.FWHM/drilling.FWHM/2-
			t*t/2/drilling.pulseWidth/drilling.pulseWidth);
	}
	void setEnvironment(double T, double Tmax, double Tprecison, double P){
		temprature = T;
		maxTemprature = Tmax;
		tempraturePrecision = Tprecison;
		pressure = P;
	}
	void initMaterial(material mat, std::string DOS){
		workPiece = mat;

		//matData = new CalcMaterial("free");
		//matData->lambda = drilling.lambda;
		//matData->mat = workPiece;
		//matData->freeDOS();
		//matData->calculate(temprature, maxTemprature, tempraturePrecision);

		
		double Vratio = 1;
		CalcMaterial ionData("ion");
		ionData.DOSmode=ionData.gas_plasma;
		//ionData.bareAtomOrbit=4;
		ionData.useIon=1;
		ionData.initMaterial(workPiece, DOS);
		matData = new CalcMaterial("lattice");
		matData->nEle=11;
		matData->useIon=1;
		matData->DOSmode = matData->lattice_plasma;
		matData->lambda = drilling.lambda;
		matData->initMaterial(workPiece, DOS);
		if(0) {
			ionData.VratioDefault = Vratio;
			ionData.calculate(temprature, maxTemprature, tempraturePrecision);
			exit(0);
		}
		
		if(0){
			matData->doPlots = 0;
			for (double i = 1; i < 20000; i*=pow(2,0.1))
			{
				Vratio = i;
				matData->VratioDefault = Vratio;
				matData->calculate(temprature, maxTemprature, tempraturePrecision);
				matData->write("V_"+std::to_string(Vratio));
			}
		}else{
			matData->VratioDefault = Vratio;
			matData->calculate(temprature, maxTemprature, tempraturePrecision);
			gasData = new CalcMaterial("gas");
			gasData->DOSmode = gasData->gas_plasma;
			gasData->initMaterial(workPiece, DOS);
			gasData->VratioDefault = Vratio;
			gasData->calculate(temprature, maxTemprature, tempraturePrecision);
			ionData.VratioDefault = Vratio;
			ionData.calculate(temprature, maxTemprature, tempraturePrecision);
			Plot plot;
			plot.init(".");
			plot.plotTripple("Degree of ionization, lattice",
				"Degree of ionization, gas",
				"Degree of ionization, nuclei",
				gasData->TsInUnit(), matData->Z_l, gasData->Z_l, ionData.Z_l, gasData->TUnit, "Z_l_comb.pdf");
			plot.plotTripple("${\\mu}$/eV lattice", 
				"${\\mu}$/eV gas",
				"${\\mu}$/eV nuclei",
				gasData->TsInUnit(), matData->mus, gasData->mus, ionData.mus, gasData->TUnit, "mu_comb.pdf");
			plot.plotTripple("${C_e}$/J${\\cdot}$K${^{-1}\\cdot}$m${^{-3}}$ lattice",
				"${C_e}$/J${\\cdot}$K${^{-1}\\cdot}$m${^{-3}}$ gas",
				"${C_e}$/J${\\cdot}$K${^{-1}\\cdot}$m${^{-3}}$ nuclei",
				gasData->TsInUnit(), matData->C_e, gasData->C_e, ionData.C_e, gasData->TUnit, "C_e_comb.pdf");
			plot.plotTripple("${U}$/J${\\cdot}$m${^{-3}}$ lattice", 
				"${U}$/J${\\cdot}$m${^{-3}}$ gas",
				"${U}$/J${\\cdot}$m${^{-3}}$ nuclei",
				gasData->TsInUnit(), matData->U_e, gasData->U_e, ionData.U_e, gasData->TUnit, "U_e_comb.pdf");
		//gasData->averageZ(&ionData);
		}
		
	}
	void setLaser(double lambda, double fwhm, double pulsewidth,
	 double pulsefrequency, double power){
		drilling.lambda = lambda;
		drilling.FWHM = fwhm;
		drilling.pulseWidth = pulsewidth;
		drilling.pulseFrequency = pulsefrequency;
		drilling.pulseEnergy = power/pulsefrequency;

		double nu = c/lambda;
		double E_ph = nu*h/q_e;
		std::cout<<"Laser Property: "<<std::endl<<
		"Photon wave length= "<<lambda*1e9<<"nm"<<std::endl<<
		"Photon Frequency= "<<nu<<"Hz"<<std::endl<<
		"Photon Energy="<<E_ph<<"eV"<<std::endl<<
		"Pulse power="<<power<<"W"<<std::endl<<
		"Pulse width="<<pulsewidth*1e15<<"fs"<<std::endl<<
		"Pulse frequency="<<pulsefrequency<<"Hz"<<std::endl<<
		"Pulse energy="<<power/pulsefrequency<<"J"<<std::endl<<
		"Instant flux="<<power/pulsefrequency/fwhm/fwhm/pulsewidth*8/pow(2*PI,1.5)*1e-4<<"W/cm2"<<std::endl<<
		"Focal diameter="<<fwhm*1e6<<"um"<<std::endl;
	}
	double laserStrength(double t, double dx){
		return drilling.pulseEnergy/pow(2*PI,1.5)/pow(drilling.FWHM/2,2)/drilling.pulseWidth*2
		*exp(-pow(t/drilling.pulseWidth,2)*2)*exp(-pow(dx/drilling.FWHM,2)*2);

	}
	void setLaser(laser las){
		drilling = las;
	}

};