#ifndef CALCMATERIAL
#define CALCMATERIAL
#include <vector>
#include <string>
#include "materials.h"
#include "algorithm.h"
#include "Plot.h"
class CalcMaterial
{

	double conductivity_a,conductivity_b,conductivity_c;
	std::vector<std::vector<double>> alphas; //skin depth (T_l, T_e)
	std::vector<std::vector<double>> Rs; //reflectivity (T_l, T_e)
	std::vector<std::vector<double>> alphasInt; //skin depth (T_l, T_e) For radiation integration
	std::vector<std::vector<double>> RsInt; //reflectivity (T_l, T_e) For radiation integration
	std::vector<std::vector<double>> radP; // Radiation power (T_e, depth)
	const double startInt = 200e-9;
	const double endInt = 1e-6;
	const double intervalInt = 1e-9;

public:
	CalcMaterial(std::string name){
		mkdir(name.c_str(),0777);
		plot.init(name);
	}
	~CalcMaterial(){}
	void initMaterial(material m, std::string filename);
	std::string TUnit="T/eV";
	std::vector<double> DOS;
	std::vector<double> epsilon_l; // energy level l
	std::vector<double> a_l;  // degeneration of l energy level
	std::vector<double> intDOS; // integrated DOS
	std::vector<double> C_e; // electron specific heat ok
	std::vector<double> C_l; // lattice specific heat
	std::vector<double> U_e; // electron internal energy
	std::vector<double> Gs; // electron-lattice heat conductivity
	std::vector<double> Ts; // temprature
	std::vector<double> TsineV; // temprature
	std::vector<double> TsinEpsilon; // temprature
	std::vector<double> mus; // chemical potential
	//std::vector<double> Ks; // electron heat conductivity
	std::vector<double> ys;
	std::vector<double> xs;
	std::vector<double> Z_l;
	Plot plot;
	material mat;
	double lambda;
	double fermiEnergy;
	double energyEnd;
	double interval;
	double nEle = 11;
	double startT;
	double endT;
	double intervalT;
	double latticeResolution = 1e-9;
	bool doPlots = 1;
	int bareAtomOrbit = 0;
	bool dimentionless = 0;
	enum mode {free, bloch, lattice_plasma, gas_plasma} DOSmode;
	void readDOS(std::string filename);
	void freeDOS();
	void plotDOSNandF(double x, double y, std::string savename);
	double DOSIntI(double x, double y, double power=0, double Vratio=1);
	double DOSInt(double x, double y, double power=0, double start=0, double end=0);
	double gasPlasmaDOSInt(double x, double y, double power=0, double Vratio=1);
	double gasPlasmaZ(double x, double y, double Vratio=1);
	double latticePlasmaDOSInt(double x, double y, double power=0);
	int calculate(double sT, double eT, double iT);
	int calculateOptical();
	std::vector<double> TsInUnit(){
		if(TUnit=="T/K") return Ts;
		else if(TUnit=="y") return ys;
		else return TsineV;
	}
	double computeG(double x, double y, double lambdaOmega2_meV);
	double DOS_incidence(double x, double y, double power=0, double start=0, double end=0);

	double convert(double T){return (T-startT)/intervalT;}

	double getC_e(double T_e){return getrTC_e(convert(T_e));}
	double getrTC_e(double rT){return interpolate(rT,C_e);}

	double getG(double T_e){return getrTG(convert(T_e));}
	double getrTG(double rT){return interpolate(rT,Gs);}

	double getC_l(double T_l){return getrTC_l(convert(T_l));}
	double getrTC_l(double rT){return interpolate(rT,C_l);}

	double getMu(double T_e){return getrTMu(convert(T_e));}
	double getrTMu(double rT){return interpolate(rT,mus);}

	double getAlpha(double T_e, double T_l){return getrTAlpha(convert(T_e), convert(T_l));}
	double getrTAlpha(double rT_e, double rT_l){return interpolate(rT_e,rT_l,alphas);}

	double getAlphaInt(double T_e, double lambda){
		if(lambda<startInt) lambda = startInt;
		return interpolate(convert(T_e), (lambda-startInt)/intervalInt, alphasInt);
	}
	double getRInt(double T_e, double lambda){
		if(lambda<startInt) lambda = startInt;
		return interpolate(convert(T_e), (lambda-startInt)/intervalInt, RsInt);
	}
	double latticePlasmaZ(double x, double y);
	double getR(double T_e, double T_l){return getrTR(convert(T_e), convert(T_l));}
	double getrTR(double rT_e, double rT_l){return interpolate(rT_e,rT_l,Rs);}
	double averageZ(CalcMaterial *source);
	double getK_e(double T_l, double T_e){return (conductivity_c+T_e)/(conductivity_a*T_l+conductivity_b*T_e*T_e)*getC_e(T_e);}
};
#endif