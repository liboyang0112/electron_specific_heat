#ifndef MATERIAL
#define MATERIAL
#include <string>

typedef struct {
	double lambda;
	double FWHM;
	double pulseWidth;
	double pulseFrequency;
	double pulseEnergy;
} laser;

//https://www.materialsproject.org/rest/v2/materials/mp-989782/vasp?API_KEY=RdVMCd0C5rJQalqe
typedef struct {
	std::string name;
	double ele_number;
	double rho;
	double n;
	double conductivity;
//	double A,B; //K. Vestentoft and P. Balling, Appl. Phys. A 84, 207 (2006)
	double rhoee; //A. H. MacDonald, Phys. Rev. Lett. 44, 489 (1980).
	double debyeTemprature; //Zhicheng Wang, Relixue Tongjiwuli
	double lambdaOmega2_meV; // Z. Lin, L. V. Zhigilei, and V. Celli, Phys. Rev. B 77, 075133 (2008).
	double workFunc; //https://chem.libretexts.org/Ancillary_Materials/Reference/Reference_Tables/Bulk_Properties/B1%3A_Workfunction_Values_(Reference_Table)
//	double epsiloninfty;
	double bindingE; 
	double YoungsModulus; 
} material;

const material Cu = {
	"Cu",
	29,
	8.960e3,
	63.546,
	5.96e7,
//	1.28e7, 1.23e11,
	2.0e-15,
	340,
	29,
	5,
//	3.8
	5,
	2e-11
};
const material Al = {
	"Al",
	13,
	2.702e3,
	26,
	5.96e7,
//	1.28e7, 1.23e11,
	2.0e-15,
	340,
	29,
	5,
//	3.8
	5,
	2e-11
};
const material Ag = {
	"Ag",
	47,
	10.49e3,
	107.8682,
	6.30e7,
//	0.932e7, 1.02e11,
	2.0e-15,
	215,
	22.5,
	4.5,
//	3.8
	5,
	2e-11
};
const material Au = {
	"Au",
	79,
	19.32e3,
	196.97,
	4.10e7,
//	1.18e7, 1.25e11,
	2.5e-15,
	169,
	23,
	5.2,
//	3.8
	5,
	2e-11
};


#endif