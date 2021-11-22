#include <cmath>
#include "materials.h"

double getFermiEnergy_eV(material mat);

double getFermiT(material mat);

double getFermiV(material mat);

double getOmegaP(material mat);

double getOhmCollisionRate(double T, double n);
double getPlasmaCollisionRate(double T, double n);

void getRandAlpha(double ep1, double ep2, double lambda, double &R, double &Alpha, double &f1, double &f2);
void getEpsilon(double omegaP, double omega, double gamma, double &ep1, double &ep2);

void getRandAlphaMat(material mat, double lambda, double T_l, double T_e,
	double &R, double &Alpha, double &f1, double &f2, double &ep1, double &ep2);
void getRandAlphaMat(material mat, double N_e, double rV, double lambda, double T_l, double T_e, 
	double &R, double &Alpha);

double debye_int(double x);