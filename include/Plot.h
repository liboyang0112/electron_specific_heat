#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;
#ifndef MYPLOT
#define MYPLOT
class Plot
{
public:
	Plot(){}
	~Plot(){}
	std::map<std::string, std::string> rcparam;
	std::map<std::string, std::string> legend_option;
	std::string outputpath;
	void init(std::string name);
	void plotFreeNandF(double maxEpsilon, int nPoints, double x2, double y, std::string savename);
	void plotSingle(std::string name, std::vector<double> x, std::vector<double> y, std::string xlabel, std::string savename);
	void plotSingleLogXY(std::string name, std::vector<double> x, std::vector<double> y, std::string xlabel, std::string savename);
	void plotSingleLogX(std::string name, std::vector<double> x, std::vector<double> y, std::string xlabel, std::string savename);
	void plotSingleLogY(std::string name, std::vector<double> x, std::vector<double> y, std::string xlabel, std::string savename);
	void plotDouble(std::string name1, std::string name2, std::vector<double> x, std::vector<double> y1, std::vector<double> y2, std::string xlabel, std::string savename);
	void plotTripple(std::string name1, std::string name2, std::string name3, std::vector<double> x, std::vector<double> y1, std::vector<double> y2, std::vector<double> y3, std::string xlabel, std::string savename);
};

#endif