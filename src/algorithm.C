#include "algorithm.h"
#include <fstream>
#include <sstream>
	double interpolate(double x, std::vector<double> &v){
		int integer = (int)x;
		double decimal = x-integer;
		if(integer>=v.size()) return v.back();
		if(decimal<1e-5) return v.at(integer);
		return v.at(integer)*(1-decimal)+v.at(integer+1)*decimal;
	}
	void interpolate(std::vector<double> &v1, std::vector<double> &v2, std::vector<double> &v3, double x){
		for (int i = 0; i < v1.size(); ++i)
		{
			v3.push_back(v1[i]*(1-x)+v2[i]*x);
		}
	}
	void interpolate(std::vector<std::vector<double>*> &v1, std::vector<double> &v2, double x){
		int integer = (int)x;
		double decimal = x-integer;
		if(decimal<1e-5) {
			v2=*(v1[integer]);
			return;
		}
		for (int i = 0; i < v1.at(0)->size(); ++i)
		{
			v2.push_back(v1[integer]->at(i)*(1-decimal)+v1[integer+1]->at(i)*decimal);
		}
	}
	double interpolate(double x, double y, std::vector<std::vector<double>> &v){
		int integerx = (int)x;
		double decimalx = x-integerx;
		if(decimalx<1e-5) return interpolate(y,v[integerx]);
		int integery = (int)y;
		double decimaly = y-integery;
		if(decimaly<1e-5) {
			return v.at(integerx).at(integery)*(1-decimalx)+v.at(integerx+1).at(integery)*decimalx;
		}
		return v.at(integerx).at(integery)*(1-decimalx-decimaly)+
		v.at(integerx+1).at(integery)*decimalx+
		v.at(integerx).at(integery+1)*decimaly;
	}
	double order_in(std::vector<double> &v, double x){
		uint up = v.size()-1;
		uint low = 0;
		if(v.at(0)>x) return 0;
		if(v.at(up)<x) return up;
		while(up-low>1){
			uint test = (low+up)/2;
			if(v.at(test)>x) {
				up = test;
			}else low = test;
		}
		return low+(x-v.at(low))/(v.at(up)-v.at(low));
	}

	void readLine(std::string filename, std::string linename, std::vector<double> &v){
		std::ifstream file(filename);
		if (!file)
		{
			std::cerr << "File " << filename << " could not be opened" << std::endl;
			exit(0);
		}
		std::string inputline;
		std::string tmpname;
		double tmp;
		while(getline(file,inputline)){
			if(inputline.find(linename)!=std::string::npos){
				std::stringstream ss(inputline);
				ss>>tmpname;
				while(ss){
					ss>>tmp;
					v.push_back(tmp);
				}
				break;
			}
		}

	}

	std::vector<double>* readLine(std::string filename, std::string linename){
		std::vector<double>* v = new std::vector<double>();
		readLine(filename, linename, *v);
		return v;
	}