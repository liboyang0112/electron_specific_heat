#include "algorithm.h"
	double interpolate(double x, std::vector<double> &v){
		int integer = (int)x;
		double decimal = x-integer;
		if(integer>=v.size()) return v.back();
		if(decimal<1e-5) return v.at(integer);
		return v.at(integer)*(1-decimal)+v.at(integer+1)*decimal;
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
		if(v.at(0)>x) return x/v.at(0);
		if(v.at(up)<x) return up;
		while(up-low>1){
			uint test = (low+up)/2;
			if(v.at(test)>x) {
				up = test;
			}else low = test;
		}
		return 1+low+(x-v.at(low))/(v.at(up)-v.at(low));
	}