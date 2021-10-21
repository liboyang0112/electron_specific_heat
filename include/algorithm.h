#include <vector>
#include <iostream>
double interpolate(double x, std::vector<double> &v);
double interpolate(double x, double y, std::vector<std::vector<double>> &v);
double order_in(std::vector<double> &v, double x);