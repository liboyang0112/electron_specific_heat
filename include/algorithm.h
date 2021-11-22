#include <vector>
#include <iostream>
#include <string>
double interpolate(double x, std::vector<double> &v);
double interpolate(double x, double y, std::vector<std::vector<double>> &v);
void interpolate(std::vector<double> &v1, std::vector<double> &v2, std::vector<double> &v3, double x);
void interpolate(std::vector<std::vector<double>*> &v1, std::vector<double> &v2, double x);
double order_in(std::vector<double> &v, double x);
void readLine(std::string file, std::string linename, std::vector<double> &v);
std::vector<double>* readLine(std::string file, std::string linename);