#include "Plot.h"
#include "constants.h"
namespace plt = matplotlibcpp;
void Plot::init(std::string name){
	std::map<std::string, std::string> rcparam;
    rcparam["axes.unicode_minus"]="False";
    rcparam["font.family"]="stixgeneral";
    rcparam["mathtext.fontset"]="stix";
	plt::rcparams(rcparam);
	std::map<std::string, std::string> legend_option;
	legend_option["fontsize"]="20";
    outputpath = name;
}

void Plot::plotFreeNandF(double maxEpsilon, int nPoints, double x2, double y, std::string savename){
    plt::clf();
	std::cout<<"plot N and f: "<<savename<<std::endl;
	std::vector<double> epsilon;
	std::vector<double> f;
	std::vector<double> N;
	double Nint = 0;
	double Uint = 0;
	double dEpsilon = maxEpsilon/nPoints;
	for (int i = 0; i < nPoints; ++i)
	{
		double epsilontmp = i*dEpsilon;
		epsilon.push_back(epsilontmp);
		f.push_back(1./(exp((epsilontmp-x2)/y)+1));
		N.push_back(f.back()*sqrt(epsilontmp));
		Nint+=N.back()*dEpsilon;
		Uint+=epsilontmp*N.back()*dEpsilon;
	}
	std::cout<<"N integral = "<<Nint<<std::endl;
	std::cout<<"U integral = "<<Uint<<std::endl;

    plt::named_semilogx("能态占有率", epsilon, f);
    //plt::save("f.pdf");
    //plt::clf();
    plt::named_semilogx("能级密度", epsilon, N);
    plt::legend();
    plt::save(savename);
}

void Plot::plotSingle(std::string name, std::vector<double> x, std::vector<double> y, std::string xlabel, std::string savename){
    plt::clf();
    if(x.size()!=y.size()) printf("Plot::plotSingle : ERROR: x length %lu not equal to y length %lu\n", x.size(), y.size());
    plt::named_plot(name, x, y);
    plt::xlabel(xlabel,legend_option);
    plt::legend(legend_option);
    std::cout<<"Plot::plotSingle : INFO Saving file to "<<outputpath+"/"+savename<<std::endl;
    plt::save(outputpath+"/"+savename);
}

void Plot::plotSingleLogX(std::string name, std::vector<double> x, std::vector<double> y, std::string xlabel, std::string savename){
    plt::clf();
    plt::named_semilogx(name, x, y);
    plt::xlabel(xlabel,legend_option);
    plt::legend(legend_option);
    std::cout<<"Plot::plotSingle : INFO Saving file to "<<outputpath+"/"+savename<<std::endl;
    plt::save(outputpath+"/"+savename);
}
void Plot::plotSingleLogY(std::string name, std::vector<double> x, std::vector<double> y, std::string xlabel, std::string savename){
    plt::clf();
    plt::named_semilogy(name, x, y);
    plt::xlabel(xlabel,legend_option);
    plt::legend(legend_option);
    std::cout<<"Plot::plotSingle : INFO Saving file to "<<outputpath+"/"+savename<<std::endl;
    plt::save(outputpath+"/"+savename);
}

void Plot::plotSingleLogXY(std::string name, std::vector<double> x, std::vector<double> y, std::string xlabel, std::string savename){
    plt::clf();
    plt::named_loglog(name, x, y);
    plt::xlabel(xlabel,legend_option);
    plt::legend(legend_option);
    std::cout<<"Plot::plotSingle : INFO Saving file to "<<outputpath+"/"+savename<<std::endl;
    plt::save(outputpath+"/"+savename);
}

void Plot::plotDouble(std::string name1, std::string name2, std::vector<double> x, std::vector<double> y1, std::vector<double> y2, std::string xlabel, std::string savename){
    plt::clf();
    plt::named_semilogx(name1, x, y1);
    plt::named_semilogx(name2, x, y2);
    plt::xlabel(xlabel,legend_option);
    plt::legend(legend_option);
    std::cout<<"Plot::plotDouble : INFO Saving file to "<<outputpath+"/"+savename<<std::endl;
    plt::save(outputpath+"/"+savename);
}

void Plot::plotTripple(std::string name1, std::string name2, std::string name3, std::vector<double> x, std::vector<double> y1, std::vector<double> y2, std::vector<double> y3, std::string xlabel, std::string savename){
    plt::clf();
    plt::named_plot(name1, x, y1);
    plt::named_plot(name2, x, y2);
    plt::named_plot(name3, x, y3);
    plt::xlabel(xlabel,legend_option);
    plt::legend(legend_option);
    std::cout<<"Plot::plotTripple : INFO Saving file to "<<outputpath+"/"+savename<<std::endl;
    plt::save(outputpath+"/"+savename);
}