#include <iostream>
#include <cmath>
#include <fstream>
#include "Constants.h"
#include "thirdparty/matplotlibcpp.h"
using namespace std;
namespace plt = matplotlibcpp;

int main(){
	vector<double> x_plt(4),y_plt(4),z_plt(4);

	x_plt.at(0) = 9*9;
	x_plt.at(1) = 17*17;
	x_plt.at(2) = 33*33;
	x_plt.at(3) = 65*65;
	y_plt.at(0) = 17.734;
	y_plt.at(1) = 7.3919;
	y_plt.at(2) = 2.0171;
	y_plt.at(3) = 0.524;
	z_plt.at(0) = 3.1276;
	z_plt.at(1) = 0.99885;
	z_plt.at(2) = 0.22964;
	z_plt.at(3) = 0.0582;
	string s = "r-o";
	string s1 = "b-o";
	plt::loglog(x_plt,y_plt,s);
	plt::loglog(x_plt,z_plt,s1);
	plt::show();
}
