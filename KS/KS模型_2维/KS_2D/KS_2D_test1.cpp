#include <iostream>
#include <cmath>
#include "algebra/Algebra_kernel.h"
#include "FEMethod/FEMethod_kernel.h"
#include "Constants.h"
#include <Eigen/Eigen>
typedef typename PDWSC::FEMethod_kernel<double, int> FEKernel;
typedef typename Eigen::SparseMatrix<double> SpMat;
typedef typename Eigen::MatrixXd Matrix;
typedef typename Eigen::VectorXd Vector;


using namespace std;
namespace plt = matplotlibcpp;
double a = 14;
double D = -1/a;
double c(double x,double y)
{
    return 1;
}
double f(double x,double y)
{
    double r = -2*a*exp(a*(2*x-1))*(3*exp(a*(2*x-1))-1)/pow(1+exp(a*(2*x-1)),3);
    /* double r = -2*a*exp(a*(2*x-1))*(1+2*a)/pow(1+exp(a*(2*x+1)),2) + 8*a*a*exp(2*a*(2*x-4))/pow(1+exp(a*(2*x-1)),3); */
    return r;
}


double analytic_solution(double x,double y)
{
    return 1/(1+exp(a*(2*x-1)));
}
double analyticSolution_1_der(double x,double y)
{
    return 1;
}



int main() {
    //一致，三角剖分，线性元***************************************************************
    //生成信息矩阵P,T,Pb,Tb
    double right = 1,top = 1,left = 0,bottom = 0;
    int N1 = 8,N2 = N1;
	double h = (right-left)/N1;
    Matrix P(2,(N1+1)*(N2+1)),T(3,N1*N2*2);
    FEKernel::generate_PT_2D(P,T,N1,N2,right,top,bottom,left);
    //201使用线性元
    int basis_type_trial = 201;
    int basis_type_test = 201;
    int nlb = 3;
    Matrix Pb(2,(N1+1)*(N2+1)),Tb(nlb,N1*N2*2);
    FEKernel::generate_PbTb_2D(Pb,Tb,N1,N2,basis_type_trial,right,top,bottom,left);
    //组装A矩阵和V向量
    int N = 2*N1*N2;
    int Arows = (N1+1)*(N2+1),Acols = Arows;
    const int Arowsd = Arows - 1;
    SpMat L(Arows,Acols);
    SpMat L1(Arows,Acols);
    SpMat L2(Arows,Acols);
    SpMat F(Arows,Acols);
    SpMat A(Arows,Acols);
    Vector b(Arows);
    Vector x(Arows);
    int Nlbtest = 3;
    int Nlbtrial = 3;
    FEKernel::assembleA_2D(N,Nlbtest,Nlbtrial,basis_type_trial,basis_type_test,P,T,Tb,Tb,L1,c,1,0,1,0);
    FEKernel::assembleA_2D(N,Nlbtest,Nlbtrial,basis_type_trial,basis_type_test,P,T,Tb,Tb,L2,c,0,1,0,1);
    L = L1 + L2;
    FEKernel::assembleA_2D(N,Nlbtest,Nlbtrial,basis_type_trial,basis_type_test,P,T,Tb,Tb,F,c,1,0,0,0);
    FEKernel::assembleV_2D(N,Nlbtest,basis_type_test,0,0,P,T,Pb,Tb,b,f);
    A = -D*L + F;
    /* cout << 'A' << endl << A << 'b' << endl << b << endl; */

    for(int i = 0; i < A.cols(); ++i) {
        if(i==0)
            A.coeffRef(0,i) = 1;
        else
            A.coeffRef(0,i) = 0;
    }
    b(0) = 1;

    /* //求解线性方程组 */
    Eigen::BiCGSTAB<SpMat> solver;
    solver.compute(A);
    if(solver.info() != 0) {
        cout << "decomposition failed" <<endl;
    }
    x = solver.solve(b);
    /* cout  << x <<endl; */


    //计算误差
    //组装解析解向量
    Vector x1(Arows);
    for(int i = 0; i != x1.size(); ++i) {
        x1[i] = analytic_solution(Pb(0,i),Pb(1,i));
    }

    /* cout << "x1" << x1 <<endl; */
    //无穷范数
    double error1 = (x1 - x).norm();
    cout << error1 <<endl;
    //2范数
    /* int der_x1=0; */
    /* int der_y1=0; */
    /* double error2 = FEKernel::compute_Hs_error_2D(N,P,T,Tb,basis_type_trial,der_x1,der_y1,x,analytic_solution); */
    /* cout << error2 <<endl; */
    //H1半范
    /* int der_x2 = 0,der_y2=1; */
    /* double error3 = FEKernel::compute_Hs_error_2D(N,P,T,Tb,basis_type_trial,der_x2,der_y2,x,analyticSolution_1_der);//传入解析解的一阶偏导数 */
    /* cout << error3 <<endl; */



    /************************************画图*/

    vector<vector<double>> X,Y,Z,ZN;
    for(int i = 0; i < N1+1; ++i) {
        vector<double> xi,yi,zi,zn;
        for(int j = 0; j < N1 + 1; ++j) {
            int k = i*(N1+1)+j;
            xi.push_back(Pb(0,k));
            yi.push_back(Pb(1,k));
            zi.push_back(analytic_solution(Pb(0,k),Pb(1,k)));
			zn.push_back(x(k));
        }
        X.push_back(xi);
        Y.push_back(yi);
        Z.push_back(zi);
		ZN.push_back(zn);
    }
    plt::plot_surface(X,Y,Z);
    plt::title("analytic_solution");
    plt::xlabel("x");
    plt::ylabel("y");
    plt::show();
	string s1 = "numerical_solution,h=";
	string s2 = ",2norm=";
	stringstream ss;
	ss <<  s1 << h << s2 << error1;
	string s = ss.str();
	plt::plot_surface(X,Y,ZN);
    plt::title(s);
    plt::xlabel("x");
    plt::ylabel("y");
    plt::show();


}
