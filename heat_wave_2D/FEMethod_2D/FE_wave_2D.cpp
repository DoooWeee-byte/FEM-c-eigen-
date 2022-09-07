#include <iostream>
#include <cmath>
#include "algebra/Algebra_kernel.h"
#include "FEMethod/FEMethod_kernel.h"
#include "Constants.h"
using namespace std;
typedef typename WHYSC::Algebra_kernel<double, int> Kernel;
typedef typename PDWSC::FEMethod_kernel<double, int> FEKernel;
typedef typename Kernel::Matrix Matrix;
typedef typename Kernel::Vector Vector;
typedef typename Kernel::CSRMatrix CSRMatrix;
//波方程，热方程中u_t换乘u_tt
double c(double x,double y)
{
    return 2;
}
double cm(double x,double y){
	return 1;
}
double f(double x,double y,double t)
{
    double r = -3*exp(x+y+t);
    return r;
}

double g(double x,double y,double t)
{
    double r = 0;
    if(x==0)
        r = exp(y+t);
    else if(x==2)
        r = exp(2+y+t);
    else if(y==0)
        r = exp(x+t);
    else if(y==1)
        r = exp(x+1+t);
    return r;
}
double p(double x,double y){
	return -exp(x-1);
}
//初始条件
double u0(double x,double y){
	return exp(x+y);
}
//u00是时间t=0时，u对时间的导数值
double u00(double x,double y){
	return exp(x+y);
}
double r(double x)
{
    return 1;
}

double Neufunc(double x,double y){
	return c(x,y)*p(x,y);
}
double analytic_solution(double x,double y,double t)
{
    return exp(x+y+t);
}
double analytic_solution_1_der(double x,double y,double t)
{
    return exp(x+y+t);
}
int main() {
    //一致，三角剖分，线性元(更换元与泊松中类似)***********************************************************
    //生成信息矩阵P,T,Pb,Tb
    double right = 2,top = 1,left = 0,bottom = 0;
	int N1 = 64,N2 = 32;
    Matrix P(2,(N1+1)*(N2+1)),T(3,N1*N2*2);
    FEKernel::generate_PT_2D(P,T,N1,N2,right,top,bottom,left);
    //201使用线性元
    int basis_type_trial = 201;
    int basis_type_test = 201;
	int nlb = 3;
    Matrix Pb(2,(N1+1)*(N2+1)),Tb(nlb,N1*N2*2);
    FEKernel::generate_PbTb_2D(Pb,Tb,N1,N2,basis_type_trial,right,top,bottom,left);

    //生成边界矩阵，边界边，边界点
    Matrix boundaryedges(4,2*N1+2*N2);
    FEKernel::generate_boundaryedges_2D(boundaryedges,N1,N2);
    Matrix boundarynodes(2,2*(N1+N2));
    FEKernel::generate_boundarynodes_2D(boundarynodes,N1,N2);
    //组装A矩阵和V向量
    int N = 2*N1*N2;
    Matrix A((N1+1)*(N2+1),(N1+1)*(N2+1));
	Matrix M((N1+1)*(N2+1),(N1+1)*(N2+1));
    Vector V((N1+1)*(N2+1));
    int Nlbtest = 3;
    int Nlbtrial = 3;
    FEKernel::assembleA_2D(N,Nlbtest,Nlbtrial,basis_type_trial,basis_type_test,P,T,Tb,Tb,A,c,1,0,1,0);
    FEKernel::assembleA_2D(N,Nlbtest,Nlbtrial,basis_type_trial,basis_type_test,P,T,Tb,Tb,A,c,0,1,0,1);
    FEKernel::assembleA_2D(N,Nlbtest,Nlbtrial,basis_type_trial,basis_type_test,P,T,Tb,Tb,M,cm,0,0,0,0);

	//求解+处理边界条件*********************************************
	double Time = 1;
	//设置delta_t,theta
	double h = (right-left)/N1;
	//方案1
	/* double delta_t = h; */
	/* double theta = 0.5; */
	//方案2
	double delta_t = 4*h*h;
	double theta = 1;

	Vector x(V.size);
	int p = 0,q = 0;
	FEKernel::ODE_WAVE_Solver(Time,delta_t,theta,A,x,N,nlb,basis_type_test,p,q,f,g,u0,u00,M,P,T,Pb,Tb,boundarynodes);

    /* //计算误差 */
    Vector x1(V.size);
    for(int i = 0;i != V.size;++i){
    	x1[i] = analytic_solution(Pb[0][i],Pb[1][i],1);
    }
	/* //无穷范数 */
    double error1 = (x - x1).maxnorm();
    cout <<"无穷范数："<< error1 <<endl;

	/* //2范数 */ 
    int der_x1=0;
    int der_y1=0;
    double error2 = FEKernel::compute_Hs_error_2D(N,P,T,Tb,basis_type_trial,der_x1,der_y1,x,analytic_solution,Time);
    cout << error2 <<endl;

	/* //H1半范 */
	int der_x2=1;
    int der_y2=0;
    double error3 = FEKernel::compute_Hs_error_2D(N,P,T,Tb,basis_type_trial,der_x2,der_y2,x,analytic_solution_1_der,Time);
    cout << error3 <<endl;









}
