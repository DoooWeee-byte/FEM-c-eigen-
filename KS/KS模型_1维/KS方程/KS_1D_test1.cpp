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


double s(double x)
{
    return sin(2*x);
}
double c1(double x){
	return -1;
}
double s_der(double x){
	return 2*cos(2*x);
}
double c2(double x){
	double chi = 1;
	return s_der(x)*chi;
}
double c3(double x){
	return 1;
}
double u0(double x){
	return 1;
}
int main() {
    //生成信息矩阵P,T,Pb,Tb***********************************************************************
    double a = 0,b = 1;
	int N = 2;
    Matrix P(1,N+1),T(2,N);//网格信息矩阵
    FEKernel::generate_PT_1d(a,b,N,P,T);


	//101代表使用线性元
    int basis_type_trial = 101;
    int basis_type_test = 101;
	int Nlb = 2;
	int Nlbtest=Nlb;
    int Nlbtrial=Nlb;//trial和test用的是一个元
    int N_test = N+1;
    int N_trial = N+1;

    Matrix Pb_trial(1,N+1),Tb_trial(2,N),Pb_test(1,N+1),Tb_test(2,N);
    FEKernel::generate_PbTb_1d(a,b,N,basis_type_trial,N_trial,Nlbtrial,Pb_trial,Tb_trial);
    FEKernel::generate_PbTb_1d(a,b,N,basis_type_test,N_test,Nlbtest,Pb_test,Tb_test);

    //组装A矩阵和V向量***********************************************************************************
    Matrix A(N_test,N_test);
    Matrix B(N_test,N_test);
    FEKernel::assembleA_1d(N,Nlbtest,Nlbtrial,basis_type_trial,basis_type_test,P,T,Tb_trial,Tb_test,A,c1,1,1);
    FEKernel::assembleA_1d(N,Nlbtest,Nlbtrial,basis_type_trial,basis_type_test,P,T,Tb_trial,Tb_test,A,c2,1,0);
    FEKernel::assembleA_1d(N,Nlbtest,Nlbtrial,basis_type_trial,basis_type_test,P,T,Tb_trial,Tb_test,B,c3,0,0);
    Matrix C(A.shape[0],A.shape[1]);


	//求解+处理边界条件************************************************************
	double Time = 1;
	//设置delta_t,theta
	double h = (b-a)/N;
	//方案1
	/* double delta_t = h; */
	/* double theta = 0.5; */
	//方案2
	double delta_t = 0.25*h*h;
	cout << delta_t << endl;
	double theta = 1;
	Vector x(N_trial);
	FEKernel::ODE_KS_1D_Solver(Time,delta_t,A,B,Pb_trial,x,u0);//求差分方程

		
    /* //计算误差 */
    /* Vector x1(V.size); */
    /* for(int i = 0;i != V.size;++i){ */
    /* 	x1[i] = analytic_solution(Pb[0][i],Pb[1][i],1); */
    /* } */
	/* //无穷范数 */
    /* double error1 = (x - x1).maxnorm(); */
    /* cout <<"无穷范数："<< error1 <<endl; */
	/* //2范数 */ 
    /* int der_x1=0; */
    /* int der_y1=0; */
    /* double error2 = FEKernel::compute_Hs_error_2D(N,P,T,Tb,basis_type_trial,der_x1,der_y1,x,analytic_solution,Time); */
    /* cout <<"L2范数："<< error2 <<endl; */
	/* //H1半范 */
	/* int der_x=1; */
    /* int der_y=0; */
    /* double error3 = FEKernel::compute_Hs_error_2D(N,P,T,Tb,basis_type_trial,der_x,der_y,x,analytic_solution_1_der,Time); */
    /* cout << "H1半范："<<error3 <<endl; */


    
}
