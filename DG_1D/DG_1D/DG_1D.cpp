#include <iostream>
#include <cmath>
#include "algebra/Algebra_kernel.h"
#include "FEMethod/FEMethod_kernel.h"
#include "Constants.h"
using namespace std;
typedef typename PDWSC::FEMethod_kernel<double, int> FEKernel;
typedef typename Eigen::SparseMatrix<double> SpMat;
typedef typename Eigen::MatrixXd Matrix;
typedef typename Eigen::VectorXd Vector;
double one(double x)
{
    return 1;
}
double a = 0;
double f(double x)
{
    double r = cos(x);
    return r;
}
double g_left(double x)
{
    return 0;
}
double r(double x)
{
    return 1;
}
double analyticSolution(double x)
{
    return sin(x);
}
double analyticSolution_1_der(double x)
{
    return cos(x);
}
int main() {

	vector<int> a{4,8,16,32,64,128};
	for(int N : a){
		//生成信息矩阵P,T,Pb,Tb***********************************************************************
		double a = 0,b = 1;
		Matrix P(1,N+1),T(2,N);//网格信息矩阵
		FEKernel::generate_PT_1d(a,b,N,P,T);

		//101代表使用线性元
		int basis_type_trial = 101;
		int basis_type_test = 101;
		int Nlb = 2; 
		Matrix Pb_trial(1,N+1),Tb_trial(2,N),Pb_test(1,N+1),Tb_test(2,N);
		int Nb = N + 1;

		//换成102，可使用二次元但一些参数的准备也需要更改
		/*int basis_type_trial = 102;
		int basis_type_test = 102;
		int Nlb = 3;
		Matrix Pb_trial(1,2*N+1),Tb_trial(3,N),Pb_test(1,2*N+1),Tb_test(3,N);*/


		int Nlbtest=Nlb;
		int Nlbtrial=Nlb;//trial和test用的是一个元
		int N_test = N+1;
		int N_trial = N+1;
		FEKernel::generate_PbTb_1d(a,b,N,basis_type_trial,N_trial,Nlbtrial,Pb_trial,Tb_trial);
		FEKernel::generate_PbTb_1d(a,b,N,basis_type_test,N_test,Nlbtest,Pb_test,Tb_test);
		
		Matrix solution(2,N);
		for(int i = 0; i < N;++i){
			if(i == 0){
				Matrix vertices(1,2);
				vertices(0,0) = P(0, 0);
				vertices(0,1) = P(0, 1);
				double x = vertices(0,1);
				double uh_pre = 0;
				uh_pre *= 0;
				Matrix A(Nlb,Nlb), B(Nlb,Nlb);
			
				FEKernel::assembleA_DG_1d(Nlb, Nlb, basis_type_trial, basis_type_test, vertices, A,one,0,1);
				FEKernel::assembleA_DG_bd_1d(x, Nlbtest, Nlbtrial, basis_type_trial, basis_type_test, vertices, B, one , 0, 0);
				Vector b1(Nlb), b2(Nlb);
				FEKernel::assembleV_1d_DG(Nlb, basis_type_test, 0, vertices, b1, f);
				FEKernel::assembleV_1d_DG_BD(x,uh_pre,Nlb,basis_type_trial,basis_type_test,0,0,vertices,vertices,b2);
				A = -1*A + B;
				b1 = b1 + a*b2;
				
				Vector slo = A.colPivHouseholderQr().solve(b1);
				/* slo(0) = slo(0); */
				/* slo(1) = slo(0) + slo(1); */
				solution.middleCols(0,1) = slo;
				continue;
			}
			
				Matrix vertices(1,2),vertices_pre(1,2);
				vertices(0,0) = P(0, i);
				vertices(0,1) = P(0, i + 1);
				vertices_pre(0,0) = P(0, i-1);
				vertices_pre(0,1) = P(0, i);
				double x = vertices(0,1);
				double x1 = vertices(0,0);
				double uh_pre = solution(1,i-1);
				Matrix A(Nlb,Nlb), B(Nlb,Nlb);
				FEKernel::assembleA_DG_1d(Nlb, Nlb, basis_type_trial, basis_type_test, vertices, A,one,0,1);
				FEKernel::assembleA_DG_bd_1d(x, Nlbtest, Nlbtrial, basis_type_trial, basis_type_test, vertices, B, one , 0, 0);
				Vector b1(Nlb), b2(Nlb);
				FEKernel::assembleV_1d_DG(Nlb, basis_type_test, 0, vertices, b1, f);
				FEKernel::assembleV_1d_DG_BD(x1,uh_pre, Nlb, basis_type_trial, basis_type_test,0,0,vertices,vertices_pre,b2);
				A = -1*A + B;
				b1 = b1 + b2;
				Vector slo = A.colPivHouseholderQr().solve(b1);
				/* slo(0) = slo(0); */
				/* slo(1) = slo(0) + slo(1); */
				solution.middleCols(i,1) = slo;
		}


		//计算误差****************************************************************************
		Vector x(Nb);
		x(0) = a;
		for(int i = 1;i < Nb -1 ;i++){
			x(i) = solution(1,i-1) + solution(0,i);
		}
		x(Nb-1) = solution(1,N-1);
		x.tail(Nb-1) = solution.bottomRows(1).transpose();
		Vector x1(Nb);
		for(int i = 0; i < Nb;i++){
			x1(i) = analyticSolution(P(0,i));
		}
		//数值解与解析节的误差向量
		Vector error = x - x1;
		cout << (x - x1).norm() << endl;
		 //计算无穷范
		 /* double error1 = FEKernel::compute_Max_error_1D(x, Pb_test, analyticSolution); */
		 /* cout << error1 <<endl; */
		/* //取0和analyticSolution时可以计算2范数， */
		/* double error2 = FEKernel::compute_Hs_error(N,P,T,Tb_test,basis_type_test,0,x,analyticSolution); */
		/* cout << error2 <<endl; */

		/* //取1和analyticSolution的导数时可以计算H1半范 */
		/* double error3 = FEKernel::compute_Hs_error(N,P,T,Tb_test,basis_type_test,1,x,analyticSolution_1_der); */
		/* cout << error3 <<endl; */
	}

}
