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
//给出泊松方程的各个参数*******************************************************************
double a = 100;
double D = -1/a;
double c(double x,double y)
{
    return 1;
}
double f(double x,double y)
{
    double r = -2*a*exp(a*(2*x-1))*(3*exp(a*(2*x-1))-1)/pow(1+exp(a*(2*x-1)),3);
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
    /* //一致，三角剖分，线性元*************************************************************** */
    /* //生成信息矩阵P,T,Pb,Tb */
    double right = 1,top = 1,left = 0,bottom = 0;
    int N1 = 2,N2 = 2;
    Matrix P(2,(N1+1)*(N2+1)),T(3,N1*N2*2);
    FEKernel::generate_PT_2D(P,T,N1,N2,right,top,bottom,left);
    /* //201使用线性元 */
    int basis_type_trial = 201;
    int basis_type_test = 201;
    int nlb = 3;
    Matrix Pb(2,(N1+1)*(N2+1)),Tb(nlb,N1*N2*2);
    FEKernel::generate_PbTb_2D(Pb,Tb,N1,N2,basis_type_trial,right,top,bottom,left);


    /* //组装A矩阵和V向量 */
    int N = 2*N1*N2;
    Matrix L((N1+1)*(N2+1),(N1+1)*(N2+1));
    Matrix F((N1+1)*(N2+1),(N1+1)*(N2+1));
    Matrix A((N1+1)*(N2+1),(N1+1)*(N2+1));
    Vector V((N1+1)*(N2+1)),x(V.size);
    int Nlbtest = 3;
    int Nlbtrial = 3;
    FEKernel::assembleA_2D(N,Nlbtest,Nlbtrial,basis_type_trial,basis_type_test,P,T,Tb,Tb,L,c,1,0,1,0);
    FEKernel::assembleA_2D(N,Nlbtest,Nlbtrial,basis_type_trial,basis_type_test,P,T,Tb,Tb,L,c,0,1,0,1);
    FEKernel::assembleA_2D(N,Nlbtest,Nlbtrial,basis_type_trial,basis_type_test,P,T,Tb,Tb,F,c,1,0,0,0);
    FEKernel::assembleV_2D(N,Nlbtest,basis_type_test,0,0,P,T,Pb,Tb,V,f);
    A = -D*L + F;
	cout << A <<endl;


    /* //一致，三角剖分，线性元*************************************************************** */
    /* //生成信息矩阵P,T,Pb,Tb */
    /* /1* double right = 1,top = 1,left = -1,bottom = -1; *1/ */
    /* /1* int N1 = 64,N2 = 64; *1/ */
    /* /1* Matrix P(2,(N1+1)*(N2+1)),T(3,N1*N2*2); *1/ */
    /* /1* FEKernel::generate_PT_2D(P,T,N1,N2,right,top,bottom,left); *1/ */
    /* //201使用线性元 */
    /* /1* int basis_type_trial = 201; *1/ */
    /* /1* int basis_type_test = 201; *1/ */
    /* /1* int nlb = 3; *1/ */
    /* /1* Matrix Pb(2,(N1+1)*(N2+1)),Tb(nlb,N1*N2*2); *1/ */
    /* /1* FEKernel::generate_PbTb_2D(Pb,Tb,N1,N2,basis_type_trial,right,top,bottom,left); *1/ */

    /* //生成边界矩阵，边界边，边界点 */
    /* /1* Matrix boundaryedges(4,2*N1+2*N2); *1/ */
    /* /1* FEKernel::generate_boundaryedges_DN_2D(boundaryedges,N1,N2); *1/ */
    /* /1* Matrix boundarynodes(2,2*(N1+N2)); *1/ */
    /* /1* FEKernel::generate_boundarynodes_DN_2D(boundarynodes,N1,N2); *1/ */

    /* //组装A矩阵和V向量 */
    /* /1* int N = 2*N1*N2; *1/ */
    /* /1* Matrix A((N1+1)*(N2+1),(N1+1)*(N2+1)); *1/ */
    /* /1* Vector V((N1+1)*(N2+1)); *1/ */
    /* /1* int Nlbtest = 3; *1/ */
    /* /1* int Nlbtrial = 3; *1/ */
    /* /1* FEKernel::assembleA_2D(N,Nlbtest,Nlbtrial,basis_type_trial,basis_type_test,P,T,Tb,Tb,A,c,1,0,1,0); *1/ */
    /* /1* FEKernel::assembleA_2D(N,Nlbtest,Nlbtrial,basis_type_trial,basis_type_test,P,T,Tb,Tb,A,c,0,1,0,1); *1/ */
    /* /1* FEKernel::assembleV_2D(N,Nlbtest,basis_type_test,0,0,P,T,Pb,Tb,V,f); *1/ */

    
    /* //一致三角剖分，二次元********************************************************************* */
    /* //生成信息矩阵P,T,Pb,Tb */
    /* double right = 1,top = 1,left = -1,bottom = -1; */
    /* int N1 = 64,N2 = 64; */
    /* int N1_basis = 2*N1; //基函数节点个数N1方向 */
    /* int N2_basis = 2*N2; //基函数节点个数N2方向 */
    /* int nbn = 2*(N1_basis+N2_basis);//边界上基函数节点的个数 */
    /* Matrix P(2,(N1+1)*(N2+1)),T(3,N1*N2*2); */
    /* FEKernel::generate_PT_2D(P,T,N1,N2,right,top,bottom,left); */
    /* //使用二次元 */
    /* int basis_type_trial = 202; */
    /* int basis_type_test = 202; */
	/* int nlb = 6; */
    /* int Nb = (N1+1)*(N2+1)+N1*(N2+1)+N2*(N1+1)+N1*N2; //总基函数节点个数 */
    /* int N = 2*N1*N2; */
    /* int Nlbtest = 6; */
    /* int Nlbtrial = 6; */
    /* Matrix Pb(2,Nb),Tb(6,N); */
    /* FEKernel::generate_PbTb_2D(Pb,Tb,N1,N2,basis_type_trial,right,top,bottom,left); */
    /* //生成边界矩阵，边界边，边界点 */
    /* Matrix boundaryedges(4,2*N1+2*N2); */
    /* FEKernel::generate_boundaryedges_DN_2D(boundaryedges,N1,N2); */
    /* Matrix boundarynodes(2,nbn); */
    /* FEKernel::generate_boundarynodes_DN_2D(boundarynodes,N1_basis,N2_basis); */
    /* //组装A矩阵和V向量 */ 
    /* Matrix A(Nb,Nb); */
    /* Vector V(Nb); */
    /* FEKernel::assembleA_2D(N,Nlbtest,Nlbtrial,basis_type_trial,basis_type_test,P,T,Tb,Tb,A,c,1,0,1,0); */
    /* FEKernel::assembleA_2D(N,Nlbtest,Nlbtrial,basis_type_trial,basis_type_test,P,T,Tb,Tb,A,c,0,1,0,1); */
    /* FEKernel::assembleV_2D(N,Nlbtest,basis_type_test,0,0,P,T,Pb,Tb,V,f); */



	/* //边界处理,处理过程与元的选取没关系******************************************************************************** */
	/* //先处理Nuemann条件再处理Dirichlet条件 */
    /* Vector b1(V.size); */
    /* FEKernel::treat_boundary_Neumann_2d(b1,boundaryedges,Pb,Tb,P,T,Neufunc,nlb,basis_type_test,0,0,4); */
    /* V = V + b1; */
    /* FEKernel::treat_boundary_Dirichlet_2d(boundarynodes,A,Pb,V,g); */
	/* //求解线性方程组 */
    /* Vector x(V.size); */
    /*  CSRMatrix R(A); */
    /* R.SOR(V,x,1.5,100000); */
    /* //计算误差 */
	/* //组装解析节向量 */
    /* Vector x1(V.size); */
    /* for(int i = 0; i != V.size; ++i) { */
    /*     x1[i] = analytic_solution(Pb[0][i],Pb[1][i]); */
    /* } */
	/* //无穷范数 */
    /* /1* double error1 = (x - x1).maxnorm(); *1/ */
    /* /1* cout << error1 <<endl; *1/ */
	/* //2范数 */        
	/* int der_x1=0; */
    /* int der_y1=0; */
    /* double error2 = FEKernel::compute_Hs_error_2D(N,P,T,Tb,basis_type_trial,der_x1,der_y1,x,analytic_solution); */
    /* cout << error2 <<endl; */
	/* //H1半范 */
	/* /1* int der_x2 = 0,der_y2=1; *1/ */
	/* /1* double error3 = FEKernel::compute_Hs_error_2D(N,P,T,Tb,basis_type_trial,der_x2,der_y2,x,analyticSolution_1_der);//传入解析解的一阶偏导数 *1/ */
    /* /1* cout << error3 <<endl; *1/ */


















}
