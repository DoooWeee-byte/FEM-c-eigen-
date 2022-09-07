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

double c(double x,double y)
{
    return 1;
}
double f(double x,double y)
{
    double r = -y*(1-y)*(1-x-x*x/2)*exp(x+y)-x*(1-x/2)*(-3*y-y*y)*exp(x+y);
    return r;
}

double g(double x,double y)
{
    double r = 0;
    if(x==-1)
        r = -1.5*y*(1-y)*exp(-1+y);
    else if(x==1)
        r = 0.5*y*(1-y)*exp(1+y);
    else if(y==-1)
        r = -2.0*x*(1-x/2.0)*exp(x - 1.0);
    else if(y==1)
        r = 0;
    return r;
}
double p(double x,double y) {
    return -exp(x-1);
}

double r(double x)
{
    return 1;
}

double analytic_solution(double x,double y)
{
    return x*y*(1-x/2)*(1-y)*exp(x+y);
}
double analyticSolution_1_der(double x,double y)
{
    return (x-0.5*x*x)*(1-y-y*y)*exp(x+y);
}
int main() {
    //一致，三角剖分，线性元***********************************************************
    //生成信息矩阵P,T,Pb,Tb
    double right = 1,top = 1,left = -1,bottom = -1;
    int N1 = 32,N2 = 32;
    Matrix P(2,(N1+1)*(N2+1)),T(3,N1*N2*2);
    FEKernel::generate_PT_2D(P,T,N1,N2,right,top,bottom,left);
    //201使用线性元
    int basis_type_trial = 201;
    int basis_type_test = 201;
    int nlb = 3;
	int N = 2*N1*N2;
    SpMat A((N1+1)*(N2+1),(N1+1)*(N2+1)),A1(A.rows(),A.cols()),A2(A.rows(),A.cols());
    Vector V((N1+1)*(N2+1)),x((N1+1)*(N2+1));
    int Nlbtest = nlb;
    int Nlbtrial = nlb;

    Matrix Pb(2,(N1+1)*(N2+1)),Tb(nlb,N1*N2*2);
    FEKernel::generate_PbTb_2D(Pb,Tb,N1,N2,basis_type_trial,right,top,bottom,left);

    //生成边界矩阵，边界边，边界点
    Matrix boundaryedges(4,2*N1+2*N2);
    FEKernel::generate_boundaryedges_2D(boundaryedges,N1,N2);
    Matrix boundarynodes(2,2*(N1+N2));
    FEKernel::generate_boundarynodes_2D(boundarynodes,N1,N2);
    //组装A矩阵和V向量
    FEKernel::assembleA_2D(N,Nlbtest,Nlbtrial,basis_type_trial,basis_type_test,P,T,Tb,Tb,A1,c,1,0,1,0);
    FEKernel::assembleA_2D(N,Nlbtest,Nlbtrial,basis_type_trial,basis_type_test,P,T,Tb,Tb,A2,c,0,1,0,1);
    FEKernel::assembleV_2D(N,Nlbtest,basis_type_test,0,0,P,T,Pb,Tb,V,f);
	A = A1 + A2;
    



    //一致三角剖分，二次元**********************************************************************************
    //生成信息矩阵P,T,Pb,Tb
    /* double right = 1,top = 1,left = -1,bottom = -1; */
	/* int N1 = 3,N2 = 3; */
    /* int N1_basis = 2*N1; //基函数节点个数N1方向 */
    /* int N2_basis = 2*N2; //基函数节点个数N2方向 */
    /* int nbn = 2*(N1_basis+N2_basis); */
    /* Matrix P(2,(N1+1)*(N2+1)),T(3,N1*N2*2); */
    /* FEKernel::generate_PT_2D(P,T,N1,N2,right,top,bottom,left); */
    /* //202使用二次元 */
    /* int basis_type_trial = 202; */
    /* int basis_type_test = 202; */
	/* int nlb = 6; */
    /* int Nb = (N1+1)*(N2+1)+N1*(N2+1)+N2*(N1+1)+N1*N2; //总基函数节点个数 */
    /* int N = 2*N1*N2; */
    /* int Nlbtest = nlb; */
    /* int Nlbtrial = nlb; */
    /* Matrix Pb(2,Nb),Tb(6,N); */
    /* FEKernel::generate_PbTb_2D(Pb,Tb,N1,N2,basis_type_trial,right,top,bottom,left); */
    /* //生成边界矩阵，边界边，边界点 */
    /* Matrix boundaryedges(4,2*N1+2*N2); */
    /* FEKernel::generate_boundaryedges_2D(boundaryedges,N1,N2); */
    /* Matrix boundarynodes(2,nbn); */
    /* FEKernel::generate_boundarynodes_2D(boundarynodes,N1_basis,N2_basis); */
    /* //组装A矩阵和V向量 */
    /* SpMat A(Nb,Nb); */
    /* Vector V(Nb),x(Nb); */
    /* FEKernel::assembleA_2D(N,Nlbtest,Nlbtrial,basis_type_trial,basis_type_test,P,T,Tb,Tb,A,c,1,0,1,0); */
    /* FEKernel::assembleA_2D(N,Nlbtest,Nlbtrial,basis_type_trial,basis_type_test,P,T,Tb,Tb,A,c,0,1,0,1); */
    /* FEKernel::assembleV_2D(N,Nlbtest,basis_type_test,0,0,P,T,Pb,Tb,V,f); */

    //处理边界条件****************************************************
    FEKernel::treat_boundary_Dirichlet_2d(boundarynodes,A,Pb,V,g);
	/* cout << "A with Eigen" << A <<endl; */
	/* cout << "b with Eigen" <<endl<< V <<endl; */

   /* /1* //求解线性方程组 *1/ */
    Eigen:: SparseLU<SpMat> solver;
	solver.analyzePattern(A);
    solver.compute(A);
    if(solver.info() != 0) {
        cout << "decomposition failed" <<endl;
    }
	/* cout << A.rows() << A.cols() <<endl; */
	/* cout << "AA" << A <<endl; */
	/* cout << "b" << V << endl; */
    x = solver.solve(V);
	/* cout << "x with Eigen" << endl<<x << endl; */

    /* /1* //计算误差 *1/ */
    /* /1* //组装解析解向量 *1/ */
    Vector x1(V.size());
    for(int i = 0; i != V.size(); ++i) {
        x1[i] = analytic_solution(Pb(0,i),Pb(1,i));
    }
	/* cout << "analytic_solution with Eigen" << endl <<  x1 <<endl; */
	/* /1* //计算无穷范数 *1/ */
    /* double error = (x - x1).norm(); */
    /* cout << error <<endl; */

    /* //2范数:der_x = der_y = 0, */
    int der_x=0;
    int der_y=0;
    double error = FEKernel::compute_Hs_error_2D(N,P,T,Tb,basis_type_trial,der_x,der_y,x,analytic_solution);
    cout << "error"<<error <<endl;

	/* //H1半范 */
	/* int der_x = 0,der_y=1; */
	/* double error = FEKernel::compute_Hs_error_2D(N,P,T,Tb,basis_type_trial,der_x,der_y,x,analyticSolution_1_der);//传入解析解的一阶偏导数 */
    /* cout << error <<endl; */























}
