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
double pi =  WHYSC::Constants::pi;
double func_lambda(double x,double y) {
    double lambda = 1;
    return lambda;
}
double func_mu(double x,double y) {
    double mu = 2;
    return mu;
}
double g1(double x,double y)
{
    double r = 0;
    return r;
}
double g2(double x,double y) {
    return 0;
}

double f1(double x,double y) {
    double lambda = 1;
    double mu = 2;
    return -(lambda + 2*mu)*(-pi*pi*sin(pi*x)*sin(pi*y))-(lambda+mu)*((2*x-1)*(2*y-1))-mu*(-pi*pi*sin(pi*x)*sin(pi*y));
}
double f2(double x,double y) {
    double lambda = 1;
    double mu = 2;
    return -(lambda + 2*mu)*(2*x*(x-1))-(lambda + mu)*(pi*pi*cos(pi*x)*cos(pi*y))-mu*(2*y*(y-1));
}
double u1(double x,double y)
{
    return sin(pi*x)*sin(pi*y);
}
double u2(double x,double y)
{
    return x*(x-1)*y*(y-1);
}

double u1_1_der(double x,double y) {
    return pi*cos(pi*y)*sin(pi*x);
}

double u2_1_der(double x,double y) {
    return x*(x-1)*(2*y-1);
}

int main() {
    //一致，三角剖分，线性元***********************************************************
    //生成信息矩阵P,T,Pb,Tb
    /* double right = 1,top = 1,left = 0,bottom = 0; */
    /* int N1 = 64,N2 = N1; */
    /* Matrix P(2,(N1+1)*(N2+1)),T(3,N1*N2*2); */
    /* FEKernel::generate_PT_2D(P,T,N1,N2,right,top,bottom,left); */
    //201使用线性元
    /* int basis_type_trial = 201; */
    /* int basis_type_test = 201; */
    /* int nlb = 3; */
    /* int N = 2*N1*N2; */
    /* SpMat A1((N1+1)*(N2+1),(N1+1)*(N2+1)),A2(A1.rows(),A1.cols()),A3(A1.rows(),A1.cols()),A4(A1.rows(),A1.cols()),A5(A1.rows(),A1.cols()),A6(A1.rows(),A1.cols()),A7(A1.rows(),A1.cols()),A8(A1.rows(),A1.cols()); */
    /* Matrix A(2*A1.rows(),2*A1.cols()); */
    /* Vector V1((N1+1)*(N2+1)),V2(V1.size()),V(2*V1.size()),x(V.size()),x1(V1.size()),x2(V1.size()); */
    /* int Nlbtest = nlb; */
    /* int Nlbtrial = nlb; */

    /* Matrix Pb(2,(N1+1)*(N2+1)),Tb(nlb,N1*N2*2); */
    /* FEKernel::generate_PbTb_2D(Pb,Tb,N1,N2,basis_type_trial,right,top,bottom,left); */

    //生成边界矩阵，边界边，边界点
    /* Matrix boundaryedges(4,2*N1+2*N2); */
    /* FEKernel::generate_boundaryedges_2D(boundaryedges,N1,N2); */
    /* Matrix boundarynodes(2,2*(N1+N2)); */
    /* FEKernel::generate_boundarynodes_2D(boundarynodes,N1,N2); */
    //组装A矩阵和V向量
    /* FEKernel::assembleA_2D(N,Nlbtest,Nlbtrial,basis_type_trial,basis_type_test,P,T,Tb,Tb,A1,func_lambda,1,0,1,0); */
    /* FEKernel::assembleA_2D(N,Nlbtest,Nlbtrial,basis_type_trial,basis_type_test,P,T,Tb,Tb,A2,func_mu,1,0,1,0); */
    /* FEKernel::assembleA_2D(N,Nlbtest,Nlbtrial,basis_type_trial,basis_type_test,P,T,Tb,Tb,A3,func_mu,0,1,0,1); */
    /* FEKernel::assembleA_2D(N,Nlbtest,Nlbtrial,basis_type_trial,basis_type_test,P,T,Tb,Tb,A4,func_lambda,0,1,1,0); */
    /* FEKernel::assembleA_2D(N,Nlbtest,Nlbtrial,basis_type_trial,basis_type_test,P,T,Tb,Tb,A5,func_mu,1,0,0,1); */
    /* FEKernel::assembleA_2D(N,Nlbtest,Nlbtrial,basis_type_trial,basis_type_test,P,T,Tb,Tb,A6,func_lambda,1,0,0,1); */
    /* FEKernel::assembleA_2D(N,Nlbtest,Nlbtrial,basis_type_trial,basis_type_test,P,T,Tb,Tb,A7,func_mu,0,1,1,0); */
    /* FEKernel::assembleA_2D(N,Nlbtest,Nlbtrial,basis_type_trial,basis_type_test,P,T,Tb,Tb,A8,func_lambda,0,1,0,1); */
    /* A.topLeftCorner(A1.rows(),A1.cols()) = A.topLeftCorner(A1.rows(),A1.cols()) + A1+2*A2+A3; */
    /* A.topRightCorner(A1.rows(),A1.cols()) = A.topRightCorner(A1.rows(),A1.cols()) + A4 + A5; */
    /* A.bottomLeftCorner(A1.rows(),A1.cols()) = A.bottomLeftCorner(A1.rows(),A1.cols()) + A6 + A7; */
    /* A.bottomRightCorner(A1.rows(),A1.cols()) = A.bottomRightCorner(A1.rows(),A1.cols()) + A8 + 2*A3 + A2; */

    /* FEKernel::assembleV_2D(N,Nlbtest,basis_type_test,0,0,P,T,Pb,Tb,V1,f1); */
    /* FEKernel::assembleV_2D(N,Nlbtest,basis_type_test,0,0,P,T,Pb,Tb,V2,f2); */
    /* V.head(V1.size()) = V1; */
    /* V.tail(V1.size()) = V2; */


    //一致三角剖分，二次元**********************************************************************************
    //生成信息矩阵P,T,Pb,Tb
    double right = 1,top = 1,left = 0,bottom = 0;
    int N1 = 16,N2 = N1;
    int N1_basis = 2*N1; //基函数节点个数N1方向
    int N2_basis = 2*N2; //基函数节点个数N2方向
    int nbn = 2*(N1_basis+N2_basis);
    Matrix P(2,(N1+1)*(N2+1)),T(3,N1*N2*2);
    FEKernel::generate_PT_2D(P,T,N1,N2,right,top,bottom,left);
    /* //202使用二次元 */
    int basis_type_trial = 202;
    int basis_type_test = 202;
    int nlb = 6;
    int Nb = (N1+1)*(N2+1)+N1*(N2+1)+N2*(N1+1)+N1*N2; //总基函数节点个数
    int N = 2*N1*N2;
    int Nlbtest = nlb;
    int Nlbtrial = nlb;
    Matrix Pb(2,Nb),Tb(6,N);
    FEKernel::generate_PbTb_2D(Pb,Tb,N1,N2,basis_type_trial,right,top,bottom,left);
    SpMat A1(Nb,Nb),A2(A1.rows(),A1.cols()),A3(A1.rows(),A1.cols()),A4(A1.rows(),A1.cols()),A5(A1.rows(),A1.cols()),A6(A1.rows(),A1.cols()),A7(A1.rows(),A1.cols()),A8(A1.rows(),A1.cols());
    Matrix A(2*A1.rows(),2*A1.cols());
    Vector V1(Nb),V2(V1.size()),V(2*V1.size()),x(V.size()),x1(V1.size()),x2(V1.size());
    /* //生成边界矩阵，边界边，边界点 */
    Matrix boundaryedges(4,2*N1+2*N2);
    FEKernel::generate_boundaryedges_2D(boundaryedges,N1,N2);
    Matrix boundarynodes(2,nbn);
    FEKernel::generate_boundarynodes_2D(boundarynodes,N1_basis,N2_basis);
    /* //组装A矩阵和V向量 */
    FEKernel::assembleA_2D(N,Nlbtest,Nlbtrial,basis_type_trial,basis_type_test,P,T,Tb,Tb,A1,func_lambda,1,0,1,0);
    FEKernel::assembleA_2D(N,Nlbtest,Nlbtrial,basis_type_trial,basis_type_test,P,T,Tb,Tb,A2,func_mu,1,0,1,0);
    FEKernel::assembleA_2D(N,Nlbtest,Nlbtrial,basis_type_trial,basis_type_test,P,T,Tb,Tb,A3,func_mu,0,1,0,1);
    FEKernel::assembleA_2D(N,Nlbtest,Nlbtrial,basis_type_trial,basis_type_test,P,T,Tb,Tb,A4,func_lambda,0,1,1,0);
    FEKernel::assembleA_2D(N,Nlbtest,Nlbtrial,basis_type_trial,basis_type_test,P,T,Tb,Tb,A5,func_mu,1,0,0,1);
    FEKernel::assembleA_2D(N,Nlbtest,Nlbtrial,basis_type_trial,basis_type_test,P,T,Tb,Tb,A6,func_lambda,1,0,0,1);
    FEKernel::assembleA_2D(N,Nlbtest,Nlbtrial,basis_type_trial,basis_type_test,P,T,Tb,Tb,A7,func_mu,0,1,1,0);
    FEKernel::assembleA_2D(N,Nlbtest,Nlbtrial,basis_type_trial,basis_type_test,P,T,Tb,Tb,A8,func_lambda,0,1,0,1);
    A.topLeftCorner(A1.rows(),A1.cols()) = A.topLeftCorner(A1.rows(),A1.cols()) + A1+2*A2+A3;
    A.topRightCorner(A1.rows(),A1.cols()) = A.topRightCorner(A1.rows(),A1.cols()) + A4 + A5;
    A.bottomLeftCorner(A1.rows(),A1.cols()) = A.bottomLeftCorner(A1.rows(),A1.cols()) + A6 + A7;
    A.bottomRightCorner(A1.rows(),A1.cols()) = A.bottomRightCorner(A1.rows(),A1.cols()) + A8 + 2*A3 + A2;

    FEKernel::assembleV_2D(N,Nlbtest,basis_type_test,0,0,P,T,Pb,Tb,V1,f1);
    FEKernel::assembleV_2D(N,Nlbtest,basis_type_test,0,0,P,T,Pb,Tb,V2,f2);
    V.head(V1.size()) = V1;
    V.tail(V1.size()) = V2;

    //处理边界条件****************************************************
    FEKernel::treat_boundary_Dirichlet_2d(boundarynodes,A,Pb,V,g1,g2);
    /* /1* cout << "A with Eigen" << A <<endl; *1/ */
    /* /1* cout << "b with Eigen" <<endl<< V <<endl; *1/ */

    /* /1* /2* //求解线性方程组 *2/ *1/ */
    x = A.colPivHouseholderQr().solve(V);
    V1 = x.head(V1.size());
    V2 = x.tail(V1.size());

    /* /1* cout << A.rows() << A.cols() <<endl; *1/ */
    /* /1* cout << "AA" << A <<endl; *1/ */
    /* /1* cout << "b" << V << endl; *1/ */
    /* x = solver.solve(V); */
    /* /1* cout << "x with Eigen" << endl<<x << endl; *1/ */

    /* /1* /2* //计算误差 *2/ *1/ */
    /* /1* /2* //组装解析解向量 *2/ *1/ */
    /* for(int i = 0; i != V1.size(); ++i) { */
    /*     x1[i] = u1(Pb(0,i),Pb(1,i)); */
    /*     x2[i] = u2(Pb(0,i),Pb(1,i)); */
    /* } */
    /* /1* cout << "analytic_solution with Eigen" << endl <<  x1 <<endl; *1/ */
    /* /1* /2* //计算无穷范数 *2/ *1/ */
    /*pass*/

    double error_max1 = FEKernel::compute_Max_error_2D(V1,Pb,u1);
    double error_max2 = FEKernel::compute_Max_error_2D(V2,Pb,u2);
    cout <<"max"<< max(error_max1,error_max2) << endl;

    /* /1* //2范数:der_x = der_y = 0, *1/ */
    int der_x=0;
    int der_y=0;
    double error1 = FEKernel::compute_Hs_error_2D(N,P,T,Tb,basis_type_trial,der_x,der_y,V1,u1);
    double error2 = FEKernel::compute_Hs_error_2D(N,P,T,Tb,basis_type_trial,der_x,der_y,V2,u2);
    double error = sqrt(error1*error1+error2*error2);
    cout << "errorl2:"<<error <<endl;

    /* /1* //H1半范 *1/ */
    int der_x1 = 0,der_y1=1;
    double error3 = FEKernel::compute_Hs_error_2D(N,P,T,Tb,basis_type_trial,der_x1,der_y1,V1,u1_1_der);//传入解析解的一阶偏导数
    double error4 = FEKernel::compute_Hs_error_2D(N,P,T,Tb,basis_type_trial,der_x1,der_y1,V2,u2_1_der);//传入解析解的一阶偏导数
    double error5 = sqrt(error3*error3+error4*error4);
    cout << "errorH1:" <<error5 <<endl;
}
