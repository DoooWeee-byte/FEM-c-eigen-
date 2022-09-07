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
const double nv = 2;
double func_nv(double x,double y) {
    return nv;
}
double func_negtiveone(double x,double y) {
    return -1;
}
double g1(double x,double y)
{
    double r = 0;
    if(x == 0) {
        r = exp(-y);
    } else if(x == 1) {
        r = y*y + exp(-y);
    } else if(y == -0.25) {
        r = 1.0/16*x*x + exp(0.25);
    } else if(y == 0) {
        r = 1;
    }
    return r;
}
double g2(double x,double y) {
    double r = 0;
    if(x == 0) {
        r = 2;
    } else if(x == 1) {
        r = -2.0/3*y*y*y+2;
    } else if(y == -0.25) {
        r = 1.0/96*x + 2 -pi*sin(pi*x);
    } else if(y == 0) {
        r = 2 - pi*sin(pi*x);
    }
    return r;
}

double f1(double x,double y) {
    return -2*nv*x*x-2*nv*y*y-nv*exp(-y)+pi*pi*cos(pi*x)*cos(2*pi*y);
}
double f2(double x,double y) {
    return 4*nv*x*y - nv*pi*pi*pi*sin(pi*x) + 2*pi*(2-pi*sin(pi*x))*sin(2*pi*y);
}
double u1(double x,double y)
{
    return x*x*y*y+exp(-y);
}
double u2(double x,double y)
{
    return -2.0/3*x*y*y*y + 2 - pi*sin(pi*x);
}

double p(double x, double y) {
    return -(2 - pi*sin(pi*x))*cos(2*pi*y);
}


double u1_1_der_x(double x,double y) {
    return 2*y*y*x;
}
double u2_1_der_x(double x,double y) {
	return -2.0/3*y*y*y - pi*pi*cos(pi*x);
}
double u1_1_der_y(double x,double y) {
    return 2*x*x*y-exp(-y);
}
double u2_1_der_y(double x,double y) {
	return -2*x*y*y;
}
double p_1_der_x(double x,double y){
	return pi*pi*cos(2*pi*y)*cos(pi*x);
}
double p_1_der_y(double x,double y){
	return 2*pi*(2 - pi*sin(pi*x))*sin(2*pi*y);
}


int main() {
    //一致，三角剖分，线性元***********************************************************
    //生成信息矩阵P,T,Pb,Tb
    double right = 1,top = 0,left = 0,bottom = -0.25;
	int N1 =8 ,N2 = N1/4;
    int N = 2*N1*N2;
    int Nb_u = (N1+1)*(N2+1)+N1*(N2+1)+N2*(N1+1)+N1*N2; //总基函数节点个数 */
    int Nb_p = (N1+1)*(N2+1);
	int N1_basis_u = 2*N1 + 1;
	int N2_basis_u = 2*N2 + 1;
	int N1_basis_p = N1 + 1;
	int N2_basis_p = N2 + 1;
    
    
    Matrix P(2,Nb_p),T(3,N);
    FEKernel::generate_PT_2D(P,T,N1,N2,right,top,bottom,left);
    //201使用线性元,202使用二次元，对速度u使用二次元，压力p使用线性元
    int basis_type_u = 202;
    int basis_type_p = 201;
    int nlb_p = 3;
    int nlb_u = 6;
 
 
    Matrix Pb_u(2,Nb_u),Tb_u(nlb_u,N);
    FEKernel::generate_PbTb_2D(Pb_u,Tb_u,N1,N2,basis_type_u,right,top,bottom,left);


    Matrix Pb_p(2,Nb_p),Tb_p(nlb_p,N);
    FEKernel::generate_PbTb_2D(Pb_p,Tb_p,N1,N2,basis_type_p,right,top,bottom,left);

    //生成边界矩阵，边界点
    int N1_basis = 2*N1; //基函数节点个数N1方向 */
    int N2_basis = 2*N2; //基函数节点个数N2方向
    int nbn = 2*(N1_basis+N2_basis);
    Matrix boundarynodes(2,nbn);
    FEKernel::generate_boundarynodes_2D(boundarynodes,N1_basis,N2_basis);




    Matrix A(Nb_u+Nb_u+Nb_p,Nb_u+Nb_u+Nb_p);
    A = A*0;
    Vector V1(Nb_u),V2(Nb_u),V3(Nb_p),V(Nb_u+Nb_u+Nb_p),x(V.size());



    //组装A矩阵和V向量
    if(1) {
        /* SpMat A1(Nb_u,Nb_u),A2(Nb_u,Nb_u),A3(Nb_u,Nb_u),A4(Nb_u,Nb_u),A5(Nb_u,Nb_p),A6(Nb_u,Nb_p),A7(Nb_p,Nb_u),A8(Nb_p,Nb_u); */
        /* FEKernel::assembleA_2D(N,nlb_u,nlb_u,basis_type_u,basis_type_u,P,T,Tb_u,Tb_u,A1,func_nv,1,0,1,0); */
        /* FEKernel::assembleA_2D(N,nlb_u,nlb_u,basis_type_u,basis_type_u,P,T,Tb_u,Tb_u,A2,func_nv,0,1,0,1); */
        /* FEKernel::assembleA_2D(N,nlb_u,nlb_u,basis_type_u,basis_type_u,P,T,Tb_u,Tb_u,A3,func_nv,0,1,1,0); */
        /* FEKernel::assembleA_2D(N,nlb_u,nlb_u,basis_type_u,basis_type_u,P,T,Tb_u,Tb_u,A4,func_nv,1,0,0,1); */
        /* FEKernel::assembleA_2D(N,nlb_u,nlb_p,basis_type_u,basis_type_p,P,T,Tb_u,Tb_p,A5,func_negtiveone,1,0,0,0); */
        /* FEKernel::assembleA_2D(N,nlb_u,nlb_p,basis_type_u,basis_type_p,P,T,Tb_u,Tb_p,A6,func_negtiveone,0,1,0,0); */
        /* FEKernel::assembleA_2D(N,nlb_p,nlb_u,basis_type_p,basis_type_u,P,T,Tb_p,Tb_u,A7,func_negtiveone,0,0,1,0); */
        /* FEKernel::assembleA_2D(N,nlb_p,nlb_u,basis_type_p,basis_type_u,P,T,Tb_p,Tb_u,A8,func_negtiveone,0,0,0,1); */


        /* A.topLeftCorner(A1.rows(),A1.cols()) = 2*A1+ A2; */
        /* A.topRightCorner(A5.rows(),A5.cols()) =  A5; */
        /* A.bottomLeftCorner(A7.rows(),A7.cols()) = A7; */
        /* A.block(0,A1.cols(),A3.rows(),A3.cols()) = A3; */
        /* A.block(A1.rows(),0,A4.rows(),A4.cols()) = A4; */
        /* A.block(A1.rows()+A4.rows(),A1.cols(),A8.rows(),A8.cols()) = A8; */
        /* A.block(A1.rows(),A1.cols(),A2.rows(),A2.cols()) = 2*A2 + A1; */
        /* A.block(A1.rows(),A1.cols()+A1.cols(),A6.rows(),A6.cols()) = A6; */
		SpMat A1(Nb_u,Nb_u),A2(Nb_u,Nb_u),A3(Nb_u,Nb_u),A5(Nb_u,Nb_p),A6(Nb_u,Nb_p);
        FEKernel::assembleA_2D(N,nlb_u,nlb_u,basis_type_u,basis_type_u,P,T,Tb_u,Tb_u,A1,func_nv,1,0,1,0);
        FEKernel::assembleA_2D(N,nlb_u,nlb_u,basis_type_u,basis_type_u,P,T,Tb_u,Tb_u,A2,func_nv,0,1,0,1);
        FEKernel::assembleA_2D(N,nlb_u,nlb_u,basis_type_u,basis_type_u,P,T,Tb_u,Tb_u,A3,func_nv,0,1,1,0);
        FEKernel::assembleA_2D(N,nlb_u,nlb_p,basis_type_u,basis_type_p,P,T,Tb_u,Tb_p,A5,func_negtiveone,1,0,0,0);
        FEKernel::assembleA_2D(N,nlb_u,nlb_p,basis_type_u,basis_type_p,P,T,Tb_u,Tb_p,A6,func_negtiveone,0,1,0,0);


        A.topLeftCorner(A1.rows(),A1.cols()) = 2*A1+ A2;
        A.topRightCorner(A5.rows(),A5.cols()) =  A5;
        A.bottomLeftCorner(A5.cols(),A5.rows()) = A5.transpose();
        A.block(0,A1.cols(),A3.rows(),A3.cols()) = A3;
        A.block(A1.rows(),0,A3.cols(),A3.rows()) = A3.transpose();
        A.block(A1.rows()+A2.rows(),A1.cols(),A6.cols(),A6.rows()) = A6.transpose();
        A.block(A1.rows(),A1.cols(),A2.rows(),A2.cols()) = 2*A2 + A1;
        A.block(A1.rows(),A1.cols()+A1.cols(),A6.rows(),A6.cols()) = A6;
    }
    
    V1 = V1*0;
    V2 = V2*0;
    FEKernel::assembleV_2D(N,nlb_u,basis_type_u,0,0,P,T,Pb_u,Tb_u,V1,f1);
    FEKernel::assembleV_2D(N,nlb_u,basis_type_u,0,0,P,T,Pb_u,Tb_u,V2,f2);
    V = V*0;
    V.head(V1.size()) = V1;
    V.segment(V1.size(),V2.size()) = V2;


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

    //处理边界条件:把p一个值固定****************************************************
    FEKernel::treat_boundary_Dirichlet_2d(boundarynodes,A,Pb_u,V,g1,g2);
   

    if(1) {
        Vector x0(A.cols());
        x0 = x0*0;
		int i = 2*Nb_u;
		int j = 0;
        x0(2*Nb_u + j ) = 1;
        A.row(i) = x0.transpose();
        V(i) = p(Pb_p(0,j),Pb_p(1,j));
    }


    /* /1* cout << "A with Eigen" << A <<endl; *1/ */
    /* /1* cout << "b with Eigen" <<endl<< V <<endl; *1/ */

    /* /1* //求解线性方程组 *1/ */
    x = A.colPivHouseholderQr().solve(V);
    V1 = x.head(V1.size());
    V2 = x.segment(V1.size(),V2.size());
    V3 = x.tail(V3.size());

	//最大模范数
	double error_max_1 = FEKernel::compute_Max_error_2D(V1,Pb_u,u1);
	double error_max_2 = FEKernel::compute_Max_error_2D(V2,Pb_u,u2);
	double error_max_p = FEKernel::compute_Max_error_2D(V3,Pb_p,p);
	double error_max_u = max(error_max_1,error_max_2);
	cout << "error_max_u:" << error_max_u << endl;
	cout << "error_max_p:" << error_max_p << endl;
    //2范数  
    int der_y=0,der_x=0;
    double error1 = FEKernel::compute_Hs_error_2D(N,P,T,Tb_u,basis_type_u,der_x,der_y,V1,u1);
    double error2 = FEKernel::compute_Hs_error_2D(N,P,T,Tb_u,basis_type_u,der_x,der_y,V2,u2);
    double error30 = FEKernel::compute_Hs_error_2D(N,P,T,Tb_p,basis_type_p,der_x,der_y,V3,p);
    double error = sqrt(error1*error1+error2*error2);
    cout << "error_L2_u:"<<error <<endl;
    cout << "error_L2_p:" << error30 << endl;
    //H1半范
    double error3 = FEKernel::compute_Hs_error_2D(N,P,T,Tb_u,basis_type_u,1,0,V1,u1_1_der_x);//传入解析解的一阶偏导数
    double error31 = FEKernel::compute_Hs_error_2D(N,P,T,Tb_u,basis_type_u,0,1,V1,u1_1_der_y);//传入解析解的一阶偏导数
    double error4 = FEKernel::compute_Hs_error_2D(N,P,T,Tb_u,basis_type_u,1,0,V2,u2_1_der_x);//传入解析解的一阶偏导数
    double error41 = FEKernel::compute_Hs_error_2D(N,P,T,Tb_u,basis_type_u,0,1,V2,u2_1_der_y);//传入解析解的一阶偏导数
    double error5 = FEKernel::compute_Hs_error_2D(N,P,T,Tb_p,basis_type_p,1,0,V3,p_1_der_x);//传入解析解的一阶偏导数
    double error51 = FEKernel::compute_Hs_error_2D(N,P,T,Tb_p,basis_type_p,0,1,V3,p_1_der_y);//传入解析解的一阶偏导数
	double error_H_u = sqrt(error3*error3+error31*error31+error4*error4+error41*error41);
	double error_H_p = sqrt(error5*error5+error51*error51);
    cout << "errorH1_u:" <<error_H_u <<endl;
    cout << "errorH1_p:" <<error_H_p <<endl;

}
