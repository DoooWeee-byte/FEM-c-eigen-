#include <iostream>
#include <fstream>
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
double nv = 2;
double func_nv(double x,double y) {
    return nv;
}
double func_negtiveone(double x,double y) {
    return -1;
}
double c(double x,double y) {
    return 1;
}
double u10(double x,double y) {
    return x*x*y*y + exp(-y);
}
double u20(double x,double y) {
    return -2.0/3*x*y*y*y + 2 - pi*sin(pi*x);
}
double p0(double x,double y) {
    return -(2-pi*sin(pi*x))*cos(2*pi*y);
}
double g1(double x,double y,double t)
{
    double r = 0;
    if(x == 0) {
        r = exp(-y)*cos(2*pi*t);
    } else if(x == 1) {
        r = (y*y + exp(-y))*cos(2*pi*t);
    } else if(y == -0.25) {
        r = (1.0/16*x*x + exp(0.25))*cos(2*pi*t);
    } else if(y == 0) {
        r = cos(2*pi*t);
    }
    return r;
}
double g2(double x,double y,double t) {
    double r = 0;
    if(x == 0) {
        r = 2*cos(2*pi*t);
    } else if(x == 1) {
        r = (-2.0/3*y*y*y+2)*cos(2*pi*t);
    } else if(y == -0.25) {
        r = (1.0/96*x + 2 -pi*sin(pi*x))*cos(2*pi*t);
    } else if(y == 0) {
        r = (2 - pi*sin(pi*x))*cos(2*pi*t);
    }
    return r;
}

double f1(double x,double y,double t) {
    return -2*pi*(x*x*y*y + exp(-y))*sin(2*pi*t)+(-2*nv*x*x-2*nv*y*y-nv*exp(-y)+pi*pi*cos(pi*x)*cos(2*pi*y))*cos(2*pi*t);
}
double f2(double x,double y,double t) {
    return -2*pi*(-2.0/3*x*y*y*y+2-pi*sin(pi*x))*sin(2*pi*t)+(4*nv*x*y - nv*pi*pi*pi*sin(pi*x) + 2*pi*(2-pi*sin(pi*x))*sin(2*pi*y))*cos(2*pi*t);
}
double u1(double x,double y)
{
    double t = 1;
    return (x*x*y*y+exp(-y))*cos(2*pi*t);
}
double u2(double x,double y)
{
    double t = 1;
    return (-2.0/3*x*y*y*y + 2 - pi*sin(pi*x))*cos(2*pi*t);
}

double p(double x, double y) {
    double t = 1;
    return -(2 - pi*sin(pi*x))*cos(2*pi*y)*cos(2*pi*t);
}

double u1_1_der_x(double x,double y) {
    return 2*x*y*y;
}

double u1_1_der_y(double x,double y){
	return 2*x*x*y-exp(-y);
}

double u2_1_der_x(double x,double y) {
	return -2.0/3*y*y*y - pi*pi*cos(pi*x);
}

double u2_1_der_y(double x,double y){
	return -2.0*y*y*x;
}

double p_1_der_x(double x,double y){
	return cos(2*pi*y)*pi*pi*cos(pi*x);
}

double p_1_der_y(double x,double y){
	return (2 - pi*sin(pi*x))*2*pi*sin(2*pi*y);
}
int main() {
    //?????????????????????????????????***********************************************************
    //??????????????????P,T,Pb,Tb
	ofstream fout("result.txt");
    double right = 1,top = 0,left = 0,bottom = -0.25;

	int N1 = 64,N2 = N1/4;
    double h = (right - left)/N1;
	fout << "StepH:" <<h << endl;
    Matrix P(2,(N1+1)*(N2+1)),T(3,N1*N2*2);
    FEKernel::generate_PT_2D(P,T,N1,N2,right,top,bottom,left);
    //201???????????????,202???????????????????????????u????????????????????????p???????????????
    int basis_type_u = 202;
    int basis_type_p = 201;
    int nlb_p = 3;
    int nlb_u = 6;
    int N = 2*N1*N2;
    int Nb_u = (N1+1)*(N2+1)+N1*(N2+1)+N2*(N1+1)+N1*N2; //???????????????????????? */
    int Nb_p = (N1+1)*(N2+1);
    Matrix Pb_u(2,Nb_u),Tb_u(nlb_u,N);
    FEKernel::generate_PbTb_2D(Pb_u,Tb_u,N1,N2,basis_type_u,right,top,bottom,left);


    Matrix Pb_p(2,(N1+1)*(N2+1)),Tb_p(nlb_p,N);
    FEKernel::generate_PbTb_2D(Pb_p,Tb_p,N1,N2,basis_type_p,right,top,bottom,left);

    //??????????????????????????????
    int N1_basis = 2*N1; //?????????????????????N1?????? */
    int N2_basis = 2*N2; //?????????????????????N2??????
    int nbn = 2*(N1_basis+N2_basis);
    Matrix boundarynodes(2,nbn);
    FEKernel::generate_boundarynodes_2D(boundarynodes,N1_basis,N2_basis);



    Matrix A(Nb_u+Nb_u+Nb_p,Nb_u+Nb_u+Nb_p),M(A.rows(),A.cols());
	A = 0*A,M = 0*M; // ????????????0
    Vector V1(Nb_u),V2(Nb_u),x(A.rows()),V3(Nb_p);
	V1 = 0*V1,V2 = 0*V2,x = 0*x,V3 = 0*V3;
	


    //??????A?????????V??????
    if(1) {
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


    if(1) {
        SpMat A1(Nb_u,Nb_u);
        FEKernel::assembleA_2D(N,nlb_u,nlb_u,basis_type_u,basis_type_u,P,T,Tb_u,Tb_u,A1,c,0,0,0,0);
        M.topLeftCorner(A1.rows(),A1.cols()) = A1;
        M.block(A1.rows(),A1.cols(),A1.rows(),A1.cols()) = A1;
    }
    Vector x0(x.size()); //???????????????
	x0 = 0*x0;
	


/**************??????????????????*************************************/
    for(int i = 0; i < Nb_u; ++i) {
        x0(i) = u10(Pb_u(0,i),Pb_u(1,i));
        x0(i + Nb_u) = u20(Pb_u(0,i),Pb_u(1,i));
    }

    for(int i = 0; i < Nb_p; ++i) {
        x0(i+2*Nb_u) = p0(Pb_p(0,i),Pb_p(1,i));
    }
	
	
	/******************??????ODE***********************************/
	int  N_t = 512;
	double dt = 1.0/N_t;
	fout << "StepT: " << dt << endl;
	double theta = 0.5;
	fout << "Theta: " << theta << endl;
	double Tmax = 1;
	double Tmin = 0;
    FEKernel::ODE_Unsteady_Stokes(Tmax,Tmin,N_t,N,nlb_u,Nb_u,basis_type_u,theta,M,A,P,T,Pb_u,Tb_u,boundarynodes,x,x0,f1,f2,g1,g2);
   
   
	V1 = x.head(V1.size());
    V2 = x.segment(V1.size(),V2.size());
    V3 = x.tail(V3.size());


	/*******************????????????*************************/
	//???????????????
	double error_max_1 = FEKernel::compute_Max_error_2D(N,P,T,Tb_u,basis_type_u,V1,u1);
	double error_max_2 = FEKernel::compute_Max_error_2D(N,P,T,Tb_u,basis_type_u,V2,u2);
	double error_max_p = FEKernel::compute_Max_error_2D(N,P,T,Tb_p,basis_type_p,V3,p);
	double error_max_u = max(error_max_1,error_max_2);
	fout << "error_max_u:" << error_max_u << endl;
	fout << "error_max_p:" << error_max_p << endl;
    //2??????  
    int der_y=0,der_x = 0;
    double error1 = FEKernel::compute_Hs_error_2D(N,P,T,Tb_u,basis_type_u,der_x,der_y,V1,u1);
    double error2 = FEKernel::compute_Hs_error_2D(N,P,T,Tb_u,basis_type_u,der_x,der_y,V2,u2);
    double error30 = FEKernel::compute_Hs_error_2D(N,P,T,Tb_p,basis_type_p,der_x,der_y,V3,p);
    double error = sqrt(error1*error1+error2*error2);
    fout << "error_L2_u:"<<error <<endl;
    fout << "error_L2_p:" << error30 << endl;
    //H1??????
    double error3 = FEKernel::compute_Hs_error_2D(N,P,T,Tb_u,basis_type_u,1,0,V1,u1_1_der_x);//?????????????????????????????????
    double error31 = FEKernel::compute_Hs_error_2D(N,P,T,Tb_u,basis_type_u,0,1,V1,u1_1_der_y);//?????????????????????????????????
    double error4 = FEKernel::compute_Hs_error_2D(N,P,T,Tb_u,basis_type_u,1,0,V2,u2_1_der_x);//?????????????????????????????????
    double error41 = FEKernel::compute_Hs_error_2D(N,P,T,Tb_u,basis_type_u,0,1,V2,u2_1_der_y);//?????????????????????????????????
    double error5 = FEKernel::compute_Hs_error_2D(N,P,T,Tb_p,basis_type_p,1,0,V3,p_1_der_x);//?????????????????????????????????
    double error51 = FEKernel::compute_Hs_error_2D(N,P,T,Tb_p,basis_type_p,0,1,V3,p_1_der_y);//?????????????????????????????????
	double error_H_u = sqrt(error3*error3+error31*error31+error4*error4+error41*error41);
	double error_H_p = sqrt(error5*error5+error51*error51);
    fout << "errorH1_u:" <<error_H_u <<endl;
    fout << "errorH1_p:" <<error_H_p <<endl;
	fout.close();



}
