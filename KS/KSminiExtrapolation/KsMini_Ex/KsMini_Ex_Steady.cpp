#include <iostream>
#include <cmath>
#include <fstream>
#include "algebra/Algebra_kernel.h"
#include "FEMethod/FEMethod_kernel.h"
#include "Constants.h"
using namespace std;
typedef typename PDWSC::FEMethod_kernel<double, int> FEKernel;
typedef typename Eigen::SparseMatrix<double> SpMat;
typedef typename Eigen::MatrixXd Matrix;
typedef typename Eigen::VectorXd Vector;

double pi =  WHYSC::Constants::pi;
double chi = 1;
double Tol = pow(10,-6);
double epsilon = -pow(10,-5);
int maxstep = 10000;
double dt = 100;


double coefunction(double x,double y){
	return chi;
}

double u(double x,double y) {

    return 0.025*(cos(2.0*pi*(y - 0.5)) + 1)*(cos(2.0*pi*(x - 0.5))+1);
}

double c(double x,double y) {
    return 0.5*u(x,y);
}
double v(double x,double y){
	return c(x,y);
}

double u0(double x,double y) {
    return (1 - 0.1*0.1)*u(x,y);
	/* return u(x,y); */
}

double c0(double x,double y) {
    return (1 - 0.1*0.1)*c(x,y);
	/* return c(x,y); */
}

double function_one(double x,double y) {
    return 1;
}

double ux(double x,double y){
	return -pi/20*(cos(2*pi*y) - 1)*sin(2*pi*x);
}

double uxx(double x,double y){
	return -(pi*pi)/10*(cos(2*pi*y) - 1)*cos(2*pi*x);
}

double uy(double x,double y){
	return -pi/20*(cos(2*pi*x) - 1)*sin(2*pi*y);
}

double uyy(double x,double y){
	return -(pi*pi)/10*(cos(2*pi*x) - 1)*cos(2*pi*y);
}

double vx(double x,double y){
	return 0.5*ux(x,y);
}

double vxx(double x,double y){
	return 0.5*uxx(x,y);
}

double vy(double x,double y){
	return 0.5*uy(x,y);
}

double vyy(double x,double y){
	return 0.5*uyy(x,y);
}


double f1(double x,double y){
	return -(uxx(x,y) + uyy(x,y) - chi*(ux(x,y)*vx(x,y)) -chi*(uy(x,y)*vy(x,y))- chi*u(x,y)*(vxx(x,y) + vyy(x,y)));
}

double f2(double x,double y){
   return -(vxx(x,y) + vyy(x,y) - v(x,y) + u(x,y));	
}

int main() {
	ofstream fout("result1.txt");
    //一致，三角剖分，线性元***********************************************************
    //生成信息矩阵P,T,Pb,Tb
    double right = 1,top = 1,left = 0,bottom = 0;
    int N1 = 64,N2 = N1;
    Matrix P(2,(N1+1)*(N2+1)),T(3,N1*N2*2);
    FEKernel::generate_PT_2D(P,T,N1,N2,right,top,bottom,left);
    //201使用线性元
    int basis_type = 201;
    int Nlb = 3;
    int N = 2*N1*N2;
    int Nb = (N1+1)*(N2+1);
    //组装Pb,Tb
    Matrix Pb(2,(N1+1)*(N2+1)),Tb(Nlb,N);
    FEKernel::generate_PbTb_2D(Pb,Tb,N1,N2,basis_type,right,top,bottom,left);

	//组装出精确解向量
	Vector analytic_solution(2*Nb);
	for(int i = 0;i < Nb;++i){
		analytic_solution(i) = u(Pb(0,i),Pb(1,i));
		analytic_solution(i+Nb) = c(Pb(0,i),Pb(1,i));
	}

    //组装M，L
    SpMat M(Nb,Nb),L(Nb,Nb);
    FEKernel::assembleA_2D(N,Nlb,Nlb,basis_type,basis_type,P,T,Tb,Tb,M,function_one,0,0,0,0);
    if(1) {
        SpMat L1(Nb,Nb),L2(Nb,Nb);
        FEKernel::assembleA_2D(N,Nlb,Nlb,basis_type,basis_type,P,T,Tb,Tb,L1,function_one,1,0,1,0);
        FEKernel::assembleA_2D(N,Nlb,Nlb,basis_type,basis_type,P,T,Tb,Tb,L2,function_one,0,1,0,1);
        L = L1 + L2;
    }

	//组装右端额外项
	Vector bextal(2*Nb);
	if(1){
		Vector tempv(Nb),tempv2(Nb);
		FEKernel::assembleV_2D(N,Nlb,basis_type,0,0,P,T,Pb,Tb,tempv,f1);
		bextal.head(Nb) = tempv;
		FEKernel::assembleV_2D(N,Nlb,basis_type,0,0,P,T,Pb,Tb,tempv2,f2);
		bextal.tail(Nb) = tempv2;
	}

	Matrix A(2*Nb,2*Nb);
	A *= 0;
	//初始化迭代向量
	Vector xn(2*Nb);
	for(int i = 0;i < Nb;++i){
		xn(i) = u0(Pb(0,i),Pb(1,i));
		xn(i + Nb) = c0(Pb(0,i),Pb(1,i));
	}
	Vector xn1 = xn;
		
	cout << "initial_L2:" << (xn - analytic_solution).norm() << endl;
	Vector  x(2*Nb);


	//迭代
	for(int i = 0;i < maxstep;++i){

		//取出vh这个向量
		Vector xlin = (2*xn - xn1).tail(Nb);
		//组装出K矩阵
		SpMat K1(Nb,Nb),K2(Nb,Nb);
		FEKernel::assembleFE_coe_2D(1,0,basis_type,xlin,N,Nlb,Nlb,basis_type,basis_type,P,T,Tb,Tb,K1,coefunction,1,0,0,0);
		FEKernel::assembleFE_coe_2D(0,1,basis_type,xlin,N,Nlb,Nlb,basis_type,basis_type,P,T,Tb,Tb,K2,coefunction,0,1,0,0);
		SpMat K = K2 + K1;
		
		//计算右端项
		Matrix M1(Nb*2,Nb*2);
		M1 *= 0;
		M1.topLeftCorner(Nb,Nb) = M;
		M1.bottomRightCorner(Nb,Nb) = M;
		Vector b = M1*xn;
		b = b + dt*bextal;
		
		//组装矩阵
		A.topLeftCorner(Nb,Nb) = M + dt*(L - K);
		A.bottomLeftCorner(Nb,Nb) = -dt*M;
		A.bottomRightCorner(Nb,Nb) = M + dt*(L + M);

		//求解代数系统
		x = A.colPivHouseholderQr().solve(b);
		if((x - xn).norm() < Tol){
			fout << "i:" << i << endl;
			cout << " L2norm: "<< (x - analytic_solution).norm() << endl;
			fout.close();
			break;
		}
		cout << "i:" << i << endl;
		cout << " L2norm: "<< (x - analytic_solution).norm() << endl;
		xn1 = xn;
		xn = x;
		A *= 0;
    }
	fout.close();

    

    
}
