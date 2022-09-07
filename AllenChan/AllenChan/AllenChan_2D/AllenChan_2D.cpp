#include <iostream>
#include <cmath>
#include <fstream>
#include "algebra/Algebra_kernel.h"
#include "FEMethod/FEMethod_kernel.h"
#include "Constants.h"
#include "thirdparty/matplotlibcpp.h"
using namespace std;
namespace plt = matplotlibcpp;
typedef typename PDWSC::FEMethod_kernel<double, int> FEKernel;
typedef typename Eigen::SparseMatrix<double> SpMat;
typedef typename Eigen::MatrixXd Matrix;
typedef typename Eigen::VectorXd Vector;

double pi =  WHYSC::Constants::pi;
double chi = 1;
double e = 0.02;
double tend=0.001;
double coefunction(double x,double y) {
    return chi;
}

double u0(double x,double y) {
    return x*x*(x-1)*(x-1)*y*y*(y-1)*(y-1);
}


double function_one(double x,double y) {
    return 1;
}

double u(double x,double y, double t) {
    return exp(t)*u0(x,y);
}

double u_ana(double x,double y) {
    return u(x,y,tend);
}

double ux(double x,double y,double t) {
    return exp(t)*(4*x*x*x + 2*x - 6*x*x)*(pow(y,4)+y*y - 2*y*y*y);
}

double ux_ana(double x,double y) {
    return ux(x,y,tend);
}
double uxx(double x,double y,double t) {
    return exp(t)*(12*x*x + 2 - 12*x)*(pow(y,4)+ y*y - 2*y*y*y);
}

double uy(double x,double y,double t) {
    return exp(t)*(4*y*y*y +2*y-6*y*y)*(pow(x,4) + x*x - 2*x*x*x);
}

double uy_ana(double x,double y) {
    return uy(x,y,tend);
}

double uyy(double x,double y,double t) {
    return exp(t)*(12*y*y+2-12*y)*(pow(x,4) + x*x - 2*x*x*x);
}

double ut(double x,double y,double t) {
    return exp(t)*x*x*(x-1)*(x-1)*y*y*(y-1)*(y-1);
}


double f(double x,double y,double t) {
    return ut(x,y,t) - uxx(x,y,t) - uyy(x,y,t) + pow(e,-2)*pow(u(x,y,t),3) - pow(e,-2)*u(x,y,t);
}

int main() {
    ofstream fout("result_Allen.txt");
    //一致，三角剖分，线性元***********************************************************
    //生成信息矩阵P,T,Pb,Tb
    double right = 1,top = 1,left = 0,bottom = 0;
    int N1 = 10,N2 = N1;
	int nt = 10000;
    double dt = tend/nt;
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
 
    //组装M，L
    SpMat M(Nb,Nb);
    FEKernel::assembleA_2D(N,Nlb,Nlb,basis_type,basis_type,P,T,Tb,Tb,M,function_one,0,0,0,0);
    SpMat L1(Nb,Nb),L2(Nb,Nb);
    FEKernel::assembleA_2D(N,Nlb,Nlb,basis_type,basis_type,P,T,Tb,Tb,L1,function_one,1,0,1,0);
    FEKernel::assembleA_2D(N,Nlb,Nlb,basis_type,basis_type,P,T,Tb,Tb,L2,function_one,0,1,0,1);
    SpMat L = L1 + L2;
    //初始化迭代向量
    Vector xn(Nb);
    for(int i = 0; i < Nb; ++i) {
        xn(i) = u0(Pb(0,i),Pb(1,i));
    }
    Vector xn1 = xn;

    vector<vector<double>> X,Y,Z;
    //数值解图像
    for(int i = 0; i < N1+1; i++) {
        vector<double> x_mark,y_mark,z_mark;
        for(int j = 0; j < N1+1; j++) {
            int k = (N1 + 1)*i + j;
            x_mark.push_back(Pb(0,k));
            y_mark.push_back(Pb(1,k));
        }
        X.push_back(x_mark);
        Y.push_back(y_mark);
    }

	Vector F(Nb);
	Vector u_old(Nb);
    //时间迭代
    for(int i = 1; i < nt+1; ++i) {
        double t = i*dt;

		FEKernel::assembleV_2D_t(t,N,Nlb,basis_type,0,0,P,T,Pb,Tb,F,f);

		Vector b = 1.0/dt*M*xn + F;
		u_old = xn;
		//牛顿迭代
		for(int j=0;j<100;++j){
			SpMat G(Nb, Nb);
			FEKernel::assembleFE_coe_pow2_2D(0,0,basis_type,u_old,N,Nlb,Nlb,basis_type,basis_type,P,T,Tb,Tb,G,coefunction,0,0,0,0);
			SpMat A = 1.0/dt*M + L - pow(e,-2)*M + pow(e,-2)*G;
			SpMat J = 1.0/dt*M + L - pow(e,-2)*M + 3*pow(e,-2)*G;
			Vector res = b - A*u_old;
			Eigen::SparseLU<SpMat> solver;
			solver.compute(J);
			xn1 = solver.solve(res) + u_old;
			cout << "j:" << j  << "norm:" << (xn1 - u_old).norm() << endl;
			if((xn1 - u_old).norm() < pow(10,-5)){
				cout << "牛顿迭代次数:" << j+1 << endl << "norm:" << (xn1 - u_old).norm() << endl;
				break;
			}
			u_old = xn1;	
		}
		xn = xn1;
		cout << "t" << t << endl; 
    }
    
    //计算终止时刻的误差
    cout << "结束时间为: " << nt*dt << endl;



    /*******************计算误差*************************/
    //最大模范数
    double error_max_u = FEKernel::compute_Max_error_2D(N,P,T,Tb,basis_type,xn1,u_ana);
    fout << "error_max_u:" << error_max_u << endl;
    //2范数
    int der_y=0,der_x = 0;
    double erroru = FEKernel::compute_Hs_error_2D(N,P,T,Tb,basis_type,der_x,der_y,xn1,u_ana);
    fout << "error_L2_u:"<<erroru <<endl;
    //H1半范
    double error3 = FEKernel::compute_Hs_error_2D(N,P,T,Tb,basis_type,1,0,xn1,ux_ana);//传入解析解的一阶偏导数
    double error31 = FEKernel::compute_Hs_error_2D(N,P,T,Tb,basis_type,0,1,xn1,uy_ana);//传入解析解的一阶偏导数
    double error_H_u = sqrt(error3*error3);
    fout << "errorH1_u:" <<error_H_u <<endl;

    fout.close();

}
