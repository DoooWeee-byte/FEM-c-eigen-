#include <iostream>
#include <cmath>
#include <fstream>
#include "algebra/Algebra_kernel.h"
#include "FEMethod/FEMethod_kernel.h"
#include "Constants.h"
#include<Eigen/SparseLU>
using namespace std;
typedef typename PDWSC::FEMethod_kernel<double, int> FEKernel;
typedef typename Eigen::SparseMatrix<double> SpMat;
typedef typename Eigen::MatrixXd Matrix;
typedef typename Eigen::VectorXd Vector;

double pi =  WHYSC::Constants::pi;
double chi = 1;
int maxstep = 51;
double dt = pow(10,-4);
namespace plt = matplotlibcpp;
double Cu = 60;
double Cv = Cu;
double endtime = 0.005;
double coefunction(double x,double y) {
    return chi;
}

double u0(double x,double y) {
    return Cu*exp(-Cu*(x*x + y*y));
}

double c0(double x,double y) {
    return Cv*exp(-Cv*(x*x + (y - 0.5)*(y - 0.5)));
}

double function_one(double x,double y) {
    return 1;
}

int main() {
    ofstream fout("poscheckerror.txt");
    //一致，三角剖分，线性元***********************************************************
    //生成信息矩阵P,T,Pb,Tb
    double right = 0.5,top = 0.5,left = -0.5,bottom = -0.5;
    int N1 =128,N2 = N1;
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

    //组装L
    SpMat L(Nb,Nb);

    if(1) {
        SpMat L1(Nb,Nb),L2(Nb,Nb);
        FEKernel::assembleA_2D(N,Nlb,Nlb,basis_type,basis_type,P,T,Tb,Tb,L1,function_one,1,0,1,0);
        FEKernel::assembleA_2D(N,Nlb,Nlb,basis_type,basis_type,P,T,Tb,Tb,L2,function_one,0,1,0,1);
        L = L1 + L2;
    }

	//组装M
	Vector Mtemp(Nb);
	FEKernel::assembleV_2D(N,Nlb,basis_type,0,0,P,T,Pb,Tb,Mtemp,function_one);
	SpMat M(Nb,Nb);
	typedef Eigen::Triplet<double> Tk;
    std::vector<Tk> tripletList;
	for(int i = 0;i < Nb;++i){
		tripletList.push_back(Tk(i,i,Mtemp(i)));
	}
	M.setFromTriplets(tripletList.begin(), tripletList.end());
	//画图向量
	std::vector<double> x_plt(maxstep),y_plt(maxstep),z_plt(maxstep);

    //初始化迭代向量
    Vector xn(2*Nb);
    for(int i = 0; i < Nb; ++i) {
        xn(i) = u0(Pb(0,i),Pb(1,i));
        xn(i + Nb) = c0(Pb(0,i),Pb(1,i));
    }
    Vector xn1 = xn;
	
	
    fout << "time:" << 0 <<"  "<<  "min:" << xn.minCoeff() << endl;
 
	x_plt.at(0) = 0;
	y_plt.at(0) = xn.head(Nb).minCoeff();
	z_plt.at(0) = xn.tail(Nb).minCoeff();
	//解向量初始化
    Vector  x(2*Nb);
	x *= 0;

	//定义两个中间向量，接受U和v的解
	Vector xu(Nb);
	Vector xv(Nb);

	    //迭代
	Eigen::SparseLU<SpMat> solver;
    for(int i = 1; i < maxstep; ++i) {
        //取出vh这个向量
        /* Vector xlin = (2*xn - xn1).tail(Nb); */
		Vector xlin = xn.tail(Nb);

		//组装出K矩阵
        SpMat K1(Nb,Nb),K2(Nb,Nb);
        FEKernel::assembleFE_coe_2D(1,0,basis_type,xlin,N,Nlb,Nlb,basis_type,basis_type,P,T,Tb,Tb,K1,coefunction,1,0,0,0);
        FEKernel::assembleFE_coe_2D(0,1,basis_type,xlin,N,Nlb,Nlb,basis_type,basis_type,P,T,Tb,Tb,K2,coefunction,0,1,0,0);
        SpMat K = K2 + K1;

        
        //组装矩阵
		SpMat A11(2*Nb,2*Nb);
		SpMat M1 = M + dt*(L - K);
		SpMat M2 = -dt*M;
		SpMat M3 = M + dt*(L + M);
		FEKernel::addSpMat(A11,M1,0,0);
		FEKernel::addSpMat(A11,M2,Nb,0);
		FEKernel::addSpMat(A11, M3, Nb, Nb);

		//计算右端项
		xu = xn.head(Nb);
		xv = xn.tail(Nb);
        Vector b(2*Nb);
		b.head(Nb) = M*xu;
		b.tail(Nb) = M*xv;

		
        //求解代数系统
		solver.compute(A11);
		x = solver.solve(b);
		
		//为下一次迭代作准备
        xn1 = xn;
        xn = x;

		//画图
		x_plt.at(i)	= i;
		y_plt.at(i) = x.head(Nb).minCoeff();
		z_plt.at(i) = x.tail(Nb).minCoeff();
    }

	//计算终止时刻的误差
    double t_end = (maxstep-1)*dt;
    cout << "结束时间为: " << t_end << endl;

	//画图
	/* string s = "o"; */
	/* plt::plot(x_plt,y_plt,s); */
	/* plt::plot(x_plt,z_plt,s); */
	plt::named_plot("min(u)",x_plt,y_plt);
	plt::named_plot("min(v)",x_plt,z_plt);
	plt::title("Cu = 60,maxstep = 50");
	plt::xlim(0,60);
	plt::legend();

	const char* filename = "./Cu60.png";
	std::cout << "Saving result to" << filename << std::endl;
	plt::save(filename);	

	//文件写入结束，关闭文件
    fout.close();
}
