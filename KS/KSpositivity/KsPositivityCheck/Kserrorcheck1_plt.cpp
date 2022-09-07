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
/* double dt = pow(10,-4); */

double Cu = 70;
double Cv = Cu;

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

double u(double x,double y, double t) {
    return Cu*exp(-Cu*(x*x + y*y))*cos(2*pi*t);
}

double v(double x,double y, double t) {
    return Cv*exp(-Cv*(x*x + (y - 0.5)*(y - 0.5)))*cos(2*pi*t);
}

double u_ana(double x,double y) {
    return u(x,y,0.5);
}

double v_ana(double x,double y) {
    return v(x,y,0.5);
}
double ux(double x,double y,double t) {
    return -2*Cu*x*u(x,y,t);
}

double ux_ana(double x,double y) {
    return ux(x,y,0.5);
}
double uxx(double x,double y,double t) {
    return -2*Cu*u(x,y,t) - 2*Cu*x*ux(x,y,t);
}

double uy(double x,double y,double t) {
    return -2*Cu*y*u(x,y,t);
}

double uy_ana(double x,double y) {
    return uy(x,y,0.5);
}

double uyy(double x,double y,double t) {
    return -2*Cu*u(x,y,t) -2*Cu*y*uy(x,y,t);
}

double vx(double x,double y,double t) {
    return -2*Cv*x*v(x,y,t);
}

double vx_ana(double x,double y) {
    return vx(x,y,0.5);
}

double vxx(double x,double y,double t) {
    return -2*Cv*v(x,y,t) -2*Cv*x*vx(x,y,t);
}

double vy(double x,double y,double t) {
    return -2*Cv*(y - 0.5)*v(x,y,t);
}

double vy_ana(double x,double y) {
    return vy(x,y,0.5);
}

double vyy(double x,double y,double t) {
    return -2*Cv*v(x,y,t) - 2*Cv*(y - 0.5)*vy(x,y,t);
}

double ut(double x,double y,double t) {
    return Cu*exp(-Cu*(x*x + y*y))*(-2*pi)*sin(2*pi*t);
}

double vt(double x,double y,double t) {
    return Cv*exp(-Cv*(x*x + (y - 0.5)*(y - 0.5)))*(-2*pi)*sin(2*pi*t);
}

double f(double x,double y,double t) {
    return ut(x,y,t) - uxx(x,y,t) - uyy(x,y,t) + ux(x,y,t)*vx(x,y,t) + u(x,y,t)*vxx(x,y,t) + uy(x,y,t)*vy(x,y,t) + u(x,y,t)*vyy(x,y,t);
}

double g(double x,double y,double t) {
    return vt(x,y,t) - vxx(x,y,t) - vyy(x,y,t) - u(x,y,t) + v(x,y,t);
}
double test(double x,double y){
    return (-1 + 1.0/16*8*pi*pi)*sin(2*pi*x)*sin(2*pi*y);
}
int main() {
    ofstream fout("result1_plt.txt");
    //一致，三角剖分，线性元***********************************************************
    //生成信息矩阵P,T,Pb,Tb
    double right = 0.5,top = 0.5,left = -0.5,bottom = -0.5;
    int N1 = 2,N2 = N1;
    double h = 1.0/N1;
    cout<< "h" << h << endl;
    double dt = pow(h,2);
    int maxstep = 0.5/dt + 1;
    cout << "dt" << dt << endl;
    fout << "C:" << Cu << endl;
    fout << "N:" << N1 << endl;
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

    cout << g(-0.5,-0.5,0.05) << "," << g(-0.5,0,0.05) << "," << g(-0.5,0.5,0.05) << endl;
    cout << g(0,-0.5,0.05)<< "," << g(0,0,0.05) << "," << g(0,0.5,0.05) <<endl;
    cout << g(0.5,-0.5,0.05) << "," << g(0.5,0,0.05) << "," << g(0.5,0.5,0.05) << endl; 
    //组装M，L
    SpMat M(Nb,Nb);
    FEKernel::assembleA_2D(N,Nlb,Nlb,basis_type,basis_type,P,T,Tb,Tb,M,function_one,0,0,0,0);
    SpMat L1(Nb,Nb),L2(Nb,Nb);
    FEKernel::assembleA_2D(N,Nlb,Nlb,basis_type,basis_type,P,T,Tb,Tb,L1,function_one,1,0,1,0);
    FEKernel::assembleA_2D(N,Nlb,Nlb,basis_type,basis_type,P,T,Tb,Tb,L2,function_one,0,1,0,1);
    SpMat L = L1 + L2;
    //初始化迭代向量
    Vector xn(2*Nb);
    for(int i = 0; i < Nb; ++i) {
        xn(i) = u0(Pb(0,i),Pb(1,i));
        xn(i + Nb) = c0(Pb(0,i),Pb(1,i));
    }
    Vector xn1 = xn;

    Vector  x(2*Nb);

    Vector btop(Nb);
    Vector bbottom(Nb);
    Vector bright(2*Nb);
    Vector b(2*Nb);
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
    //迭代
    dt = 0.025;
    Eigen::SparseLU<SpMat> solver;
    for(int i = 1; i < 3; ++i) {
        double t = i*dt;
        cout << "t" << t << endl;
        //取出vh这个向量
        Vector xlin = (2*xn - xn1).tail(Nb);
        //组装出K矩阵
        SpMat K1(Nb,Nb),K2(Nb,Nb);
        FEKernel::assembleFE_coe_2D(1,0,basis_type,xlin,N,Nlb,Nlb,basis_type,basis_type,P,T,Tb,Tb,K1,coefunction,1,0,0,0);
        FEKernel::assembleFE_coe_2D(0,1,basis_type,xlin,N,Nlb,Nlb,basis_type,basis_type,P,T,Tb,Tb,K2,coefunction,0,1,0,0);
        SpMat K = K2 + K1;

        //组装总矩阵
        SpMat A11 = M + dt*(L - K);
        SpMat A21 = -dt*M;
        SpMat A22 = M + dt*(L + M);
        SpMat A(2*Nb,2*Nb);
        FEKernel::addSpMat(A,A11,0,0);
        FEKernel::addSpMat(A,A21,Nb,0);
        FEKernel::addSpMat(A,A22, Nb, Nb);

        
        //计算右端项
        b.head(Nb) = M*xn.head(Nb);
        b.tail(Nb) = M*xn.tail(Nb);
        
        FEKernel::assembleV_2D_t(t,N,Nlb,basis_type,0,0,P,T,Pb,Tb,btop,f);
        FEKernel::assembleV_2D_t(t,N,Nlb,basis_type,0,0,P,T,Pb,Tb,bbottom,g);
        
        //cout  << "-----------" << endl << bbottom << endl;
        //cout <<"f:" << f(0,0,0.025) << endl;
        //cout <<"g"<< g(0,0,0.025)<<endl;

        bright.head(Nb) = btop*dt;
        bright.tail(Nb) = bbottom*dt;
        b += bright;
        
        //求解代数系统
        A.makeCompressed();
        solver.compute(A);
        x = solver.solve(b);
        xn1 = xn;
        xn = x;
        
        cout << x.tail(Nb) << endl << "-----------------------------"<< endl;
    }
    
    //计算终止时刻的误差
    double t_end = (maxstep-1)*dt;
    cout << "结束时间为: " << t_end << endl;



    btop = x.head(Nb);
    bbottom = x.tail(Nb);

    /*******************计算误差*************************/
    //最大模范数
    double error_max_u = FEKernel::compute_Max_error_2D(N,P,T,Tb,basis_type,btop,u_ana);
    double error_max_v = FEKernel::compute_Max_error_2D(N,P,T,Tb,basis_type,bbottom,v_ana);
    fout << "error_max_u:" << error_max_u << endl;
    fout << "error_max_v:" << error_max_v << endl;
    //2范数
    int der_y=0,der_x = 0;
    double erroru = FEKernel::compute_Hs_error_2D(N,P,T,Tb,basis_type,der_x,der_y,btop,u_ana);
    double errorv = FEKernel::compute_Hs_error_2D(N,P,T,Tb,basis_type,der_x,der_y,bbottom,v_ana);
    fout << "error_L2_u:"<<erroru <<endl;
    fout << "error_L2_v:" << errorv << endl;
    //H1半范
    double error3 = FEKernel::compute_Hs_error_2D(N,P,T,Tb,basis_type,1,0,btop,ux_ana);//传入解析解的一阶偏导数
    double error31 = FEKernel::compute_Hs_error_2D(N,P,T,Tb,basis_type,0,1,btop,uy_ana);//传入解析解的一阶偏导数
    double error4 = FEKernel::compute_Hs_error_2D(N,P,T,Tb,basis_type,1,0,bbottom,vx_ana);//传入解析解的一阶偏导数
    double error41 = FEKernel::compute_Hs_error_2D(N,P,T,Tb,basis_type,0,1,bbottom,vy_ana);//传入解析解的一阶偏导数
    double error_H_u = sqrt(error3*error3+error31*error31);
    double error_H_v = sqrt(error4*error4+error41*error41);
    fout << "errorH1_u:" <<error_H_u <<endl;
    fout << "errorH1_v:" <<error_H_v <<endl;

    //作图
    //数值解图像
    for(int i = 0; i < N1+1; i++) {
        vector<double> z_mark;
        for(int j = 0; j < N1+1; j++) {
            int k = (N1 + 1)*i + j;
            z_mark.push_back(btop(k));
        }
        Z.push_back(z_mark);
    }

    plt::plot_surface(X,Y,Z);
    plt::title("numerical solution u,t=0.5");
    plt::xlabel("x");
    plt::ylabel("y");
    plt::set_zlabel("numerical_u");
	const char* filename1="./1.png";
	plt::save(filename1);
    Z.clear();

    for(int i = 0; i < N1+1; i++) {
        vector<double> x_mark,y_mark,z_mark;
        for(int j = 0; j < N1+1; j++) {
            int k = (N1 + 1)*i + j;
            z_mark.push_back(bbottom(k));
        }
        Z.push_back(z_mark);
    }
    plt::plot_surface(X,Y,Z);
    plt::title("numerical solution v,t=0.5");
    plt::xlabel("x");
    plt::ylabel("y");
    plt::set_zlabel("numerical_v");
	const char* filename2="./2.png";
	plt::save(filename2);
    Z.clear();

    for(int i = 0; i < N1+1; i++) {
        vector<double> x_mark,y_mark,z_mark;
        for(int j = 0; j < N1+1; j++) {
            int k = (N1 + 1)*i + j;
            z_mark.push_back(abs(btop(k) - u_ana(Pb(0,k),Pb(1,k))));
        }
        Z.push_back(z_mark);
    }
    plt::plot_surface(X,Y,Z);
    plt::title("solution error u,t=0.5");
    plt::xlabel("x");
    plt::ylabel("y");
    plt::set_zlabel("error_u");
	const char* filename3="./3.png";
	plt::save(filename3);
    Z.clear();

    for(int i = 0; i < N1+1; i++) {
        vector<double> x_mark,y_mark,z_mark;
        for(int j = 0; j < N1+1; j++) {
            int k = (N1 + 1)*i + j;
            z_mark.push_back(abs(bbottom(k) - v_ana(Pb(0,k),Pb(1,k))));
        }
        Z.push_back(z_mark);
    }
    plt::plot_surface(X,Y,Z);
    plt::title("solution error v,t=0.5");
    plt::xlabel("x");
    plt::ylabel("y");
    plt::set_zlabel("error_v");
	const char* filename4="./4.png";
	plt::save(filename4);
    Z.clear();

    for(int i = 0; i < N1+1; i++) {
        vector<double> x_mark,y_mark,z_mark;
        for(int j = 0; j < N1+1; j++) {
            int k = (N1 + 1)*i + j;
            z_mark.push_back( u_ana(Pb(0,k),Pb(1,k)));
        }
        Z.push_back(z_mark);
    }
    plt::plot_surface(X,Y,Z);
    plt::title("exact solution u,t=0.5");
    plt::xlabel("x");
    plt::ylabel("y");
    plt::set_zlabel("u");
	const char* filename5="./5.png";
	plt::save(filename5);
    Z.clear();


    for(int i = 0; i < N1+1; i++) {
        vector<double> x_mark,y_mark,z_mark;
        for(int j = 0; j < N1+1; j++) {
            int k = (N1 + 1)*i + j;
            z_mark.push_back(v_ana(Pb(0,k),Pb(1,k)));
        }
        Z.push_back(z_mark);
    }
    plt::plot_surface(X,Y,Z);
    plt::title("exact solution v,t=0.5");
    plt::xlabel("x");
    plt::ylabel("y");
    plt::set_zlabel("v");
	const char* filename6="./6.png";
	plt::save(filename6);
    fout.close();

}
