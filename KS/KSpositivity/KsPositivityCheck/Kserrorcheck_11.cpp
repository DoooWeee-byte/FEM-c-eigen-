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
double dt = pow(10,-3);

double Cu = 70;
double Cv = Cu;
double endtime = 0.5;
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
    return u(x,y,endtime);
}

double v_ana(double x,double y) {
    return v(x,y,endtime);
}
double ux(double x,double y,double t) {
    return -2*Cu*x*u(x,y,t);
}

double ux_ana(double x,double y) {
    return ux(x,y,endtime);
}
double uxx(double x,double y,double t) {
    return -2*Cu*u(x,y,t) - 2*Cu*x*ux(x,y,t);
}

double uy(double x,double y,double t) {
    return -2*Cu*y*u(x,y,t);
}

double uy_ana(double x,double y) {
    return uy(x,y,endtime);
}

double uyy(double x,double y,double t) {
    return -2*Cu*u(x,y,t) -2*Cu*y*uy(x,y,t);
}

double vx(double x,double y,double t) {
    return -2*Cv*x*v(x,y,t);
}

double vx_ana(double x,double y) {
    return vx(x,y,endtime);
}

double vxx(double x,double y,double t) {
    return -2*Cv*v(x,y,t) -2*Cv*x*vx(x,y,t);
}

double vy(double x,double y,double t) {
    return -2*Cv*(y - 0.5)*v(x,y,t);
}

double vy_ana(double x,double y) {
    return vy(x,y,endtime);
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

int main() {
    ofstream fout("checkerrorI.txt");
    //一致，三角剖分，线性元***********************************************************
    //生成信息矩阵P,T,Pb,Tb
    double right = 0.5,top = 0.5,left = -0.5,bottom = -0.5;
    int N1 = 16,N2 = N1;
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
    for(int i = 0; i < Nb; ++i) {
        tripletList.push_back(Tk(i,i,Mtemp(i)));
    }
    M.setFromTriplets(tripletList.begin(), tripletList.end());


    //初始化迭代向量
    Vector xn(2*Nb);
    for(int i = 0; i < Nb; ++i) {
        xn(i) = u0(Pb(0,i),Pb(1,i));
        xn(i + Nb) = c0(Pb(0,i),Pb(1,i));
    }
    Vector xn1 = xn;


    fout << "time:" << 0 <<"  "<<  "min:" << xn.minCoeff() << endl;


    //解向量初始化
    Vector  x(2*Nb);
    x *= 0;
    //定义两个中间向量，接受U和v的解
    Vector xu(Nb);
    Vector xv(Nb);
    Vector b(2*Nb),bright(2*Nb);
    //迭代
    Eigen::SparseLU<SpMat> solver;
    for(int i = 1; i < maxstep; ++i) {
        double t = i*dt;
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
        b.head(Nb) = M*xn.head(Nb);
        b.tail(Nb) = M*xn.tail(Nb);
        FEKernel::assembleV_2D_t(t,N,Nlb,basis_type,0,0,P,T,Pb,Tb,xu,f);
        FEKernel::assembleV_2D_t(t,N,Nlb,basis_type,0,0,P,T,Pb,Tb,xv,g);
        bright.head(Nb) = xu*dt;
        bright.tail(Nb) = xv*dt;
        b += bright;

        //求解代数系统

        solver.compute(A11);
        x = solver.solve(b);

        //为下一次迭代作准备
        xn1 = xn;
        xn = x;
    }

//计算终止时刻的误差
    double t_end = (maxstep-1)*dt;
    cout << "结束时间为: " << t_end << endl;



    xu = x.head(Nb);
    xv = x.tail(Nb);

    /*******************计算误差*************************/
    //最大模范数
    double error_max_u = FEKernel::compute_Max_error_2D(N,P,T,Tb,basis_type,xu,u_ana);
    double error_max_v = FEKernel::compute_Max_error_2D(N,P,T,Tb,basis_type,xv,v_ana);
    fout << "error_max_u:" << error_max_u << endl;
    fout << "error_max_v:" << error_max_v << endl;
    //2范数
    int der_y=0,der_x = 0;
    double erroru = FEKernel::compute_Hs_error_2D(N,P,T,Tb,basis_type,der_x,der_y,xu,u_ana);
    double errorv = FEKernel::compute_Hs_error_2D(N,P,T,Tb,basis_type,der_x,der_y,xv,v_ana);
    fout << "error_L2_u:"<<erroru <<endl;
    fout << "error_L2_v:" << errorv << endl;
    //H1半范
    double error3 = FEKernel::compute_Hs_error_2D(N,P,T,Tb,basis_type,1,0,xu,ux_ana);//传入解析解的一阶偏导数
    double error31 = FEKernel::compute_Hs_error_2D(N,P,T,Tb,basis_type,0,1,xu,uy_ana);//传入解析解的一阶偏导数
    double error4 = FEKernel::compute_Hs_error_2D(N,P,T,Tb,basis_type,1,0,xv,vx_ana);//传入解析解的一阶偏导数
    double error41 = FEKernel::compute_Hs_error_2D(N,P,T,Tb,basis_type,0,1,xv,vy_ana);//传入解析解的一阶偏导数
    double error_H_u = sqrt(error3*error3+error31*error31);
    double error_H_v = sqrt(error4*error4+error41*error41);
    fout << "errorH1_u:" <<error_H_u <<endl;
    fout << "errorH1_v:" <<error_H_v <<endl;

    //画图
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

    //数值解图像
    for(int i = 0; i < N1+1; i++) {
        vector<double> z_mark;
        for(int j = 0; j < N1+1; j++) {
            int k = (N1 + 1)*i + j;
            z_mark.push_back(xu(k));
        }
        Z.push_back(z_mark);
    }

    plt::plot_surface(X,Y,Z);
    plt::title("numerical solution u,t=0.5");
    plt::xlabel("x");
    plt::ylabel("y");
    plt::set_zlabel("numerical_u");
    const char* filename1="./11.png";
    plt::save(filename1);
    Z.clear();

    for(int i = 0; i < N1+1; i++) {
        vector<double> x_mark,y_mark,z_mark;
        for(int j = 0; j < N1+1; j++) {
            int k = (N1 + 1)*i + j;
            z_mark.push_back(xv(k));
        }
        Z.push_back(z_mark);
    }
    plt::plot_surface(X,Y,Z);
    plt::title("numerical solution v,t=0.5");
    plt::xlabel("x");
    plt::ylabel("y");
    plt::set_zlabel("numerical_v");
    const char* filename2="./22.png";
    plt::save(filename2);
    Z.clear();

    for(int i = 0; i < N1+1; i++) {
        vector<double> x_mark,y_mark,z_mark;
        for(int j = 0; j < N1+1; j++) {
            int k = (N1 + 1)*i + j;
            z_mark.push_back(abs(xu(k) - u_ana(Pb(0,k),Pb(1,k))));
        }
        Z.push_back(z_mark);
    }
    plt::plot_surface(X,Y,Z);
    plt::title("solution error u,t=0.5");
    plt::xlabel("x");
    plt::ylabel("y");
    plt::set_zlabel("error_u");
    const char* filename3="./33.png";
    plt::save(filename3);
    Z.clear();

    for(int i = 0; i < N1+1; i++) {
        vector<double> x_mark,y_mark,z_mark;
        for(int j = 0; j < N1+1; j++) {
            int k = (N1 + 1)*i + j;
            z_mark.push_back(abs(xv(k) - v_ana(Pb(0,k),Pb(1,k))));
        }
        Z.push_back(z_mark);
    }
    plt::plot_surface(X,Y,Z);
    plt::title("solution error v,t=0.5");
    plt::xlabel("x");
    plt::ylabel("y");
    plt::set_zlabel("error_v");
    const char* filename4="./44.png";
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
    const char* filename5="./55.png";
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
    const char* filename6="./66.png";
    plt::save(filename6);
    plt::show();

    /* //文件写入结束，关闭文件 */
    fout.close();







}
