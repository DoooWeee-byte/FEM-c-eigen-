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
double nv = 1;
double rho = 1000;
double mu = 0.001;
double func_nv(double x,double y) {
    return nv;
}
double f(double x,double y){
	return 1;
}
double func_negtiveone(double x,double y) {
    return -1;
}
double c(double x,double y) {
    double r = 1;
    return  r;
}
double g1(double x,double y)
{
    double r = 0;
    if(x == 0) {
        r = 0;
    } else if(x == 3) {
        r = 0;
    } else if(y == -1) {
        r = 0.5*sin(pi*x);
    } else if(y == 1) {
        r = 0.5*sin(pi*x);
    }
    return r;
}
double g2(double x,double y) {
    double r = 0;
    if(x == 0) {
        r = -3*sin(pi/3*y);
    } else if(x == 3) {
        r = 3*sin(pi/3*y);
    } else if(y == -1) {
        r = 3*sqrt(3)/2*cos(pi*x);
    } else if(y == 1) {
        r = -3*sqrt(3)/2*cos(pi*x);
    }
    return r;
}


double u1(double x,double y)
{
    return sin(pi*x)*cos(pi/3*y);
}
double u2(double x,double y)
{
    return -3*cos(pi*x)*sin(pi/3*y);
}

double u1_1_der_x(double x,double y) {
    return pi*cos(pi*x)*cos(pi/3*y);
}


double u1_1_der_y(double x,double y) {
    return -pi/3*sin(pi*x)*sin(pi/3*y);
}

double u1_2_der_y(double x,double y){
	return -pi*pi/9*sin(pi*x)*cos(pi/3*y);
}

double f1(double x,double y){
	return rho*(u1(x,y)*u1_1_der_x(x,y) + u1_1_der_y(x,y)*u2(x,y)) - mu*u1_2_der_y(x,y);
}

double u2_1_der_x(double x, double y){
	return 3*pi*sin(pi*x)*sin(pi/3*y);
}

double u2_1_der_y(double x, double y){
	return -pi*cos(pi*x)*cos(pi/3*y);
}

int main() {

	ofstream fout("result_burger.txt");
   //一致，三角剖分，线性元***********************************************************
    //生成信息矩阵P,T,Pb,Tb
    double right = 3,top = 1,left = 0,bottom = -1;
	int N1 = 33,N2 = 22;
	/* cout << "StepH: " <<  1.0/N1 << endl; */
	fout << "StepH: " <<  1.0/N1 << endl;
    int N = 2*N1*N2;
    int Nb_u = (N1+1)*(N2+1)+N1*(N2+1)+N2*(N1+1)+N1*N2; //总基函数节点个数 */
	int N1_basis_u = 2*N1 + 1;
	int N2_basis_u = 2*N2 + 1;
    
    
    Matrix P(2,Nb_u),T(3,N);
    FEKernel::generate_PT_2D(P,T,N1,N2,right,top,bottom,left);

    int basis_type_u = 202;
    int nlb_u = 6;
 
    Matrix Pb_u(2,Nb_u),Tb_u(nlb_u,N);
    FEKernel::generate_PbTb_2D(Pb_u,Tb_u,N1,N2,basis_type_u,right,top,bottom,left);
    //生成边界矩阵，边界点
    int N1_basis = 2*N1; //基函数节点个数N1方向 */
    int N2_basis = 2*N2; //基函数节点个数N2方向
    int nbn = 2*(N1_basis+N2_basis);
    Matrix boundarynodes(2,nbn);
    FEKernel::generate_boundarynodes_2D(boundarynodes,N1_basis,N2_basis);
    
    Vector V1(Nb_u),V2(Nb_u),VB1(Nb_u),VB2(Nb_u);
	V1 *= 0;
	V2 *= 0;
	
    //组装A矩阵和V向量
	SpMat A4(Nb_u,Nb_u),A6(Nb_u,Nb_u),A7(Nb_u,Nb_u);
	FEKernel::assembleA_2D(N,nlb_u,nlb_u,basis_type_u,basis_type_u,P,T,Tb_u,Tb_u,A4,func_nv,0,1,0,1);
	FEKernel::assembleA_2D(N,nlb_u,nlb_u,basis_type_u,basis_type_u,P,T,Tb_u,Tb_u,A6,func_nv,0,0,1,0);
	FEKernel::assembleA_2D(N,nlb_u,nlb_u,basis_type_u,basis_type_u,P,T,Tb_u,Tb_u,A7,func_nv,0,0,0,1);

	//source
	Vector Vc(Nb_u);
	Vc *= 0;
    FEKernel::assembleV_2D(N,nlb_u,basis_type_u,0,0,P,T,Pb_u,Tb_u,Vc,f1);
    //牛顿迭代
    int L = 50; //牛顿迭代次数
    double xi = 0.000001;

    for(int i = 0; i < L; ++i) {
        SpMat AN(2*Nb_u,2*Nb_u);
        Vector VN(2*Nb_u);
        if(1) {
            int basis_type_coe = 202;
            SpMat A1(Nb_u,Nb_u),A2(Nb_u,Nb_u),A3(Nb_u,Nb_u),A5(Nb_u,Nb_u);
			
			//前两个是test,后两个是trial
            FEKernel::assembleFE_coe_2D(1,0,basis_type_coe,V1,N,nlb_u,nlb_u,basis_type_u,basis_type_u,P,T,Tb_u,Tb_u,A1,c,0,0,0,0);
            FEKernel::assembleFE_coe_2D(0,0,basis_type_coe,V1,N,nlb_u,nlb_u,basis_type_u,basis_type_u,P,T,Tb_u,Tb_u,A2,c,0,0,1,0);
            FEKernel::assembleFE_coe_2D(0,0,basis_type_coe,V1,N,nlb_u,nlb_u,basis_type_u,basis_type_u,P,T,Tb_u,Tb_u,A3,c,0,0,0,1);
            FEKernel::assembleFE_coe_2D(0,1,basis_type_coe,V2,N,nlb_u,nlb_u,basis_type_u,basis_type_u,P,T,Tb_u,Tb_u,A5,c,0,0,0,0);
		   
			AN = AN*0;
			SpMat Atl = rho*A1 + rho*A2 + rho*A3 + A4;
			FEKernel::addSpMat(AN,Atl ,0,0);
			SpMat Atr = rho*A5;
			FEKernel::addSpMat(AN, Atr, 0, Nb_u);
			FEKernel::addSpMat(AN, A6, Nb_u, 0);
			FEKernel::addSpMat(AN, A7, Nb_u, Nb_u);
        }
        if(1) {
            int basis_type_coe_1 = 202,basis_type_coe_2 = 202;
            Vector VN1(V1.size()),VN2(V1.size());

			//V1 和 V2 不要反了
            FEKernel::assembleFE_coe_2D_V(1,0,basis_type_coe_1,V1,0,0,basis_type_coe_2,V1,N,nlb_u,basis_type_u,0,0,P,T,Pb_u,Tb_u,VN1,f);
            FEKernel::assembleFE_coe_2D_V(0,1,basis_type_coe_1,V1,0,0,basis_type_coe_2,V2,N,nlb_u,basis_type_u,0,0,P,T,Pb_u,Tb_u,VN2,f);
            
			VN = VN*0;
			VN.head(V1.size()) = rho*VN1 + rho*VN2  + Vc;
        }
        
        /*******处理边界条件*************/
		SpMat AN_temp(AN.rows(),AN.cols());
        FEKernel::treat_boundary_Dirichlet_2d_SpMat_sp(boundarynodes,AN,AN_temp,Pb_u,VN,g1,g2); 
        /**********求解，并更新解***********************/
		Eigen::SparseLU<SpMat> solver;
		solver.compute(AN_temp);
		if(solver.info() != Eigen::Success){
			cout << "decomposition failed!" << endl;
		}
		Vector x;
		x = solver.solve(VN);
		cout << "solve complete" << endl;
        VB1 = V1;
        VB2 = V2;
        V1 = x.head(V1.size());
        V2 = x.segment(V1.size(),V2.size());
        double error = sqrt(pow((V1-VB1).norm(),2)+pow((V2-VB2).norm(),2));
		cout << "i:" << i << endl << "error:" << error << endl;
        if(error < xi) {
            /* cout << "NewtonIterationNumber: " << i <<endl; */
            /* cout << "IterationError: " << error << endl; */
			fout << "NewtonIterationNumber: " << i <<endl;
            fout << "IterationError: " << error << endl;
            break;
        }
        if(i == L -1) {
			fout << "MaxNewtonIterationNumber: " << i <<endl;
            fout << "IterationError: " << error << endl;
            cout << error << endl;
        }
    }
    
	//max_norm
	double error_max_1 = FEKernel::compute_Max_error_2D(V1,Pb_u,u1);
	double error_max_2 = FEKernel::compute_Max_error_2D(V2,Pb_u,u2);
	double error_max_u = max(error_max_1,error_max_2);
	fout << "error_max_u:" << error_max_u << endl;

    //2范数:der_x = der_y = 0,
    int der_x=0;
    int der_y=0;
    double error1 = FEKernel::compute_Hs_error_2D(N,P,T,Tb_u,basis_type_u,der_x,der_y,V1,u1);
    double error2 = FEKernel::compute_Hs_error_2D(N,P,T,Tb_u,basis_type_u,der_x,der_y,V2,u2);
    double error = sqrt(error1*error1+error2*error2);
    fout << "error_L2_u:"<<error <<endl;

    //H1半范
    double error3 = FEKernel::compute_Hs_error_2D(N,P,T,Tb_u,basis_type_u,1,0,V1,u1_1_der_x);//传入解析解的一阶偏导数
    double error31 = FEKernel::compute_Hs_error_2D(N,P,T,Tb_u,basis_type_u,0,1,V1,u1_1_der_y);//传入解析解的一阶偏导数
    double error4 = FEKernel::compute_Hs_error_2D(N,P,T,Tb_u,basis_type_u,1,0,V2,u2_1_der_x);//传入解析解的一阶偏导数
    double error41 = FEKernel::compute_Hs_error_2D(N,P,T,Tb_u,basis_type_u,0,1,V2,u2_1_der_y);//传入解析解的一阶偏导数
	double error_H_u = sqrt(error3*error3+error31*error31+error4*error4+error41*error41);
    fout << "errorH1_u:" <<error_H_u <<endl;
	fout.close();
}
