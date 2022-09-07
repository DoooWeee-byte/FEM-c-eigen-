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
double chi = 1;
double sigma_u = 2;
double rand_num = 0.5;
double epsilon = -pow(10,-5);
double Tol = pow(10,-6);
int maxstep_time = 10;
double dt = 0.01;

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
    return (1 - 0.01)*u(x,y);
}

double c0(double x,double y) {
    return (1 - 0.01)*c(x,y);
}

double function_one(double x,double y) {
    return 1;
}

double ux(double x,double y){
	return pi/20*(-cos(2*pi*y) + 1)*sin(2*pi*x);
}

double uxx(double x,double y){
	return pi*pi/10*(-cos(2*pi*y) + 1)*cos(2*pi*x);
}

double uy(double x,double y){
	return pi/20*(-cos(2*pi*x) + 1)*sin(2*pi*y);
}

double uyy(double x,double y){
	return pi*pi/10*(-cos(2*pi*x) + 1)*cos(2*pi*y);
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
	/* return -pi*pi/10*((cos(2*pi*y)-1)*cos(2*pi*x) + (cos(2*pi*x) - 1)*cos(2*pi*y)) - chi*pi*pi/800*(pow((cos(2*pi*y)-1)*sin(2*pi*x),2)+pow((cos(2*pi*x)-1)*sin(2*pi*y),2))+chi*pi*pi/20*((cos(2*pi*y)-1)*cos(2*pi*x) + (cos(2*pi*x) - 1)*cos(2*pi*y))*u(x,y); */
	return -(uxx(x,y) + uyy(x,y) - chi*(ux(x,y)*vx(x,y)) -chi*(uy(x,y)*vy(x,y))- chi*u(x,y)*(vxx(x,y) + vyy(x,y)));
}

double f2(double x,double y){
	/* return -pi*pi/20*((cos(2*pi*y)-1)*cos(2*pi*x) + (cos(2*pi*x) - 1)*cos(2*pi*y)) - c(x,y) + u(x,y); */
   return -(vxx(x,y) + vyy(x,y) - v(x,y) + u(x,y));	
}

int main() {
    //一致，三角剖分，线性元***********************************************************
    //生成信息矩阵P,T,Pb,Tb
    double right = 1,top = 1,left = 0,bottom = 0;
    int N1 = 8,N2 = N1;
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

    //组装M，L，S
    SpMat M(Nb,Nb),L(Nb,Nb);
    FEKernel::assembleA_2D(N,Nlb,Nlb,basis_type,basis_type,P,T,Tb,Tb,M,function_one,0,0,0,0);
    if(1) {
        SpMat L1(Nb,Nb),L2(Nb,Nb);
        FEKernel::assembleA_2D(N,Nlb,Nlb,basis_type,basis_type,P,T,Tb,Tb,L1,function_one,1,0,1,0);
        FEKernel::assembleA_2D(N,Nlb,Nlb,basis_type,basis_type,P,T,Tb,Tb,L1,function_one,0,1,0,1);
        L = L1 + L2;
    }

	//组装右端额外项
	Vector bextal(2*Nb);
	if(1){
		Vector tempv(Nb);
		FEKernel::assembleV_2D(N,Nlb,basis_type,0,0,P,T,Pb,Tb,tempv,f1);
		bextal.head(Nb) = tempv;
		FEKernel::assembleV_2D(N,Nlb,basis_type,0,0,P,T,Pb,Tb,tempv,f2);
		bextal.tail(Nb) = tempv;
	}

	Matrix A(2*Nb,2*Nb);
	A.bottomLeftCorner(Nb,Nb) = -dt*M;
	A.bottomRightCorner(Nb,Nb) = M + dt*(L+M);
	//初始化迭代向量
	Vector xn(2*Nb),xn1(2*Nb);
	if(1){
		Vector xn_temp1(Nb),xn_temp2(Nb);
		for(int i = 0; i < Nb;++i){
			xn_temp1(i) = u0(Pb(0,i),Pb(1,i));
			xn_temp2(i) = c0(Pb(0,i),Pb(1,i));
		}
		xn.head(Nb) = xn_temp1;
		xn.tail(Nb) = xn_temp2;
		xn1 = xn;
	}

	Vector bn(2*Nb);
	//时间迭代
	for(int i = 0;i < maxstep_time;++i){
		Vector xtemp1 = xn.head(Nb);
		Vector xtemp2 = xn.tail(Nb);
		bn.head(Nb) = M*xtemp1;
		bn.tail(Nb) = M*xtemp2;
		bn = bn +  bextal;
		/* Vector xlin = 2*xn - xn1; */
		Vector xlin = xn;
		//组装出K矩阵
		SpMat K(Nb,Nb),K1(Nb,Nb),K2(Nb,Nb);
		Vector clin = xlin.tail(Nb);
		FEKernel::assembleFE_coe_2D(1,0,basis_type,clin,N,Nlb,Nlb,basis_type,basis_type,P,T,Tb,Tb,K1,function_one,1,0,0,0);
		FEKernel::assembleFE_coe_2D(0,1,basis_type,clin,N,Nlb,Nlb,basis_type,basis_type,P,T,Tb,Tb,K2,function_one,0,1,0,0);
		K = (K2 + K1)*chi;
		A.topLeftCorner(Nb,Nb) = M + dt*(L - K);
		xn1 = xn;
		xn = A.colPivHouseholderQr().solve(bn);
		if((xn - analytic_solution).norm() < Tol){
			cout << "i:" << i << endl;
			break;
		}
		cout << "i:" << i << endl;
		cout << " L2norm: "<< (xn - analytic_solution).norm() << endl; 
      }

    

    //组装A矩阵和V向量
    /*if(1) {*/
    /*    SpMat A1(Nb_u,Nb_u),A2(Nb_u,Nb_u),A3(Nb_u,Nb_u),A4(Nb_u,Nb_u),A5(Nb_u,Nb_p),A6(Nb_u,Nb_p),A7(Nb_p,Nb_u),A8(Nb_p,Nb_u);*/
    /*    FEKernel::assembleA_2D(N,nlb_u,nlb_u,basis_type_u,basis_type_u,P,T,Tb_u,Tb_u,A1,func_nv,1,0,1,0);*/
    /*    FEKernel::assembleA_2D(N,nlb_u,nlb_u,basis_type_u,basis_type_u,P,T,Tb_u,Tb_u,A2,func_nv,0,1,0,1);*/
    /*    FEKernel::assembleA_2D(N,nlb_u,nlb_u,basis_type_u,basis_type_u,P,T,Tb_u,Tb_u,A3,func_nv,0,1,1,0);*/
    /*    FEKernel::assembleA_2D(N,nlb_u,nlb_u,basis_type_u,basis_type_u,P,T,Tb_u,Tb_u,A4,func_nv,1,0,0,1);*/
    /*    FEKernel::assembleA_2D(N,nlb_u,nlb_p,basis_type_u,basis_type_p,P,T,Tb_u,Tb_p,A5,func_negtivenoe,1,0,0,0);*/
    /*    FEKernel::assembleA_2D(N,nlb_u,nlb_p,basis_type_u,basis_type_p,P,T,Tb_u,Tb_p,A6,func_negtivenoe,0,1,0,0);*/
    /*    FEKernel::assembleA_2D(N,nlb_p,nlb_u,basis_type_p,basis_type_u,P,T,Tb_p,Tb_u,A7,func_negtivenoe,0,0,1,0);*/
    /*    FEKernel::assembleA_2D(N,nlb_p,nlb_u,basis_type_p,basis_type_u,P,T,Tb_p,Tb_u,A8,func_negtivenoe,0,0,0,1);*/


    /*    A.topLeftCorner(A1.rows(),A1.cols()) = A.topLeftCorner(A1.rows(),A1.cols()) + 2*A1+ A2;*/
    /*    A.topRightCorner(A5.rows(),A5.cols()) =  A5;*/
    /*    A.bottomLeftCorner(A7.rows(),A7.cols()) = A7;*/
    /*    A.block(0,A1.cols(),A3.rows(),A3.cols()) = A3;*/
    /*    A.block(A1.rows(),0,A4.rows(),A4.cols()) = A4;*/
    /*    A.block(A1.rows()+A4.rows(),A1.cols(),A8.rows(),A8.cols()) = A8;*/
    /*    A.block(A1.rows(),A1.cols(),A2.rows(),A2.cols()) += 2*A2 + A1;*/
    /*    A.block(A1.rows(),A1.cols()+A3.cols(),A6.rows(),A6.cols()) = A6;*/
    /*}*/
    /*FEKernel::assembleV_2D(N,nlb_u,basis_type_u,0,0,P,T,Pb_u,Tb_u,V1,f1);*/
    /*FEKernel::assembleV_2D(N,nlb_u,basis_type_u,0,0,P,T,Pb_u,Tb_u,V2,f2);*/
    /*V = V*0;*/
    /*V.head(V1.size()) = V1;*/
    /*V.segment(V1.size(),V2.size()) = V2;*/
    /*V1 = V1*0;*/
    /*V2 = V2*0;*/
    /*Vector VB1(V1.size()),VB2(V2.size()),VB3(V3.size());*/

    /*//牛顿迭代*/
    /*int L = 30; //牛顿迭代次数*/
    /*double xi = 0.000000000001;*/
    /*for(int i = 0; i < L; ++i) {*/
    /*    Matrix AN(A.rows(),A.cols());*/
    /*    AN = AN*0;*/
    /*    if(1) {*/
    /*        int basis_type_coe = 202;*/
    /*        SpMat A1(Nb_u,Nb_u),A2(Nb_u,Nb_u),A3(Nb_u,Nb_u),A4(Nb_u,Nb_u),A5(Nb_u,Nb_u),A6(Nb_u,Nb_u);*/
    /*        FEKernel::assembleFE_coe_2D(1,0,basis_type_coe,V1,N,nlb_u,nlb_u,basis_type_u,basis_type_u,P,T,Tb_u,Tb_u,A1,c,0,0,0,0);*/
    /*        FEKernel::assembleFE_coe_2D(0,0,basis_type_coe,V1,N,nlb_u,nlb_u,basis_type_u,basis_type_u,P,T,Tb_u,Tb_u,A2,c,0,0,1,0);*/
    /*        FEKernel::assembleFE_coe_2D(0,0,basis_type_coe,V2,N,nlb_u,nlb_u,basis_type_u,basis_type_u,P,T,Tb_u,Tb_u,A3,c,0,0,0,1);*/
    /*        FEKernel::assembleFE_coe_2D(0,1,basis_type_coe,V1,N,nlb_u,nlb_u,basis_type_u,basis_type_u,P,T,Tb_u,Tb_u,A4,c,0,0,0,0);*/
    /*        FEKernel::assembleFE_coe_2D(1,0,basis_type_coe,V2,N,nlb_u,nlb_u,basis_type_u,basis_type_u,P,T,Tb_u,Tb_u,A5,c,0,0,0,0);*/
    /*        FEKernel::assembleFE_coe_2D(0,1,basis_type_coe,V2,N,nlb_u,nlb_u,basis_type_u,basis_type_u,P,T,Tb_u,Tb_u,A6,c,0,0,0,0);*/
    /*        AN.topLeftCorner(A1.rows(),A1.cols()) = A1+ A2 +A3;*/
    /*        AN.block(0,A1.cols(),A4.rows(),A4.cols()) = A4;*/
    /*        AN.block(A1.rows(),0,A5.rows(),A5.cols()) = A5;*/
    /*        AN.block(A1.rows(),A1.cols(),A2.rows(),A2.cols()) = A2+A3+A6;*/
    /*    }*/
    /*    Vector VN(V.size());*/
    /*    VN = VN*0;*/
    /*    if(1) {*/
    /*        int basis_type_coe_1 = 202,basis_type_coe_2 = 202;*/
    /*        Vector VN1(V1.size()),VN2(V1.size()),VN3(V2.size()),VN4(V2.size());*/
    /*        FEKernel::assembleFE_coe_2D_V(0,0,basis_type_coe_1,V1,1,0,basis_type_coe_2,V1,N,nlb_u,basis_type_u,0,0,P,T,Pb_u,Tb_u,VN1,f);*/
    /*        FEKernel::assembleFE_coe_2D_V(0,0,basis_type_coe_1,V2,0,1,basis_type_coe_2,V1,N,nlb_u,basis_type_u,0,0,P,T,Pb_u,Tb_u,VN2,f);*/
    /*        FEKernel::assembleFE_coe_2D_V(0,0,basis_type_coe_1,V1,1,0,basis_type_coe_2,V2,N,nlb_u,basis_type_u,0,0,P,T,Pb_u,Tb_u,VN3,f);*/
    /*        FEKernel::assembleFE_coe_2D_V(0,0,basis_type_coe_1,V2,0,1,basis_type_coe_2,V2,N,nlb_u,basis_type_u,0,0,P,T,Pb_u,Tb_u,VN4,f);*/
    /*        VN.head(V1.size()) = VN1 + VN2;*/
    /*        VN.segment(V1.size(),V2.size()) = VN3 + VN4;*/
    /*    }*/
    /*    AN += A;*/
    /*    VN += V;*/
        /* *******处理边界条件************ */
    /*    FEKernel::treat_boundary_Dirichlet_2d(boundarynodes,AN,Pb_u,VN,g1,g2);*/
    /*    if(1) {*/
    /*        Vector x0(AN.cols());*/
    /*        x0 = x0*0;*/
    /*        x0(AN.cols()-1) = 1;*/
    /*        AN.row(AN.rows()-1) = x0.transpose();*/
    /*        VN(VN.size()-1) = analytic_solution_p(Pb_p(0,Pb_p.cols()-1),Pb_p(1,Pb_p.cols()-1));*/
    /*    }*/
    /*   ***********求解，并更新解************************/
    /*    x = AN.colPivHouseholderQr().solve(VN);*/
    /*    VB1 = V1;*/
    /*    VB2 = V2;*/
    /*    VB3 = V3;*/
    /*    V1 = x.head(V1.size());*/
    /*    V2 = x.segment(V1.size(),V2.size());*/
    /*    V3 = x.tail(V3.size());*/
    /*    double error = sqrt(pow((V1-VB1).norm(),2)+pow((V2-VB2).norm(),2));*/
    /*    if(error < xi) {*/
    /*        cout << "i:" << i <<endl;*/
    /*        cout << error << endl;*/
    /*        break;*/
    /*    }*/
    /*    if(i == L -1)*/
    /*        cout << error << endl;*/
    /*}*/




    /* /1* /2* //计算误差 *2/ *1/ */
    /* /1* /2* //组装解析解向量 *2/ *1/ */
    /* for(int i = 0; i != V1.size(); ++i) { */
    /*     x1[i] = u1(Pb(0,i),Pb(1,i)); */
    /*     x2[i] = u2(Pb(0,i),Pb(1,i)); */
    /* } */
    /* /1* cout << "analytic_solution with Eigen" << endl <<  x1 <<endl; *1/ */
    /* /1* /2* //计算无穷范数 *2/ *1/ */
    /*pass*/

    //2范数:der_x = der_y = 0,
    /* int der_x=0; */
    /* int der_y=0; */
    /* double error1 = FEKernel::compute_Hs_error_2D(N,P,T,Tb_u,basis_type_u,der_x,der_y,V1,u1); */
    /* double error2 = FEKernel::compute_Hs_error_2D(N,P,T,Tb_u,basis_type_u,der_x,der_y,V2,u2); */
    /* double error3 = FEKernel::compute_Hs_error_2D(N,P,T,Tb_p,basis_type_p,der_x,der_y,V3,analytic_solution_p); */
    /* double error = sqrt(error1*error1+error2*error2); */
    /* cout << "error_U_L2:"<<error <<endl; */
    /* cout << "error_V_L2:" << error3 << endl; */
    /* /1* //H1半范 *1/ */
    /* int der_x1 = 0,der_y1=1; */
    /* double error3 = FEKernel::compute_Hs_error_2D(N,P,T,Tb,basis_type_trial,der_x1,der_y1,V1,u1_1_der);//传入解析解的一阶偏导数 */
    /* double error4 = FEKernel::compute_Hs_error_2D(N,P,T,Tb,basis_type_trial,der_x1,der_y1,V2,u2_1_der);//传入解析解的一阶偏导数 */
    /* double error5 = sqrt(error3*error3+error4*error4); */
    /* cout << "errorH1:" <<error5 <<endl; */


}
