#ifndef GAUSSQUADRATRUE_H
#define GAUSSQUADRATRUE_H

#include <vector>
#include <cmath>
#include "basefunction.h"
using namespace std;
namespace PDWSC {
namespace FEMethod {

template<typename Vector,typename Matrix>
void generate_Gauss_reference_triangle(Vector &Gc,Matrix &Gp)
{
    int gpn = Gc.size;
    if(gpn==4) {
        Vector v = {(1-1/sqrt(3))/8,(1-1/sqrt(3))/8,(1+1/sqrt(3))/8,(1+1/sqrt(3))/8};
        Matrix m = {{(1/sqrt(3)+1)/2,(1/sqrt(3)+1)/2,(-1/sqrt(3)+1)/2,(-1/sqrt(3)+1)/2},{(1-1/sqrt(3))*(1+1/sqrt(3))/4,(1-1/sqrt(3))*(1-1/sqrt(3))/4,(1+1/sqrt(3))*(1+1/sqrt(3))/4,(1+1/sqrt(3))*(1-1/sqrt(3))/4}};
        Gc = v;
        Gp = m;
    }
    else if(gpn==9) {
        Vector v = {64/81*(1-0)/8,100/324*(1-sqrt(3/5))/8,100/324*(1-sqrt(3/5))/8,100/324*(1+sqrt(3/5))/8,100/324*(1+sqrt(3/5))/8,40/81*(1-0)/8,40/81*(1-0)/8,40/81*(1-sqrt(3/5))/8,40/81*(1+sqrt(3/5))/8};
        Matrix m = {{(1.0+0)/2,(1+sqrt(3/5))/2,(1+sqrt(3/5))/2,(1-sqrt(3/5))/2,(1-sqrt(3/5))/2,(1.0+0)/2,(1.0+0)/2,(1+sqrt(3/5))/2,(1-sqrt(3/5))/2},{(1.0-0)*(1+0)/4,(1-sqrt(3/5))*(1+sqrt(3/5))/4,(1-sqrt(3/5))*(1-sqrt(3/5))/4,(1+sqrt(3/5))*(1+sqrt(3/5))/4,(1+sqrt(3/5))*(1-sqrt(3/5))/4,(1-0)*(1+sqrt(3/5))/4,(1-0)*(1-sqrt(3/5))/4,(1-sqrt(3/5))*(1+0)/4,(1+sqrt(3/5))*(1+0)/4}};
        Gc = v;
        Gp = m;
    }
    else if(gpn==3) {
        Vector v = {1.0/6,1.0/6,1.0/6};
        Matrix m = {{1.0/2,1.0/2,0},{0,1.0/2,1.0/2}};
        Gc = v;
        Gp = m;
    }
    else {
        cout << "Wrong gpn!" <<endl;
    }
}

template<typename F,typename I,typename Matrix,typename Vector>
void generate_Gauss_local_triangle(Matrix &vertices,Vector &Glc,Matrix &Glp,Vector &Gc,Matrix &Gp)
{
    generate_Gauss_reference_triangle<Vector,Matrix>(Gc,Gp);
    F x1=vertices[0][0];
    F y1=vertices[1][0];
    F x2=vertices[0][1];
    F y2=vertices[1][1];
    F x3=vertices[0][2];
    F y3=vertices[1][2];
    F Jacobi=abs((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1));
    Glc = Gc*Jacobi;
    for(I i = 0; i != Glp.shape[1]; ++i) {
        Glp[0][i] = x1+(x2-x1)*Gp[0][i]+(x3-x1)*Gp[1][i];
        Glp[1][i] = y1+(y2-y1)*Gp[0][i]+(y3-y1)*Gp[1][i];
    }
}

template<typename F,typename I,typename Matrix,typename Vector>
F Gauss_quadratrue_M_2D(Matrix &vertices,I basis_type_trial,I basis_type_test,I alfa,I beta,I r,I s,I p,I q,I Gpn,F (*c)(F,F)) {
// alfa对应trial,beta对应test
    Vector Glc(Gpn);
    Matrix Glp(2,Gpn);
    Vector Gc(Gpn);
    Matrix Gp(2,Gpn);
    generate_Gauss_local_triangle<F,I,Matrix,Vector>(vertices,Glc,Glp,Gc,Gp);
    F ret = 0;
    for(I i = 0; i != Gpn; ++i) {
        //高斯权重*某个系数函数*test函数的某阶导数*trial函数的某阶导数
        ret += Glc[i]*c(Glp[0][i],Glp[1][i])*local_basis_2D(Glp[0][i],Glp[1][i],vertices,basis_type_trial,alfa,r,s)*local_basis_2D(Glp[0][i],Glp[1][i],vertices,basis_type_test,beta,p,q);
    }
    return ret;
}

template<typename F,typename I,typename Matrix,typename Vector>
F Gauss_quadratrue_V_2D(Matrix &vertices,I basis_type_test,I index,I Gpn,I p,I q,F (*f)(F,F,F),F t) {
    Vector Glc(Gpn);
    Matrix Glp(2,Gpn);
    Vector Gc(Gpn);
    Matrix Gp(2,Gpn);
    generate_Gauss_local_triangle<F,I,Matrix,Vector>(vertices,Glc,Glp,Gc,Gp);
    F int_value = 0;
    for(I i = 0; i != Gpn; ++i) {

        int_value += Glc[i]*f(Glp[0][i],Glp[1][i],t)*local_basis_2D(Glp[0][i],Glp[1][i],vertices,basis_type_test,index,p,q);
    }
    return int_value;
}

template<typename F,typename I,typename Matrix,typename Vector>
F Gauss_int_error_2D(Matrix &vertices,I basis_type,I der_x,I der_y,Vector &uh_local_vec,F (*analytic_solution)(F,F,F),F t,I Gpn) {
    Vector Glc(Gpn);
    Matrix Glp(2,Gpn);
    Vector Gc(Gpn);
    Matrix Gp(2,Gpn);
    generate_Gauss_local_triangle<F,I,Matrix,Vector>(vertices,Glc,Glp,Gc,Gp);
    F int_value = 0;
    for(auto i = 0; i != Gpn; ++i) {
        F zz = analytic_solution(Glp[0][i],Glp[1][i],t)- local_FE_function_2D<F,I,Matrix,Vector>(Glp[0][i],Glp[1][i],vertices,basis_type,der_x,der_y,uh_local_vec);
        int_value += Glc[i]*zz*zz;
    }
    return int_value;
}

//以下为处理边界条件的GS积分

template<typename F,typename I>
void generate_Gauss_reference_1D(I Gauss_point_number,vector<F> &Gc,vector<F> &Gp)
{
    if(Gauss_point_number==4) {
        Gc = {0.3478548451,0.3478548451,0.6521451549,0.6521451549};
        Gp = {0.8611363116,-0.8611363116,0.3399810436,-0.3399810436};
    }
    else if(Gauss_point_number==8) {
        Gc = {0.1012285363,0.1012285363,0.2223810345,0.2223810345,0.3137066459,0.3137066459,0.3626837834,0.3626837834};
        Gp = {0.9602898565,-0.9602898565,0.7966664774,-0.7966664774,0.5255324099,-0.5255324099,0.1834346425,-0.1834346425};
    }
    else {
        Gc = {1,1};
        Gp = {-1/sqrt(3),1/sqrt(3)};
    }
}

template<typename F,typename I>
void generate_Gauss_local_1D(F upper_bound,F lower_bound,vector<F> &Gc,vector<F> &Gp,vector<F> &Glc,vector<F> &Glp)
{
    for(auto i = 0; i != Gc.size(); ++i) {
        F lc = (upper_bound-lower_bound)*Gc[i]/2;
        F lp = (upper_bound-lower_bound)*Gp[i]/2+(upper_bound+lower_bound)/2;
        Glc.push_back(lc);
        Glp.push_back(lp);
    }
}


template<typename F,typename I,typename Matrix,typename Vector>
F Gauss_quadratrue_line_test_2D(Vector &end_point_1,Vector &end_point_2,Matrix &vertices,I basis_type,I basis_index,I der_x,I der_y,F (*f)(F,F),I Gpn) {
	F int_value = 0;
    if(end_point_1[1]==end_point_2[1]) {
        //水平线
        F lower_bound = min(end_point_1[0],end_point_2[0]);
        F upper_bound = max(end_point_1[0],end_point_2[0]);
        vector<F> Gc,Gp,Glc,Glp;
        generate_Gauss_reference_1D<F,I>(Gpn,Gc,Gp);
        generate_Gauss_local_1D<F,I>(upper_bound,lower_bound,Gc,Gp,Glc,Glp);
        for(I i = 0; i != Gpn; ++i) {
            int_value += Glc[i]*f(Glp[i],end_point_1[1])*local_basis_2D(Glp[i],end_point_1[1],vertices,basis_type,basis_index,der_x,der_y);
        }
    }
    else if(end_point_1[0]==end_point_2[0]) {
        //垂直线
        F lower_bound = min(end_point_1[1],end_point_2[1]);
        F upper_bound = max(end_point_1[1],end_point_2[1]);
        vector<F> Gc,Gp,Glc,Glp;
        generate_Gauss_reference_1D<F,I>(Gpn,Gc,Gp);
        generate_Gauss_local_1D<F,I>(upper_bound,lower_bound,Gc,Gp,Glc,Glp);
		for(I i=0;i != Gpn;++i){
			int_value += Glc[i]*f(end_point_1[0],Glp[i])*local_basis_2D(end_point_1[0],Glp[i],vertices,basis_type,basis_index,der_x,der_y);

		}
    }
	else{
		//直线斜率存在且不为0
		F lower_bound = min(end_point_1[0],end_point_2[0]);
        F upper_bound = max(end_point_1[0],end_point_2[0]);
        vector<F> Gc,Gp,Glc,Glp;
        generate_Gauss_reference_1D<F,I>(Gpn,Gc,Gp);
        generate_Gauss_local_1D<F,I>(upper_bound,lower_bound,Gc,Gp,Glc,Glp);
		F slope = (end_point_2[1]-end_point_1[1])/(end_point_2[0]-end_point_1[0]);
		F Jacobi = sqrt(slope*slope+1);
		for(I i=0;i != Gpn;++i){
			F x = Glp[i];
			F y = slope*(x - end_point_1[0]) + end_point_1[1];
			int_value += Glc[i]*Jacobi*f(x,y)*local_basis_2D(x,y,vertices,basis_type,basis_index,der_x,der_y);	
		}
	}
	return int_value; 
}

template<typename F,typename I,typename Matrix,typename Vector>
F Gauss_quadratrue_line_trial_test_2D(Vector &end_point_1,Vector &end_point_2,Matrix &vertices,I basis_type_trial,I basis_index_trial,I basis_type_test,I basis_index_test,I der_x_trial,I der_y_trial,I der_x_test,I der_y_test, F (*f)(F,F),I Gpn){
	F int_value = 0;
    if(end_point_1[1]==end_point_2[1]) {
        //水平线
        F lower_bound = min(end_point_1[0],end_point_2[0]);
        F upper_bound = max(end_point_1[0],end_point_2[0]);
        vector<F> Gc,Gp,Glc,Glp;
        generate_Gauss_reference_1D<F,I>(Gpn,Gc,Gp);
        generate_Gauss_local_1D<F,I>(upper_bound,lower_bound,Gc,Gp,Glc,Glp);
        for(I i = 0; i != Gpn; ++i) {
            int_value += Glc[i]*f(Glp[i],end_point_1[1])*local_basis_2D(Glp[i],end_point_1[1],vertices,basis_type_trial,basis_index_trial,der_x_trial,der_y_trial)*local_basis_2D(Glp[i],end_point_1[1],vertices,basis_type_test,basis_index_test,der_x_test,der_y_test);
        }
    }
    else if(end_point_1[0]==end_point_2[0]) {
        //垂直线
        F lower_bound = min(end_point_1[1],end_point_2[1]);
        F upper_bound = max(end_point_1[1],end_point_2[1]);
        vector<F> Gc,Gp,Glc,Glp;
        generate_Gauss_reference_1D<F,I>(Gpn,Gc,Gp);
        generate_Gauss_local_1D<F,I>(upper_bound,lower_bound,Gc,Gp,Glc,Glp);
		for(I i=0;i != Gpn;++i){
			int_value += Glc[i]*f(end_point_1[0],Glp[i])*local_basis_2D(end_point_1[0],Glp[i],vertices,basis_type_trial,basis_index_trial,der_x_trial,der_y_trial)*local_basis_2D(end_point_1[0],Glp[i],vertices,basis_type_test,basis_index_test,der_x_test,der_y_test);
		}
    }
	else{
		//直线斜率存在且不为0
		F lower_bound = min(end_point_1[0],end_point_2[0]);
        F upper_bound = max(end_point_1[0],end_point_2[0]);
        vector<F> Gc,Gp,Glc,Glp;
        generate_Gauss_reference_1D<F,I>(Gpn,Gc,Gp);
        generate_Gauss_local_1D<F,I>(upper_bound,lower_bound,Gc,Gp,Glc,Glp);
		F slope = (end_point_2[1]-end_point_1[1])/(end_point_2[0]-end_point_1[0]);
		F Jacobi = sqrt(slope*slope+1);
		for(I i=0;i != Gpn;++i){
			F x = Glp[i];
			F y = slope*(x - end_point_1[0]) + end_point_1[1];
			int_value += Glc[i]*Jacobi*f(x,y)*local_basis_2D(x,y,vertices,basis_type_trial,basis_index_trial,der_x_trial,der_y_trial)*local_basis_2D(x,y,vertices,basis_type_test,basis_index_test,der_x_test,der_y_test);	
		}
	}
	return int_value; 	
} 


} // end of namespace FEMethod
} // end of namespace PDWSC
#endif
