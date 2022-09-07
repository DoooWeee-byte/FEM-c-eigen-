#ifndef COMPUTEERROR_H
#define COMPUTEERROR_H

#include "GaussQuadratrue.h"
namespace PDWSC {
namespace FEMethod {


template<typename Matrix,typename Vector>
double compute_Hs_error_2D(int N,Matrix &P,Matrix &T,Matrix &Tb,int basis_type,int der_x,int der_y,Vector &solution_vec,double (*analytic_solution)(double,double),int Gpn){
	double error = 0;
	for(int n = 1;n != N+1;++n){
//取得端点坐标放到vertices中
		Matrix vertices(P.rows(),T.rows());
		for(int j = 0; j != T.rows();++j){
			int j1 = T(j,n-1);
			for(int i = 0;i != P.rows();++i){
				vertices(i,j) = P(i,j1-1);
			}
		}
//uh_local_vec代表第n个单元上的局部解向量
//取出Tb的第n列的值就是局部解向量在全局解向量中的位置
		Vector uh_local_vec(Tb.rows());
		int row = Tb.rows();
		for(int i=0;i < row;++i){
			int j = Tb(i,n-1) - 1;
			uh_local_vec(i) = solution_vec(Tb(i,n-1) - 1);
		}
		double int_value = Gauss_int_error_2D<double,int,Matrix,Vector>(vertices,basis_type,der_x,der_y,uh_local_vec,analytic_solution,Gpn);
		error += int_value;
	}
	error = sqrt(error);
	return error;
}

template<typename Matrix,typename Vector,typename F,typename I>
F compute_Max_error_2D(const Vector &solution_vec,const Matrix &Pb,F (*analytic_solution)(F,F)){
	Vector analytic_solution_vec(solution_vec.size());
	analytic_solution_vec *= 0;
	for(I i = 0;i < Pb.cols();i++){
		analytic_solution_vec(i) = analytic_solution(Pb(0,i),Pb(1,i));
	}
	F result = 0;
	analytic_solution_vec -= solution_vec;
	for(I i = 0;i < solution_vec.size();++i){
		if(abs(analytic_solution_vec(i)) > result)
			result = abs(analytic_solution_vec(i));	
	}
	return result;
}

}  // end of namespace FEMethod
}  // end of namespace PDWSC
#endif
