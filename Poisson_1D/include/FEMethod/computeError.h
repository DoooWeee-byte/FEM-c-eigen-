#ifndef COMPUTEERROR_H
#define COMPUTEERROR_H

#include "GaussQuadratrue.h"
namespace PDWSC {
namespace FEMethod {


template<typename Matrix,typename Vector>
double compute_Hs_error(int N,Matrix &P,Matrix &T,Matrix &Tb,int basis_type,int der_x,Vector &solution_vec,double (*analytic_solution)(double)){
	double error = 0;
	for(int n = 1;n != N+1;++n){
//取得端点坐标放到vertices中
		Matrix vertices(P.shape[0],T.shape[0]);
		for(int j = 0; j != T.shape[0];++j){
			int j1 = T[j][n-1];
			for(int i = 0;i != P.shape[0];++i){
				vertices[i][j] = P[i][j1-1];
			}
		}
//uh_local_vec代表第n个单元上的局部解向量
//取出Tb的第n列的值就是局部解向量在全局解向量中的位置
		Vector uh_local_vec(Tb.shape[0]);
		for(int i=0;i != Tb.shape[0];++i){
			int j = Tb[i][n-1] - 1;
			uh_local_vec[i] = solution_vec[j];
		}
		double int_value = Gauss_int_error_1D(vertices,basis_type,der_x,uh_local_vec,analytic_solution);
		error += int_value;
	}
	error = sqrt(error);
	return error;
}

}  // end of namespace FEMethod
}  // end of namespace PDWSC
#endif
