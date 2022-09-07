#ifndef ASSEMBLEA_H
#define ASSEMBLEA_H

#include "GaussQuadratrue.h"
#include "basefunction.h"
#include <Eigen/Sparse>
#include <vector>
namespace PDWSC{
namespace FEMethod{

template<typename F,typename I,typename Matrix>
void assembleA_1d(I N,I Nlbtest,I Nlbtrial,I basis_type_trial,I basis_type_test,const Matrix &P,const Matrix &T,const Matrix &Tb_trial,const Matrix &Tb_test,Matrix &A,F (*c)(F),I der_trial,I der_test,I Gpn = 4)
{
	for(I n = 1;n != N+1;++n){
		Matrix vertices(P.rows(),T.rows());
		for(I j = 0; j != T.rows();++j){
			I j1 = T(j,n-1);
			for(I i = 0;i != P.rows();++i){
				vertices(i,j) = P(i,j1-1);
			}
		}
		for(I alfa = 1;alfa != Nlbtrial+1; ++alfa){
			for(I beta = 1;beta != Nlbtest+1; ++beta){
				F r = 0;
				r = Gauss_quadratrue_M_1d<F,I,Matrix>(vertices,basis_type_trial,basis_type_test,alfa,beta, der_trial,der_test,Gpn,c);
				I i = Tb_test(beta-1,n-1),j = Tb_trial(alfa-1,n-1);
				A(i-1,j-1) += r;
			}	
		}
	}
}

template<typename F,typename I,typename Matrix>
void assembleA_DG_1d(I Nlbtest,I Nlbtrial,I basis_type_trial,I basis_type_test,Matrix &vertices, Matrix &A,F (*c)(F),I der_trial,I der_test,I Gpn = 4)
{
        A *= 0;
		for(I alfa = 1;alfa != Nlbtrial+1; ++alfa){
			for(I beta = 1;beta != Nlbtest+1; ++beta){
				F r = 0;
				r = Gauss_quadratrue_M_1d<F,I,Matrix>(vertices,basis_type_trial,basis_type_test,alfa,beta, der_trial,der_test,Gpn,c);
				A(beta-1,alfa-1) += r;
			}	
		}
}

template<typename F,typename I,typename Matrix>
void assembleA_DG_bd_1d(F x,I Nlbtest,I Nlbtrial,I basis_type_trial,I basis_type_test,Matrix &vertices,Matrix &A,F (*c)(F),I der_trial,I der_test,I Gpn = 4)
{
        A *= 0;
		for(I alfa = 1;alfa != Nlbtrial+1; ++alfa){
			for(I beta = 1;beta != Nlbtest+1; ++beta){
				F r = 0;
				r = c(x)*basis_local_function_1d(x, vertices, basis_type_trial, alfa, der_trial)*basis_local_function_1d(x, vertices, basis_type_test, beta, der_test);
				A(beta-1,alfa-1) += r;
			}	
		}
}


}  //end of namespace FEMethod
}  //end of namespace PDWSC
#endif
