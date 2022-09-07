#ifndef ASSEMBLEA_H
#define ASSEMBLEA_H

#include "GaussQuadratrue.h"
namespace PDWSC{
namespace FEMethod{

template<typename F,typename I,typename Matrix>
void assembleA_1d(I N,I Nlbtest,I Nlbtrial,I basis_type_trial,I basis_type_test,const Matrix &P,const Matrix &T,const Matrix &Tb_trial,const Matrix &Tb_test,Matrix &A,F (*c)(F),I der_trial,I der_test,I Gpn = 4)
{
	for(I n = 1;n != N+1;++n){
		Matrix vertices(P.shape[0],T.shape[0]);
		for(I j = 0; j != T.shape[0];++j){
			I j1 = T[j][n-1];
			for(I i = 0;i != P.shape[0];++i){
				vertices[i][j] = P[i][j1-1];
			}
		}
		for(I alfa = 1;alfa != Nlbtrial+1; ++alfa){
			for(I beta = 1;beta != Nlbtest+1; ++beta){
				F r = 0;
				r = Gauss_quadratrue_M_1d<F,I,Matrix>(vertices,basis_type_trial,basis_type_test,alfa,beta, der_trial,der_test,Gpn,c);
				I i = Tb_test[beta-1][n-1],j = Tb_trial[alfa-1][n-1];
				A[i-1][j-1] += r;
			}	
		}
	}
}



}  //end of namespace FEMethod
}  //end of namespace PDWSC
#endif
