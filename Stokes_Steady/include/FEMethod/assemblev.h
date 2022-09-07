#ifndef ASSEMBLEV_H
#define ASSEMBLEV_H

#include "GaussQuadratrue.h"

namespace PDWSC{
namespace FEMethod{

template<typename F,typename I,typename Matrix,typename Vector>
void assembleV_2D(I N,I Nlb,I basis_type_test,I p,I q,const Matrix &P,const Matrix &T,const Matrix &Pb,const Matrix &Tb,Vector &v,F (*f)(F,F),I Gpn = 4)
{
	v = 0*v;
	for(I n = 1;n != N+1;++n){
	// vertices = P(:,T(:,n))
		Matrix vertices(P.rows(),T.rows());
		for(I j = 0; j != T.rows();++j){
			I j1 = T(j,n-1);
			for(I i = 0;i != P.rows();++i){
				vertices(i,j) = P(i,j1-1);
			}
		}
		for(I beta = 1;beta != Nlb + 1;++beta){
			F temp=0;
			temp = Gauss_quadratrue_V_2D<F,I,Matrix,Vector>(vertices,basis_type_test,beta,Gpn,p,q,f);
			I tb = Tb(beta-1,n-1);
			v(tb-1) += temp; 
		}
	}	
}




} //end of namespace FEMethod
} //end of namespace PDWSC

#endif  // end of ASSEMBLEV_H
