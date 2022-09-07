#ifndef ASSEMBLEV_H
#define ASSEMBLEV_H

#include "GaussQuadratrue.h"

namespace PDWSC{
namespace FEMethod{

template<typename F,typename I,typename Matrix,typename Vector>
void assembleV_1d(I N,I Nlb,I basis_type_test,I der_x,const Matrix &P,const Matrix &T,const Matrix &Pb,const Matrix &Tb,Vector &v,F (*f)(F),I Gpn = 4)
{
	
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
			F r=0;
			r = Gauss_quadratrue_V_1d(vertices,basis_type_test,beta,der_x,Gpn,f);
			I tb = Tb(beta-1,n-1);
			v(tb-1) += r; 
		}
	}	
}

template<typename F,typename I,typename Matrix,typename Vector>
void assembleV_1d_DG(I Nlb,I basis_type_test,I der_x,Matrix &vertices,Vector &v,F (*f)(F),I Gpn = 4)
{   
        v *= 0;	
		for(I beta = 1;beta != Nlb + 1;++beta){
			F r=0;
			r = Gauss_quadratrue_V_1d(vertices,basis_type_test,beta,der_x,Gpn,f);
			v(beta-1) += r; 
		}
}

template<typename F,typename I,typename Matrix,typename Vector>
void assembleV_1d_DG_BD(F x,F uh_pre,I Nlb,I basis_type_trial,I basis_type_test,I der_x,I der_x_pre,Matrix &vertices,Matrix &vertices_pre,Vector &v,I Gpn = 4)
{   
        v *= 0;
        if((vertices - vertices_pre).norm() == 0){
          for(I beta = 1;beta != Nlb + 1;++beta){
			F r=0;
			r = basis_local_function_1d(x, vertices, basis_type_test,beta,der_x);
			v(beta-1) += r; 
		    }
        }else{	
		for(I beta = 1;beta != Nlb + 1;++beta){
			F r=0;
			r = uh_pre*basis_local_function_1d(x, vertices, basis_type_test,beta,der_x);
			v(beta-1) += r; 
		    }
        }
}

} //end of namespace FEMethod
} //end of namespace PDWSC

#endif  // end of ASSEMBLEV_H
