#ifndef ASSEMBLEV_H
#define ASSEMBLEV_H

#include "GaussQuadratrue.h"

namespace PDWSC {
namespace FEMethod {

template<typename F,typename I,typename Matrix,typename Vector>
void assembleV_2D(I N,I Nlb,I basis_type_test,I p,I q,const Matrix &P,const Matrix &T,const Matrix &Pb,const Matrix &Tb,Vector &v,F (*f)(F,F,F),F t,I Gpn = 4)
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
			temp = Gauss_quadratrue_V_2D<F,I,Matrix,Vector>(vertices,basis_type_test,beta,Gpn,p,q,f,t);
			I tb = Tb(beta-1,n-1);
			v(tb-1) += temp; 
		}
	}	
}

template<typename F,typename I,typename Matrix,typename Vector>
void assembleFE_coe_2D_V(I der_x_coe_1,I der_y_coe_1,I basis_type_coe_1,Vector &solution_vec_1,I der_x_coe_2,I der_y_coe_2,I basis_type_coe_2,Vector &solution_vec_2,I N,I Nlb,I basis_type_test,I p,I q,const Matrix &P,const Matrix &T,const Matrix &Pb,const Matrix &Tb,Vector &v,F (*f)(F,F),I Gpn = 4)
{
    v = v*0;
    for(I n = 1; n != N+1; ++n) {
        // vertices = P(:,T(:,n))
        Matrix vertices(P.rows(),T.rows());
        for(I j = 0; j != T.rows(); ++j) {
            I j1 = T(j,n-1);
            for(I i = 0; i != P.rows(); ++i) {
                vertices(i,j) = P(i,j1-1);
            }
        }
        /****取出局部解向量*/
		Vector uh_local_vec_1(Tb.rows()),uh_local_vec_2(Tb.rows());
		int row = Tb.rows();
		for(int i=0;i < row;++i){
			int j = Tb(i,n-1) - 1;
			uh_local_vec_1(i) = solution_vec_1(Tb(i,n-1) - 1);
			uh_local_vec_2(i) = solution_vec_2(Tb(i,n-1) - 1);
		}	
        for(I beta = 1; beta != Nlb + 1; ++beta) {
            F temp=0;
            temp = Gauss_quadratrueFE_coe_2D_V<F,I,Matrix,Vector>(vertices,basis_type_test,beta,Gpn,p,q,f,uh_local_vec_1,der_x_coe_1,der_y_coe_1,basis_type_coe_1,uh_local_vec_2,der_x_coe_2,der_y_coe_2,basis_type_coe_2);
            I tb = Tb(beta-1,n-1);
            v(tb-1) += temp;
        }
    }
}



} //end of namespace FEMethod
} //end of namespace PDWSC

#endif  // end of ASSEMBLEV_H
