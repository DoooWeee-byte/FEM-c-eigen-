#ifndef ODESOLVER_H
#define ODESOLVER_H

#include "assemblev.h"
#include "boundary.h"

namespace PDWSC {
namespace FEMethod {

template<typename Matrix,typename Vector,typename I>
void fiexd_pressure_at_one_point(Matrix &Atilde,Vector &b,Matrix &Pb,I Nb){
	for(I i = 0;i < Atilde.cols();i++){
		Atilde(2*Nb,i) = 0;
	}
	Atilde(2*Nb,2*Nb) = 1;
	b(2*Nb) = 0;
}

template<typename F,typename I,typename Matrix,typename Vector>
void ODE_Unsteady_Stokes(F Tmax,F Tmin,I N_t,I N,I Nlb,I Nb,I basis_type_test,F theta,Matrix &M,Matrix &A,Matrix &P,Matrix &T,Matrix &Pb,Matrix &Tb,Matrix &boundarynodes,Vector &x,Vector &x0,F (*f1)(F,F,F),F (*f2)(F,F,F),F (*g1)(F,F,F),F (*g2)(F,F,F)) {
    F  dt = (Tmax - Tmin)/N_t;
    Matrix Atilde = M/dt + theta*A;
    Matrix A_temp = M/dt - (1-theta)*A;
    Vector xm(x.size());
    for(I m = 0; m < N_t; m++) {
        if(m == 0) {
            xm = x0;
        }
        F t_m1 = (m+1)*dt;
        F t_m = m*dt;
        Vector b_m(x.size()),b_m1(x.size()),b1_m(Nb),b1_m1(Nb);
		b_m = 0*b_m;
		b_m1 = 0*b_m1;
        assembleV_2D(N,Nlb,basis_type_test,0,0,P,T,Pb,Tb,b1_m,f1,t_m);
        assembleV_2D(N,Nlb,basis_type_test,0,0,P,T,Pb,Tb,b1_m1,f2,t_m);
        b_m.head(Nb) = b1_m;
        b_m.segment(Nb,Nb) = b1_m1;
        assembleV_2D(N,Nlb,basis_type_test,0,0,P,T,Pb,Tb,b1_m,f1,t_m1);
        assembleV_2D(N,Nlb,basis_type_test,0,0,P,T,Pb,Tb,b1_m1,f2,t_m1);
        b_m1.head(Nb) = b1_m;
        b_m1.segment(Nb,Nb) = b1_m1;
        Vector b = theta*b_m1 + (1 - theta)*b_m + A_temp*xm;
		treat_boundary_Dirichlet_2d<F,I,Matrix,Vector>(boundarynodes,Atilde,Pb,b,g1,g2,t_m1);
		fiexd_pressure_at_one_point<Matrix,Vector,I>(Atilde,b,Pb,Nb); //case by case
		x = Atilde.completeOrthogonalDecomposition().solve(b);
        xm = x;
    }
}



}// end of namespace FEMethod
}// end of namespace PDWSC
#endif
