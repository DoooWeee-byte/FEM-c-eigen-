#ifndef ODESOLVER_H
#define ODESOLVER_H



#include <iostream>
#include <cmath>
#include "algebra/Algebra_kernel.h"
#include "Constants.h"

#include "assembleA.h"
#include "GaussQuadratrue.h"
#include "basefunction.h"
#include "algebra/Algebra_kernel.h"
#include "assemblev.h"
#include "boundary.h"

using namespace std;
namespace PDWSC {
namespace FEMethod {

template<typename F,typename I,typename Vector,typename Matrix>
void ODE_HEAT_Solver(F Time,F delta_t,F theta,Matrix &A,Vector &x,I N,I Nlb,I basis_type_test,I p,I q,F (*f)(F,F,F),F (*g)(F,F,F),F (*u0)(F,F),Matrix &M,Matrix &P,Matrix &T,Matrix &Pb,Matrix &Tb,Matrix &boundarynodes)
{
    I Mm = ceil(Time/delta_t);
    F temp1 = 1/delta_t;
    Matrix AQuta(A.shape[0],A.shape[1]);

    AQuta = temp1*M + theta*A;
    Matrix ATemp(A.shape[0],A.shape[1]);
    ATemp = temp1*M - (1-theta)*A;
    for(I j = 0; j != Pb.shape[1]; ++j) {
        x[j] = u0(Pb[0][j],Pb[1][j]);
    }
    for(I m = 0; m != Mm; ++m) {
        Vector v_m1(x.size);
        Vector v_m(x.size);

        Vector b_m1(x.size);
        F t_m1 = (m+1)*delta_t;
        F t_m = m*delta_t;
        assembleV_2D(N,Nlb,basis_type_test,p,q,P,T,Pb,Tb,v_m1,f,t_m1);//t_m1
        assembleV_2D(N,Nlb,basis_type_test,p,q,P,T,Pb,Tb,v_m	,f,t_m);//t_m
        b_m1 = theta*v_m1 + (1-theta)*v_m + ATemp*x;
        treat_boundary_Dirichlet_2d<F,I,Matrix,Vector>(boundarynodes,AQuta,Pb,b_m1,g,t_m1);

        typename WHYSC::Algebra_kernel<F,I>::CSRMatrix A_quta(AQuta);
        A_quta.SOR(b_m1,x,1.5,10000);
    }
}
template<typename F,typename I,typename Vector,typename Matrix>
void ODE_WAVE_Solver(F Time,F delta_t,F theta,Matrix &A,Vector &x,I N,I Nlb,I basis_type_test,I p,I q,F (*f)(F,F,F),F (*g)(F,F,F),F (*u0)(F,F),F (*u00)(F,F),Matrix &M,Matrix &P,Matrix &T,Matrix &Pb,Matrix &Tb,Matrix &boundarynodes)
{
    I Mm = ceil(Time/delta_t);
    F temp1 = 1/(delta_t*delta_t);
    Matrix ATemp1(A.shape[0],A.shape[1]);
    ATemp1 = 2*temp1*M-0.5*A;
    Matrix ATemp2(A.shape[0],A.shape[1]);
    ATemp2 = temp1*M+0.25*A;
    Vector x0(x.size);
    Matrix AQuta(A.shape[0],A.shape[1]);
    AQuta = temp1*M + 0.25*A;
  
    for(I j = 0; j != Pb.shape[1]; ++j) {
        x[j] = u0(Pb[0][j],Pb[1][j]);
    }
    for(I j = 0; j != Pb.shape[1]; ++j) {
        x0[j] = -u00(Pb[0][j],Pb[1][j])*delta_t+u0(Pb[0][j],Pb[1][j]);
    }
    for(I m = 0; m != Mm; ++m) {
        Vector v_m(x.size);

        Vector b_m1(x.size);
        F t_m = m*delta_t;
        F t_m1 = (m+1)*delta_t;
        assembleV_2D(N,Nlb,basis_type_test,p,q,P,T,Pb,Tb,v_m,f,t_m);//t_m

        b_m1 = v_m + ATemp1*x - ATemp2*x0;
	x0 = x;
        treat_boundary_Dirichlet_2d<F,I,Matrix,Vector>(boundarynodes,AQuta,Pb,b_m1,g,t_m1);

        typename WHYSC::Algebra_kernel<F,I>::CSRMatrix A_quta(AQuta);
        A_quta.SOR(b_m1,x,1.5,10000);
    }
}


} // end of namespace FEMethod
} // end of namespace PDWSC
#endif
