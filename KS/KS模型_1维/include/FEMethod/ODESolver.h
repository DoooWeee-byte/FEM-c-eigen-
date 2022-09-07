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
void ODE_KS_1D_Solver(F Time,F delta_t,Matrix &A,Matrix &B,Matrix &Pb,Vector &x,F (*u0)(F))
{
    I Mm = ceil(Time/delta_t);
    Matrix AQuta(A.shape[0],A.shape[1]);

    AQuta = B - delta_t*A;
 /*   cout << "AQuta" << AQuta << endl;
    cout << "B"  << B <<endl;*/
    typename WHYSC::Algebra_kernel<F,I>::CSRMatrix A_quta(AQuta);
    Vector x0(x.size);
    for(I j = 0; j != Pb.shape[1]; ++j) {
        x0[j] = u0(Pb[0][j]);
    }
    cout << "A"<< A << endl;
    cout << "delta_t*A" << delta_t*A <<endl;
    cout << delta_t*A <<endl;
/*    cout << "x0" << x0 <<endl;*/
    Vector b1(x.size);
    b1 = B*x0;
/*    cout  <<"b1" <<  b1 <<endl;*/
    for(I m = 0; m != Mm; ++m) {
	Vector b(x.size);
	b = B*x0;
        A_quta.SOR(b,x,1.5,10000);
        x0 = x;
    }
}


} // end of namespace FEMethod
} // end of namespace PDWSC
#endif
