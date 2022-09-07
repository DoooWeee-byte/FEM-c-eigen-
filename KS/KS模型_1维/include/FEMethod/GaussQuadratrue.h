#ifndef GAUSSQUADRATRUE_H
#define GAUSSQUADRATRUE_H

#include <vector>
#include <cmath>
#include "basefunction.h"
using namespace std;
namespace PDWSC {
namespace FEMethod {

template<typename F,typename I>
void generate_Gauss_reference_1D(I Gauss_point_number,vector<F> &Gc,vector<F> &Gp)
{
    if(Gauss_point_number==4) {
        Gc = {0.3478548451,0.3478548451,0.6521451549,0.6521451549};
        Gp = {0.8611363116,-0.8611363116,0.3399810436,-0.3399810436};
    }
    else if(Gauss_point_number==8) {
        Gc = {0.1012285363,0.1012285363,0.2223810345,0.2223810345,0.3137066459,0.3137066459,0.3626837834,0.3626837834};
        Gp = {0.9602898565,-0.9602898565,0.7966664774,-0.7966664774,0.5255324099,-0.5255324099,0.1834346425,-0.1834346425};
    }
    else {
        Gc = {1,1};
        Gp = {-1/sqrt(3),1/sqrt(3)};
    }
}

template<typename F,typename I>
void generate_Gauss_local_1D(F upper_bound,F lower_bound,vector<F> &Gc,vector<F> &Gp,vector<F> &Glc,vector<F> &Glp)
{
    for(auto i = 0; i != Gc.size(); ++i) {
        F lc = (upper_bound-lower_bound)*Gc[i]/2;
        F lp = (upper_bound-lower_bound)*Gp[i]/2+(upper_bound+lower_bound)/2;
        Glc.push_back(lc);
        Glp.push_back(lp);
    }
}

template<typename F,typename I,typename Matrix>
F Gauss_quadratrue_M_1d(Matrix &vertices,I basis_type_trial,I basis_type_test,I alfa,I beta,I der_trial,I der_test,I pn,F (*c)(F)) {
// alfa对应trial,beta对应test
    F lower_bound = vertices[0][0];
    F upper_bound = vertices[0][1];
    vector<F> Gc,Gp,Glc,Glp;
    generate_Gauss_reference_1D<F,I>(pn,Gc,Gp);
    generate_Gauss_local_1D<F,I>(upper_bound,lower_bound,Gc,Gp,Glc,Glp);
    F ret = 0;
    for(auto i = 0; i != Glc.size(); ++i) {
        //高斯权重*某个系数函数*test函数的某阶导数*trial函数的某阶导数
        ret += Glc[i]*c(Glp[i])*basis_local_function_1d(Glp[i],vertices,basis_type_trial,alfa,der_trial)*basis_local_function_1d(Glp[i],vertices,basis_type_test,beta,der_test);
    }
    return ret;
}

template<typename F,typename I,typename Matrix>
F Gauss_quadratrue_V_1d(Matrix &vertices,I basis_type,I basis_index,I der_x,I pn,F (*f)(F)) {
    F upper_bound = vertices[0][1];
    F lower_bound = vertices[0][0];
    vector<F> Gc,Gp,Glc,Glp;
    generate_Gauss_reference_1D<F,I>(pn,Gc,Gp);
    generate_Gauss_local_1D<F,I>(upper_bound,lower_bound,Gc,Gp,Glc,Glp);
    F int_value = 0;
    for(auto i = 0; i != Glc.size(); ++i) {
        int_value += Glc[i]*f(Glp[i])*basis_local_function_1d(Glp[i],vertices,basis_type,basis_index,der_x);
    }
    return int_value;
}

template<typename F,typename I,typename Matrix,typename Vector>
F Gauss_int_error_1D(Matrix &vertices,I basis_type,I der_x,Vector &uh_local_vec,F (*analytic_solution)(F),I pn = 4) {
    F upper_bound = vertices[0][1];
    F lower_bound = vertices[0][0];
    vector<F> Gc,Gp,Glc,Glp;
    generate_Gauss_reference_1D<F,I>(pn,Gc,Gp);
    generate_Gauss_local_1D<F,I>(upper_bound,lower_bound,Gc,Gp,Glc,Glp);
    F int_value = 0;
    for(auto i = 0; i != Glc.size(); ++i) {
		F zz = analytic_solution(Glp[i])- local_FE_function_1D(Glp[i],vertices,basis_type,der_x,uh_local_vec);
        int_value += Glc[i]*zz*zz;
    }
    return int_value;

}




} // end of namespace FEMethod
} // end of namespace PDWSC
#endif
