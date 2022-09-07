#ifndef EFMETHOD_KERNEL_H
#define EFMETHOD_KERNEL_H

#include "generate_PT.h"
#include "assembleA.h"
#include "GaussQuadratrue.h"
#include "basefunction.h"
#include "algebra/Algebra_kernel.h"
#include "assemblev.h"
#include "boundary.h"
#include "computeError.h"

namespace PDWSC {

template<typename F = double,typename I = int>
class FEMethod_kernel
{
public:
    typedef typename WHYSC::AlgebraObject::Matrix<F, I> Matrix;
    typedef typename WHYSC::AlgebraObject::Vector<F, I> Vector;
    typedef typename WHYSC::AlgebraObject::CSRMatrix<F, I> CSRMatrix;
public:
    static void generate_PT_1d(F a,F b,I N,Matrix &P,Matrix &T)
    {
        return FEMethod::generate_PT_1d<F,I,Matrix>(a,b,N,P,T);
    }
    static void generate_PbTb_1d(F a,F b,I N,I basis_type,I &Nb,I &Nlb,Matrix &Pb,Matrix &Tb) 
    {
        return FEMethod::generate_PbTb_1d<F,I,Matrix>(a,b,N,basis_type,Nb,Nlb,Pb,Tb);
    }
	static void generate_boundraynodes(Matrix &boundarynodes,Matrix &Tb,int boundary_type = 101) {
        return FEMethod::generate_boundraynodes<F,I,Matrix>(boundarynodes,Tb,boundary_type) ;
    }

   static void assembleA_1d(I N,I Nlbtest,I Nlbtrial,I basis_type_trial,I basis_type_test,const Matrix &P,const Matrix &T,const Matrix &Tb_trial,const Matrix &Tb_test,Matrix &A,F (*c)(F),I der_trial,I der_test,I Gpn = 4){
		return FEMethod::assembleA_1d<F,I,Matrix>(N,Nlbtest,Nlbtrial,basis_type_trial,basis_type_test,P,T,Tb_trial,Tb_test,A,c,der_trial,der_test,Gpn); 
   }
  static void assembleV_1d(I N,I Nlb,I basis_type_test,I der_x,const Matrix &P,const Matrix &T,const Matrix &Pb,const Matrix &Tb,Vector &v,F (*f)(F),I Gpn = 4){
		return FEMethod::assembleV_1d<F,I,Matrix,Vector>(N,Nlb,basis_type_test,der_x,P,T,Pb,Tb,v,f,Gpn);
  } 
    static void treat_boundary_1d(Matrix &boundarynodes,Matrix &A,Matrix &Pb,Vector &b,F (*g)(F)) {
        return FEMethod::treat_boundary_1d<F,I,Matrix,Vector>(boundarynodes,A,Pb,b,g);
    }

	static F compute_Hs_error(I N,Matrix &P,Matrix &T,Matrix &Tb,I basis_type,I der_x,Vector &solution_vec,F (*analytic_solution)(F)){
		return FEMethod::compute_Hs_error(N,P,T,Tb,basis_type,der_x,solution_vec,analytic_solution);	
	}



};
} // end of namespace PDWSC
#endif // end of EFMETHOD_KERNEL_H
