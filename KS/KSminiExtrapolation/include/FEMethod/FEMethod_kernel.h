#ifndef EFMETHOD_KERNEL_H
#define EFMETHOD_KERNEL_H


#include "Eigen/Eigen"
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
    typedef typename Eigen::MatrixXd Matrix;
    typedef typename Eigen::VectorXd Vector;
    typedef typename Eigen::SparseMatrix<F> SpMat;
    typedef typename WHYSC::AlgebraObject::CSRMatrix<F, I> CSRMatrix;
public:
    static void generate_PT_2D(Matrix &P,Matrix &T,I N1,I N2,F right,F top,F bottom = 0,F left = 0)
    {
        return FEMethod::generate_PT_2D<F,I,Matrix>(P,T,N1,N2,right,top,bottom,left);
    }

    static void generate_boundaryedges_2D(Matrix &boundaryedges,I N1,I N2) {
        return FEMethod::generate_boundaryedges_2D<F,I,Matrix>(boundaryedges,N1,N2);
    }


    static F compute_Hs_error_2D(I N,Matrix &P,Matrix &T,Matrix &Tb,I basis_type,I der_x,I der_y,Vector &solution_vec,F (*analytic_solution)(F,F),I Gpn=4) {
        return FEMethod::compute_Hs_error_2D<Matrix,Vector>(N,P,T,Tb,basis_type,der_x,der_y,solution_vec,analytic_solution,Gpn);
    }

    static void generate_PbTb_2D(Matrix &Pb,Matrix &Tb,I N1,I N2,I basis_type,F right,F top,F bottom = 0,F left = 0) {
        return FEMethod::generate_PbTb_2D(Pb,Tb,N1,N2,basis_type,right,top,bottom,left);
    }


    static void assembleA_2D(I N,I Nlbtest,I Nlbtrial,I basis_type_test,I basis_type_trial,const Matrix &P,const Matrix &T,const Matrix &Tb_test,const Matrix &Tb_trial,SpMat &A,F (*c)(F,F),I r,I s,I p,I q,I Gpn = 4)
    {
        return FEMethod::assembleA_2D<F,I,Matrix,Vector,SpMat>(N,Nlbtest,Nlbtrial,basis_type_test,basis_type_trial,P,T,Tb_test,Tb_trial,A,c,r,s,p,q,Gpn);
    }
    static void assembleV_2D(I N,I Nlb,I basis_type_test,I p,I q,const Matrix &P,const Matrix &T,const Matrix &Pb,const Matrix &Tb,Vector &v,F (*f)(F,F),I Gpn = 4)
    {
        return FEMethod::assembleV_2D<F,I,Matrix,Vector>(N,Nlb,basis_type_test,p,q,P,T,Pb,Tb,v,f,Gpn);
    }
    static void generate_boundaryedges_DN_2D(Matrix &boundaryedges,I N1,I N2) {
        return FEMethod::generate_boundaryedges_DN_2D<F,I,Matrix>(boundaryedges,N1,N2);
    }
    static void generate_boundarynodes_2D(Matrix &boundarynodes,I N1,I N2)
    {
        return FEMethod::generate_boundarynodes_2D<F,I,Matrix>(boundarynodes,N1,N2);
    }
    static void treat_boundary_Dirichlet_2d(Matrix &boundarynodes,Matrix &A,Matrix &Pb,Vector &b,F (*g1)(F,F),F (*g2)(F,F))
    {
        return FEMethod::treat_boundary_Dirichlet_2d<F,I,Matrix,Vector>(boundarynodes,A,Pb,b,g1,g2);
    }
    static void treat_boundary_Neumann_2d(Vector &b1,Matrix &boundaryedges,Matrix &Pb,Matrix &Tb,Matrix &P,Matrix &T,F (*Neufunc)(F,F),I nlb,I basis_type,I der_x,I der_y,I Gpn = 4)
    {
        return FEMethod::treat_boundary_Neumann_2d<F,I,Matrix,Vector>(b1,boundaryedges,Pb,Tb,P,T,Neufunc,nlb,basis_type,der_x,der_y,Gpn);
    }
    static void treat_boundary_Robin_2d(SpMat &A1,Vector &b1,Matrix &boundaryedges,Matrix &T,Matrix &P,Matrix &Tb_trial,Matrix &Tb_test,I basis_type_trial,I basis_type_test,I der_x_trial,I der_y_trial,I der_x_test,I der_y_test,I nlb_trial,I nlb_test,F (*Robinfunc)(F,F),F (*Neufunc)(F,F),I Gpn=4) {

        return FEMethod::treat_boundary_Robin_2d<F,I,Matrix,Vector,SpMat>(A1,b1,boundaryedges,T,P,Tb_trial,Tb_test,basis_type_trial,basis_type_test,der_x_trial,der_y_trial,der_x_test,der_y_test,nlb_trial,nlb_test,Robinfunc,Neufunc,Gpn);
    }
    static void generate_boundaryedges_DR_2D(Matrix &boundaryedges,I N1,I N2) {

        return FEMethod::generate_boundaryedges_DR_2D<F,I,Matrix>(boundaryedges,N1,N2);
    }
    static void generate_boundarynodes_DR_2D(Matrix &boundarynodes,I N1,I N2) {
        return FEMethod::generate_boundarynodes_DR_2D<F,I,Matrix>(boundarynodes,N1,N2);
    }

    static void generate_boundarynodes_DN_2D(Matrix &boundarynodes,I N1,I N2) {
        return FEMethod::generate_boundarynodes_DN_2D<F,I,Matrix>(boundarynodes,N1,N2);
    }


    static void assembleFE_coe_2D(I der_x_coe,I der_y_coe,I basis_type_coe,Vector &solution_vec,I N,I Nlbtest,I Nlbtrial,I basis_type_test,I basis_type_trial,const Matrix &P,const Matrix &T,const Matrix &Tb_test,const Matrix &Tb_trial,SpMat &A,F (*c)(F,F),I r,I s,I p,I q,I Gpn = 4) {
        return FEMethod::assembleFE_coe_2D<F,I,Matrix,Vector, SpMat>(der_x_coe,der_y_coe,basis_type_coe,solution_vec,N,Nlbtest,Nlbtrial,basis_type_test,basis_type_trial,P,T,Tb_test,Tb_trial,A,c,r,s,p,q,Gpn);
    }
	
	static void assembleFE_coe_2D_V(I der_x_coe_1,I der_y_coe_1,I basis_type_coe_1,Vector &solution_vec_1,I der_x_coe_2,I der_y_coe_2,I basis_type_coe_2,Vector &solution_vec_2,I N,I Nlb,I basis_type_test,I p,I q,const Matrix &P,const Matrix &T,const Matrix &Pb,const Matrix &Tb,Vector &v,F (*f)(F,F),I Gpn = 4){
		return FEMethod::assembleFE_coe_2D_V<F,I,Matrix,Vector>(der_x_coe_1,der_y_coe_1,basis_type_coe_1,solution_vec_1,der_x_coe_2,der_y_coe_2,basis_type_coe_2,solution_vec_2,N,Nlb,basis_type_test,p,q,P,T,Pb,Tb,v,f,Gpn);
	}
		static F compute_Max_error_2D(const Vector &solution_vec,const Matrix &Pb,F (*analytic_solution)(F,F)){
		return FEMethod::compute_Max_error_2D<Matrix,Vector,F,I>(solution_vec,Pb,analytic_solution);
	}
	static F local_FE_function_2D(F x,F y,const Matrix &vertices,I basis_type,I der_x,I der_y,Vector &uh_local_vec){
		return FEMethod::local_FE_function_2D<F,I,Matrix,Vector>(x,y,vertices,basis_type,der_x,der_y,uh_local_vec);
	}

};
} // end of namespace PDWSC
#endif // end of EFMETHOD_KERNEL_H
