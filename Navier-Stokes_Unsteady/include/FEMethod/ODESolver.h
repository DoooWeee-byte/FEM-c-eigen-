#ifndef ODESOLVER_H
#define ODESOLVER_H

#include "assemblev.h"
#include "assembleA.h"
#include "boundary.h"
#include "addSpMat.h"
namespace PDWSC {
namespace FEMethod {

template<typename Matrix,typename Vector,typename I>
void fiexd_pressure_at_one_point(Matrix &Atilde,Vector &b,Matrix &Pb,I Nb) {
    for(I i = 0; i < Atilde.cols(); i++) {
        Atilde(2*Nb,i) = 0;
    }
    Atilde(2*Nb,2*Nb) = 1;
    b(2*Nb) = 0;
}

template<typename Matrix,typename Vector,typename I,typename SpMat>
void fiexd_pressure_at_one_point_SpMat(SpMat &Atilde,SpMat &Am1l,Vector &b,Matrix &Pb,I Nb) {
	typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
	tripletList.push_back(T(2*Nb,2*Nb,1));
    for(I k = 0; k < Atilde.outerSize(); ++k) {
        for(Eigen::SparseMatrix<double>::InnerIterator it(Atilde,k); it; ++it) {
			I row = it.row();
			if(row != 2*Nb){
				tripletList.push_back(T(it.row(),it.col(),it.value()));
			}
        }
    }
	Am1l.setFromTriplets(tripletList.begin(),tripletList.end());
}

double c(double x,double y) {
    return 1;
}

double f(double x,double y) {
    return 1;
}

template<typename F,typename I,typename Matrix,typename Vector,typename SpMat>
void ODE_Unsteady_NaStokes(F Tmax,F Tmin,I N_t,I N,I Nlb,I Nb,I basis_type_test,F theta,SpMat &M,SpMat &A,Matrix &P,Matrix &T,Matrix &Pb,Matrix &Tb,Matrix &boundarynodes,Vector &x,Vector &x0,F (*f1)(F,F,F),F (*f2)(F,F,F),F (*g1)(F,F,F),F (*g2)(F,F,F)) {
    F  dt = (Tmax - Tmin)/N_t;
    Vector xm(x.size());
    I nlb_u = Nlb;
    I basis_type_u = 202;

    for(I m = 0; m < N_t; m++) {
        if(m == 0) {
            xm = x0;//上一个时间层的解
        }
        F t_m1 = (m+1)*dt;
        F t_m = m*dt;
        //组装b_m1
        Vector b_m1(x.size()),b1_m(Nb),b1_m1(Nb);
        b_m1 = 0*b_m1;
        assembleV_2D_t(t_m1,N,Nlb,basis_type_test,0,0,P,T,Pb,Tb,b1_m,f1);
        assembleV_2D_t(t_m1,N,Nlb,basis_type_test,0,0,P,T,Pb,Tb,b1_m1,f2);
        b_m1.head(Nb) = b1_m;
        b_m1.segment(Nb,Nb) = b1_m1;

        //装填初始牛顿迭代向量
        Vector V1(Nb),V2(Nb);
        if(m == 0) {
            V1 = x0.head(Nb);
            V2 = x0.segment(Nb,Nb);
        }
        else {
            V1 = x.head(Nb);
            V2 = x.segment(Nb,Nb);
        }


        //牛顿迭代
        I L = 30;
        F Tol = pow(10,-6);
        Vector x_old(x.size());//储存上一个牛顿迭代的解
        for(I i = 1; i < L+1; ++i) {
            //组装牛顿迭代矩阵
            SpMat AN(A.rows(),A.cols());
            Vector VN(x.size());
            if(1) {
                int basis_type_coe = 202;
                SpMat A1(Nb,Nb),A2(Nb,Nb),A3(Nb,Nb),A4(Nb,Nb),A5(Nb,Nb),A6(Nb,Nb);

                assembleFE_coe_2D(1,0,basis_type_coe,V1,N,nlb_u,nlb_u,basis_type_u,basis_type_u,P,T,Tb,Tb,A1,c,0,0,0,0);
                assembleFE_coe_2D(0,0,basis_type_coe,V1,N,nlb_u,nlb_u,basis_type_u,basis_type_u,P,T,Tb,Tb,A2,c,0,0,1,0);
                assembleFE_coe_2D(0,0,basis_type_coe,V2,N,nlb_u,nlb_u,basis_type_u,basis_type_u,P,T,Tb,Tb,A3,c,0,0,0,1);
                assembleFE_coe_2D(0,1,basis_type_coe,V1,N,nlb_u,nlb_u,basis_type_u,basis_type_u,P,T,Tb,Tb,A4,c,0,0,0,0);
                assembleFE_coe_2D(1,0,basis_type_coe,V2,N,nlb_u,nlb_u,basis_type_u,basis_type_u,P,T,Tb,Tb,A5,c,0,0,0,0);
                assembleFE_coe_2D(0,1,basis_type_coe,V2,N,nlb_u,nlb_u,basis_type_u,basis_type_u,P,T,Tb,Tb,A6,c,0,0,0,0);

                SpMat A123 = A1 + A2 + A3;
                addSpMat(AN,A123,0,0);
                addSpMat<I,SpMat>(AN,A4,0,A1.cols());
                addSpMat<I,SpMat>(AN,A5,A1.rows(),0);
                SpMat A236 = A2 + A3 + A6;
                addSpMat<I,SpMat>(AN,A236,A1.rows(),A1.cols());
            }
            if(1) {
                int basis_type_coe_1 = 202,basis_type_coe_2 = 202;
                Vector VN1(V1.size()),VN2(V1.size()),VN3(V2.size()),VN4(V2.size());

                assembleFE_coe_2D_V(0,0,basis_type_coe_1,V1,1,0,basis_type_coe_2,V1,N,nlb_u,basis_type_u,0,0,P,T,Pb,Tb,VN1,f);
                assembleFE_coe_2D_V(0,0,basis_type_coe_1,V2,0,1,basis_type_coe_2,V1,N,nlb_u,basis_type_u,0,0,P,T,Pb,Tb,VN2,f);
                assembleFE_coe_2D_V(0,0,basis_type_coe_1,V1,1,0,basis_type_coe_2,V2,N,nlb_u,basis_type_u,0,0,P,T,Pb,Tb,VN3,f);
                assembleFE_coe_2D_V(0,0,basis_type_coe_1,V2,0,1,basis_type_coe_2,V2,N,nlb_u,basis_type_u,0,0,P,T,Pb,Tb,VN4,f);

                VN = VN*0;
                VN.head(V1.size()) = VN1 + VN2;
                VN.segment(V1.size(),V2.size()) = VN3 + VN4;
            }

            //生成迭代解矩阵和向量
            /* Matrix Am1l = M/dt + A + AN; */
            SpMat Am1l = M/dt + A + AN;
			SpMat Atemp(Am1l.rows(),Am1l.cols());
            Vector bm1l = b_m1 + (M/dt)*xm + VN;
            //处理Dirichlet边界条件
            treat_boundary_Dirichlet_2d_SpMat<F,I,Matrix,Vector,SpMat>(boundarynodes,Am1l,Atemp,Pb,bm1l,g1,g2,t_m1);
            /* treat_boundary_Dirichlet_2d<F,I,Matrix,Vector>(boundarynodes,Am1l,Pb,bm1l,g1,g2,t_m1); */
            fiexd_pressure_at_one_point_SpMat(Atemp,Am1l,bm1l,Pb,Nb);
            Am1l.makeCompressed();
            //求解线性系统
            /* Eigen:: BiCGSTAB<SpMat>  solver; */
			/* solver.compute(Am1l); */
			/* x = solver.solve(bm1l); */

			Eigen::SparseLU<SpMat> solver;
            solver.analyzePattern(Am1l);
            solver.factorize(Am1l);
			if(solver.info() != Eigen::Success){
				cout << "solve fail!" << endl;
			}
			x = solver.solve(bm1l);

            /* fiexd_pressure_at_one_point(Am1l,bm1l,Pb,Nb); */
            /* x = Am1l.colPivHouseholderQr().solve(bm1l); */
            V1 = x.head(Nb);
            V2 = x.segment(Nb,Nb);
            //判断牛顿迭代是否停止
            if((x - x_old).norm() < Tol)
            {
                cout << "NewtonEndTimes:" << i << endl;
                break;
            }
            x_old = x;
        }
        xm = x;
    }
}





}// end of namespace FEMethod
}// end of namespace PDWSC
#endif
