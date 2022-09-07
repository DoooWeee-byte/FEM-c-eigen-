#ifndef ASSEMBLEA_H
#define ASSEMBLEA_H

#include "GaussQuadratrue.h"
#include <Eigen/Sparse>
#include <vector>

namespace PDWSC {
namespace FEMethod {

template<typename F,typename I,typename Matrix,typename Vector,typename SpMat>
void assembleA_2D(I N,I Nlbtest,I Nlbtrial,I basis_type_test,I basis_type_trial,const Matrix &P,const Matrix &T,const Matrix &Tb_test,const Matrix &Tb_trial,SpMat &A,F (*c)(F,F),I r,I s,I p,I q,I Gpn = 4)
{

    typedef Eigen::Triplet<double> Tk;
    std::vector<Tk> tripletList;

    for(I n = 1; n != N+1; ++n) {
        Matrix vertices(P.rows(),T.rows());
        for(I j = 0; j != T.rows(); ++j) {
            I j1 = T(j,n-1);
            for(I i = 0; i != P.rows(); ++i) {
                vertices(i,j) = P(i,j1-1);
            }
        }
        for(I alfa = 1; alfa != Nlbtrial+1; ++alfa) {
            for(I beta = 1; beta != Nlbtest+1; ++beta) {
                F temp = 0;
                temp =Gauss_quadratrue_M_2D<F,I,Matrix,Vector>(vertices,basis_type_test,basis_type_trial,alfa,beta,r,s,p,q,Gpn,c) ;
                I i = Tb_test(beta-1,n-1),j = Tb_trial(alfa-1,n-1);
                tripletList.push_back(Tk(i-1,j-1,temp));
            }
        }
    }
    /**	for(auto i : tripletList){
    		        cout <<"i"<< i.row() << endl;
            cout <<"j"<< i.col() << endl;
            cout << "value"<<i.value() << endl;
    		cout << "*********************"<<endl;	}*/
    A.setFromTriplets(tripletList.begin(),tripletList.end());
}


template<typename F,typename I,typename Matrix,typename Vector,typename SpMat>
void assembleFE_coe_2D(I der_x_coe,I der_y_coe,I basis_type_coe,Vector &solution_vec,I N,I Nlbtest,I Nlbtrial,I basis_type_test,I basis_type_trial,const Matrix &P,const Matrix &T,const Matrix &Tb_test,const Matrix &Tb_trial,SpMat &A,F (*c)(F,F),I r,I s,I p,I q,I Gpn = 4)
{
    typedef Eigen::Triplet<double> Tk;
    std::vector<Tk> tripletList;

    for(I n = 1; n != N+1; ++n) {
        Matrix vertices(P.rows(),T.rows());
        for(I j = 0; j != T.rows(); ++j) {
            I j1 = T(j,n-1);
            for(I i = 0; i != P.rows(); ++i) {
                vertices(i,j) = P(i,j1-1);
            }
        }
        /****取出局部解向量*/
        I row = Tb_trial.rows();
        Vector uh_local_vec(row);
        for(I i=0; i < row; ++i){
            int j = Tb_trial(i,n-1) - 1;
            uh_local_vec(i) = solution_vec(j);
        }

        /*********************************************************/
        for(I alfa = 1; alfa != Nlbtrial+1; ++alfa) {
            for(I beta = 1; beta != Nlbtest+1; ++beta) {
                F temp = 0;
                temp =Gauss_quadratrueFE_coe_2D<F,I,Matrix,Vector>(vertices,basis_type_test,basis_type_trial,alfa,beta,r,s,p,q,Gpn,c,uh_local_vec,der_x_coe,der_y_coe,basis_type_coe) ;
                I i = Tb_test(beta-1,n-1),j = Tb_trial(alfa-1,n-1);
                tripletList.push_back(Tk(i-1,j-1,temp));
            }
        }
    }
    A.setFromTriplets(tripletList.begin(),tripletList.end());
}

}  //end of namespace FEMethod
}  //end of namespace PDWSC
#endif
