#ifndef BOUNDARY_H
#define BOUNDARY_H

namespace PDWSC {
namespace FEMethod {



template<typename F,typename I,typename Matrix,typename Vector,typename SpMat>
void treat_boundary_Dirichlet_2d(Matrix &boundarynodes,SpMat &A,Matrix &Pb,Vector &b,F (*g)(F,F))
{
    for(I k = 0; k != boundarynodes.cols(); ++k) {
        if(boundarynodes(0,k) == -1) {
            I i = boundarynodes(1,k);
            for(I j = 0; j != A.cols(); ++j) {
            	A.coeffRef(i-1,j) = 0;
		A.makeCompressed();
            }
            A.coeffRef(i-1,i-1) = 1;
	    A.makeCompressed();
            b(i-1) = g(Pb(0,i-1),Pb(1,i-1));
        }
    }
}

template<typename F,typename I,typename Matrix,typename Vector>
void treat_boundary_Neumann_2d(Vector &b1,Matrix &boundaryedges,Matrix &Pb,Matrix &Tb,Matrix &P,Matrix &T,F (*Neufunc)(F,F),I nlb,I basis_type,I der_x,I der_y,I Gpn =4)
{   //用test的Pb,Tb.
    I nbn = boundaryedges.cols();
    Vector end_ponit_1(2);
    Vector end_point_2(2);
    for(I k = 0; k != nbn; ++k) {
        if(boundaryedges(0,k) == -2) {
            //n:代表第几个单元
            I n = boundaryedges(1,k);
            I n_end_1 = boundaryedges(2,k);
            I n_end_2 = boundaryedges(3,k);
            end_ponit_1(0) = P(0,n_end_1-1);
            end_ponit_1(1) = P(1,n_end_1-1);
            end_point_2(0) = P(0,n_end_2-1);
            end_point_2(1) = P(1,n_end_2-1);
            Matrix vertices(P.rows(),T.rows());
	    for(I j = 0; j != T.rows();++j){
		I j1 = T(j,n-1);
		for(I i = 0;i != P.rows();++i){
			vertices(i,j) = P(i,j1-1);
		}
        	}
            for(I alfa=1; alfa != nlb+1; ++alfa) {
                F temp = 0;
                temp = Gauss_quadratrue_line_test_2D(end_ponit_1,end_point_2,vertices,basis_type,alfa,der_x,der_y,Neufunc,Gpn);
                I location = Tb(alfa-1,n-1);
                b1(location-1) += temp;
            }
        }
    }
}

template<typename F,typename I,typename Matrix,typename Vector,typename SpMat>
void treat_boundary_Robin_2d(SpMat &A1,Vector &b1,Matrix &boundaryedges,Matrix &T,Matrix &P,Matrix &Tb_trial,Matrix &Tb_test,I basis_type_trial,I basis_type_test,I der_x_trial,I der_y_trial,I der_x_test,I der_y_test,I nlb_trial,I nlb_test,F (*Robinfunc)(F,F),F (*Neufunc)(F,F),I Gpn)
{
    I nbn = boundaryedges.cols();
    Vector end_point_1(2);
    Vector end_point_2(2);
    for(I k = 0; k != nbn; ++k) {
        if(boundaryedges(0,k) == -3) {
            //n:代表第n个单元
            I n = boundaryedges(1,k);
            I n_end_1 = boundaryedges(2,k);
            I n_end_2 = boundaryedges(3,k);
            end_point_1(0) = P(0,n_end_1-1);
            end_point_1(1) = P(1,n_end_1-1);
            end_point_2(0) = P(0,n_end_2-1);
            end_point_2(1) = P(1,n_end_2-1);
            Matrix vertices(P.rows(),T.rows());
            for(I j = 0; j != T.rows(); ++j) {
                I j1 = T(j,n-1);
                for(I i = 0; i != P.rows(); ++i) {
                    vertices(i,j) = P(i,j1-1);
                }
            }
			for(I alfa = 1;alfa != nlb_trial +1;++alfa){
				for(I beta = 1;beta != nlb_test+1;++beta){
					F temp = 0;
					temp = Gauss_quadratrue_line_trial_test_2D(end_point_1,end_point_2,vertices,basis_type_trial,alfa,basis_type_test,beta,der_x_trial,der_y_trial,der_x_test,der_y_test,Robinfunc,Gpn);
					I location_i = Tb_trial(alfa-1,n-1);
					I location_j = Tb_test(beta-1,n-1);
					A1.coeffRef(location_i-1,location_j-1) += temp;
				}
			}
            for(I alfa=1; alfa != nlb_test+1; ++alfa) {
                F temp = 0;
                temp = Gauss_quadratrue_line_test_2D(end_point_1,end_point_2,vertices,basis_type_test,alfa,der_x_test,der_y_test,Neufunc,Gpn);
                I location = Tb_test(alfa-1,n-1);
                b1(location-1) += temp;
            }
        }
    }
}



} //end of namespace FEMethod
} //end of namespace PDWSC
#endif // end of BOUNDARY_H
