#ifndef ADDSPMAT_H
#define ADDSPMAT_H

#include "GaussQuadratrue.h"
#include <Eigen/Sparse>
#include <vector>

namespace PDWSC {
namespace FEMethod {

template<typename I,typename SpMat>
void addSpMat(SpMat &A,SpMat &M,I rowstart,I colstart)
{

    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    SpMat A1(A.rows(),A.cols());
    for (I k=0; k<M.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(M,k); it; ++it) {
            if(it.value() != 0) {
                tripletList.push_back(T(it.row()+rowstart,it.col()+colstart,it.value()));
            }
        }
    }
	A1.setFromTriplets(tripletList.begin(),tripletList.end());
	A += A1;
}




}  //end of namespace FEMethod
}  //end of namespace PDWSC
#endif
