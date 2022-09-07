#ifndef BOUNDARY_H
#define BOUNDARY_H

namespace PDWSC{
namespace FEMethod{

template<typename F,typename I,typename Matrix,typename Vector>
void treat_boundary_1d(Matrix &boundarynodes,Matrix &A,Matrix &Pb,Vector &b,F (*g)(F))
{
	for(I k = 0;k != boundarynodes.cols();++k){
		if(boundarynodes(0,k) == -1){
			I i = boundarynodes(1,k);
			for(I j = 0;j != A.cols();++j){
				A(i-1,j) = 0;
			}
			A(i-1,i-1) = 1;
			b(i-1) = g(Pb(0,i-1));
		}
	}
}


} //end of namespace FEMethod
} //end of namespace PDWSC
#endif // end of BOUNDARY_H
