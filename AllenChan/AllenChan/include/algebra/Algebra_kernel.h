#ifndef AlgebraKernel_h
#define AlgebraKernel_h

#include "Matrix.h"
#include "CSRMatrix.h"
#include "Vector.h"
#include "linalg.h"
#include "MatrixFactory.h"

namespace WHYSC {

template<typename F=double, typename I=int>
class Algebra_kernel
{
public:
    typedef typename AlgebraObject::Matrix<F, I> Matrix;
    typedef typename AlgebraObject::Vector<F, I> Vector;

    typedef typename AlgebraObject::CSRMatrix<F, I> CSRMatrix;
    //typedef AlgebraAlgrithom::lu<Matrix> lu;
public:
    static void lu(Matrix & A, Matrix & L, Matrix & U)
    {
        return AlgebraAlgrithom::lu<Matrix>(A, L, U);
    }

    static void qr_gs(Matrix & A, Matrix & Q, Matrix & R)
    {
        return AlgebraAlgrithom::qr_gs<Matrix>(A, Q, R);
    }
	 static void qr_mgs(Matrix & A, Matrix & Q, Matrix & R)
    {
        return AlgebraAlgrithom::qr_mgs<Matrix>(A, Q, R);
    }
     static void qr_householder(Matrix & A, Matrix & Q, Matrix & R)
    {
        return AlgebraAlgrithom::qr_householder<Matrix>(A, Q, R);
    }
     static void qr_givens(Matrix & A, Matrix & Q, Matrix & R)
    {
        return AlgebraAlgrithom::qr_givens<Matrix>(A, Q, R);
    }
    static void laplace_1d(const I n, CSRMatrix & m)
    {
        return AlgebraObject::laplace_1d(n, m);
    }

    static void laplace_2d(const I N, CSRMatrix & m)
    {
        return AlgebraObject::laplace_2d(N, m);
    }

    static void laplace_3d(const I M, CSRMatrix & m)
    {
        return AlgebraObject::laplace_3d(M, m);
    }
};

} // end of namespace WHYSC
#endif // end of AlgebraKernel_h
