#ifndef CSRMatrix_h
#define CSRMatrix_h

#include <iostream>
#include <cmath>
#include <initializer_list>
#include <cassert>
#include <string>
#include <algorithm>
#include "thirdparty/matplotlibcpp.h"

using namespace std;
namespace plt = matplotlibcpp;


namespace WHYSC {
namespace AlgebraObject {

template<typename F=double, typename I=int>
struct CSRMatrix
{
    typedef F Float;
    typedef I Int;

    F * data; // 非零元数组
    I * indices; // 非零元对应列指标数组
    I * indptr; // 非零元行起始位置数组
    I shape[2]; // 矩阵的阶数
    I nnz; // 非零元的个数

    static std::string format; // format = "csr"

    /*
     *
     * Notes
     * -----
     *  生成一个空稀疏矩阵.
     *
     */
    CSRMatrix()
    {
        data = NULL;
        indices = NULL;
        indptr = NULL;
        shape[0] = 0;
        shape[1] = 0;
        nnz = 0;
    }

    /*
     *
     * Notes
     * -----
     *  生成一个有 n 个非零元的稀疏矩阵, 形状为 (nrow, ncol)
     *  
     */

    CSRMatrix(
            const I nrow, 
            const I ncol, 
            const I nnz_, 
            const I row[], 
            const I col[], 
            const F data_[])
    {
        from_coo(nrow, ncol, nnz_, row, col, data_);
    }

    void from_coo(
            const I nrow,
            const I ncol,
            const I nnz_,
            const I row[],
            const I col[],
            const F data_[])
    {
        clear();
        nnz = nnz_;
        shape[0] = nrow;
        shape[1] = ncol;

        // 分配内存空间
        data = new F[nnz];
        indices = new I[nnz];
        indptr = new I[shape[0] + 1];

        std::fill(indptr, indptr + shape[0], 0);

        // 计算每一行非零元的个数 
        for(I n = 0; n < nnz; n++)
        {
            indptr[row[n]]++;
        }


        for(I i=0, cumsum =0; i < shape[0]; i++)
        {
            I temp = indptr[i];
            indptr[i] = cumsum;
            cumsum += temp;
        }
        indptr[shape[0]] = nnz;

        for(I n = 0; n < nnz; n++)
        {
            I i = row[n];
            I k = indptr[i];

            indices[k] = col[n];
            data[k] = data_[n];

            indptr[i]++;
        }

        for(I i = 0, last = 0; i <= shape[0]; i++)
        {
            I temp = indptr[i];
            indptr[i] = last;
            last = temp;
        }
    }


    template<typename Matrix> 
    CSRMatrix(Matrix & m)
    {
        shape[0] = m.shape[0];
        shape[1] = m.shape[1];
        nnz = 0;
        for(auto i=0; i < shape[0]; i++)
        {
            for(auto j=0; j < shape[1]; j++)
            {
                if( m[i][j] != 0)
                    nnz++;
            }
        }

        data = new F[nnz];
        indices = new I[nnz];
        indptr = new I[shape[0] + 1];
        indptr[shape[0]] = nnz;
        auto c = 0;
        for(auto i=0; i < shape[0]; i++)
        {
            indptr[i] = c;
            for(auto j=0; j < shape[1]; j++)
            {
                if( m[i][j] != 0)
                {
                    data[c] = m[i][j];
                    indices[c] = j;
                    c++;
                }
            }
        }
    }

    const F operator() (const I i, const I j) const
    {
        for(auto k = indptr[i]; k < indptr[i+1]; k++)
        {
            if(indices[k] == j)
            {
                return data[k];
            }
        }
        return 0.0;
    }

    void clear()
    {
        shape[0] = 0;
        shape[1] = 0;
        nnz = 0;
        if(data != NULL)
            delete [] data;

        if(indices != NULL)
            delete [] indices;

        if(indptr != NULL)
            delete [] indptr;
    }

    ~CSRMatrix()
    {
        if(data != NULL)
            delete [] data;

        if(indices != NULL)
            delete [] indices;

        if(indptr != NULL)
            delete [] indptr;
    }

    template<typename Vector>
    void jacobi(Vector & b, Vector & x, I maxit=100, F tol=1e-8)
    {
        F d = 0.0;
        F bnorm = b.norm();
        I count = 0;
        Vector r0 = b - (*this)*x;
        F r0norm = r0.norm();
        Vector x1(x.size);
        Vector rab(maxit);
        
       
        for(I m=0; m < maxit; m++)
        {
            for(I i=0; i < x.size; i++)
            {
                x1[i] = b[i];
                for(I k=indptr[i]; k < indptr[i+1]; k++)
                {
                    if(indices[k] == i)
                        d = data[k];
                    else
                        x1[i] -= data[k]*x[indices[k]];
                }
                x1[i] /= d;
            }

            for(I i=0; i< x.size; i++)
            {
                x[i] = x1[i];
            }

            Vector r = b - (*this)*x;
            F rnorm = r.norm();
            rab[m] = rnorm;
            count += 1;
            std::cout << m << ": " << rnorm << std::endl;
            if( rnorm < tol*bnorm)
                break;
        }
        std::vector<F> xx(count);
        std::vector<F> yy(count);
        xx.at(0) = 0;
        yy.at(0) = r0norm;
        
        for(I i =1 ;i < count;i++)
        {
        	xx.at(i) = i;
        	yy.at(i) = rab[i-1];
        }
        plt::named_plot("jacobi",xx,yy,"r--");
    }

    template<typename Vector>
    void gauss_seidel_forward(Vector & b, Vector & x, I maxit=100, F tol=1e-8)
    {
        F d = 0.0;
        F bnorm = b.norm();
        I count = 0;
        Vector r0 = b - (*this)*x;
        F r0norm = r0.norm();
        Vector rab(maxit);

        for(I m=0; m < maxit; m++)
        {
             for(I i=0; i < x.size; i++)
             {
                x[i] = b[i];
                for(I k = indptr[i]; k < indptr[i+1]; k++)
                {
                    if (indices[k]==i)
                        d = data[k];
                    else
                       x[i] -= data[k]*x[indices[k]];
                }

                x[i] /= d;
            }
            Vector r = b - (*this)*x;
            F rnorm = r.norm();
            rab[m] = rnorm;
            count += 1;
            std::cout << m << ": " << rnorm << std::endl;
            if( rnorm < tol*bnorm)
                break;
        }
        std::vector<F> xx(count);
        std::vector<F> yy(count);
        xx.at(0) = 0;
        yy.at(0) = r0norm;
        for(I i =1 ;i < count;i++)
        {
        	xx.at(i) = i;
        	yy.at(i) = rab[i-1];
        }
        plt::named_plot("gauss",xx,yy,"b--");
    }

    template<typename Vector>
    void gauss_seidel_backward(Vector & b, Vector & x, I maxit=100, F tol=1e-8)
    {
        F d = 0.0;
        F bnorm = b.norm();
        for(I m=0; m < maxit; m++)
        {
             for(I i = x.size-1; i > -1; i--)
             {
                x[i] = b[i];
                for(I k = indptr[i]; k < indptr[i+1]; k++)
                {
                    if (indices[k]==i)
                        d = data[k];
                    else
                       x[i] -= data[k]*x[indices[k]];
                }

                x[i] /= d;
            }


            Vector r = b - (*this)*x;
            F rnorm = r.norm();
            std::cout << m << ": " << rnorm << std::endl;
            if( rnorm < tol*bnorm)
                break;
        }
        
    }

    template<typename Vector>
    void SOR(Vector & b,Vector & x,F w,I maxit=1000,F tol=1e-8)
    {
    	F d = 0.0;
	F bnorm = b.norm();
	Vector x1(b.size);
	F Rnorm = 0.0;
	for(I k =0; k < maxit ;k++)
	{
		for(I i = 0;i < x.size ;i++)
		{
			x1[i] = b[i]; 
			for(I j = indptr[i];j < indptr[i+1];j++)
			{
				if(indices[j]==i)
				{
					d = data[j];
				}
				else
				{
					x1[i] -= x1[indices[j]]*data[j];
				}
			}
			x1[i] = (1-w)*x[i] + w*x1[i]/d; 	
		}
		for(I i = 0;i < x.size ;i++)
		{
			x[i] = x1[i];
		}
		Vector r = b - (*this)*x;
		F rnorm = r.norm();
		Rnorm = rnorm;
		if(rnorm < bnorm*tol){
			std::cout << k <<": "<<Rnorm<<std::endl;
			break;
		}
		else if(rnorm >= bnorm*tol && k==maxit-1){
			std::cout << k<<":" << Rnorm <<std::endl; 
		}
      }        	
    }
    
    template<typename Vector>
    void CG(Vector & b ,Vector &x,I maxit = 100,F tol=1e-8)
    {
	Vector r = b - (*this)*x;
	Vector p = b - (*this)*x;
	
	F alfa = 0.0;
	F beta = 0.0;
	F bnorm = b.norm();
	I count = 0;
        Vector r0 = b - (*this)*x;
        F r0norm = r0.norm();
        Vector rab(maxit);
	
	for(I j = 0;j < maxit;j++)
	{
		F pap = 0.0;
		F rj = r.norm()*r.norm();
		Vector q = (*this)*p;
		for(I i = 0;i < q.size;i++)
		{
			pap += p[i]*q[i];
		}
		alfa = rj/pap;
		for(I i = 0;i < q.size;i++)
		{
			x[i] = x[i] + alfa*p[i];
		}
		Vector AP = (*this)*p;
		for(I i = 0;i < r.size;i++)
		{
			r[i] -= alfa*AP[i];
		}
		F rnorm = r.norm();
		std::cout << j <<": "<<rnorm<<std::endl;
		rab[j] = rnorm;
		count += 1;
		if(r.norm()< bnorm*tol)
		{
			break;
		}
		beta = r.norm()*r.norm()/rj;
		for(I i = 0;i < p.size;i++)
		{
			p[i] = r[i] + beta*r[i];
		} 
		     	
	} 
	std::vector<F> xx(count);
        std::vector<F> yy(count);
        xx.at(0) = 0;
        yy.at(0) = r0norm;
        for(I i =1 ;i < count;i++)
        {
        	xx.at(i) = i;
        	yy.at(i) = rab[i-1];
        }
        plt::named_plot("CG",xx,yy,"g--");		  	
    }
};

template<typename F, typename I>
std::string CSRMatrix<F, I>::format = "csr";


template<typename F, typename I, typename Vector>
inline Vector operator * (const CSRMatrix<F, I> & m, 
        const Vector & v0)
{
    assert( m.shape[1] == v0.size );
    Vector v1(m.shape[0], 0.0);
    for(auto i = 0; i < m.shape[0]; i++)
    {
        for(auto k = m.indptr[i]; k < m.indptr[i+1]; k++)
        {
            auto j = m.indices[k]; // 列指标
            v1[i] += m.data[k]*v0[j];
        }
    }
    return v1;
}

template<typename F, typename I>
std::ostream& operator << (std::ostream & os, 
        const CSRMatrix<F, I> & m)
{
    std::cout << "CSRMatrix("<< m.shape[0] << ","
        << m.shape[1] << "):" << std::endl;

    for(auto i = 0; i < m.shape[0]; i++)
    {
        for(auto j = m.indptr[i]; j < m.indptr[i+1]; j++)
        {
            os << "(" << i << ", " << m.indices[j] << 
                ") " << m.data[j] << std::endl;
        }
    }
    os << "\n";
    return os;
}

} // end of namespace AlgebraObject
} // end of namespace WHYSC

#endif // end of CSRMatrix_h
