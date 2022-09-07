#ifndef linalg_h
#define linalg_h

#include <iostream>
#include <initializer_list>

namespace WHYSC {
namespace AlgebraAlgrithom {

template<typename Matrix>
void lu(Matrix & A, Matrix & L, Matrix & U)
{
    typedef typename Matrix::Int  I;
    I n = A.shape[0];
    for(I k=0; k<n; k++)
	{
		for (I i=k+1; i<n; i++)
		{
			L[i][k] = A[i][k]/A[k][k];
		}

		for(I j=k; j<n; j++)
		{
			U[k][j] = A[k][j];
		}

		for(I i=k+1; i<n ;i++)
		{
			for(I j=k+1;j<n;j++)
			{
				A[i][j] -= L[i][k]*U[k][j];
			}
		}

	}

    for(I i=0;i<n;i++)
        L[i][i] =1 ;

    return;
}

template<typename Matrix>
void lu_kij(Matrix & A)
{
    typedef typename Matrix::Int  I;
    I n = A.shape[0];
    for(I k=0; k<n; k++)
	{
		for (I i=k+1; i<n; i++)
		{
			A[i][k] = A[i][k]/A[k][k];
		}

		for(I i=k+1; i<n ;i++)
		{
			for(I j=k+1;j<n;j++)
			{
				A[i][j] -= A[i][k]*A[k][j];
			}
		}

	}

    return;
}

template<typename Matrix>
void lu_ikj(Matrix & A)
{
    typedef typename Matrix::Int  I;
    I n = A.shape[0];
    for(I i = 0;i < n; i++)
    	{
   		for(I k = 0;k < i; k++)
   		{
   			A[i][k]=A[i][k]/A[k][k];
   			for(I j = k + 1; j < n; j++)
   			{
   				A[i][j] -= A[i][k]*A[k][j];
   			}
   		}
    	}
    return;
}

template<typename Matrix>
void lu_kij(Matrix & A, Matrix & L, Matrix & U)
{
    typedef typename Matrix::Int  I;
    I n = A.shape[0];
    for(I k=0; k<n; k++)
	{
		for (I i=k+1; i<n; i++)
		{
			L[i][k] = A[i][k]/A[k][k];
		}

		for(I j=k; j<n; j++)
		{
			U[k][j] = A[k][j];
		}

		for(I i=k+1; i<n ;i++)
		{
			for(I j=k+1;j<n;j++)
			{
				A[i][j] -= L[i][k]*U[k][j];
			}
		}

	}

    for(I i=0;i<n;i++)
        L[i][i] =1 ;

    return;
}


template<typename Matrix>
void qr_gs(Matrix & A, Matrix & Q, Matrix & R)
{
    typedef typename Matrix::Int  I;
    I m = A.shape[0];
    I n = A.shape[1];

    R[0][0] = A.col_norm_l2(0);

    for(I i = 0; i < m; i++)
        Q[i][0] = A[i][0]/R[0][0];

    for(I j = 1; j < n; j++)
    {
        for(I i = 0; i < m; i++)
            Q[i][j] = A[i][j];

        for(I k =0; k < j; k++)
        {
            for(I i = 0; i < m; i++)
            {
                R[k][j] += Q[i][k]*Q[i][j];
            }

            for(I i = 0; i < m; i++)
                Q[i][j] -= R[k][j]*Q[i][k]; 
        }

        R[j][j] = Q.col_norm_l2(j);

        for(I i = 0; i < m; i++)
            Q[i][j] /= R[j][j];
    }
}  

template<typename Matrix>
void qr_gs(Matrix & A)
{
}

template<typename Matrix>
void qr_mgs(Matrix & A, Matrix & Q, Matrix & R)
{
    typedef typename Matrix::Float  F;	
    typedef typename Matrix::Int  I;
    typedef typename  WHYSC::AlgebraObject::Vector<F,I>  Vector;
	Vector a1(A.shape[0]);
	for(I i = 0;i < A.shape[0];i++)
	{
		a1[i] = A[i][0];
	}
	double a1_norm = a1.norm();
	if(a1_norm  == 0)
	{
		for(I i = 0;i < Q.shape[0];i++)
		{
			Q[i][0] = 0;
		}
	}
	else
	{
		R[0][0] = a1_norm;
		for(I i = 0;i < Q.shape[0];i++)
		{
			Q[i][0] = a1[i]/a1_norm;
		}
	}
	I n = Q.shape[1];
	for(I k = 1;k < n;k++)
	{
		for(I i = 0;i < Q.shape[0];i++)
		{
			Q[i][k] = A[i][k];
		}
		for(I i = 0;i < k;i++)
		{
			for(I j = 0;j< Q.shape[0];j++)
			{
				R[i][k] += Q[j][i]*Q[j][k];
			}
			for(I m = 0;m < Q.shape[0];m++)
			{
				Q[m][k] = Q[m][k] - R[i][k]*Q[m][i];
			}
		}
		Vector qk(Q.shape[0]);
		for(I i = 0;i < Q.shape[0];i++)
		{
			qk[i] = Q[i][k];
		}
		double qk_norm = qk.norm();
		if(qk_norm != 0)
		{
			R[k][k] = qk_norm;
			for(I i = 0;i < Q.shape[0];i++)
			{
				Q[i][k] /= R[k][k];
			}
		}
	}
}

template<typename Matrix>
void qr_mgs(Matrix & A)
{

}

template<typename Matrix>
void qr_householder(Matrix & A, Matrix & Q, Matrix & R)
{
    typedef typename Matrix::Float  F;
    typedef typename Matrix::Int  I;
    typedef typename  WHYSC::AlgebraObject::Vector<F,I>  Vector;
    I m = A.shape[0];
    I n = A.shape[1];
    for(I i = 0;i < m;i++)
    {
    	for(I j = 0;j < m;j++)
    	{
    		if(i == j)
    			Q[i][j] = 1;
    		else
    			Q[i][j] = 0;
    	}
    }
    for(I k = 0;k < n;k++)
    {
    	Vector x(m-k);
    	Vector v(m-k);
    	for(I i = 0;i < m-k;i++)
    	{
    		x[i] = A[i+k][k];
    	}
    	F beta = x.House(v);
    	Matrix A1(m-k,n-k);
    	for(I i = 0;i < m-k;i++)
    	{
    		for(I j = 0; j < n-k;j++)
    		{
    			A1[i][j] = A[i+k][j+k];
    		}
    	}
    	Matrix Q1(m,m-k);
    	for(I i = 0;i < m;i++)
    	{
    		for(I j = 0;j < m-k;j++)
    		{
    			Q1[i][j] = Q[i][j+k];
    		}
    	}
    	Matrix A2 = A1-beta*v*(v*A1);
    	Matrix Q2 = Q1 - beta*(Q1*v)*v;
    	for(I i = 0;i < m-k;i++)
    	{
    		for(I j = 0; j < n-k;j++)
    		{
    			A[i+k][j+k] = A2[i][j];
    		}
    	}
    	for(I i = 0;i < m;i++)
    	{
    		for(I j = 0;j < m-k;j++)
    		{
    			Q[i][j+k] = Q2[i][j];
    		}
    	}
    }
}

template<typename Matrix>
void qr_givens(Matrix & A, Matrix & Q, Matrix & R)
{
    typedef typename Matrix::Float  F;
    typedef typename Matrix::Int  I;
    typedef typename  WHYSC::AlgebraObject::Vector<F,I>  Vector;
	I m = A.shape[0];
    I n = A.shape[1];
    for(I i = 0;i < m;i++)
    {
    	for(I j = 0;j < m;j++)
    	{
    		if(i == j)
    			Q[i][j] = 1;
    		else
    			Q[i][j] = 0;
    	}
    }
    for(I k = 0;k < n;k++)
    {
    	for(I i = k+1;i < m;i++)
    	{
 			Vector v(2);
 			v.givens(A[k][k],A[i][k]);
 			F c = v[0];
 			F s = v[1];
 			Matrix G(2,2);
 			G[0][0] = c;
 			G[0][1] = s;
 			G[1][0] = -s;
 			G[1][1] = c;
 			Matrix GT(2,2);
 			GT[0][0] = c;
 			GT[0][1] = -s;
 			GT[1][0] = s;
 			GT[1][1] = c;
 			Matrix A1(2,n-k);
 			Matrix Q1(m,2);
 			for(I r = 0 ;r < 2;r++)
 			{
 				for(I l = 0;l < n-k;l++)
 				{
 					if(r == 0)
 						A1[r][l] = A[k][l+k];
 					else
 						A1[r][l] = A[i][l+k];
 				}
 			}
 			for(I l = 0 ;l < 2;l++)
 			{
 				for(I r = 0;r < m;r++)
 				{
 					if(l == 0)
 						Q1[r][l] = Q[r][k];
 					else
 						Q1[r][l] = Q[r][i];
 				}
 			}
 			Matrix A2 = G*A1;
 			Matrix Q2 = Q1*GT;
 			for(I r = 0 ;r < 2;r++)
 			{
 				for(I l = 0;l < n-k;l++)
 				{
 					if(r == 0)
 						A[k][l+k] = A2[r][l];
 					else
 						A[i][l+k] = A2[r][l];
 				}
 			}
 			for(I l = 0 ;l < 2;l++)
 			{
 				for(I r = 0;r < m;r++)
 				{
 					if(l == 0)
 						Q[r][k] = Q2[r][l];
 					else
 						Q[r][i] = Q2[r][l];
 				}
 			}	
    	}
    }
    A[A.shape[0]-1][A.shape[1]-1] *= -1;
    for(I r = 0;r < Q.shape[0];r++)
    {
    	Q[r][Q.shape[1]-1] *= -1;
    }
}

} // end of AlgebraAlgrithom
} // end of WHYSC
#endif // end of linalg_h
