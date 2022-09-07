#ifndef Vector_h
#define Vector_h

#include <iostream>
#include <cmath>
#include <initializer_list>

namespace WHYSC {
namespace AlgebraObject {

template<typename F=double, typename I=int>
struct Vector 
{
    typedef F Float;
    typedef I Int;

    F * data;
    I size;
    bool T = false;

    /*
     * 默认构造函数
     */
    Vector()
    {
        data = NULL;
        size = 0;
    }

    Vector(I s, F val=0.0)
    {
        size = s;
        init(val);
    }

    Vector(const std::initializer_list<F> &l)
    {
        size = l.size();
        data = new F[size];
        I i = 0;
        for(auto & val: l)
        {
            data[i] = val;
            i++;
        }
    }

    void init(F val=0.0)
    {
        data = new F[size];
        for(I i=0; i < size; i++)
            data[i] = val;
    }

    ~Vector()
    {
        if(data != NULL)
            delete [] data;
    }

    F norm()
    {
        F sum = 0.0;
        for(I i=0; i < size; i++)
            sum += data[i]*data[i];
        return std::sqrt(sum);
    }

    F maxnorm()
    {
        F r = 0.0;
        for(I i=0; i < size; i++)
            if( r < std::abs(data[i]))
                r = std::abs(data[i]);
        return r;
    }

    F & operator[](const I i) 
    {
        return data[i];
    }

    const F & operator[](const I i) const
    {
        return data[i];
    }
    
    F House(Vector &v)
    {
    	
    	I n = size;
    	F sigma = 0;
    	F beta = 0;
    	F xnorm = this -> norm();
    	if(v.size == 1)
    	{
    		beta = 2;
    		v[0] = 1;
    	}
    	else
    	{
    	data[0] /= xnorm;
    	for(I i = 1;i < n;i++)
    	{
    		data[i] /= xnorm; 
    		sigma += data[i]*data[i];
    	}
    	for(I i = 0;i < n;i++)
    	{
    		v[i] = data[i];
    	}
    	if(sigma == 0)
    	{
    		if(data[0] < 0)
    		{
    			v[0] = 2*data[0];
    			beta = 2/v[0]*v[0];	
    		}
    		else
    		{
    			v[0] = 0;
    			beta = 0;
    		}
    	}
    	else
    	{
    		F alfa = this->norm();
    		if(data[0] < 0)
    		{
    			v[0] = data[0] - alfa;
    		}
    		else
    		{
    			v[0] = -sigma/(data[0]+alfa);
    		}
    		beta = 2/(v[0]*v[0] + sigma);
    	}
    	}
    	return beta;
    }
	void givens(F a,F b)
	{
		F c = 0;
		F s = 0;
		if(size == 2)
		{
			if( b == 0)
			{
				if(a >= 0)
				{
					c = 1;
					s = 0;	
				}
				else
				{
					c = -1;
					s = 0;
				}
			}
			else
			{
				if(abs(b)>abs(a))
				{
					F t = a/b;
					if(b>0)
					{
						s = 1/sqrt(1+t*t);
						c = s*t;
					}
					else
					{
						s = -1/sqrt(1+t*t);
						c = s*t;
					}
				}
				else
				{
					F t = b/a;
					if(a>=0)
					{
						c = 1/sqrt(1+t*t);
						s = c*t;
					}
					else
					{
						c = -1/sqrt(1+t*t);
						s = c*t;
					}
				}
			}
		}
		else
		{
			std::cout << "wrong input"<<std::endl;
		}
		data[0] = c;
		data[1] = s;
	}
	void operator = (const Vector &v1){
		auto n1 = v1.size;
		for(int i = 0;i < n1;++i){
			data[i] = v1[i];
		}
	}
};

template<typename F, typename I>
inline Vector<F, I> operator - (const Vector<F, I> & v0, 
        const Vector<F, I> & v1)
{
    Vector<F, I> r(v0.size);
    for(auto i=0; i < v0.size; i++)
    {
        r[i] = v0[i] - v1[i];
    }
    return r;
}
template<typename F, typename I>
inline Vector<F, I> operator + (const Vector<F, I> & v0, 
        const Vector<F, I> & v1)
{
    Vector<F, I> r(v0.size);
    for(auto i=0; i < v0.size; i++)
    {
        r[i] = v0[i] + v1[i];
    }
    return r;
}
template<typename F, typename I>
inline Vector<F, I> operator * (const Vector<F, I> & v0, 
        const Matrix<F, I> & m0)
{
    Vector<F, I> r(v0.size);
    for(auto j=0; j < m0.shape[1]; j++)
    {
    	for(auto i = 0;i < m0.shape[0];i++)
    	{
    		r[j] += v0[i]*m0[i][j];
    	}
    }
    return r;
}
template<typename F, typename I>
inline Vector<F, I> operator * (const Matrix<F, I> & m0, 
        const Vector<F, I> & v0)
{
    Vector<F, I> r(m0.shape[0]);
    for(auto i=0; i < m0.shape[0]; i++)
    {
    	for(auto j = 0;j < m0.shape[1];j++)
    	{
    		r[i] += v0[j]*m0[i][j];
    	}
    }
    return r;
}

template<typename F, typename I>
inline Matrix<F, I> operator * (const Vector<F, I> & v0, 
        const Vector<F, I> & v1)
{
    	Matrix<F,I> r(v0.size,v1.size);
    	for(I i = 0;i < v0.size;i++)
    	{
    		for(I j = 0; j < v1.size;j++)
    		{
    			r[i][j] = v0[i]*v1[j];
    		}
    	}
    	return r;
}
template<typename F, typename I>
inline Vector<F, I> operator * (const F a, 
        const Vector<F, I> & v1)
{
    	Vector<F,I> r(v1.size);
    	for(I i = 0;i < v1.size;i++)
    	{
    		r[i] = v1[i]*a;
    		
    	}
    	return r;
}
template<typename F, typename I>
inline Vector<F, I> operator * (const Vector<F, I> & v1,const F a)
{
    	Vector<F,I> r(v1.size);
    	for(I i = 0;i < v1.size;i++)
    	{
    		r[i] = v1[i]*a;
    		
    	}
    	return r;
}
template<typename F, typename I>
std::ostream& operator << (std::ostream & os, const Vector<F, I> & v)
{
    std::cout << "Vector("<< v.size <<")" << std::endl;

    for(I i = 0; i < v.size; i++)
    {
        os << v[i] << " ";
    }
    os << std::endl;
    return os;
}

} // end of namespace AlgebraObject

} // end of namespace WHYSC
#endif // end of Vector_h
