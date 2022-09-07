#include <iostream>
#include <cmath>
#include <string>
#include <algebra/Possion1d_Kernel.h>
#include "thirdparty/matplotlibcpp.h"
using namespace std;
namespace plt = matplotlibcpp;


//如Afunc是一个常函数，然后k0,kL为0，又是均匀剖分，则A是一个奇异矩阵
double Afunc(double x)
{
	return x;
}

double func(double x)
{
	return -6*pow(x,2)+24*x;
}

int main()
{
	/*  -(au`)` = f  , x 属于 I=[0,L]
	 *   au`(0) = k0(u(0)-g0)
	 *  -au`(L) = kL(u(L)-gL)
	 *  本例使用题目：
	 *  -xT'' = -6x^2 + 24x;
	 *  T'(0) = 0;
	 *  T'(8) = 0; 0<x<8;
	 *  真解为T = X^3 -12X^2
	
	  n：区间分割的段数
	  L：区间终点的坐标（起点为0）
	  k0,kL,g0,gL都是给定的参数
	*/
	int N = 10;
	double L = 8;
	double k0 = 0;
	double kL = 0;
	double g0 = 0;
	double gL = 0;
	
	//假设是一个均匀剖分
	double h = L/N;
	
	//生成刚度矩阵和重载向量的容器
//	SourceAssembler1D I(n+1,n+1);
//	StiffnessAssembler1D A(n+1,n+1);

	//组装刚度矩阵和重载向量
/*	I.Assemble(h,func);
	A.Assemble(h,Afunc);
	A[0][0] += k0;
	A[n][n] += kL;
	I[0]  += k0*g0;
        I[n]  += kL*gL;

	//输出装载好的刚度矩阵和重载向量
	cout << A <<endl;
	cout << I <<endl;
		
	//利用GaussSeidel迭代求方程组
	Vector x(n+1);
	GaussSeidel(A,I,x,10000);
	//输出方程组的求解结果*/
	// A*I;
/*	Matrix A(2,3);
	A[0][0] = 1;
	A[0][1] = 2;
	A[1][1] = 1;
	Vector I(3);
	I[0] = 1;
	I[1] = 1;
	I[2] = 1;
//	cout << I <<endl;
	Vector I1 = A*I;
 	cout <<I1 <<endl;	 */
/*	Matrix P(1,N+1);
	for(int i=0;i < N+1;i++)
	{
		double H = h*i;
		P[0][i] = H;
	}
	Matrix T(2,N);
	for(int j=0;j<N;j++)
	{
		for(int i=0;i<2;i++)
		{
			T[i][j] = j+i; 
		}
	}*/
/*	cout << P <<endl;
	cout << T <<endl;*/
/*	Vector x(4);
	x[0]= 1;
	x[1]= 1;
	x[2]= 1;
	x[3]= 1;
	Vector x1(4);*/
/*	Assembler1D A(N+1,N+1);
	A.Assemble(P,T,Afunc);
	cout << A <<endl;
	SourceAssembler1D I(N+1,N+1);
	I.Assemble(h,func);
	cout << I <<endl;
	Vector x(I.size);
	A.SOR(I,x,1.5,1000);
	cout<<"x\n"<<x<<endl;*/

/*	A.Gauss_Seidel(x,x1,1000);
	A.SOR(x,x1,1.5);
	cout << x1<<endl;
	cout << P <<endl;
	cout << T <<endl;*/
//	FE_solver_1D_Poisson();
//	Matrix A(2,5);
//	Matrix C(1,5);
//	C.generate_Pb(0,10,101);
//	A.generate_Tb(101);
//	cout << C<<endl;
//	cout << A <<endl;
	Matrix P(1,5);
	Matrix T(2,4);
	generate_PT_1D(P,T,0,1,101);
	cout << P<<endl;
	cout <<T <<endl;	
}
