#include <iostream>
#include <cmath>
#include "algebra/Algebra_kernel.h"
#include "FEMethod/FEMethod_kernel.h"
#include "Constants.h"
#include <Eigen/Eigen>
typedef typename WHYSC::Algebra_kernel<double, int> Kernel;
typedef typename PDWSC::FEMethod_kernel<double, int> FEKernel;
typedef typename Kernel::Matrix Matrix;
typedef typename Kernel::Vector Vector;
typedef typename Kernel::CSRMatrix CSRMatrix;
double a = 100;
double D = -1/a;
double c(double x,double y)
{
    return 1;
}
double f(double x,double y)
{
    double r = -2*a*exp(a*(2*x-1))*(3*exp(a*(2*x-1))-1)/pow(1+exp(a*(2*x-1)),3);
    return r;
}


double analytic_solution(double x,double y)
{
    return 1/(1+exp(a*(2*x-1)));
}
double analyticSolution_1_der(double x,double y)
{
    return 1;
}
int main() {
	using namespace Eigen;
	using namespace std;
/*****************************************************/	
	/* MatrixXd m(2,2); */
	/* m(0,0) = 3; */
	/* m(1,0) = 2.5; */
	/* m(0,1) = -1; */
	/* m(1,1) = m(0,1)+m(1,0); */
	/* cout << m << endl; */

/*******************************************************/	
	/* int N1 = 2,N2 = 2; */
	/* MatrixXd P(2,(N1+1)*(N2+1)),T(3,N1*N2*2); */
	/* cout<<P <<endl; */
	/* cout << T <<endl; */



/*******************************************************/
	/* Matrix3d m = Matrix3d::Random(); */
	/* m = (m + Matrix3d::Constant(1.2))*50; */
	/* cout << "m = " << endl << m << endl; */
	/* Vector3d v(1,2,3) */
	/* VectorXd v1(3); */
	/* v << 1,2,3; */
	/* cout << "m*v="<<endl << m*v <<endl; */


/********************************************************/
	/* Matrix3d a; */
	/* a(0,0) = 1; */
	/* a << 1,2,3, */
	/*      4,5,61, */
	/* 	 7,8,9; */
	/* cout << a(0,0) << endl; */
	/* VectorXd v(5); */

	/* cout << v[4] << endl; */

/*********************************************************/
	/* MatrixXd m(2,5); */
	/* m.resize(1,2); */
	/* cout << "m行数" << m.rows() <<endl; */
	/* cout << "m列数" << m.cols() <<endl; */

/*********************************************************/
	/* MatrixXd m1(2,2); */
	/* MatrixXd m2(3,3); */
	/* m1 = m2; */

/***转置和共轭*****************************************************/
	/* Matrix2i b; */
	/* b << 1,2,3,4; */
	/* cout << "b" << b <<endl; */
	/* cout << "转置" << b.transpose() << endl; */
	/* cout << "共轭" << b.conjugate() << endl; */
   /* cout << "转置共轭" << b.adjoint() << endl; */
	/* b.transposeInPlace(); */
	/* cout << "b" << b <<endl; */ 	

/************************点积和叉积******************************/
	Vector3d v(1,2,3);
	Vector3d w(0,1,2);


	cout << "点积" << v.dot(w) <<endl;
	double dp = v.adjoint()*w;
	cout << "通过矩阵乘积的点积" << dp <<endl;
	cout << "叉积：\n" <<v.cross(w) << endl;

/******************************************
 *相关函数：
 mat.sum();    //求所有元素和
 mat.prod();   //求所有元素积
 mat.mean();   //求所有元素平均值
 mat.minCoeff();  //求最小值
 mat.maxCoeff();  // 求最大值
 mat.trace()      // 求迹
 * *****************************************/	

/*******************稀疏矩阵*******************************************/
	typedef Eigen::SparseMatrix<double> SpMat;
	SpMat U(2,2);
	/* cout << U <<endl; */
	/* VectorXi h = VectorXi::Constant(5,6); */
	/* cout << h <<endl; */
	/* U.reserve(h); */
	/* cout << U <<endl; */	
	/* U.insert(0,0) = 1; */
	/* cout << U <<endl; */
	/* U.makeCompressed(); */
	/* cout << U <<endl; */
	/* U.insert(0,1) = 2; */
	/* U.coeffRef(0,0) += 3; */
	/* U.coeffRef(0,1) += 2; */
	/* U.coeffRef(0,0) += 4; */
	/* cout << U <<endl; */
	/* U.makeCompressed(); */
	/* cout << U <<endl; */
	for(int i = 0; i < 2;++i){
		for(int j = 0;j < 2;++j){
			if(j+i != 0)
				U.coeffRef(i,j) += j+i;
			U.makeCompressed();
			cout << U <<endl;
		}
	}
	Vector2d b,x;
	b << 2,5;
	cout << "b" << b <<endl;
	LeastSquaresConjugateGradient<SpMat> solver;
	solver.compute(U);
	cout <<"info1"<< solver.info() << endl;
	/* Matrix2d I; */
	/* I << 1,0,0,1; */
	/* cout << "UI" << U*I <<endl; */
	x = solver.solve(b);
	cout << "info2" << solver.info() <<endl;
	cout << "x"<<x <<endl;
	cout << "Ux" << U*x <<endl;


