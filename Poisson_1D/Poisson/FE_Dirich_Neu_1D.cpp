#include <iostream>
#include <cmath>
#include "algebra/Algebra_kernel.h"
#include "FEMethod/FEMethod_kernel.h"
#include "Constants.h"
using namespace std;
typedef typename WHYSC::Algebra_kernel<double, int> Kernel;
typedef typename PDWSC::FEMethod_kernel<double, int> FEKernel;
typedef typename Kernel::Matrix Matrix;
typedef typename Kernel::Vector Vector;
typedef typename Kernel::CSRMatrix CSRMatrix;
// 求解：-d(c(x)du(x)/dx)/dx = f(x) (a <= x <= b)  D边界条件:ga gb N边界条件:ra,rb   R边界条件：u'(b) + qb*u(b) = pb
double c(double x)
{
    return exp(x);
}
double f(double x)
{
    double r = -exp(x)*(cos(x) - 2*sin(x) - x*cos(x) - x*sin(x));
    return r;
}
double g(double x)
{
    return cos(x);
}
double r(double x)
{
    return 1;
}
double analyticSolution(double x)
{
    return x*cos(x);
}
double analyticSolution_1_der(double x)
{
    return cos(x) - x*sin(x);
}
int main() {
    //生成信息矩阵P,T,Pb,Tb***********************************************************************
    double a = 0,b = 1;
	int N = 128;
    Matrix P(1,N+1),T(2,N);//网格信息矩阵
    FEKernel::generate_PT_1d(a,b,N,P,T);


	//101代表使用线性元
    /* int basis_type_trial = 101; */
    /* int basis_type_test = 101; */
	/* int Nlb = 2; */
    /* Matrix Pb_trial(1,N+1),Tb_trial(2,N),Pb_test(1,N+1),Tb_test(2,N); */

	//换成102，可使用二次元但一些参数的准备也需要更改
    int basis_type_trial = 102;
    int basis_type_test = 102;
	int Nlb = 3;
    Matrix Pb_trial(1,2*N+1),Tb_trial(3,N),Pb_test(1,2*N+1),Tb_test(3,N);


    int Nlbtest=Nlb;
    int Nlbtrial=Nlb;//trial和test用的是一个元
    int N_test = N+1;
    int N_trial = N+1;
    FEKernel::generate_PbTb_1d(a,b,N,basis_type_trial,N_trial,Nlbtrial,Pb_trial,Tb_trial);
    FEKernel::generate_PbTb_1d(a,b,N,basis_type_test,N_test,Nlbtest,Pb_test,Tb_test);

    //生成边界信息矩阵。矩阵的第一行储存边界条件类型的信息：-1表示Dirichlet，-2表示Neumann,-3表示Robin.第二行存全局编号
    Matrix boundarynodes(2,2);
	
	//只有两个边界点，使用直接赋值也很简单
    FEKernel::generate_boundraynodes(boundarynodes,Tb_test);
    boundarynodes[0][0] = -1;
    boundarynodes[0][1] = -1;
    boundarynodes[1][0] = 1;
    if(basis_type_test == 101)//这里的边界点是指基函数的边界点而不是网格的
        boundarynodes[1][1] = N + 1;
    if(basis_type_test == 102)
        boundarynodes[1][1] = 2*N + 1;

    //组装A*******************************************************************************
    Matrix A(N_test,N_test);
    FEKernel::assembleA_1d(N,Nlbtest,Nlbtrial,basis_type_trial,basis_type_test,P,T,Tb_trial,Tb_test,A,c,1,1);


    //组装v*******************************************************************************
    Vector v(N_test);
    FEKernel::assembleV_1d(N,Nlbtest,basis_type_test,0,P,T,Pb_test,Tb_test,v,f);


    //处理边界条件************************************************************************
    FEKernel::treat_boundary_1d(boundarynodes,A,Pb_trial,v,g);//处理Dirichle边界条件，遍历边界矩阵，擦除，对角线赋值为1，修改b中相应位置的值
    v[0] = 0;//本例的要求



	//处理Neumann边界条件(本算例只有Drichlet边界条件，以下处理Neumann边界条件的代码无效，如果有Neumann边界条件，在边界点的信息矩阵中把边界点类型写为-2,那么以下代码就可以处理Neumann边界条件)****************************************************************
    if(boundarynodes[0][0] == -2)//在右端的b向量中相应的位置加上相应的项即可
        v[0] += (-r(a)*c(a));
    if(boundarynodes[0][1] == -2)
        v[N_test-1] += r(b)*c(b);


    //求解***********************************************************************************
    Vector x(N_test);
    CSRMatrix R(A);
    R.SOR(v,x,1,150000);



    //计算误差****************************************************************************
    //装填解析解向量
     Vector r(x.size); 
	 for(int i=0;i != r.size;++i){ 
     r[i] = analyticSolution(Pb_trial[0][i]); 
     }

	 //计算无穷范
     double error1 = (x - r).maxnorm();
	 cout << error1 <<endl;
    //取0和analyticSolution时可以计算2范数，
	double error2 = FEKernel::compute_Hs_error(N,P,T,Tb_test,basis_type_test,0,x,analyticSolution);
    cout << error2 <<endl;

	//取1和analyticSolution的导数时可以计算H1半范
    double error3 = FEKernel::compute_Hs_error(N,P,T,Tb_test,basis_type_test,1,x,analyticSolution_1_der);
    cout << error3 <<endl;


}
