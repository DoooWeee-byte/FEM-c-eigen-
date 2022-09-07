#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <vector>
#include <iostream>

typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double
typedef Eigen::Triplet<double> T;
using namespace std;

int main() {
/*    std::vector<T> tripletList;
    for(int i = 0; i < 4; ++i) {
        tripletList.push_back(T(i,i+1,i+2));
    }
    tripletList.push_back(T(0,1,5));
    tripletList.push_back(T(1,2,-5));*/
	/* tripletList.push_back(T(1,2,0)); */
/*	tripletList.push_back(T(0,0,0.5));

	tripletList.push_back(T(0,0,0.5));
    SpMat mat(5,5),mat1(5,5);*/
    /* for(auto i : tripletList) { */

    /*     cout <<"i"<< i.row() << endl; */
    /*     cout <<"j"<< i.col() << endl; */
    /*     cout << "value"<<i.value() << endl; */
		/* cout << "*********************"<<endl; */
    /* } */
   /* mat.setFromTriplets(tripletList.begin(),tripletList.end());
    mat1 = mat;
    std::cout << mat1 << std::endl;*/

	/*****************关于向量分块操作的测试
	Eigen::ArrayXf v(6);
	v << 1,2,3,4,5,6;
	cout << "v.head(3) =" << endl << v.head(3) << endl << endl;
	cout << "v.tail<3>()=" << endl << v.tail<3>() <<endl << endl;
	v.segment(1,4) *= 2;
	cout << "after 'v.segment(1,4) *= 2',v=" << endl << v << endl;
*********************************************************/

	/****************关于向量分块赋值的操作
	Eigen::VectorXf v(6);
	Eigen::VectorXf v1(5);
    v1 << 1,2,3,4,5;
	v.head(5) = v1;
	cout << v << endl;
****************************************/


	/******************关于稀疏矩阵的分块赋值操作*/
	//先生成一个4*4的矩阵
	std::vector<T> tripletList;
	for(int i = 0;i < 4;++i){
		tripletList.push_back(T(i,i,5));
	}
	SpMat mat(4,4),mat2(4,4);
	mat.setFromTriplets(tripletList.begin(),tripletList.end());
	/* cout << mat << endl; */
	//生成一个5*5的空矩阵,让它最后一行变为1
	SpMat mat1(5,5);
	std::vector<T> tripletList1;
	for(int i = 0;i < mat1.cols();++i){
		tripletList1.push_back(T(mat1.rows()-1,i,1));
	}
	mat1.setFromTriplets(tripletList1.begin(),tripletList1.end());
	/* cout << mat1 << endl; */
	//分块操作，让5*5矩阵的上4*4被赋值
	Eigen::MatrixXd M(5,5);
	for(int i = 0;i < 5;++i){
		M(i,i) = 1;
	}

	cout << M*mat1 << endl;
}
