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
    /**	std::vector<T> tripletList;
    	for(int i = 0;i < 4;++i){
    		tripletList.push_back(T(i,i,5));
    	}
    	SpMat mat(4,4),mat2(4,4);
    	mat.setFromTriplets(tripletList.begin(),tripletList.end()); **/
    /* cout << mat << endl; */
    //生成一个5*5的空矩阵,让它最后一行变为1
    /**	SpMat mat1(5,5);
    	std::vector<T> tripletList1;
    	for(int i = 0;i < mat1.cols();++i){
    		tripletList1.push_back(T(i,i,i +1));
    	}
    	mat1.setFromTriplets(tripletList1.begin(),tripletList1.end());
    	cout << mat1 << endl;  **/
    //分块操作，让4*4矩阵被赋上5*5右下4*4的值
    /**	mat = mat1.block(1,1,4,4);
    	cout << mat << endl;  **/
    /* SpMat mat1(5,5); */
    /* std::vector<T> tripletList1; */
    /* for(int i = 0; i < mat1.cols(); ++i) { */
    /* 	if(i == 0){ */
    /* 		tripletList1.push_back(T(i,i+1,5)); */
    /* 		tripletList1.push_back(T(i,i+2,5)); */
    /* 	} */
    /* tripletList1.push_back(T(i,i,i +2)); */
    /* 	if(i+1 < mat1.cols()) */
    /* 		tripletList1.push_back(T(i+1,i,i+2)); */
    /* } */
    /* mat1.setFromTriplets(tripletList1.begin(),tripletList1.end()); */
    /* cout << mat1 << endl; */
    /* int cols = mat1.cols(); */
    /* Eigen::VectorXf v(cols); */
    /* v = v*0; */
    /* v(0) = 1; */
    /* cout << v << endl; */
    /* /1* mat1.topLeftCorner(1,cols) = v.transpose(); *1/ */
    /* for(int i = 0;i < cols;++i){ */
    /* 	if(i==0) */
    /* 		mat1.coeffRef(0,i) = 1; */
    /* 	else */
    /* 		mat1.coeffRef(0,i) = 0; */
    /* } */
    /* cout << mat1<< endl; */
    /* Eigen::MatrixXd M(5,5); */
    /* cout << M << endl; */
    /* Eigen::MatrixXd A = M*1+M*5; */
    /* cout << A << endl; */
    /* delete A; */
    /* Eigen::VectorXd x(20); */
    /* Eigen::VectorXd a(3); */
    /* a<<1,2,3; */
    /* x.head(3) = a; */
    /* cout << x << endl << a << endl; */
    /* x(0) = 15; */
    /* cout << x << endl << a << endl; */
    /* Eigen::VectorXd p(5); */
    /* p<<-1,1,5,-6,6; */
    /* cout<<p.minCoeff() << endl; */
    std::vector<T> tripletList,tripletList1;
    for(int i = 0; i < 4; ++i) {
        for(int j = 0; j < 4; ++j) {
            tripletList.push_back(T(i,j,i+j));
        }
    }
    SpMat mat(4,4),mat1(4,4);
    mat.setFromTriplets(tripletList.begin(),tripletList.end());
    cout << mat << endl;
    for(int i = 0; i < 4; ++i) {
        for(int j = 0; j < 4; ++j) {
            tripletList1.push_back(T(i,j,i+j+100));
        }
    }
    mat1.setFromTriplets(tripletList1.begin(),tripletList1.end());
    cout << mat1 << endl;
    SpMat A(8,8);
    std::vector<T> tripletList2,tripletList3;
    int rowstart = 0,colstart = 0;
    for (int k=0; k<mat.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(mat,k); it; ++it) {
            /* cout << "row:" << it.row()<<endl; */
            /* cout << "col:" << it.col()<<endl; */
            /* cout << "value:"<< it.value()<<endl; */
            if(it.value() != 0) {
                tripletList2.push_back(T(it.row()+rowstart,it.col()+colstart,it.value()));
            }
        }
    }
    A.setFromTriplets(tripletList2.begin(),tripletList2.end());
    cout << A << endl;
	SpMat A1(8,8);
	rowstart = 4;
	colstart = 4;
    for (int k=0; k<mat.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(mat1,k); it; ++it) {
            /* cout << "row:" << it.row()<<endl; */
            /* cout << "col:" << it.col()<<endl; */
            /* cout << "value:"<< it.value()<<endl; */
            if(it.value() != 0) {
                tripletList3.push_back(T(it.row()+rowstart,it.col()+colstart,it.value()));
            }
        }
    }
    A1.setFromTriplets(tripletList3.begin(),tripletList3.end());
    cout << A1 << endl;
	cout << A1+A << endl;

}
