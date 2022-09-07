#ifndef BASEFUNCTION
#define BASEFUNCTION
namespace PDWSC {
namespace FEMethod {


template<typename F,typename I,typename Matrix>
F basis_local_function_1d(F x,const Matrix &vertices,I basis_type,I basis_index,I der_x) {
// x:自变量；lower_bound:当前遍历区间的左端点坐标；upper_bound:当前遍历区间的右端点坐标；basis_type:基函数的类型，101代表线性元；der_x:x的导数阶
    F ret = 0;
    if(basis_type==101) {
        F upper_bound = vertices(0,1);
        F lower_bound = vertices(0,0);
        F h = vertices(0,1) - vertices(0,0);
        if(basis_index == 1) {
            if(der_x == 0) {
                ret = (upper_bound - x)/h;
            }
            else if(der_x == 1) {
                ret = -1/h;
            }
            else if(der_x >= 2) {
                ret = 0;
            }
            else {
				cout << " Wrong der_x!" << endl;
            }
        }
        else if(basis_index == 2) {
            if(der_x == 0) {
                ret = (x - lower_bound)/h;
            }
            else if(der_x == 1) {
                ret = 1/h;
            }
            else if(der_x >= 2) {
                ret = 0;
            }
            else {
                cout << " Wrong der_x!" << endl;
            }
        }
        else {
            cout << "Wrong basis_index!" << endl;
        }
    }
    else if(basis_type == 102) {
        F upper_bound = vertices(0,1);
        F lower_bound = vertices(0,0);
        F h = upper_bound - lower_bound;
		F y = (x - lower_bound)/h;
        if(basis_index == 1) {
            if(der_x == 0) {
                ret = 2*y*y-3*y+1;
            }
            else if(der_x == 1) {
                ret = 4*y/h-3/h;
            }
            else if(der_x == 2) {
                ret = 4/(h*h);
            }
            else if(der_x >= 3) {
                ret = 0;
            }
            else{
                cout << " Wrong der_x!" <<endl;
            }
        }
        else if(basis_index == 2) {
            if(der_x == 0) {
                ret = 2*y*y - y;
            }
            else if(der_x == 1) {
                ret = 4*y/h - 1/h;
            }
            else if(der_x == 2) {
                ret = 4/(h*h);
            }
            else if(der_x >= 3) {
                ret = 0;
            }
            else {
                cout << " Wrong der_x!" <<endl;
            }
        }
        else if(basis_index == 3) {
            if(der_x == 0) {
                ret = -4*y*y +4*y;
            }
            else if(der_x == 1) {
                ret = -8*y/h + 4/h;
            }
            else if(der_x == 2) {
                ret = -8/(h*h);
            }
            else if(der_x >= 3) {
                ret = 0;
            }
            else {
                cout << " Wrong der_x!" <<endl;
            }
        }
        else {
            cout << "Wrong basis_index!" << endl;
        }
    }
    else if(basis_type==111) {
        F upper_bound = vertices(0,1);
        F lower_bound = vertices(0,0);
        F h = vertices(0,1) - vertices(0,0);
        if(basis_index == 1) {
            if(der_x == 0) {
                ret = 1;
            }
            else if(der_x >= 1) {
                ret = 0;
            }
            else {
				cout << " Wrong der_x!" << endl;
            }
        }
        else if(basis_index == 2) {
            if(der_x == 0) {
                ret = (x - lower_bound)/h;
            }
            else if(der_x == 1) {
                ret = 1/h;
            }
            else if(der_x >= 2) {
                ret = 0;
            }
            else {
                cout << " Wrong der_x!" << endl;
            }
        }
        else {
            cout << "Wrong basis_index!" << endl;
        }
    }
    return ret;
}


template<typename F,typename I,typename Matrix,typename Vector>
F local_FE_function_1D(F x,const Matrix &vertices,I basis_type,I der_x,Vector &uh_local_vec){
	F ret = 0;
	I nlb = uh_local_vec.size();
	for(I k = 0;k != nlb;++k){
		ret += uh_local_vec(k)*basis_local_function_1d(x,vertices,basis_type,k+1,der_x); 
	}
	return ret;
}


}	// end of namespace FEMethod
}	// end of namespace PDWSC

#endif
