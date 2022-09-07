#ifndef BASEFUNCTION
#define BASEFUNCTION
namespace PDWSC {
namespace FEMethod {

template<typename F,typename I,typename Matrix>
F reference_basis_2D(F xh,F yh,I basis_type,I basis_index,I der_x,I der_y) {
	F ret=0;
    if(basis_type == 201) {
		if(der_x==0&&der_y==0){
			if(basis_index==1)
				ret = 1-xh-yh;
			else if(basis_index==2)
				ret = xh;
			else if(basis_index==3)
				ret = yh;
		}
		else if(der_x==1&&der_y==0){
			if(basis_index==1)
				ret = -1;
			else if(basis_index==2)
				ret = 1;
			else if(basis_index==3)
				ret = 0;	
		}
		else if(der_x==0&&der_y==1){
			if(basis_index==1)
				ret =-1;
			else if(basis_index==2)
				ret =0;
			else if(basis_index==3)
				ret =1;
		}
    }
    else if(basis_type == 202) {
        if(basis_index == 1) {
            if(der_x == 0&&der_y ==0)
                return 2*xh*xh + 2*yh*yh + 4*xh*yh - 3*yh - 3*xh +1;
            else if(der_x == 1&&der_y ==0)
                return 4*xh + 4*yh - 3;
            else if(der_x == 2&&der_y ==0)
                return 4;
            else if(der_x == 0&&der_y==1)
                return 4*yh + 4*xh - 3;
            else if(der_x == 1&&der_y == 1)
                return 4;
            else if(der_x == 0&&der_y == 2)
                return 4;
            else if((der_x+der_y) > 2)
                return 0;
            else {
                cout << "Wrong der!" <<endl;
                return 0;
            }
        }
        else if(basis_index == 2) {
            if(der_x == 0&&der_y ==0)
                return 2*xh*xh - xh;
            else if(der_x == 1&&der_y ==0)
                return 4*xh-1;
            else if(der_x == 2&&der_y ==0)
                return 4;
            else if(der_x == 0&&der_y==1)
                return 0;
            else if(der_x == 1&&der_y == 1)
                return 0;
            else if(der_x == 0&&der_y == 2)
                return 0;
            else if((der_x+der_y) > 2)
                return 0;
            else {
                cout << "Wrong der!" <<endl;
                return 0;
            }
        }
        else if(basis_index == 3) {
            if(der_x == 0&&der_y ==0)
                return 2*yh*yh - yh;
            else if(der_x == 1&&der_y ==0)
                return 0;
            else if(der_x == 2&&der_y ==0)
                return 0;
            else if(der_x == 0&&der_y==1)
                return 4*yh - 1;
            else if(der_x == 1&&der_y == 1)
                return 0;
            else if(der_x == 0&&der_y == 2)
                return 4;
            else if((der_x+der_y) > 2)
                return 0;
            else {
                cout << "Wrong der!" <<endl;
                return 0;
            }
        }
        else if(basis_index == 4) {
            if(der_x == 0&&der_y ==0)
                return -4*xh*xh-4*xh*yh + 4*xh;
            else if(der_x == 1&&der_y ==0)
                return -8*xh - 4*yh + 4;
            else if(der_x == 2&&der_y ==0)
                return -8;
            else if(der_x == 0&&der_y==1)
                return -4*xh;
            else if(der_x == 1&&der_y == 1)
                return -4;
            else if(der_x == 0&&der_y == 2)
                return 0;
            else if((der_x+der_y) > 2)
                return 0;
            else {
                cout << "Wrong der!" <<endl;
                return 0;
            }
        }
        else if(basis_index == 5) {
            if(der_x == 0&&der_y ==0)
                return 4*xh*yh;
            else if(der_x == 1&&der_y ==0)
                return 4*yh;
            else if(der_x == 2&&der_y ==0)
                return 0;
            else if(der_x == 0&&der_y==1)
                return 4*xh;
            else if(der_x == 1&&der_y == 1)
                return 4;
            else if(der_x == 0&&der_y == 2)
                return 0;
            else if((der_x+der_y) > 2)
                return 0;
            else {
                cout << "Wrong der!" <<endl;
                return 0;
            }
        }
        else if(basis_index == 6) {
            if(der_x == 0&&der_y ==0)
                return -4*yh*yh - 4*xh*yh + 4*yh;
            else if(der_x == 1&&der_y ==0)
                return -4*yh;
            else if(der_x == 2&&der_y ==0)
                return 0;
            else if(der_x == 0&&der_y==1)
                return -8*yh -4*xh + 4;
            else if(der_x == 1&&der_y == 1)
                return -4;
            else if(der_x == 0&&der_y == 2)
                return -8;
            else if((der_x+der_y) > 2)
                return 0;
            else {
                cout << "Wrong der!" <<endl;
                return 0;
            }
        } 
		else {
            cout << "Wrong basis_index!" << endl;
            return 0;
        }
	}
	else{
		cout << "Wrong basis_type!" << endl;
		return 0;
	}
	return ret;
}

template<typename F,typename I,typename Matrix>
F local_basis_2D(F x,F y,Matrix &vertices,I basis_type,I basis_index,I der_x,I der_y) {
    F xn1 = vertices(0,0);
    F xn2 = vertices(0,1);
    F xn3 = vertices(0,2);
    F yn1 = vertices(1,0);
    F yn2 = vertices(1,1);
    F yn3 = vertices(1,2);
	F J_11 = xn2 - xn1;
	F J_12 = xn3 - xn1;
	F J_21 = yn2 - yn1;
	F J_22 = yn3 - yn1;
    F J_det = J_11*J_22 - J_12*J_21;

    F xh = (J_22*(x - xn1) - J_12*(y - yn1))/J_det;
    F yh = (-J_21*(x - xn1) + J_11*(y - yn1))/J_det;
	
	F r = 0;
    if(basis_type == 201) {
        if((der_x==0)&&(der_y==0))
            r= reference_basis_2D<F,I,Matrix>(xh,yh,basis_type,basis_index,0,0);
        else if((der_x == 1)&&(der_y == 0))
            r = (reference_basis_2D<F,I,Matrix>(xh,yh,basis_type,basis_index,1,0)*J_22 + reference_basis_2D<F,I,Matrix>(xh,yh,basis_type,basis_index,0,1)*(-J_21))/J_det;
        else if((der_x == 0)&&(der_y ==1))
            r = (reference_basis_2D<F,I,Matrix>(xh,yh,basis_type,basis_index,1,0)*(-J_12)+reference_basis_2D<F,I,Matrix>(xh,yh,basis_type,basis_index,0,1)*J_11)/J_det;		
		return r;
        }
    else if(basis_type == 202) {
        if((der_x==0)&&(der_y==0))
            return reference_basis_2D<F,I,Matrix>(xh,yh,basis_type,basis_index,0,0);
        else if((der_x == 1)&&(der_y == 0))
            return (reference_basis_2D<F,I,Matrix>(xh,yh,basis_type,basis_index,1,0)*J_22 + reference_basis_2D<F,I,Matrix>(xh,yh,basis_type,basis_index,0,1)*(-J_21))/J_det;
        else if((der_x == 0)&&(der_y ==1))
            return (reference_basis_2D<F,I,Matrix>(xh,yh,basis_type,basis_index,1,0)*(-J_12)+reference_basis_2D<F,I,Matrix>(xh,yh,basis_type,basis_index,0,1)*J_11)/J_det;	
        else if((der_x ==2)&&(der_y == 0)) {
            J_det = J_det*J_det;
            return reference_basis_2D<F,I,Matrix>(xh,yh,basis_type,basis_index,2,0)*((yn3 - yn1)*(yn3 - yn1)/J_det) +2* reference_basis_2D<F,I,Matrix>(xh,yh,basis_type,basis_index,1,1)*((yn3 - yn1)*(yn1 - yn2)/J_det) + reference_basis_2D<F,I,Matrix>(xh,yh,basis_type,basis_index,0,2)*((yn1 - yn2)*(yn1 - yn2)/J_det);
        }
        else if((der_x ==0)&&(der_y == 2)) {
            J_det = J_det*J_det;
            return reference_basis_2D<F,I,Matrix>(xh,yh,basis_type,basis_index,2,0)*((xn1 - xn3)*(xn1 - xn3)/J_det) +2* reference_basis_2D<F,I,Matrix>(xh,yh,basis_type,basis_index,1,1)*((xn1 - xn3)*(xn2 - xn1)/J_det) + reference_basis_2D<F,I,Matrix>(xh,yh,basis_type,basis_index,0,2)*((xn2 - xn1)*(xn2 - xn1)/J_det);
        }
        else if((der_x ==1)&&(der_y == 1)) {
            J_det = J_det*J_det;
            return reference_basis_2D<F,I,Matrix>(xh,yh,basis_type,basis_index,2,0)*((xn1 - xn3)*(yn3 - yn1)/J_det) + reference_basis_2D<F,I,Matrix>(xh,yh,basis_type,basis_index,1,1)*((xn1 - xn3)*(yn1 - yn2)/J_det) +  reference_basis_2D<F,I,Matrix>(xh,yh,basis_type,basis_index,1,1)*((xn2 - xn1)*(yn3 - yn1)/J_det) +reference_basis_2D<F,I,Matrix>(xh,yh,basis_type,basis_index,0,2)*((xn2 - xn1)*(yn1 - yn2)/J_det);
        }
        else if(der_x + der_y > 2)
            return 0;
        else {
            cout << "Wrong der!" << endl;
            return 0;
        }
    }
    else {
        cout << "Wrong basis_type!" <<endl;
        return 0;

    }
}
template<typename F,typename I,typename Matrix,typename Vector>
F local_FE_function_2D(F x,F y,const Matrix &vertices,I basis_type,I der_x,I der_y,Vector &uh_local_vec) {
    F ret = 0;
    I nlb = uh_local_vec.size();
    for(I k = 0; k != nlb; ++k) {
        ret += uh_local_vec(k)*local_basis_2D(x,y,vertices,basis_type,k+1,der_x,der_y) ;
    }
    return ret;
}



}	// end of namespace FEMethod
}	// end of namespace PDWSC

#endif
