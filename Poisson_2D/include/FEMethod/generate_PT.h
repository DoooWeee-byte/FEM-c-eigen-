#ifndef GENERAGTE_PT_H
#define GENERAGTE_PT_H

namespace PDWSC {
namespace FEMethod {
template<typename F,typename I,typename Matrix>
void generate_PT_2D(Matrix &P,Matrix &T,I N1,I N2,F right,F top,F bottom = 0,F left = 0) {
    F h1 = (right - left)/N1;
    F h2 = (top - bottom)/N2;
    for(I rn = 1; rn != N2+2; ++rn) {
        for(I cn = 1; cn != N1+2; ++cn) {
            F x = left + (cn - 1)*h1;
            F y = bottom + (rn - 1)*h2;
            I n = (cn - 1)*(N2+1) + rn;
            P[0][n-1] = x;
            P[1][n-1] = y;
        }
    }
    for(I re = 1; re != N2+1; ++re) {
        for(I ce = 1; ce != N1+1; ++ce) {
            I n = (ce - 1)*N2 + re;
            I n1 = 2*(n-1)+1;
            I n2 = n1 + 1;
            I rn1 = re;
            I cn1 = ce;
            I j1 = (cn1 -1)*(N2+1) + rn1;
            I rn2 = rn1;
            I cn2 = cn1 + 1;
            I j2 = (cn2 -1)*(N2+1) + rn2;
            I rn3 = rn1 + 1;
            I cn3 = cn1;
            I j3 = (cn3 -1)*(N2+1) + rn3;
            I rn4 = rn1 + 1;
            I cn4 = cn1 + 1;
            I j4 = (cn4 -1)*(N2+1) + rn4;
            T[0][n1-1] = j1;
            T[1][n1-1] = j2;
            T[2][n1-1] = j3;
            T[0][n2-1] = j3;
            T[1][n2-1] = j2;
            T[2][n2-1] = j4;
        }
    }
}
template<typename F,typename I,typename Matrix>
void generate_PbTb_2D(Matrix &Pb,Matrix &Tb,I N1,I N2,I basis_type,F right,F top,F bottom = 0,F left = 0) {
    if(basis_type == 201) {
        F h1 = (right - left)/N1;
        F h2 = (top - bottom)/N2;
        for(I rn = 1; rn != N2+2; ++rn) {
            for(I cn = 1; cn != N1+2; ++cn) {
                F x = left + (cn - 1)*h1;
                F y = bottom + (rn - 1)*h2;
                I n = (cn - 1)*(N2+1) + rn;
                Pb[0][n-1] = x;
                Pb[1][n-1] = y;
            }
        }
        for(I re = 1; re != N2+1; ++re) {
            for(I ce = 1; ce != N1+1; ++ce) {
                I n = (ce - 1)*N2 + re;
                I n1 = 2*(n-1)+1;
                I n2 = n1 + 1;
                I rn1 = re;
                I cn1 = ce;
                I j1 = (cn1 -1)*(N2+1) + rn1;
                I rn2 = rn1;
                I cn2 = cn1 + 1;
                I j2 = (cn2 -1)*(N2+1) + rn2;
                I rn3 = rn1 + 1;
                I cn3 = cn1;
                I j3 = (cn3 -1)*(N2+1) + rn3;
                I rn4 = rn1 + 1;
                I cn4 = cn1 + 1;
                I j4 = (cn4 -1)*(N2+1) + rn4;
                Tb[0][n1-1] = j1;
                Tb[1][n1-1] = j2;
                Tb[2][n1-1] = j3;
                Tb[0][n2-1] = j3;
                Tb[1][n2-1] = j2;
                Tb[2][n2-1] = j4;
            }
        }

    } else if(basis_type == 202) {
        F h1 = (right - left)/(2*N1);
        F h2 = (top - bottom)/(2*N2);
        for(I rn = 1; rn != 2*N2 + 2; ++rn) {
            for(I cn = 1; cn != 2*N1 +2; ++cn) {
                F x = left + (cn - 1)*h1;
                F y = bottom + (rn - 1)*h2;
                I n = (cn - 1)*(2*N2+1) + rn;
                Pb[0][n-1] = x;
                Pb[1][n-1] = y;
            }
        }
        for(I re = 1; re != N2+1; ++re) {
            for(I ce = 1; ce != N1+1; ++ce) {
                I n = (ce - 1)*N2 + re;
                I n1 = 2*(n-1)+1;
                I n2 = n1 + 1;
                I j = (re - 1)*2 + 1;
                I i = (ce - 1)*2 + 1;
                Tb[0][n1-1] = (i - 1)*(2*N2 +1) + j;
                Tb[1][n1-1] = (i + 1)*(2*N2 +1) + j;
                Tb[2][n1-1] = (i - 1)*(2*N2 +1) + j + 2;
                Tb[3][n1-1] = i*(2*N2 +1) + j;
                Tb[4][n1-1] = i*(2*N2 +1) + j + 1;
                Tb[5][n1-1] = (i - 1)*(2*N2 +1) + j + 1;
                Tb[0][n2-1] = (i - 1)*(2*N2 +1) + j + 2;
                Tb[1][n2-1] = (i + 1)*(2*N2 +1) + j;
                Tb[2][n2-1] = (i + 1)*(2*N2 +1) + j + 2;
                Tb[3][n2-1] = i*(2*N2 +1) + j + 1;
                Tb[4][n2-1] = (i + 1)*(2*N2 +1) + j + 1;
                Tb[5][n2-1] = i*(2*N2 +1) + j + 2;
            }
        }
    }
}

//全diriclet条件
template<typename F,typename I,typename Matrix>
void generate_boundaryedges_2D(Matrix &boundaryedges,I N1,I N2) {
    for(I j=1; j != N1 + 1; ++j) {
        boundaryedges[0][j-1] = -1;
        boundaryedges[1][j-1] = (j - 1)*(2*N2) + 1;
        boundaryedges[2][j-1] = (j - 1)*(N2 + 1) + 1;
        boundaryedges[3][j-1] = j*(N2 + 1) + 1;
    }
    for(I j=N1+1; j != N1+N2+1; ++j) {
        boundaryedges[0][j-1] = -1;
        boundaryedges[1][j-1] = (2*N2)*(N1-1)+(j-N1)*2;
        boundaryedges[2][j-1] = N1*(N2+1)+j-N1;
        boundaryedges[3][j-1] = N1*(N2+1)+j-N1+1;
    }
    for(I j=N1+N2+1; j != 2*N1+N2+1; ++j) {
		I k = j - N1 - N2;
        boundaryedges[0][j-1] = -1;
        boundaryedges[1][j-1] = (N1 - k+1)*(2*N2);
        boundaryedges[2][j-1] = (N1 - k+2)*(1+N2);
        boundaryedges[3][j-1] = (N1 - k+1)*(1+N2);
    }
	for(I j=2*N1+N2+1;j !=2*N1+2*N2+1;++j){
		I k = j - 2*N1-N2;
		boundaryedges[0][j-1] = -1;
        boundaryedges[1][j-1] = 2*N2 - (k-1)*2 - 1;
        boundaryedges[2][j-1] = N2-k+2;
        boundaryedges[3][j-1] = N2-k+1;
	}
}

//第一条边为Neumann条件，其余为Dirichlet条件
template<typename F,typename I,typename Matrix>
void generate_boundaryedges_DN_2D(Matrix &boundaryedges,I N1,I N2) {
    for(I j=1; j != N1 + 1; ++j) {
        boundaryedges[0][j-1] = -2;
        boundaryedges[1][j-1] = (j - 1)*(2*N2) + 1;
        boundaryedges[2][j-1] = (j - 1)*(N2 + 1) + 1;
        boundaryedges[3][j-1] = j*(N2 + 1) + 1;
    }
    for(I j=N1+1; j != N1+N2+1; ++j) {
        boundaryedges[0][j-1] = -1;
        boundaryedges[1][j-1] = (2*N2)*(N1-1)+(j-N1)*2;
        boundaryedges[2][j-1] = N1*(N2+1)+j-N1;
        boundaryedges[3][j-1] = N1*(N2+1)+j-N1+1;
    }
    for(I j=N1+N2+1; j != 2*N1+N2+1; ++j) {
		I k = j - N1 - N2;
        boundaryedges[0][j-1] = -1;
        boundaryedges[1][j-1] = (N1 - k+1)*(2*N2);
        boundaryedges[2][j-1] = (N1 - k+2)*(1+N2);
        boundaryedges[3][j-1] = (N1 - k+1)*(1+N2);
    }
	for(I j=2*N1+N2+1;j !=2*N1+2*N2+1;++j){
		I k = j - 2*N1-N2;
		boundaryedges[0][j-1] = -1;
        boundaryedges[1][j-1] = 2*N2 - (k-1)*2 - 1;
        boundaryedges[2][j-1] = N2-k+2;
        boundaryedges[3][j-1] = N2-k+1;
	}
}

//第一条边为Robin条件，其余为Dirichlet条件
template<typename F,typename I,typename Matrix>
void generate_boundaryedges_DR_2D(Matrix &boundaryedges,I N1,I N2) {
    for(I j=1; j != N1 + 1; ++j) {
        boundaryedges[0][j-1] = -3;
        boundaryedges[1][j-1] = (j - 1)*(2*N2) + 1;
        boundaryedges[2][j-1] = (j - 1)*(N2 + 1) + 1;
        boundaryedges[3][j-1] = j*(N2 + 1) + 1;
    }
    for(I j=N1+1; j != N1+N2+1; ++j) {
        boundaryedges[0][j-1] = -1;
        boundaryedges[1][j-1] = (2*N2)*(N1-1)+(j-N1)*2;
        boundaryedges[2][j-1] = N1*(N2+1)+j-N1;
        boundaryedges[3][j-1] = N1*(N2+1)+j-N1+1;
    }
    for(I j=N1+N2+1; j != 2*N1+N2+1; ++j) {
		I k = j - N1 - N2;
        boundaryedges[0][j-1] = -1;
        boundaryedges[1][j-1] = (N1 - k+1)*(2*N2);
        boundaryedges[2][j-1] = (N1 - k+2)*(1+N2);
        boundaryedges[3][j-1] = (N1 - k+1)*(1+N2);
    }
	for(I j=2*N1+N2+1;j !=2*N1+2*N2+1;++j){
		I k = j - 2*N1-N2;
		boundaryedges[0][j-1] = -1;
        boundaryedges[1][j-1] = 2*N2 - (k-1)*2 - 1;
        boundaryedges[2][j-1] = N2-k+2;
        boundaryedges[3][j-1] = N2-k+1;
	}
}

template<typename F,typename I,typename Matrix>
void generate_boundarynodes_2D(Matrix &boundarynodes,I N1,I N2)
{
	for(I i = 0;i < boundarynodes.shape[1];++i){
		boundarynodes[0][i] = -1;
	}
	for(I i = 1;i < 2*(N1+N2)+1;++i){
		if(i < N1 + 1){
			boundarynodes[1][i-1] = (i-1)*(N2 + 1) + 1;
		}
		if(i> N1&&(i< N1+N2+1)){
			boundarynodes[1][i-1] = N1*(N2+1) + i - N1;
		}
		if((i > N1+N2)&&(i < 2*N1+N2+1)){
			I k = N1 + 2 - (i - N1 - N2);
			boundarynodes[1][i-1] = k*(N2+1);
		}
		if(i > 2*N1 + N2){
			I k = N2 + 2 - (i - 2*N1 - N2);
			boundarynodes[1][i-1] = k;
		}
	}	
}

template<typename F,typename I,typename Matrix>
void generate_boundarynodes_DN_2D(Matrix &boundarynodes,I N1,I N2)
{
	for(I i = 0;i < boundarynodes.shape[1];++i){
		boundarynodes[0][i] = -1;
	}
	for(I i = 1;i < 2*(N1+N2)+1;++i){
		if(i < N1 + 1){
			boundarynodes[0][i-1] = -2;
			boundarynodes[1][i-1] = (i-1)*(N2 + 1) + 1;
		}
		boundarynodes[0][0] = -1;
		if(i> N1&&(i< N1+N2+1)){
			boundarynodes[1][i-1] = N1*(N2+1) + i - N1;
		}
		if((i > N1+N2)&&(i < 2*N1+N2+1)){
			I k = N1 + 2 - (i - N1 - N2);
			boundarynodes[1][i-1] = k*(N2+1);
		}
		if(i > 2*N1 + N2){
			I k = N2 + 2 - (i - 2*N1 - N2);
			boundarynodes[1][i-1] = k;
		}
	}	
}
template<typename F,typename I,typename Matrix>
void generate_boundarynodes_DR_2D(Matrix &boundarynodes,I N1,I N2)
{
	for(I i = 0;i < boundarynodes.shape[1];++i){
		boundarynodes[0][i] = -1;
	}
	for(I i = 1;i < 2*(N1+N2)+1;++i){
		if(i < N1 + 1){
			boundarynodes[0][i-1] = -3;
			boundarynodes[1][i-1] = (i-1)*(N2 + 1) + 1;
		}
		boundarynodes[0][0] = -1;
		if(i> N1&&(i< N1+N2+1)){
			boundarynodes[1][i-1] = N1*(N2+1) + i - N1;
		}
		if((i > N1+N2)&&(i < 2*N1+N2+1)){
			I k = N1 + 2 - (i - N1 - N2);
			boundarynodes[1][i-1] = k*(N2+1);
		}
		if(i > 2*N1 + N2){
			I k = N2 + 2 - (i - 2*N1 - N2);
			boundarynodes[1][i-1] = k;
		}
	}	
}
}  // end of namespace FEMethod
}  // end of namespace PDWSC
#endif
