#ifndef GENERAGTE_PT_H
#define GENERAGTE_PT_H

namespace PDWSC {
namespace FEMethod {
template<typename F,typename I,typename Matrix>
void generate_PT_1d(F a,F b,I N,Matrix &P,Matrix &T) {
    F h = (b - a)/N;
    for(I i = 0; i < N+1; ++i) {
        P[0][i] = a + h*i;
    }
    for(I i = 0; i < 2 ; ++i) {
        for(I j = 0; j < N; ++j) {
            T[i][j] = i + j + 1;
        }
    }
}
template<typename F,typename I,typename Matrix>
void generate_PbTb_1d(F a,F b,I N,I basis_type,I &Nb,I &Nlb,Matrix &Pb,Matrix &Tb) {
    if(basis_type == 101) {
        F h = (b - a)/N;
        for(I i = 0; i < N+1; ++i) {
            Pb[0][i] = a + h*i;
        }
        for(I i = 0; i < 2 ; ++i) {
            for(I j = 0; j < N; ++j) {
                Tb[i][j] = i + j + 1;
            }
        }
		Nb = N + 1;
		Nlb = 2;
    } else if(basis_type == 102) {
        F h = (b - a)/(2*N);
        for(I i = 0; i < 2*N+1; ++i) {
            Pb[0][i] = a + i*h;
        }
        for(I i = 0; i < 3; ++i) {
            for(I j = 0; j < N; ++j) {
                if(i == 0) {
                    Tb[i][j] = 2*j+1;
                }
                else if(i == 1) {
                    Tb[i][j] = 2*j+3;
                }
                else {
                    Tb[i][j] = 2*j+2;
                }
            }
        }
		Nb = 2*N + 1;
	    Nlb = 3; 	
    }
}


template<typename F,typename I,typename Matrix>
void generate_boundraynodes(Matrix &boundarynodes,Matrix &Tb,int boundary_type = 101)
{
    
}
}  // end of namespace FEMethod
}  // end of namespace PDWSC
#endif
