#include<iostream>
#include<armadillo>
#include<cstdlib>
using namespace arma
int main ()
{
 int M = 201;
	sp_cx_mat A = sp_cx_mat((M-2)*(M-2),(M-2)*(M-2))
	sp_cx_mat B = sp_cx_mat((M-2)*(M-2),(M-2)*(M-2))
	cx_mat C = cx_mat(M-2,M-2, fill::zeros)
	cx_mat D = cx_mat(M-2,M-2, fill::zeros)
	cx_mat E = cx_mat(M-2,M-2, fill::zeros)
	cx_mat F = cx_mat(M-2,M-2, fill::zeros)
  cx_mat G = cx_mat(M-2,M-2, fill::zeros)

	srand(time(NULL));

	for(int jj=0; jj<M; jj++){
		for(int ii=0; ii<M; ii++ ){
			C(ii,jj) = (double) 1.0/rand()%10;
			D(ii,jj) = (double) 1.0/rand()%10;
			E(ii,jj) = (double) 1.0/rand()%10;
			F(ii,jj) = (double) 1.0/rand()%10;
			G(ii,jj) = (double) 1.0/rand()%10;
		}
	}


  for(int jj=0; jj<(M-2)*(M-2); jj++){
		A(jj, jj) = jj;
		A(jj, jj+1) = 1;
		B(jj+1, jj) = 1;
		B(jj, jj) 
	}
	}
