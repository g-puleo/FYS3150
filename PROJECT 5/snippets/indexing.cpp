#include<armadillo>
#include<iostream>


int main () {

  arma::mat A = arma::mat(3,3);

  for (int jj = 0; jj<3; jj++){
    for (int ii = 0; ii<3 ; ii++){
      A(jj,ii) = ii+3*jj;
    }
  }
  std::cout << A << "\n" <<  std::endl;

  std::cout << A(4) << std::endl;
}
