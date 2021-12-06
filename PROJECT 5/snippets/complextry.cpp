#include<complex>
#include<armadillo>
const arma::cx_double i(0.0,1.0);
arma::cx_double testfunction (arma::cx_double a_in){
  return std::exp(i*a_in);
}
int main () {

   arma::cx_double a = arma::cx_double(2.0, 3.0);

   arma::cx_double b = 2.0+3.14*1.0j;
   arma::cx_double c = std::exp(b);
   std::cout << "exp(2+i*pi)= " << std::exp(b) << std::endl;
   std::cout << "exp(2+3i)=" << std::exp(a)<< std::endl;
   std::cout << "exp(i*pi)=" << testfunction(3.14) << std::endl;
}
