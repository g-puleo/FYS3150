#include"Lattice.hpp"
//initialize a global object "distribution" (want to sample random numbers in the interval [0,1])
std::mt19937 generator; //uses the Mersenne twister RNG
std::uniform_real_distribution<double>unif_01(0.0,1.0);

//constructor : L_in is the length of one row (or one column)
// order is a bool: True if you want to start from the ground state.
//                  False if you want to initialize each spin randomly.
//seed_ is the seed for the random number generator.
Lattice::Lattice( int L_in , bool order, int seed_)
{

  generator.seed(seed_);
  //store the length of one row:
  L_ = L_in;
  //initialize the spin grid to be a matrix L_in+1 x L_in+1
  spin_grid_ = arma::Mat<short int>(L_in+1, L_in+1, arma::fill::zeros);
  int r;//this is going to be a random int
  if (!order){
  for(int i=0; i<L_; i++)
  {
    for(int j=0; j<L_; j++)
    {
      r = generator()%2; //r  is randomly chosen between 0 and 1
      spin_grid_(i,j)=2*r-1; //transform 0 --> -1 and 1 --> 1
    }
   }
  }
  else
  {
    for(int i=0; i<L_; i++)
    {
      for(int j=0; j<L_; j++)
      {
        r = 1;//(i+j)%2; //initialize the grid with the most improbable situation (max energy)
        //picture this as a chessboard, where white squares are +1 and black squares are -1
        //each neighbouring pair will have opposite sign i.e. a positive contribution to the
        //total energy.
        spin_grid_(i,j)=2*r-1; //this serves to convert (0,1) to (-1,1)
      }
     }
  }

      //copies the first row (and column) to the last row (and column)
      spin_grid_(arma::span(L_in,L_in), arma::span(0,L_in-1)) = spin_grid_(arma::span(0,0), arma::span(0,L_in-1) );
      spin_grid_( arma::span(0,L_in-1), arma::span(L_in,L_in)) = spin_grid_(arma::span(0,L_in-1), arma::span(0,0) );

      //store energy and magnetization
      E_ = tot_energy();
      M_ = tot_magnetization();
}


//this function calculates the total energy:
double Lattice::tot_energy( void )
{
  double Etot = 0; //initialize total energy
  for(int i=0 ; i<L_; i++)
    {
      for(int j=0; j<L_; j++)
      {
        Etot+= spin_grid_(i,j)*(spin_grid_(i,j+1)+spin_grid_(i+1,j));
      }

    }
    return -Etot; //note the minus sign here
}

//this function calculates total magnetization
double Lattice::tot_magnetization( void )
{
  double Mtot = 0;
  for(int i=0; i<L_; i++)
    {
      for(int j=0; j<L_; j++)
      {
        Mtot+= spin_grid_(i,j);
      }
    }
    return Mtot;
}

///THIS IS THE IMPLEMENTATION OF THE METROPOLIS ALGORITHM TO GENERATE A CHAIN OF SAMPLES.
void Lattice::evolve_Metropolis( const double* delta )
{
  int N = L_*L_;
   for( int flips_in_a_cycle = 0; flips_in_a_cycle<N; flips_in_a_cycle++)
   {
        //start generating two random numbers (aka a position in the grid)

        int k = generator()%L_;
        int j = generator()%L_;

        //then try to flip the (k,j) spin
        spin_grid_(k,j)*=-1;
        //find the contribution to the total energy due to the final configuration around the flipped spin
        int final_E_count = 0;

          for(int h1 = -1 ; h1<2; h1+=2)
          {
            int a=k+h1; //if the random generated position in the grid belongs to the first row (or column)
            int b=j+h1; //then the upper (or left hand) neighbour is found in the last row (or column)
              if(a<0)   //therefore, the two if's serve to avoid "index-out-of-bound" error
                a=L_-1;

              if(b<0)
                b=L_-1;
            //use -= because energy is defined with a minus sign in front of the sum
              final_E_count-= spin_grid_(k,b)+spin_grid_(a, j);
            }

        final_E_count*=spin_grid_(k,j);
        //final_E_count has now 5 different possible values {-4 -2 0 2 4}
        //each of them corresponds to energy changes DeltaE = {-8 -4 0 4 8}
        //the exp(-beta*DeltaE) values are contained in the double array "delta" , which is defined
        //in main() and is passed to the evolve_MCMC function as argument. This is done to avoid
        //multiple calls to the exp(function).

        //note that proper array indexing is achieved by means of the funciton f(x)=2+ x/2, such that:
        //x = -4 ----> f(x) = 0
        //x = -2 ----> f(x) = 1
        //x = 0  ----> f(x) = 2
        //x = 2  ----> f(x) = 3
        //x = 4  ----> f(x) = 4
        double A = delta[ 2+ final_E_count/2]; //ratio p(s-new)/p(s-old)
        double r = unif_01(generator);  //random number between 0 and 1
        // std::cout << "random pick in (0,1): " << r << std::endl;
        // std::cout << "X = 2+ final_E_count/2 = "  << 2+ final_E_count/2 << "\n" << std::endl;
        // std::cout << "A = p_ratio[X] = " << A << std::endl;
        if (r>A) //if r>A, rejects flip,
        {
          spin_grid_(k,j)*=-1;
        }
        else  //else , properly update the lattice to ensure periodic boundary condition
        //note: the very last row and column in the lattice are "fictitious": they are not really
        //part of the lattice, but they only serve to calculate the energy in an easy way when considering
        //the periodic boundary condition.
        {
          if(k==0)
          { //when flipping a spin in the 1st row, update very last  row to ensure periodic boundary condition
            spin_grid_(L_,j)*=-1;
          }
          if(j==0)
          { //when flipping a spin in the 1st column, update the very last column to ensure periodic bc
            spin_grid_(k, L_)*=-1;
          }
          //in case of flipped spin: update energy and magnetization
          E_ += final_E_count*2;
          M_ += 2*spin_grid_(k,j); //if final spin is -1 (was +1), M decreased by 2
                                  //if final spin is +1 (was -1) , M increased by 2
        } //end of else

   } //end of for loop: each iteration is a Monte Carlo cycle
}
