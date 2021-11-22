# CODES AND PLOTS FOR PROJECT 4<br>
This repository hosts files related to the simulation of a 2D spin lattice using the Ising Model and<br>
a Markov chain Monte Carlo algorithm.<br><br>

BUILD COMMAND IN UNIX TERMINAL: g++ main.cpp src/Lattice.cpp -I include/ -larmadillo -fopenmp -O3 -o main <br>

The code is parallelized using OpenMP, in order to run a series of simulations for different values of the temperature.<br>
This is properly set up during the program execution. To adjust the temperature interval and the number of tested temperatures,<br>
edit the variables  T_start, T_end and N_T in the main.cpp file<br>
The program takes then two command line arguments: run it as<br>

          main <Ncycles> <L> <br>
 
where <Ncycles> is the number of Monte Carlo cycles you want, and <L> is the lattice size that you want.<br><br>
We also provide a little bash script, called run_MC, which runs (and times) the code for different sizes of the grid, using 1MLN MC cycles.<br>







