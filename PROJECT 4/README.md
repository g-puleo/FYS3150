# CODES AND PLOTS FOR PROJECT 4
This repository hosts files related to the simulation of a 2D spin lattice using the Ising Model and<br>
a Markov chain Monte Carlo algorithm.<br><br>

BUILD COMMAND IN UNIX TERMINAL: 

            g++ main.cpp src/Lattice.cpp -I include/ -larmadillo -fopenmp -O3 -o main

The code is parallelized using OpenMP, in order to run a series of simulations for different values of the temperature.<br>
This is properly set up during the program execution. To adjust the temperature interval and the number of tested temperatures,<br>
edit the variables  T_start, T_end and N_T in the main.cpp file<br>
The program takes then two command line arguments: run it as<br>

            main <Ncycles> <L>
            
where <Ncycles> is the number of Monte Carlo cycles you want, and <L> is the lattice size that you want.<br>
As it is, the outputs one file for each temperature, storing the final estimates of the thermodynamical variables which we studied.<br>
These files are ready to be merged with one simple command in the terminal, such as:<br>
            
            cat LatticeData_40_*.csv > tempdata40.csv

(this is just an example for the files output when L=40)<br>
If you want to see how numerical estimates evolve during the simulation, see the comments at the end of the main.cpp file, and edit it properly.

We also provide a little bash script, called run_MC, which runs (and times) the code for different sizes of the grid, using 1MLN MC cycles.<br>
<br>
Python script run: python3 temp_analysis.py <br>
Reads five files: tempdata40.csv, tempdata50.csv, tempdata60.csv, tempdata80.csv, tempdata100.csv <br>
Performs ordinary least squares regression to find optimal polynomial fits for the data. <br>
The script yields six plots. Four of the plots are of average energy per spin, average magnetization per spin, <br>
specific heat, and magnetic susceptibility, all against temperature. The last two are a fitted line against the <br>
maxima of the specific heat and magnetic susceptibility against the inverse of the lattice size. 




