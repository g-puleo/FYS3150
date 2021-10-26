# Codes and plots for PROJECT 3<br>


build command for simulation WITHOUT oscillating electric field:<br>
g++ main9.cpp src/Particle.cpp src/PenningTrap9.cpp -I include/ -larmadillo -o prob9.exe<br>
to compile the code for simulation WITH oscillating electric field:<br>
g++ main10.cpp src/Particle.cpp src/PenningTrap10.cpp -I include/ -larmadillo -o prob10.exe<br>

execute the program with the following command line arguments:<br>

./prob9.exe algo N c_int <br>
./prob10.exe algo N c_int osc <br>

where <br>

	algo is a string, and it must be either "RK4" or "FE" depending on the algorithm you want to choose.
	N is an integer and it is the number of particles you want in the Penning trap
	c_int is a string, and must be either "on" or "off", depending on whether you want coulomb interactions or not.
	osc is a string , and must be either "on" or "off", depending on whether you want the oscillating E field or not.
	
The program output are some .csv files organized as follows:<br><br>

	"PenningData_algo_X_N.csv": contains data related to the trajectory of particle N. The label X is a letter between "ABCDE", 
	each of the labels corresponds to a certain value of h (the first you set will have the label A and so on)
	output has 7 columns: elapsed_time; x ; y; z; v_x; v_y; v_z ;
	"r_ana_N.csv": this output is turned on only if you set the flag analytical_prediction to be true in main.cpp,
	and run the code with the coulomb interactions turned off. Contains analytical predictions for each of the particles. 
	output has 6 columns: elapsed_time; x; y; z; relative_error; absolute_error <
	"PD_algo_AMP_amplitudeindex.csv": contains the number of particles remaining in trap after simulation run for each 
	frequency. Headers: Frequency; N;, where the Frequency column holds frequencies of the applied perturbation in the 
	electric field, and the N column holds the number of particles remaining in trap. 
	
In order to make plots we have 4 python scripts: <br><br>
	
 	comparison.py
 	allows to make a comparison between the analytical solution and the numerical one
 	you should organize your outputs (for 5 values of h) into 2 folders:
 	rel_err_FE and rel_err_RK4
 	each must contain both numerical and analytical solution .csv files, output by ./prob9.exe
 	paste the python script in a folder containing those two folders, and run it.
 	
 	view_data.csv
 	allows to compare results for 2 particles, when coulomb interactions are switched on/off
 	paste your script into a folder, then organize your output as follows:
 	paste the data with c-int turned off into yourfolder/N/experiment_1
 	paste the data with c-int turned on into yourfolder/N/experiment_2
	N must be an integer, allows to organize multiple simulation data.
 	
 	simple_zplot.py
 	a very short script, run it in the same folder where you have your (one) output file. Plots
 	only the motion along the z-axis.
 	
 	plot_particles_left.py
 	this script needs to be in the same folder as output files. Script where you need to 
 	modify the arguments for the filenames to suit your needs. Plots the fraction of a 
 	simulation run with 100 particles against frequency. Can also plot a chosen number of
 	trajectories of these simulation runs. 
 	
 	
 	
 	
 		 	
