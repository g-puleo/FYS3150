# Codes and plots for PROBLEM 9<br>

compile the code with the following command on UNIX terminal:<br>

g++ main.cpp Particle.cpp PenningTrap.cpp -larmadillo -o prob9.exe<br>

execute the program with the following command line arguments:<br>

./prob9.exe algo N c_int <br>

where <br>

	algo is a string, and it must be either "RK4" or "FE" depending on the algorithm you want to choose.<br>
	N is an integer and it is the number of particles you want in the Penning trap<br>
	c_int is a string, and must be either "on" or "off", depending on whether you want coulomb interactions or not.<br>

The program output are some .csv files organized as follows:<br><br>

	"PenningData_algo_X_N.csv": contains data related to the trajectory of particle N. The label X is a letter between "ABCDE", <br>
	each of the labels corresponds to a certain value of h (the first you set will have the label A and so on)<br>
	output has 7 columns: elapsed_time; x ; y; z; v_x; v_y; v_z ;<br><br>
	"r_ana_X.csv": this output is turned on only if you set the flag analytical_prediction to be true in main.cpp,<br>
	and run the code with the coulomb interactions turned off. Contains analytical predictions for each of the particles. <br>
	output has 6 columns: elapsed_time; x; y; z; relative_error; absolute_error <br><br><br>
	
In order to make plots we have 3 python scripts: <br><br>
	
 	comparison.py<br>
 	allows to make a comparison between the analytical solution and the numerical one<br>
 	you should organize your outputs (for 5 values of h) into 2 folders:<br>
 	rel_err_FE and rel_err_RK4<br>
 	each must contain both numerical and analytical solution .csv files, output by ./prob9.exe<br>
 	paste the python script in a folder containing those two folders, and run it.<br><br>
 	
 	view_data.csv<br>
 	allows to compare results for 2 particles, when coulomb interactions are switched on/off<br>
 	paste your script into a folder, then organize your output as follows:<br>
 	paste the data with c-int turned off into yourfolder/N/experiment_1<br>
 	paste the data with c-int turned on into yourfolder/N/experiment_2<br>
	N must be an integer, allows to organize multiple simulation data.<br><br>
 	
 	simple_zplot.py<br>
 	a very short script, run it in the same folder where you have your (one) output file. Plots<br>
 	only the motion along the z-axis<br><br><br>
 	
 	
 	
 	
 		 	
