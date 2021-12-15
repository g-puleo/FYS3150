# CODE FOR THE SOLUTION OF THE 2D SCHROEDINGER EQATION IN A SQUARED BOX

build command in UNIX terminal:<br>

		g++ main.cpp src/Box.cpp src/utils.cpp -I include/ -larmadillo  -o main

run:<br>

		main initial_conditions.csv PotentialFiles/filename.dat
		
initial_conditions.csv is a .csv file containing the parameters of the problem.<br>
filename.dat is the name of the file storing the potential matrix.<br>
PotentialFiles is a folder containing various potential matrices in .dat format and the python script we used to create those files. The file names are self explaining. <br><br>
Once you have run the program, you can visualize results by just running:<br>

		python3 sol_animation.py

This python script saves the animation of the simulation to a .gif file. <br>
Some examples are available in the Animations folder in this repository.<br> 
		

