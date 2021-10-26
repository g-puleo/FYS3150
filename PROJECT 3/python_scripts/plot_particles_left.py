"""
Created Oct 22 2021

@author: jebetten
"""
import os
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

PROJECT_ROOT_DIR = "Results"
FIGURE_ID = "Results/Plots"

if not os.path.exists(PROJECT_ROOT_DIR):
    os.mkdir(PROJECT_ROOT_DIR)

if not os.path.exists(FIGURE_ID):
    os.makedirs(FIGURE_ID)

def image_path(fig_id):
    return os.path.join(FIGURE_ID, fig_id)

def save_fig(fig_id):
    plt.savefig(image_path(fig_id) + ".pdf", format='pdf')

MEDIUM = 18
#algo = "FE_"
algo = "RK4_"
"""
MAKING FRACTIONAL PLOTS AGAINST FREQUENCY
This script has been reused again and again, with different datapoints to produce 
plots for project 3, (problem 10)
"""
# Reading data
f1 = "PD_RK4_AMP_0.csv" # modify here to read proper file
data1 = pd.read_csv(f1, delimiter=";")
f2 = "PD_RK4_AMP_1.csv" # modify here to read proper file
data2 = pd.read_csv(f2, delimiter=";")
f3 = "PD_RK4_AMP_2.csv" # modify here to read proper file
data3 = pd.read_csv(f3, delimiter=";")

# Making plots for problem 10 
plt.figure(0)
plt.plot(data1["Frequency"], data1["N"]/100.0, label=r'$f=0.1$') 
plt.plot(data2["Frequency"], data2["N"]/100.0, label=r'$f=0.4$') 
plt.plot(data3["Frequency"], data3["N"]/100.0, label=r'$f=0.7$') 
plt.ylabel(r'$\frac{N}{N_0}$', fontsize=MEDIUM) 
plt.xlabel(r'$\omega_v$ (MHz)', fontsize=MEDIUM) 
plt.xticks([0.35, 0.40, 0.45, 0.50, 0.55])
plt.xticks(fontsize=MEDIUM); plt.yticks(fontsize=MEDIUM)
plt.title(r'Fraction left in trap vs applied $\omega_v$', fontsize=MEDIUM)
plt.grid(True)
plt.legend(fontsize=MEDIUM)
plt.tight_layout()
save_fig("N100range_0.35_0.55") 
plt.show()
"""
MAKING PLOTS IN SPECIFIC INDEX RANGE
OF PARTICLE TRAJECTORIES
"""
"""
indexes = [i for i in range(3000)]
N = 5 # number of particles to make plots for. 
for jj in range(N):
    particle = jj+1	
    filename = "PenningData_" + algo + str(particle) + "_f_0_AMP_0off.csv" # modify filename to read proper files
    data = pd.read_csv(filename, delimiter=";")
    #d = data.to_numpy()  ;
    
    plt.figure(1)
    plt.plot(data["Time"][indexes],data["pos_z"][indexes], label="Particle {}".format(particle))
    
    plt.figure(2)
    plt.plot(data["pos_x"][indexes],data["pos_y"][indexes], label="Particle {}".format(particle), marker="*" , markerfacecolor="r", markeredgecolor="r", markersize=10.0, markevery=[0])
    
    plt.figure(3) 
    plt.plot(data["pos_x"][indexes], data["vel_x"][indexes], label = "Particle {}".format(particle))
plt.figure(1); plt.title(r'$z(t)$ over $t\in[0, 300]\mu s$ with $\omega_v\approx 2\omega_z$ and $f=0.4$', fontsize=MEDIUM); plt.xticks(fontsize=MEDIUM); plt.yticks(fontsize=MEDIUM); plt.xlabel(r'Time ($\mu s$)', fontsize=MEDIUM); plt.ylabel(r'z(t) ($\mu m$)', fontsize=MEDIUM); save_fig("zprob10_2zCOFF")
plt.figure(2); plt.xlabel(r'x in $\mu m$', fontsize=MEDIUM); plt.ylabel(r'y in $\mu m$', fontsize=MEDIUM); plt.xticks(fontsize=MEDIUM); plt.yticks(fontsize=MEDIUM); plt.title(r'X,Y-plane with $t\in[0, 500]\mu s$, $\omega_v\approx 2\omega_z$ and $f=0.4$.', fontsize=MEDIUM); plt.grid(True); save_fig("xyprob10_2zCOFF");
plt.figure(3); plt.legend(); plt.show()
"""

