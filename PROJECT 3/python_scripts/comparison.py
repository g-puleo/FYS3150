#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 21 17:56:17 2021

@author: gianmarco
"""

#this code compares the analytical predictions with the numerical solution

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
doublecheck = False ; #set this to True if you want to plot also a comparison
                      #between numerical and analytical trajectory in xy plane

export = True; 
                    
#os.chdir("./rel_err_FE");
#algo = "RK4"
algos = ["FE", "RK4"]

for algo in algos:
    os.chdir("./rel_err_"+ algo)
    N= 5 #number of different stepsizes
    h = np.array([0.11, 0.05, 0.023, 0.01, 0.0045])
    file_label = ['A', 'B', 'C', 'D', 'E']
    plt.figure(1, figsize=(7,7.2))
    
    #plt.title("Relative error as function of time (" + algo + ")")
    plt.xlabel("t ($\mu$s)", fontsize=14);
    plt.ylabel("relative error", fontsize=14);
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    # if algo == "RK4":
    #     plt.ylim(1e-20, 1e-4)
    DeltaMax = []
    for j in range(N):
        ana_data  = pd.read_csv("r_ana_"+ file_label[j] + "_1.csv", delimiter=";")
        an_array = ana_data.to_numpy() ; 
        num_data = pd.read_csv("PenningData_" + algo + "_" + file_label[j] + "_1.csv")
        #print(len(an_array[:,0]))
         #plt.plot()
        # plt.plot(an_array[:,0], an_array[:,3] , label="analytical result")
        plt.plot(an_array[:,0], an_array[:,4], label=f"h={h[j]} $\mu$s")
        
        #this is to calculate the error convergence rate
        DeltaMax.append(np.max(an_array[:,5]))
        
    plt.grid()
    ax = plt.gca()
    ax.set_yscale('log')
    #ax.set_ylim(1e-20, 7e-4)
    plt.legend(loc=4, fontsize=14)
    if export:
        plt.savefig("REL_ERR"+algo+".pdf", dpi=300, bbox_inches='tight')
    plt.show()
    
    #calculate the error convergence rate (ECR)
    ECR  = 0;
    if algo ==  "RK4":
        nmax = 5
    else:
        nmax = 5
        
    for k in range(1,nmax):
       ECR = ECR+ np.log(DeltaMax[k]/DeltaMax[k-1])/np.log(h[k]/h[k-1])
    ECR = 0.25*ECR;
    print("ECR for " + algo +  ":" , ECR )
    
    
    
    if doublecheck:
        for j in range(N):
            plt.figure(j+2)
            ana_data  = pd.read_csv("r_ana_"+ file_label[j] + "_1.csv", delimiter=";")
            an_array = ana_data.to_numpy() ; 
            num_data = pd.read_csv("PenningData_" + algo + "_" + file_label[j] + "_1.csv", delimiter=";")
               #plt.plot()
            # plt.plot(an_array[:,0], an_array[:,3] , label="analytical result")
            num_array = num_data.to_numpy()
            #plot xy curve
            plt.plot(an_array[:,1], an_array[:,2], label="analytical solution")
            plt.plot(num_array[:,1], num_array[:,2], label="numerical results")
            plt.legend()
            plt.title(f"xy trajectory of particle (h={h[j]} $\mu$s)")
            plt.xlabel("x ($\mu$m)", fontsize=12)
            plt.ylabel("y ($\mu$m)", fontsize=12)
            plt.grid()
            plt.show()
    
    os.chdir("../")
    
    
    
    
    