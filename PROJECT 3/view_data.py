#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 14 11:00:15 2021

@author: gianmarco
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


#algo = "FE_"
algo = "RK4_"
fig_xymotion = plt.figure(1);
fig_zmotion = plt.figure(2);
N = 2;
for jj in range(N):
    filename = "PenningData_" + algo + str(jj+1) + ".csv" 
    data = pd.read_csv(filename, delimiter=";")
    d = data.to_numpy()  ;
    
    plt.figure(1)
    plt.plot(d[:,0],d[:,3])
    plt.title("Motion along z-axis")

    
    plt.figure(2)
    plt.plot(d[:,1],d[:,2])
    plt.title("Trajectory in x-y plane")
