#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 14 11:00:15 2021

@author: gianmarco
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os


#algo = "FE_"
algo = "RK4_"
N = 2; #number of particles
file_label = ['A', 'B', 'C', 'D', 'E'] #this could be useful if you ran code for more than one value of h
                                    # however we don't do this for the case of 2 particles.
d_on=[]
d_off=[]
export = False
onoff = ["on", "off"]
experiment_number = 10;
for ii in range(1,3):
    os.chdir("./"+str(experiment_number)+"/experiment_"+str(ii) ) 
    for jj in range(N):
        
        filename = "PenningData_" + algo + file_label[0] + "_" + str(jj+1)+".csv" 
        data = pd.read_csv(filename, delimiter=";", header=None)
        if ii==1:
            d_off.append( data.to_numpy()  )
        else:
            d_on.append(data.to_numpy() )
    os.chdir("../../");
print("Initial conditions:")
for kk in range(2):
    print("Particle " + str(kk+1) + ":")
    print("position:")
    print(np.transpose(d_off[kk][0,1:4]))
    print("velocity:")
    print(np.transpose(d_off[kk][0,4:7]))


#now that I have all of my data imported, start plotting stuff

#%%first plot the motion on the z axis 

for ii in range(2):
    plt.figure(ii+1)
    if ii==0:
        plt.plot(d_on[0][:,0], d_on[0][:,3], label="Particle 1")
        plt.plot(d_on[1][:,0], d_on[1][:,3], label="Particle 2")
    else:
        plt.plot(d_off[0][:,0], d_off[0][:,3],label="Particle 1")
        plt.plot(d_off[1][:,0], d_off[1][:,3],label="Particle 2")
    
    
    plt.xlabel("t ($\mu$s)", fontsize=12)
    plt.ylabel("z ($\mu$m)", fontsize=12)
    plt.legend(fontsize=12)
    plt.grid()
    
#%%then plot trajectories in xy plane

for ii in range(2):
    plt.figure(ii+3)
    if ii==0:
        plt.plot(d_on[0][:,1], d_on[0][:,2], label="Particle 1")
        plt.plot(d_on[1][:,1], d_on[1][:,2], label="Particle 2")
        plt.scatter(d_on[0][0,1], d_on[0][0,2], marker='*', label="Start Particle 1")
        plt.scatter(d_on[1][0,1], d_on[1][0,2], marker='*', label="Start Particle 2")
    else:
        plt.plot(d_off[0][:,1], d_off[0][:,2],label="Particle 1")
        plt.plot(d_off[1][:,1], d_off[1][:,2],label="Particle 2")
        plt.scatter(d_off[0][0,1], d_off[0][0,2], marker='*', label="Start Particle 1")
        plt.scatter(d_off[1][0,1], d_off[1][0,2], marker='*', label="Start Particle 2")
    plt.xlabel("x ($\mu$m)", fontsize=13)
    plt.ylabel("y ($\mu$m)", fontsize=13)
    # dx = 25
    # dy = 5
    # plt.xlim(-2000-dx, -2000+dx)
    # plt.ylim(-4000-dy, -4000+dy)
    plt.legend(fontsize=14)
    plt.xticks(fontsize=13)
    plt.yticks(fontsize=13)
    plt.grid()
    if export:
        plt.savefig( "xyplot_large" + onoff[ii]+".pdf", dpi=300,  bbox_inches='tight')

#%% phase space plots
for ll in range(3): #loops on x y and z
    labels = "xyz"
    for ii in range(2):
        fig = plt.figure(ii+5+2*ll)
        ax = fig.add_subplot()
        if ii==0:
            plt.plot(d_on[0][:,1+ll], d_on[0][:,4+ll], label="Particle 1")
            plt.plot(d_on[1][:,1+ll], d_on[1][:,4+ll], label="Particle 2")
            plt.scatter(d_on[0][0,1+ll], d_on[0][0,4+ll], marker='*', label="Start Particle 1")
            plt.scatter(d_on[1][0,1+ll], d_on[1][0,4+ll], marker='*', label="Start Particle 2")
        else:
            plt.plot(d_off[0][:,1+ll], d_off[0][:,4+ll],label="Particle 1")
            plt.plot(d_off[1][:,1+ll], d_off[1][:,4+ll],label="Particle 2")
            plt.scatter(d_off[0][0,1+ll], d_off[0][0,4+ll], marker='*' )
            plt.scatter(d_off[1][0,1+ll], d_off[1][0,4+ll], marker='*')
        plt.xlabel(labels[ll]+" ($\mu$m)", fontsize=12)
        plt.ylabel("$v_"+labels[ll]+"$ ($\mu$m)", fontsize=12)
        plt.legend( loc="best", fontsize=12, framealpha=0.5)
        
        plt.grid()

        ax.text()
#%%finally 3d plot of trajectories

for ii in range(2):
    fig = plt.figure(ii+11)
    ax = fig.add_subplot(111, projection='3d')
    
    if ii==0:
         ax.plot3D(d_on[0][:,1], d_on[0][:,2], d_on[0][:,3] ,label="Particle 1")
         ax.plot3D(d_on[1][:,1], d_on[1][:,2], d_on[1][:,3], label="Particle 2")
    else:
         ax.plot3D(d_off[0][:,1], d_off[0][:,2],d_off[0][:,3],label="Particle 1")
         ax.plot3D(d_off[1][:,1], d_off[1][:,2],d_off[1][:,3],label="Particle 2")
    
    plt.legend()
    
#%% distance (t)  plot
mm = np.size(d_on[0][:,1])
r = np.zeros((mm,1))
for j in range(mm):
    r[j] = np.linalg.norm(d_on[0][j,1:4]-d_on[1][j,1:4] )
t = d_on[0][:,0]
plt.figure(13)
plt.plot(t,r, label="Distance between particles")
plt.xlabel("t ($\mu$s)", fontsize=13)
plt.ylabel("$\Delta r$ ($\mu$m)", fontsize=13)
plt.grid()
plt.legend(fontsize=13)
plt.xticks(fontsize=13)
plt.yticks(fontsize=13)
if export:
    plt.savefig("distance(t).pdf", dpi=300, bbox_inches='tight')
