#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 23 12:12:14 2021

@author: gianmarco

This short script plots the z-motion of one particle as function of simulation time.

"""


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

data = pd.read_csv("PenningData_RK4_A_1.csv", delimiter=";");
numdata = data.to_numpy();

plt.figure()
plt.plot(numdata[:,0], numdata[:,3])
plt.xlabel("t ($\mu$s)", fontsize=14)
plt.ylabel("z ($\mu$m)", fontsize=14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.grid()
plt.savefig("zmotion.pdf",dpi=300, bbox_inches='tight')
plt.show()
