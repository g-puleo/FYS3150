import numpy as np
import matplotlib.pyplot as plt

M = 201

V = np.zeros((M,M))
# Double slit
for i in range(98, 102, 1): 
	for j in range(0, 86, 1): 
		V[i][j]=1
	for j in range(95, 106, 1): 
		V[i][j]=1
	for j in range(115, M, 1):
		V[i][j]=1

np.savetxt("double_slit.dat", V) 


V = np.zeros((M,M))
# Triple slit
for i in range(98, 102, 1): 
	for j in range(0, 76, 1): 
		V[i][j] = 1
	for j in range(85, 96, 1): 
		V[i][j] = 1
	for j in range(105, 116, 1):
		V[i][j] = 1
	for j in range(125, M, 1): 
		V[i][j] = 1
		
np.savetxt("triple_slit.dat", V) 

V = np.zeros((M,M))
# Single slit
for i in range(98, 102, 1): 
	for j in range(0, 96, 1): 
		V[i][j] = 1
	for j in range(105, M, 1): 
		V[i][j] = 1
		
np.savetxt("single_slit.dat", V) 

V = np.zeros((M,M))
# Wall
for i in range(98, 102, 1): 
	for j in range(0, M): 
		V[i][j] = 1
		
np.savetxt("no_slit.dat", V) 
		
