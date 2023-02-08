import numpy as np
import matplotlib.pyplot as plt

M = 201
h = 0.005
# V = np.zeros((M,M))
# # Double slit
# for i in range(98, 102, 1): 
# 	for j in range(0, 85, 1): 
# 		V[i][j]=1
# 	for j in range(95, 105, 1): 
# 		V[i][j]=1
# 	for j in range(115, M, 1):
# 		V[i][j]=1

# np.savetxt("double_slit.dat", V) 


# V = np.zeros((M,M))
# # Triple slit
# for i in range(98, 102, 1): 
# 	for j in range(0, 75, 1): 
# 		V[i][j] = 1
# 	for j in range(85, 95, 1): 
# 		V[i][j] = 1
# 	for j in range(105, 115, 1):
# 		V[i][j] = 1
# 	for j in range(125, M, 1): 
# 		V[i][j] = 1
		
# np.savetxt("triple_slit.dat", V) 

# V = np.zeros((M,M))
# # Single slit
# for i in range(98, 102, 1): 
# 	for j in range(0, 95, 1): 
# 		V[i][j] = 1
# 	for j in range(105, M, 1): 
# 		V[i][j] = 1
		
# np.savetxt("single_slit.dat", V) 

# V = np.zeros((M,M))
# # Wall
# for i in range(98, 102, 1): 
# 	for j in range(M): 
# 		V[i][j] = 1
		
# np.savetxt("barrier.dat", V) 

V = np.zeros((M,M))

omega=1;
m = 1;

for i in range(M):

	for j in range(M):
		
		V[i][j]  = 0.5*m*omega**2*((i*h-0.5)**2 + (j*h-0.5)**2)
		print(V[i][j], end=' ')
	print("\n")

np.savetxt("harmonic.dat", V)
		
