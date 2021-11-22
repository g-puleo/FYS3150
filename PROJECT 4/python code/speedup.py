import numpy as np
import matplotlib.pyplot as plt

times = np.zeros( (4, 6))
times[0,:] = [ 8.298, 8.313, 8.306, 8.295, 8.286, 8.294 ]
times[1,:] = [ 7.447, 7.775, 7.390, 7.867, 7.778, 7.484 ]
times[2,:] = [ 7.328, 7.222, 6.691, 6.587, 6.677, 6.560 ]
times[3,:] = [ 6.901, 6.836, 6.834, 6.788, 6.781, 6.593 ]

numthreads = [ 1, 2, 3, 4] 

speedup = np.zeros( (4,1) ) 

for jj in range(4):
	speedup[jj] = np.mean(times[0,:])/np.mean( times[jj,:] )
	print(speedup[jj], numthreads[jj])
plt.figure()
plt.plot( numthreads, speedup, marker='*', color='deepskyblue' ) 
plt.xlabel('number of threads', fontsize=14)
plt.ylabel('speed up factor' , fontsize=14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.grid(visible=True)
#plt.savefig('speedup.pdf', dpi=300, bbox_inches='tight')
#plt.show()

