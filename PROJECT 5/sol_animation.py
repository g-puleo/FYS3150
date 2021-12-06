import numpy as np 
import matplotlib 
import matplotlib.pyplot as plt
import pyarma as pa
from matplotlib.animation import FuncAnimation

M = 10; #number of points including boundary

h = 1/(M-1)  # space step 
deltaT = 2.5e-5 #timestep 

#load simulation results
psi_t = pa.cube();
psi_t.load('schrod.bin')
psi_t = np.array(psi_t)
#number of timesteps
Nt, Nrows, Ncols = psi_t.shape;
for jj in range(Nt): #transpose all matrices because conversion to numpy switches x and y
    psi_t[jj] = psi_t[jj].T
xpoints = np.arange(0,1+h, h)
ypoints = np.arange(0,1+h, h)

tpoints = np.arange(0, Nt*deltaT, deltaT)

x,y = np.meshgrid(xpoints, ypoints, sparse=True)



#settings
fsize = 12
t_min = tpoints[0];
x_min, x_max = xpoints[0], xpoints[-1]
y_min, y_max = ypoints[0], ypoints[-1]


fig = plt.figure()
axes = plt.gca()

mynorm = matplotlib.cm.colors.Normalize(vmin=0.0, vmax = np.max(psi_t[0]) )
img = axes.imshow(psi_t[0], extent=[x_min, x_max, y_min, y_max], cmap=plt.get_cmap("viridis"), norm=mynorm)

plt.xlabel("x", fontsize = fsize)
plt.ylabel("y", fontsize = fsize)
plt.xticks(fontsize= fsize)
plt.yticks(fontsize= fsize)

cbar = fig.colorbar(img, ax=axes)
cbar.set_label("z(x,y,t)", fontsize=fsize)
cbar.ax.tick_params(labelsize=fsize)

time_txt = plt.text(0.95, 0.95, "t = {:.3e}".format(t_min), color="white", 
                    horizontalalignment="right", verticalalignment="top", fontsize=fsize)

# Function that takes care of updating the z data and other things for each frame
def animation(i):
    # Normalize the colour scale to the current frame?
    mynorm = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(psi_t[i]) )
    img.set_norm(mynorm)

    # Update z data
    img.set_data(psi_t[i])

    # Update the time label
    current_time = t_min + i * deltaT
    time_txt.set_text("t = {:.3e}".format(current_time))

    return img

# Use matplotlib.animation.FuncAnimation to put it all together
anim = FuncAnimation(fig, animation, interval=1, frames=np.arange(0, Nt, 2), repeat=False, blit=0)

# Run the animation!
plt.show()
