import numpy as np 
import matplotlib 
import matplotlib.pyplot as plt
import pyarma as pa
from matplotlib.animation import FuncAnimation

M = 201; #number of points including boundary

h = 1/(M-1)  # space step 
deltaT = 2.5e-5 #timestep 

#load simulation results
psi_t2 = pa.cube();
psi_t2.load('born_prob.bin')
psi_t2 = np.array(psi_t2)
real_psi = pa.cube(); 
real_psi.load('real_wf.bin')
real_psi = np.array(real_psi)
imag_psi = pa.cube(); 
imag_psi.load('imag_wf.bin')
imag_psi = np.array(imag_psi)
#number of timesteps
Nt, Nrows, Ncols = psi_t2.shape;
print(Nt)
sum_of_probabilities = np.zeros(Nt); 
for jj in range(Nt): #transpose all matrices because conversion to numpy switches x and y
    psi_t2[jj] = psi_t2[jj].T
    real_psi[jj] = real_psi[jj].T
    imag_psi[jj] = imag_psi[jj].T
    sum_of_probabilities[jj] = np.ndarray.sum(psi_t2[jj]) 
xpoints = np.arange(0,1+h, h)
ypoints = np.arange(0,1+h, h)

tpoints = np.arange(0, Nt*deltaT, deltaT)
print(psi_t2.shape)
x,y = np.meshgrid(xpoints, ypoints, sparse=True)
fsize = 15
# Plotting 1D probability density behind slits
prob_givent = psi_t2[80]
prob_givenyt = prob_givent[100, :]
prob_givenxt = prob_givent[:, 160]
prob_givenxt = prob_givenxt/np.sum(prob_givenxt)
plt.figure()
plt.plot(y, prob_givenxt, label=r'$p(y|x=0.8, t=0.002)$')
plt.xlabel("y", fontsize=fsize)
plt.ylabel(r'$p(y|x=0.8, t=0.002)$', fontsize=fsize)
#plt.legend(fontsize=fsize)
plt.xticks(fontsize=fsize)
plt.yticks(fontsize=fsize)
plt.tight_layout()
plt.grid(True)
plt.savefig("Py_SS.pdf", format="pdf")

plt.figure()
plt.plot(y, prob_givenyt, label=r'$p(x|y=0.5, t=0.002)$')
plt.xlabel("x", fontsize=fsize)
plt.ylabel(r'$p(x|y=0.5, t=0.002)$', fontsize=fsize)
#plt.legend(fontsize=fsize)
plt.xticks(fontsize=fsize)
plt.yticks(fontsize=fsize)
plt.tight_layout()
plt.grid(True)
plt.savefig("Px.pdf", format="pdf")

#settings
t_min = tpoints[0];
x_min, x_max = xpoints[0], xpoints[-1]
y_min, y_max = ypoints[0], ypoints[-1]


fig = plt.figure()
axes = plt.gca()

mynorm = matplotlib.cm.colors.Normalize(vmin=0.0, vmax = np.max(psi_t2[0]) )
img = axes.imshow(psi_t2[0], extent=[x_min, x_max, y_min, y_max], cmap=plt.get_cmap("inferno"), norm=mynorm)

plt.xlabel("x", fontsize = fsize)
plt.ylabel("y", fontsize = fsize)
plt.xticks(fontsize= fsize)
plt.yticks(fontsize= fsize)

cbar = fig.colorbar(img, ax=axes)
cbar.set_label(r'$|u(x,y,t)|^2$', fontsize=fsize)
cbar.ax.tick_params(labelsize=fsize)

time_txt = plt.text(0.95, 0.95, "t = {:.3e}".format(t_min), color="white", 
                    horizontalalignment="right", verticalalignment="top", fontsize=fsize)

# Function that takes care of updating the z data and other things for each frame
def animation(i):
    # Normalize the colour scale to the current frame?
    mynorm = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(psi_t2[i]) )
    img.set_norm(mynorm)

    # Update z data
    img.set_data(psi_t2[i])

    # Update the time label
    current_time = t_min + i * deltaT
    time_txt.set_text("t = {:.3e}".format(current_time))

    return img

# Use matplotlib.animation.FuncAnimation to put it all together
anim = FuncAnimation(fig, animation, interval=1, frames=np.arange(0, Nt, 1), repeat=False, blit=0)
anim.save("Test.gif", writer="ffmpeg", fps=24)
# Run the animation!
plt.show()


def plot_imshow(data, extent_vec, colours, fsize, figure_name, colourbar_label): 
	fig = plt.figure()
	axes = plt.gca()
	mynorm = matplotlib.cm.colors.Normalize(vmin=np.min(data), vmax=np.max(data))
	img = axes.imshow(data, extent=extent_vec, cmap=plt.get_cmap(colours), norm=mynorm)
	plt.xlabel("x", fontsize=fsize)
	plt.ylabel("y", fontsize=fsize)
	plt.xticks(fontsize=fsize)
	plt.yticks(fontsize=fsize)
	cbar = fig.colorbar(img, ax=axes)
	cbar.set_label(colourbar_label, fontsize=fsize)
	cbar.ax.tick_params(labelsize=fsize)
	plt.savefig(figure_name, format="pdf")
	
def plot_imshow_side_by_side(data1, data2, extent_vec, colours, fsize, figure_name, colourbar_label): 
	fig = plt.figure()
	ax1 = fig.add_subplot(1,2,1) 
	ax1 = plt.gca()
	mynorm = matplotlib.cm.colors.Normalize(vmin=np.min(data1), vmax=np.max(data1))
	img = ax1.imshow(data1, extent=extent_vec, cmap=plt.get_cmap(colours), norm=mynorm)
	plt.xlabel("x", fontsize=fsize)
	plt.ylabel("y", fontsize=fsize)
	plt.xticks(fontsize=fsize)
	plt.yticks(fontsize=fsize)
	ax2 = fig.add_subplot(1,2,2)
	ax2 = plt.gca()
	img2 = ax2.imshow(data2, extent=extent_vec, cmap=plt.get_cmap(colours), norm=mynorm)
	plt.xlabel("x", fontsize=fsize)
	plt.ylabel("y", fontsize=fsize)
	plt.xticks(fontsize=fsize)
	plt.yticks(fontsize=fsize)
	cbar = fig.colorbar(img2, ax=ax1)
	cbar.set_label(colourbar_label, fontsize=fsize)
	cbar.ax.tick_params(labelsize=fsize)
	plt.savefig(figure_name, format="pdf")
	

# imshows
colour_scheme = "inferno"
# Probability at different time steps
plot_imshow(psi_t2[0], [x_min, x_max, y_min, y_max], colour_scheme, fsize, "p8_P_0.pdf", r'$|u(x,y,t=0.000)|^2$')
plot_imshow(psi_t2[40], [x_min, x_max, y_min, y_max], colour_scheme, fsize, "p8_P_1.pdf", r'$|u(x,y,t=0.001)|^2$')
plot_imshow(psi_t2[80], [x_min, x_max, y_min, y_max], colour_scheme, fsize, "p8_P_2.pdf", r'$|u(x,y,t=0.002)|^2$')

# Real part of wavefunction at different time steps
plot_imshow(real_psi[0], [x_min, x_max, y_min, y_max], colour_scheme, fsize, "p8_R_0.pdf", r'Re$(u(x,y,t=0.000))$')
plot_imshow(real_psi[40], [x_min, x_max, y_min, y_max], colour_scheme, fsize, "p8_R_1.pdf", r'Re$(u(x,y,t=0.001))$')
plot_imshow(real_psi[80], [x_min, x_max, y_min, y_max], colour_scheme, fsize, "p8_R_2.pdf", r'Re$(u(x,y,t=0.002))$')

# Real part of wavefunction at different time steps
plot_imshow(imag_psi[0], [x_min, x_max, y_min, y_max], colour_scheme, fsize, "p8_I_0.pdf", r'Im$(u(x,y,t=0.000))$')
plot_imshow(imag_psi[40], [x_min, x_max, y_min, y_max], colour_scheme, fsize, "p8_I_1.pdf", r'Im$(u(x,y,t=0.001))$')
plot_imshow(imag_psi[80], [x_min, x_max, y_min, y_max], colour_scheme, fsize, "p8_I_2.pdf", r'Im$(u(x,y,t=0.002))$')


"""
plot_imshow_side_by_side(real_psi[0], imag_psi[0], [x_min, x_max, y_min, y_max], colour_scheme, fsize, "p8_IR_0.pdf", r'$\Psi(x,y,t=0.000)$')


# Plotting total probability against time
plt.figure()
plt.plot(tpoints, abs(sum_of_probabilities-1), label=r'$\Delta=\sum_{\ell, j}|u_{\ell,j}|^2-1$')
plt.xlabel("Time, t", fontsize=fsize) 
plt.ylabel(r'$|\Delta|$', fontsize=fsize)
plt.xticks([0.000, 0.002, 0.004, 0.006, 0.008], fontsize=fsize)
plt.yticks(fontsize=fsize)
plt.legend(fontsize=fsize)
plt.grid(True)
plt.tight_layout()
plt.savefig("Dev_DS.pdf", format="pdf")
"""

