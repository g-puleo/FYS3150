import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

Tstr = 'Temperature  ($J/k_B$)'

# defining R2
def R2(y_data, y_model):
    return 1 - np.sum((y_data - y_model) ** 2) / np.sum((y_data - np.mean(y_data)) ** 2)

# defining MSE
def MSE(y_data, y_model): 
    return np.mean((y_data-y_model)**2)  
   
# function that finds the optimal polynomial fit of the data
def find_optimal_polydegree(maxDegree, x, target): 

	polynomial_degree = np.zeros(maxDegree+1)
	r2_score = np.zeros(maxDegree+1)
	for polydegree in range(1, maxDegree+1): 
		polynomial_degree[polydegree] = polydegree
		X = np.zeros((len(x), polydegree+1))
		X[:, 0] = 1
		for degree in range(polydegree+1): 
			X[:, degree] = x**(degree)
		OLS_coefficients = np.linalg.pinv(X.T @ X) @ X.T @ target
		prediction = X @ OLS_coefficients
		MSE_score[polydegree] = MSE(target, prediction)
	
	best_fit = np.argmax(MSE_score)
	return best_fit

# function that returns the coefficients, prediction and r2 score of a fit of the data
def make_polyfit(polydegree, x, target): 
	X = np.zeros((len(x), polydegree+1))
	X[:, 0] = 1
	for degree in range(polydegree+1): 
		X[:, degree] = x**(degree)
	OLS_coefficients = np.linalg.pinv(X.T @ X) @ X.T @ target
	prediction = X @ OLS_coefficients
	r2 = R2(target, prediction)
	return OLS_coefficients, prediction, r2

# fontsize
NORMAL = 14

# Making plots of data
colors = ['b', 'g', 'r', 'c', 'm']
critical_temperatures_susc = []
critical_temperatures_cv = []
for L in ['40', '50', '60' , '80', '100']:
	if L=='40': 
		col = colors[0]
	if L=='50': 
		col = colors[1]
	if L=='60': 
		col = colors[2]
	if L=='80': 
		col = colors[3]
	if L=='100': 
		col = colors[4]
	filename = 'tempdata'+ L + 'F.csv'
	data = pd.read_csv(filename)
	temps = data["T"][-50:-1]
	Llabel = 'L = ' + L
	energies = data["average energy"][-50:-1]
	
	optimal_energy_polydegree = find_optimal_polydegree(100, temps, energies)
	coeffs, polyFitEnergies, r2Energies = make_polyfit(optimal_energy_polydegree, temps, energies)
	plt.figure(1)
	plt.plot(temps, energies, linestyle ='none', marker='.',markeredgecolor=col, markerfacecolor=col, label=Llabel)
	plt.plot(temps, polyFitEnergies, linestyle="dashed", color=col)
	plt.ylabel(r'$\langle\epsilon\rangle$ ($j/$spin)', fontsize=NORMAL)
	plt.xlabel(Tstr, fontsize=NORMAL)
	plt.yticks([-1.45, -1.40, -1.35, -1.30, -1.25])
	
	mags = data["average magnetization"][-50:-1]
	optimal_mags_polydegree = find_optimal_polydegree(100, temps, mags)
	coeffs, polyFitMags, r2Mags = make_polyfit(optimal_mags_polydegree, temps, mags)
	plt.figure(2)
	plt.plot(temps, mags, linestyle='none', marker = '.',markeredgecolor=col, markerfacecolor=col, label=Llabel)
	plt.plot(temps, polyFitMags, linestyle="dashed", color = col)
	plt.ylabel(r'$\langle |m|\rangle$ (1/spin)', fontsize=NORMAL)
	plt.xlabel(Tstr, fontsize=NORMAL)
	
	cv = data["specific heat"][-50:-1]
	optimal_cv_polydegree = find_optimal_polydegree(100, temps, cv)
	coeffs, polyFitCv, r2Cv = make_polyfit(optimal_cv_polydegree, temps, cv)
	plt.figure(3)
	plt.plot(temps, cv,linestyle='none', marker  = '.', markeredgecolor=col, markerfacecolor=col, label=Llabel);
	plt.plot(temps, polyFitCv, linestyle="dashed", color = col)
	plt.xlabel(Tstr, fontsize=NORMAL)
	plt.ylabel('Specific heat $C_V$ ($k_B/$spin)', fontsize=NORMAL)

	susc = data["susceptibility"][-50:-1]
	optimal_susc_polydegree = find_optimal_polydegree(100, temps, susc)
	coeffs, polyFitSusc, r2Susc = make_polyfit(optimal_susc_polydegree, temps, susc)
	plt.figure(4)
	plt.plot(temps, susc,linestyle='none', marker  = '.', markeredgecolor=col, markerfacecolor=col, label=Llabel);
	plt.plot(temps, polyFitSusc, linestyle="dashed", color=col)
	plt.xlabel(Tstr, fontsize=NORMAL)
	plt.ylabel('Magnetic susceptibility $\chi$ (spin/$J$)', fontsize=NORMAL)
	
	print("R2 scores: ") 
	print("Energy: %.2f"%(r2Energies))
	print("Magnetization: %.2f"%(r2Mags))
	print("Specific heat: %.2f"%(r2Cv))
	print("Magnetic susceptibility: %.2f"%(r2Susc))
	
	print("Temperatures where the gradient of the fit is zero: ") 
	print("Index of max: ", np.argmax(polyFitSusc))
	critical_temperatures_susc.append(2.272 + np.argmax(polyFitSusc)*0.002)
	critical_temperatures_cv.append(2.252 + np.argmax(polyFitCv)*0.002)
print("Tc susc: ")
print(critical_temperatures_susc)
print("Tc cv: ") 
print(critical_temperatures_cv)
	
# Predicting a straight line for Tc using x=1/L
x = np.array([1.0/40.0, 1.0/50.0, 1.0/60.0, 1.0/80.0, 1.0/100.0])

# Calculating optimal fitted line for Tc
degree = 1
coefficients_susc, best_fit_susc, r2susc = make_polyfit(degree, x, critical_temperatures_susc)
coefficients_cv, best_fit_cv, r2cv = make_polyfit(degree, x, critical_temperatures_cv)


print(coefficients_susc) 
print(coefficients_cv)

# Finalizing plots of data
savestr = ["avg_nrg", "avg_mag", "spec_heat", "mag_susc"]
for ii in range(1,5):
	plt.xticks(fontsize=NORMAL)
	plt.yticks(fontsize=NORMAL)
	plt.figure(ii)
	plt.grid(visible=True)
	plt.legend(fontsize=NORMAL)
	textstr = savestr[ii-1]
	plt.tight_layout()
	plt.savefig("Temperature_" + textstr + ".pdf", format="pdf")
plt.show()




# Calculating R2-scores
R2cvfit = R2(best_fit_cv, critical_temperatures_cv)
R2suscfit = R2(best_fit_susc, critical_temperatures_susc)


# Plotting best fitted lines of the calculated Tc temperatures
plt.plot(x, best_fit_susc, linestyle="dashed", label="Fitted lattice size dependence") 
plt.plot(x, critical_temperatures_susc, linestyle="none", marker=".", label=r'$T_C$ of polynomial fits ($\chi(L)$)')
plt.xlabel(r'Inverted lattice size, ($1/$L)', fontsize=NORMAL)
plt.xticks([0.01, 0.0175 ,0.025],fontsize=NORMAL)
plt.yticks(fontsize=NORMAL)
plt.ylabel(r'Critical temperature $T_c$ ($J/k_B$)', fontsize=NORMAL)
texstr = '\n'.join(('Properties of fitted line:', 
	 r'$T_C(1/L) = \alpha(\frac{1}{L}) + T_C(0)$,', 
	 r'$T_C(0)\approx %.4f$, $\alpha \approx %.1f$,'%(coefficients_susc[0],coefficients_susc[1], ), 
	 r'$R^2=%.2f$.'%(R2suscfit, )))
props = dict(boxstyle='round', facecolor='white', alpha=1.0)

# place a text box in upper left in axes coords
plt.text(0.0172, 2.3, texstr, fontsize=NORMAL,
        verticalalignment='top', bbox=props)
plt.legend(fontsize=12)
plt.grid(visible=True)
plt.tight_layout()
plt.savefig("FitSUSC.pdf", format="pdf")
plt.show()

plt.plot(x, best_fit_cv, linestyle="dashed", label="Fitted lattice size dependence") 
plt.plot(x, critical_temperatures_cv, linestyle="none", marker=".", label=r'$T_C$s of polynomial fits ($C_V(L)$)')
plt.xlabel(r'Inverted lattice size, ($1/$L)', fontsize=NORMAL)
plt.ylabel(r'Critical temperature $T_c$ ($J/k_B$)', fontsize=NORMAL)
plt.xticks([0.01, 0.0175, 0.025],fontsize=NORMAL)
plt.yticks(fontsize=NORMAL)
texstr = '\n'.join(('Properties of fitted line:',
	 r'$T_C(\frac{1}{L}) = \alpha(\frac{1}{L}) + T_C(0)$,', 
	 r'$T_C(0)\approx %.4f$, $\alpha\approx %.1f$,'%(coefficients_cv[0],coefficients_cv[1], ), 
	 r'$R^2=%.2f$.'%(R2cvfit, )))
props = dict(boxstyle='round', facecolor='white', alpha=1.0)

# place a text box in upper left in axes coords
plt.text(0.0174, 2.28175, texstr, fontsize=NORMAL,
        verticalalignment='top', bbox=props)
plt.grid(visible=True)
plt.legend(fontsize=NORMAL)
plt.tight_layout()
plt.savefig("FitCV.pdf", format="pdf")
plt.show()

