import numpy as np
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# default units
units = {'energy': "ev",
		 'distance': "A"}

unitconv=1
## note curve fitting is not working with different unit settings
## currently curve fitting is working with unitconv=1 only
# unit converter
if unitconv == 1:
	eng_conv=96.485
	dis_conv=0.1
	units['energy']="kJ/mol"
	units['distance']="nm"
elif unitconv==2:
	eng_conv=23.061
	dis_conv=1
	units['energy']="kcal/mol"
	units['distance']="A"
else:
	eng_conv=1
	dis_conv=1

### Outputs from this script should be used carefully
### C12 and C6 obtained from curve fitting can be used in GROMACS but is little difficult in LAMMPS if used with other potentials
### Use sigma and epsilon for LAMMPS 

# O-O interactions
# A=63840.199 # ev
# rho=0.198913 # Angstrom
# C=27.899008 # ev/A^6
# REPLUSIVE=0
# # Jiawei values
# C6=3.33628E-03 # kJ/mol.nm^6
# C12=1.43228E-06 # kJ/mol.nm^12
# # curve fitting parameters
# min_fit=2.62*dis_conv
# max_fit=3.8*dis_conv
# xi=13.772 # either 13.772 or 12


## Ca-O interactions
# A=3161.6335
# rho=0.271511
# C=0.0000000
# REPLUSIVE=1
# C6=-5.62629E-03 # kJ/mol.nm^6
# C12= 4.80260E-07 # kJ/mol.nm^12
# ##curve fitting parameters
# min_fit=3.1*dis_conv
# max_fit=3.5*dis_conv
# xi=12 # either 13.772 or 12

## Ca-C interactions
A=120000000
rho=0.12
C=0.0000000
REPLUSIVE=1
C6=8.87246E-04 # kJ/mol.nm^6
C12=7.25825E-07 # kJ/mol.nm^12
## curve fitting parameters
min_fit=2.6*dis_conv
max_fit=20*dis_conv
xi=12 # either 13.772 or 12

# def LJfunc(x, a, b):
# 	return a/x**12-b/x**6

def LJfunc(x, D, R):
 	return D*(pow(R/x,12)-2*pow(R/x,6))

A=A*eng_conv
rho=rho*dis_conv
C=C*eng_conv*dis_conv**6

print ("-------------------------------")
print ("Units convertion: eng_conv=%g dis_conv=%g" %(eng_conv, dis_conv))
print ("Inputs: A=%g (%s) rho=%g (%s) C=%g (%s.%s^6)" %(A, units['energy'], rho, units['distance'], C, units['energy'], units['distance']))
x_data = np.linspace(0.5*dis_conv, 20*dis_conv, num=1e3)
y_data = A*np.exp(-x_data/rho)-C/x_data**6

D=A*(xi-6)/(6*np.exp(xi))
L1=6*np.exp(xi)/(xi-6)
R=xi*rho
L2=xi*R**6/(xi-6)
sigma=R/1.12246204830937298143

print ("Calculated variables: D=%g (%s), L1=%g (-), R=%g (%s), L2=%g (%s^6)" %(D, units['energy'], L1, R, units['distance'], L2, units['distance']))

U_B_Loose= D*(L1*np.exp(-x_data/rho)-L2/pow(x_data,6))
U_LJ_Fit= D*((R/x_data)**12-2*(R/x_data)**6)
U_LJ_Fit_Jiawei= C12/x_data**12-C6/x_data**6

if (REPLUSIVE):
	U_B_Loose= D*L1*np.exp(-x_data/rho)
	U_LJ_Fit= D*(R/x_data)**12

print ("Minmax of distance: x=%g (%s) y=%d (%s)" %(np.min(x_data), units['distance'], np.max(x_data), units['distance']))
print ("LJ parameters (TC Lim): epsilon=%g (%s), sigma=%g (%s)" % (D, units['energy'], sigma, units['distance']))
print ("                      : C_12=%g (%s.%s^12), C_6=%g (%s.%s^6)" % (D*R**12, units['energy'], units['distance'], 2*D*R**6, units['energy'], units['distance']))


# And plot it
plt.figure(figsize=(12, 8))
plt.plot(x_data, y_data, 'r-', label="Bukcingham Pot.")
plt.plot(x_data, U_B_Loose, 'g-', label="Bukcingham (Loose form) Pot.")
plt.plot(x_data, U_LJ_Fit, 'b-', label="LJ Fitted Pot. (TC. Lim) C12=%g %s.%s^12, C6=%g %s.%s^6" % (D*R**12, units['energy'], units['distance'], 2*D*R**6, units['energy'], units['distance']))
plt.plot(x_data, U_LJ_Fit_Jiawei, 'c-', label="LJ Fitted Pot. (Jiawei) C12=%g %s.%s^12, C6=%g %s.%s^6" % (C12, units['energy'], units['distance'], C6, units['energy'], units['distance']))

#sys.exit()
print ("Fitting range: min=%g (%s), max=%g (%s)" %(min_fit, units['distance'], max_fit, units['distance']))

x_data_range=x_data[(x_data>min_fit) & (x_data<max_fit)]
y_data_range=y_data[(x_data>min_fit) & (x_data<max_fit)]

#popt, pcov = curve_fit(LJfunc, x_data_range, y_data_range)
#popt, pcov = curve_fit(LJfuncREP, x_data_range, y_data_range)
# if(REPLUSIVE):
# 	b_fixed=0
# 	popt1,pcov1 = curve_fit(lambda x, a : LJfunc(x, a, b_fixed), x_data_range, y_data_range)
# 	popt1=np.append(popt1,b_fixed)

popt, pcov = curve_fit(LJfunc, x_data_range, y_data_range)
eps=popt[0]
sig=popt[1]/1.12246204830937298143
print ("LJ parameters from curve fitting: epsilon = %g (%s) , sigma = %g (%s)" %(eps, units['energy'], sig, units['distance']))
print ("                       C12 = %g (%s.%s^12) , C6 = %g (%s.%s^6)" %(4*eps*sig**12, units['energy'], units['distance'], 4*eps*sig**6, units['energy'], units['distance']))

# if(REPLUSIVE):
# 	plt.plot(x_data, LJfunc(x_data, *popt1), '+',label='LJ fit values: C12=%g %s.%s^12, C6=%g %s.%s^6' % (popt[0], units['energy'], units['distance'], b_fixed, units['energy'], units['distance']))
# 	print ("Curve 2nd fitting outpts: %g %g " %(popt1[0], popt1[1]))

plt.plot(x_data, LJfunc(x_data, *popt), '-',label='LJ fit values: C12=%g %s.%s^12, C6=%g %s.%s^6' % (4*eps*sig**12, units['energy'], units['distance'], 4*eps*sig**6, units['energy'], units['distance']))


# U_LJ_es1= 4*eps*((sig/x_data)**12-(sig/x_data)**6)
# plt.plot(x_data, U_LJ_es1, '+', label="LJ EPSSIG Pot.")
# U_LJ_es2= D*((R/x_data)**12-2*(R/x_data)**6)
# plt.plot(x_data, U_LJ_es2, 'o', label="LJ EPSSIG REP Pot.")

print ("-------------------------------")

plt.xlabel('distance (%s)' %(units['distance']))
plt.ylabel('engery (%s)' %(units['energy']))
plt.legend()
axes = plt.gca()
axes.set_ylim([-1*eng_conv,10*eng_conv])
plt.legend()
plt.show()
plt.savefig('fittingCurve.png')
