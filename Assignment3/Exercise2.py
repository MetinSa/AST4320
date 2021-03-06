import numpy as np 
import matplotlib.pyplot as plt 
from astropy import constants as const

"""
Global parameters/constants

"""

n = int(1e4)					# number of integration points
c = const.c.cgs.value			# Speed of light in cm
sigma = 6.65e-25				# Cross-section electron scattering
omega_l = 0.692					# CDM density fraction
omega_m = 0.308					# Mass density fraction
omega_r = 0						# Radiation density fraction
H_0 = 2.193548387e-18			# Hubble constant today in cm


"""
Functions

"""

def n_e(z):

	"""
	Returns the numberdensity of electrons in IGM. 
	(assuming that IGM consists of only fully ionized hydrogen)

	"""

	return 1.9e-7 * (1+z)**3


def H(z):

	"""
	Returns the Hubble parameter for a given redshift.

	"""

	H2 = H_0**2 * (omega_l*(1+z) + omega_m*(1+z)**3 + omega_r*(1+z)**4)

	return (np.sqrt(H2))


def OpticalDepth(z_min, z_max):

	"""
	Returns the Optical depth of the ionised IGM at a given redshift range.

	"""

	z, dz	= np.linspace(z_min,z_max, n, retstep = True)		# Redshift-array
	tau = np.zeros(n)									# Optical depth-array
	tau[0] = (n_e(z[0])*dz)/((1+z[0])*H(z[0]))			# Initial condition for tau

	# Computing tau
	for i in range(n-1):
		tau[i+1] = tau[i] + (n_e(z[i])*dz)/((1+z[i])*H(z[i]))

	tau = c*sigma*tau

	return z, tau


def printInfo(z, tau, z_min, z_max, z_any):

	print("Tau(z = %g) = %g \nTau(z = %g) = %g \nTau(z = %g) = %g" %(z_min, tau[0], z_max, tau[-1],\
		z_any, tau[np.where(z >= z_any)[0][0]]))

"""
Plotting and printing the results.

"""
z_min = 0
z_max = 10
z, tau = OpticalDepth(z_min,z_max)

printInfo(z, tau, z_min, z_max, 3)

plt.plot(z,tau, color = "royalblue", label = r" Optical Depth, $\tau_e(z)$")
plt.grid(linestyle="--")
plt.legend()
plt.xlabel("z")
plt.ylabel(r"$\tau_e$")
# plt.savefig("figures/tau.pdf")
plt.show()