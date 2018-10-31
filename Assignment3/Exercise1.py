import astropy.constants as const
import numpy as np 
import matplotlib.pyplot as plt 
import sys as sys 

"""
Constants/Parameters

"""
n = int(1e4)			# Integration points
G = const.G.value		# Gravitational Constant
k_B = const.k_B.value	# Stefan-Boltzmann Constant
m_p = const.m_p.value 	# Proton mass
p_c = 1e-26				# Critical density
T = 1e4					# Temperature
mu = 0.5882				# Mean molecular weight
omega_l = 0.692			# Density parameter lambda
omega_m = 0.308			# Density parameter matter
omega_b = 0.048			# Density parameter baryons
omega_r = 0				# Density parameter radiation
H_0 = 2.193548387e-18	# Hubble constant today in cm



"""
Functions

"""

def JeansLength(z):

	"""
	Returns the JeansLength for a given redshift z.

	"""

	return (np.sqrt((np.pi*k_B*T)/(G*mu*m_p*omega_b*p_c)) * (1+z)**(-3/2) )


def JeansMass(z):

	"""
	Returns the JeansMass for a given redshift z.

	"""

	return (np.pi**(5/2)/6) * ((k_B*T)/(G*mu*m_p))**(3/2) * np.sqrt(omega_b*p_c) * (1+z)**(3/2)

def WaveNumber(z):

	"""
	Returns the wave number corresponding to the Jeans length.

	"""

	return (2*np.pi/JeansLength(z))


def H(z):

	"""
	Returns the Hubble parameter for a given redshift.

	"""

	H2 = H_0**2 * (omega_m*(1+z)**3 + omega_r*(1+z)**4 + omega_l*(1+z))

	return (np.sqrt(H2))


def VelocityWidth(z):

	"""
	Returns the velocity width of a massive object.

	"""

	return H(z)*JeansLength(z)

def ThermalBroadening(T):

	"""
	Returns the thermal broadening (scale) velocity of a protongas for a given temperature.

	"""

	return np.sqrt(2*k_B*T/m_p)


def printInfo(z):

	"""
	Printing info about the above functions for a given redshift z.

	"""
	print("=================================================")
	print("Calculated quantities for exercise 1:\n")
	print("Jeans Mass: M_J(z = %g) = %g \nJeans Length: lambda_J(z = %g) = %g \nWave number: k(z = %g) = %g\nVelocity Width(z = %g): %g\nThermal Broadening velocity (T = %g) = %g\
		"%(z,JeansMass(z),z, JeansLength(z),z, WaveNumber(z),z,VelocityWidth(z), T, ThermalBroadening(T)))
	print("=================================================")


"""
Plotting and calulculations

"""

z = float(sys.argv[1])
printInfo(z)








