import astropy.constants as const
import numpy as np 
import matplotlib.pyplot as plt 

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



"""
Functions

"""

def JeansLength(z):

	"""
	Returns the JeansLength for a given redshift z.

	"""

	return (np.sqrt((np.pi*k_B*T)/(G*mu*m_p*omega_m*p_c)) * (1+z)**(-3/2) )


def JeansMass(z):

	"""
	Returns the JeansMass for a given redshift z.

	"""

	return (np.pi**(5/2)/6) * ((k_B*T)/(G*mu*m_p))**(3/2) * np.sqrt(omega_m*p_c) * (1+z)**(3/2)

def WaveNumber(z):

	"""
	Returns the wave number corresponding to the Jeans length.

	"""

	return (2*np.pi/JeansLength(z))

def printInfo(z):

	"""
	Printing info about the above functions for a given redshift z.

	"""

	print("Jeans Mass today: M_J(%g) = %g \nJeans Length today: lambda_J(%g) = %g \nWave number today: k(%g) = %g"\
	 %(z,JeansMass(z),z, JeansLength(z),z, WaveNumber(z)))



"""
Plotting and calulculations

"""

z = 4
# z = np.linspace(0,10,n)
# M_J = JeansMass(z)
# l_J = JeansLength(z)
# plt.plot(z,M_J)
# plt.plot(z,l_J)
# plt.show()
printInfo(z)








