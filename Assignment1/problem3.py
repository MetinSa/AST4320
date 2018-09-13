import numpy as np
import matplotlib.pyplot as plt

def constants():

	"""
	Calculating constants.

	"""
	T_dec = 3000			#K, kelvin
	z_dec = 1090			#Redshift
	a_dec = 1/(z_dec-1)

	B = T_dec*a_dec
	C = T_dec*a_dec**2

	return B, C


def T_gamma(a,B):

	"""
	Temperature of photons.

	"""

	return B/a


def T_gas(a,C):

	"""
	Temperature of gas.

	"""

	return C/(a**2)



def main():
	#gridpoints
	n = 1000

	#scalefactor a
	a = np.linspace(1e-4, 1, n)

	#computing the constants
	B, C = constants()

	#computing the temperatures
	T_rad = T_gamma(a,B)
	T_g = T_gas(a,C)

	# plotting
	plt.loglog(a,T_rad, label = r"$T_{rad}$")
	plt.loglog(a,T_g, label = r"$T_{gas}$")
	plt.title("Temperature of Gas and Radiation")
	plt.xlabel("$a$")
	plt.ylabel("T[k]")
	plt.legend()
	plt.grid(linestyle="--")
	# plt.savefig("Temp.pdf")
	plt.show()

main()























