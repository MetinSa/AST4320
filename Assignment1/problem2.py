import numpy as np
import matplotlib.pyplot as plt

def H(H_0, a, omega_m, omega_l):

	"""
	Hubble parameter H.
	
	"""

	return H_0 * np.sqrt(omega_m * a**(-3) + omega_l)



def a_dot(H_0, a, omega_m, omega_l):

	"""
	Derivative of scale factor a.

	"""
	return H(H_0, a, omega_m, omega_l)*a



def dvdt(H_0, a, omega_m, omega_l, u, v):

	"""
	ddelta/dt, where we sub delta = u for readability.

	"""

	h = H(H_0, a, omega_m, omega_l)			#computing H
	return ((3/2) * u * h**2) - 2*h*v



def EulerLoop(n, H_0, a_0, a_1, omega_m, omega_l, u_0):

	"""
	Integrationloop using Eulers Method.

	"""

	#making intial arrays and setting bonudary conditions
	a, da = np.linspace(a_0, a_1, n, retstep = True)
	z = (a_1/a) - 1
	u = np.zeros(n);				u[0] = u_0		
	v = np.zeros(n);				v[0] = a_dot(H_0, a[0], omega_m, omega_l)


	#integration
	for i in range(n-1):
		#computing dt
		dt = da/a_dot(H_0, a[i], omega_m, omega_l)

		v[i+1] = v[i] + dvdt(H_0, a[i], omega_m, omega_l, u[i], v[i])*dt
		u[i+1] = u[i] + v[i+1]*dt

	return a, u, z



def GrowthFactor(u, a):

	"""
	Computing the growth factor f

	"""

	return np.diff(np.log(u))/np.diff(np.log(a))



def main():

	"""
	Executes the program and displays the results.

	"""

	#constants and initial values
	H_0 = 2.269e-18
	u_0 = 1e-3
	a_0 = 1e-3
	a_1 = 1
	n = int(1e5)

	#fractional densities
	omega_list = np.array([1.0, 0.3, 0.8])

	#extracting the results
	for omega_m in omega_list:

		#starting by finding lambda fractional density
		omega_l = 1 - omega_m

		#computing a and u for the given densities
		a, u, z = EulerLoop(n, H_0, a_0, a_1, omega_m, omega_l, u_0)

		#plotting the results
		plt.loglog(a, u, label = r"$\Omega_m$ = %g , $\Omega_\Lambda$ = %g" %(omega_m, omega_l))

	plt.title("Time evolution of perturbations")
	plt.xlabel("$a$")
	plt.ylabel(r"$\delta$")
	plt.legend()
	plt.grid(linestyle="--")
	plt.savefig("2b.pdf")
	plt.show()

	#extracting the results
	for omega_m in omega_list:

		#starting by finding lambda fractional density
		omega_l = 1 - omega_m

		#computing a and u for the given densities
		a, u, z = EulerLoop(n, H_0, a_0, a_1, omega_m, omega_l, u_0)
		#computing the growth factor
		gf = GrowthFactor(u, a)

		#plotting the results
		plt.loglog(z[1:], gf , label = r"$\Omega_m$ = %g , $\Omega_\Lambda$ = %g" %(omega_m, omega_l))

	plt.title("Growth Factor for different cosmologies")
	plt.xlabel("$z$")
	plt.ylabel("$f$")
	plt.legend()
	plt.grid(linestyle="--")
	plt.savefig("growth.pdf")
	plt.show()

main()