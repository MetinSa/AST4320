import numpy as np
import matplotlib.pyplot as plt
import sys as sys

#constants
H_0 = 2.269e-18

def f(u, v, t):
	return (2*u)/(3*t**2) - (4*v)/(3*t)

#density fractions
omega_m = np.array([1.0, 0.3, 0.8])
omega_l = np.array([0.0, 0.7, 0.2])

t00 = (1/2*H_0) * (omega_m[2]/((1 - omega_m[2])**(2/3))) * ( (2/omega_m[2]) \
	* np.sqrt(1- omega_m[2]) - (1/np.cosh( (2/omega_m[2]) -1 )))


#integration definitions
n = 100000										#integration points
t_stop = 1/H_0									#end of time
t_start = t_stop*10**(-9/2)					#where our time begins
dt = (t_stop-t_start)/n 						#steplength

#defining arrays
t = np.zeros(n)									#time array
u = np.zeros(n)									#u = delta
v = np.zeros(n)									#v = u' = delta'
a = np.zeros(n)


#initial conditions
t[0] = t_start 									#intial time
u[0] = 10**(-3)								#initial perturbation delta
v[0] = (2/3)*t_stop**(-2/3)*t_start**(-1/3)		#initial change in perturbation 
a[0] = 10**(-3)


# print(t[0], t00)

#forward euler method for integration
for i in range(0,n-1):
	t[i+1] = t[i] + dt 
	a[i+1] = (t[i+1]/t_stop)**(2/3)
	u[i+1] = u[i] + v[i]*dt
	v[i+1] = v[i] + f(u[i],v[i],t[i])*dt


# print(t_stop, t[-1])
#plotting
# plt.loglog(t/t_stop,u, color = "royalblue", label = r"$\delta(t)$")
plt.loglog(a,u, color = "royalblue", label = r"$\delta(t)$")
plt.title("Time evolution of perturbations in EdS")
plt.xlabel("$a$")
plt.ylabel(r"$\delta$")
plt.legend()
plt.grid(linestyle="--")
plt.show()


