import numpy as np
import matplotlib.pyplot as plt
import sys as sys


G = 1
rho_0 = 1
H_0 = 2.269e-18

def f(u, v, t):
	return (2*u)/(3*t**2) - (4*v)/(3*t)


#density fractions
omega_m = np.array([1.0, 0.3, 0.8])
omega_l = np.array([0.0, 0.7, 0.2])



#integration points
n = 100000
t_stop = 2/(3*H_0)
t_start = 1
dt = (t_stop - t_start)/n

#defining arrays
t = np.zeros(n)
a = np.zeros(n)
u = np.zeros(n)
v = np.zeros(n)

#initial conditions

t[0] = t_start
# print(t[0])
u[0] = np.exp(-3)
v[0] = (2/3)*t[0]**(-2/3)
a[0] = np.exp(-3)

#Euler loop
for i in range(0,n-1):
	t[i+1] = t[i] + dt 
	a[i+1] = ((3/2) * H_0*t[i+1])**(2/3)
	u[i+1] = u[i] + v[i]*dt
	v[i+1] = v[i] + f(u[i],v[i],t[i])*dt

# print(a[0],a[-1])
# print(u)
plt.plot(a,u)
plt.show()
