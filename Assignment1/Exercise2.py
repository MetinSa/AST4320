import numpy as np
import matplotlib.pyplot as plt
import sys as sys


G = 1
rho_0 = 1
H_0 = 1#2.269e-18

def f(y, v, t):
	return - (4/(3*t))*v + (2/(3*t**2))*y


#density fractions
omega_m = np.array([1.0, 0.3, 0.8])
omega_l = np.array([0.0, 0.7, 0.2])



#integration points
n = 10000
dt = 1/n

#defining arrays
t = np.zeros(n)
a = np.zeros(n)
y = np.zeros(n)
v = np.zeros(n)

#initial conditions

t[0] = 2/(3*H_0)
print(t[0])
y[0] = np.exp(-3)
v[0] = 3/2 * H_0
a[0] = np.exp(-3)

#Euler loop
for i in range(0,n-1):
	t[i+1] = t[i] + dt 
	a[i+1] = ((3/2) * H_0*t[i+1])**(2/3)
	y[i+1] = y[i] + v[i]*dt
	v[i+1] = v[i] + f(y[i],v[i],t[i])*dt

# print(a[0],a[-1])
# print(y)
plt.loglog(t,y)
plt.show()
