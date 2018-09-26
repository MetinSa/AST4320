import numpy as np 
import matplotlib.pyplot as plt 

def WindowFourier(n, k, R):

	"""
	Analytic expression for the Fourier transform of the Window function W(x).

	"""

	return (2/k)*np.sin(k*R)


def FWHM(W):

	"""
	Returns the Full Width Half Maximum for a given distribution W.

	"""

	return 2*np.sqrt(2*np.log(2))*np.std(W)


"""
==========
problem 1b 
==========
"""

# Grid points
n = 1000

# Smoothing scale
R = 2

# k-array
k = np.linspace(-4*np.pi,4*np.pi,n)

# Analytic expression for W_tilde(k)
W = WindowFourier(n, k, R)

# Extracting the FWHM of W
fwhm = FWHM(W)

# Plotting the results
plt.plot(k,W, color = "royalblue", label = "FWHM = %.3f" %fwhm)
plt.title(r"Fourier Conjugate $\tilde{W}(k)$, R = 2")
plt.xlabel("k")
plt.ylabel(r"$\tilde{W}(k)$")
plt.legend()
plt.savefig("Window.pdf")
plt.show()







