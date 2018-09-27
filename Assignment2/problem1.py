import numpy as np 
import matplotlib.pyplot as plt 

def WindowFourier(n, k, R):

	"""
	Analytic expression for the Fourier transform of the Window function W(x).

	"""

	return (2/k)*np.sin(k*R)


def FWHM(W, k):

	"""
	Returns the Full Width Half Maximum for a given distribution W.

	"""

	# Finding the half max value of W
	W_half = np.max(W)/2

	# Initializing W_diff
	W_diff = np.zeros(len(W))

	# Finding the index of k which lies closest to the half maximum value
	for i in range(len(W)):

		W_diff[i] = (W[i] - W_half)

	i_min = np.argmin(np.abs(W_diff))

	return np.abs(k[i_min])

"""
==========
problem 1b 
==========
"""

# Grid points
n = 10000

# Smoothing scale
R = 2

# k-array
k = np.linspace(-4*np.pi,4*np.pi,n)

# Analytic expression for W_tilde(k)
W = WindowFourier(n, k, R)

# Extracting the FWHM of W
fwhm = FWHM(W, k)

# Plotting the results
plt.plot(k,W, color = "royalblue", label = "FWHM = %.3f" %fwhm)
plt.title(r"Fourier Conjugate $\tilde{W}(k)$, R = 2")
plt.xlabel("k")
plt.ylabel(r"$\tilde{W}(k)$")
plt.legend()
# plt.savefig("Window.pdf")
plt.show()







