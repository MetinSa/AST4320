import numpy as np 
import random as random
import matplotlib.pyplot as plt 

class RandomWalk():

	"""
	Class made to study the random walk process using random numbers
	from Gaussian random distributions. 

	"""

	def __init__(self, eps, initial_var):

		# Initializing variables and constants
		self.eps = eps
		self.k = 2*np.pi*(np.pi/initial_var)**(-1/4)


	def variance(self):

		"""
		Variance for a given smoothing scale Sc.

		"""

		return np.pi/(self.Sc**4)


	def PDF(self, delta, Sc):

		"""
		Analytic expression for the probability density function of
		a Gaussian random field with overdensity delta and smoothed 
		scale Sc.

		"""

		# Definition of the variance
		variance = (np.pi/(Sc**4))

		return (1/(np.sqrt(2*np.pi*variance))) * np.exp(- (delta**2)/(2*variance))


	def PDF_nc(self, delta, Sc):

		"""
		Analytic expression for PDF for delta smaller or equal to delta critical.

		"""

		# Definition of the variance
		variance = (np.pi/(Sc**4))

		# Defining critical value of delta
		delta_crit = 1

		return (1/(np.sqrt(2*np.pi*variance)))*( np.exp( -(delta**2)/(2*variance) ) - np.exp( -((2*delta_crit - delta)**2)/(2*variance)) )



	def resetWalk(self):

		"""
		Function that resets the the smoothing scale and delta after 
		a random walk have been taken.

		"""

		self.Sc = 2*np.pi/self.k
		self.delta = 0


	def walk(self):

		"""
		Calculating a sequence of random steps until the smoothing scale
		Sc = 1.0. This is called a single "realization" of a random walk.

		"""

		# Begin by reseting Sc and delta
		self.resetWalk()

		# Setting iteration conditions
		it = 0;		self.max_it = 1e4

		# Begin realization process
		while self.Sc >= 1 and it < self.max_it:

			# Computing the variance of the current/previous Sc 
			var = self.variance()

			# Reducing Sc by epsilon
			self.Sc -= self.eps

			# Computing the difference between the new and previous variance
			var_diff = self.variance() - var

			# Extracting random normal distributed numbers with above defined var_diff
			self.delta += np.random.normal(scale = np.sqrt(var_diff))

			it += 1

		# Returning the latest calculated delta and Sc = 1
		return self.Sc, self.delta


	def performWalk(self, n_walks):

		"""
		Utilizing the walk function to perform multiple walks, and
		extracting the data.

		"""

		# Initializing delta and Sc arrays
		deltas = np.zeros(n_walks);		Sc = np.zeros(n_walks)
		delta_crit = [];				Sc_crit = []
		# Computing multiple walks and saving the data
		for i in range(n_walks):

			S, d = self.walk()
			Sc[i] = S
			deltas[i] = d

			# Returning only critical threshold values
			if d < 1:

				Sc_crit.append(S)
				delta_crit.append(d)

		# Returning the final results. Sc is an array, but should be approx 1
		return deltas, np.array(delta_crit), Sc, np.array(Sc_crit)


	def distribution(self, n_walks):

		"""
		Function used to study the random walk distribution.

		"""

		# Performing the walks and extracting the data
		deltas, delta_crit , Sc, Sc_crit = self.performWalk(n_walks)

		# Plotting a histogram of the distribution
		plt.hist(deltas, bins=35, normed = True, color = "royalblue", label = "Random walk simulation")

		# Linspace of deltas to be used in analytic expression
		delta_list = np.linspace(np.min(deltas),np.max(deltas),len(deltas))

		# Plotting the analytic PDF for comparison with randomwalk simulation
		plt.plot(delta_list, self.PDF(delta_list, Sc), color = "black", linewidth = 1, label = "Analytic PDF")
		plt.xlabel(r"$\delta$")
		plt.ylabel("Probablity")
		plt.title("Random Walk Simulation and analytic PDF")
		plt.legend(fancybox = False, edgecolor = "black")
		# plt.savefig("randomwalk.pdf")

		plt.show()


	def critDistribution(self, n_walks):

		"""
		Function used to study delta values that never crosses the
		critical threshold of delta_crit = 1.

		"""

		# Performing the walks and extracting the data
		deltas, delta_crit, Sc, Sc_crit = self.performWalk(n_walks)

		# Plotting a histogram of the critical distribution
		plt.hist(delta_crit, bins="auto", normed = True, color = "royalblue", label = "Random walk simulation")

		# Linspace of deltas_crit to be used in analytic expression
		delta_crit_list = np.linspace(np.min(delta_crit),np.max(delta_crit),len(delta_crit))

		# Plotting the analytic PDF for comparison with randomwalk simulation
		plt.plot(delta_crit_list, self.PDF_nc(delta_crit_list, Sc_crit), color = "black", linewidth = 1, label = "Analytic PDF")
		plt.xlabel(r"$\delta$")
		plt.ylabel("Probablity")
		plt.title("Random Walk Simulation and analytic PDF")
		plt.legend(fancybox = False, edgecolor = "black")
		plt.savefig("randomcritwalk.pdf")

		plt.show()



if __name__ == '__main__':

	# Initial condition on the variance
	initial_var = 0.9e-4

	# Epsilon
	eps = 1e-1

	# Number of walks
	n_walks = int(1e4)


	RW = RandomWalk(eps, initial_var)
	# path = RW.distribution(n_walks)
	path = RW.critDistribution(n_walks)








