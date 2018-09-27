import numpy as np 
import random as random
import matplotlib.pyplot as plt 

class RandomWalk():

	"""
	Class that calculate a random walk.

	"""

	def __init__(self, eps, initial_var):

		#initializing variables
		self.eps = eps
		self.k = 2*np.pi*(np.pi/initial_var)**(-1/4)


	def variance(self):

		#Computing the assumed variance.
		return np.pi/(self.Sc**4)

	def newSc(self):

		#Subtracts epilon from current Sc
		self.Sc -= self.eps


	def resetWalk(self):

		#Resetting the process after a random walk.
		self.Sc = 2*np.pi/self.k
		self.var = self.variance()
		self.path = []
		self.delta = 0


	def walk(self):

		"""
		Calculating a single "realization" of a random walk.

		"""
		self.resetWalk()
		it = 0
		self.max_it = 1e5

		while self.Sc > 1 and it < self.max_it:

			old_var = self.variance()

			self.newSc()
			new_var = self.variance() - old_var
			self.delta = np.random.normal(scale = np.sqrt(new_var))
			self.path.append([self.Sc, self.delta])

			it += 1


		return np.array(self.path)


	def performWalk(self, n_walks):

		"""
		Performing multiple walks.

		"""

		deltas = np.zeros(n_walks)

		for i in range(n_walks):
			self.resetWalk()
			path = self.walk()
			deltas[i] = path[-1,1]

		return deltas, path[-1,0]


	def distribution(self,n_walks):
		delta , Sc = self.performWalk(n_walks)

		cont_delta = np.linspace(np.min(delta),np.max(delta),1000)

		plt.hist(delta,bins="auto", normed = True)
		plt.plot(cont_delta, self.PDF(cont_delta,Sc))
		plt.show()


	def PDF(self, delta, Sc):

		var = (np.pi/(Sc**4))

		return (1/(np.sqrt(2*np.pi)*var)) * np.exp(- (delta**2)/(2*var**2))


initial_var = 0.9e-4
eps = 1e-1
n_walks = int(1e5)
RW = RandomWalk(eps, initial_var)
path = RW.distribution(n_walks)

# print(path)
# print(path[0,-1])






