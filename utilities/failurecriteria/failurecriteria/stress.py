import numpy as np

class Stress():
	'''
	Class to store stress states
	'''

	def __init__(self, sig11=0., sig22=0., sig33=0., sig12=0., sig13=0., sig23=0.):
		self.sig11=sig11
		self.sig22=sig22
		self.sig33=sig33
		self.sig12=sig12
		self.sig13=sig13
		self.sig23=sig23


	@classmethod
	def fromVec(cls, stress_vec):
		'''
		Create stress object from a 6 component vector
		'''
		return cls(stress_vec[0], stress_vec[1], stress_vec[2], stress_vec[3], stress_vec[4], stress_vec[5])

	def __repr__(self):
		return str(self.mat)

	@property
	def mat(self):
		return np.array([[self.sig11, self.sig12, self.sig13], [self.sig12, self.sig22, self.sig23], [self.sig13, self.sig23, self.sig33]])
	
	@property
	def vec(self):
		return np.array([self.sig11, self.sig22, self.sig33, self.sig12, self.sig23, self.sig13])

	@property	
	def str_oneline(self):
		stress_dict = {
		'sig11': self.vec[0], 
		'sig22': self.vec[1], 
		'sig33': self.vec[2], 
		'sig12': self.vec[3], 
		'sig23': self.vec[4], 
		'sig13': self.vec[5]
		}
		return ', '.join(['{}={:.1f}'.format(k,v) for k,v in stress_dict.items() if v != 0.])


	def copy(self):
		return Stress.fromVec(self.vec)