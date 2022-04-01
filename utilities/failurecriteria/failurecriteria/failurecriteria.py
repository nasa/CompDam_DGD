import os
import numpy as np
import sys
import matplotlib.pyplot as plt

from .stress import Stress

class FailureCriteria():
	'''
	Base class for functionality common to all failure criteria (find failure envelop, plot, etc)
	Expects that child class creates method:
	FI(sig11, sig22, sig12, return_dict=False)  #Other stress components should be optional arguments
	'''

	def __init__(self, mat, verbose=False):
		# Verbose can be Bool or integer
		self.verbose = int(verbose)
		# Load in the material properties
		for k, v in mat.items():
			setattr(self, k, float(v))
			if self.verbose > 1:
				print('Set {0} = {1}'.format(k, float(v)))
		# Make sure alpha0 is radians
		if self.alpha0 > np.pi:
			if self.verbose:
				print('Found alpha0 = {0}. Assuming the value is degrees, converting to radians.'.format(self.alpha0))
			self.alpha0 = self.alpha0*np.pi/180.


	def failure_envelop_generator(self, x_anchors, y_anchors, theta, failure_index_fun, npts, stress_ratios=[], tol=1e-3, max_iterations=50):
		'''
		Uses the bisection method to attempt to find the failure envelop
		x_anchor: 2 element vector with -strength and +strength on x-axis; use -1 or 1 for unused quadrants
		y_anchor: 2 element vector with -strength and +strength on y-axis; use -1 or 1 for unused quadrants
		theta: 2 element vector with min_theta and max_theta in radians, measured ccw from +x axis specifying the window of failure envelop to calculate
		failure_index_fun: function to compute the failure index; assumes it takes two arguments: x,y and returns the failure index
		npts: number of points per quadrant if integer, otherwise 4 element vector with number of points specified for each quadrant

		This function uses the bisection method to iteratively find the failure envelop 
		(failure_index_fun = 0). The calculation is performed separately for each quadrant using 
		starting points that assume a linear interaction between x and y. A simple cutback
		strategy is included.
		'''

		def write_nonconvergence_data():
			print('WARNING: Stress ratio: {0} of {1} did not converge (reached max_iterations)'.format(ii, number_pts_in_quadrant))
			if True:
				file_name = 'not_converged_{}.csv'.format(ii)
				print('Writing history variables to {0}'.format(file_name))
				data_to_write = np.hstack((sigx.reshape(max_iterations,1), sigy.reshape(max_iterations,1), fi.reshape(max_iterations,1), aid.reshape(max_iterations,1)))
				np.savetxt(file_name, data_to_write, delimiter=',', header='sigx,sigy,fi,aid')

		# Error if theta out of range
		if max(theta) > 2*np.pi:
			raise ValueError('Expecting 0 <= theta < 2pi')
		if min(theta) < 0:
			raise ValueError('Expecting 0 <= theta < 2pi')
		if len(theta) != 2:
			raise ValueError('Expecting theta to have two items [min_theta, max_theta]')
		if theta[1] <= theta[0]:
			raise ValueError('Expecting theta to have two items [min_theta, max_theta]')

		# Number of quadrants
		start_quadrant = int(np.floor(theta[0]/(np.pi/2.)))+1
		end_quadrant = int(np.ceil(theta[1]/(np.pi/2.)))
		quadrants = tuple(range(start_quadrant, end_quadrant+1))
		
		# For output
		x_values = []
		y_values = []

		if isinstance(npts, int):
			npts = [npts]*4

		if self.verbose:
			print('Creating failure envelop from {0:.1f} deg to {1:.1f} deg. Quadrants: {2}. Number points: {3}'.format(theta[0]*180/np.pi, theta[1]*180/np.pi, quadrants, (np.sum(npts))))
		else:
			print('Creating failure envelop using {0} points'.format(np.sum(npts)))

		# Generate failure envelop for each quadrant separately
		for q in quadrants:

			if self.verbose > 1:
				print('Quadrant: {0}'.format(q))

			# For each stress ratio
			if stress_ratios:
				n = np.array(stress_ratios)
			else:
				n = np.linspace(0,1,npts[q-1])
			number_pts_in_quadrant = len(n)
			for ii in range(number_pts_in_quadrant):

				if self.verbose:
					print('Stress ratio = {} ({} of {})'.format(n[ii], ii, number_pts_in_quadrant))

				# Computed stress values (keep results for all iterations)
				sigx = np.zeros(max_iterations)
				sigy = np.zeros(max_iterations)

				# Start points
				if q == 1:
					sigx[0] = (1-n[ii])*x_anchors[1]
					sigy[0] = n[ii]*y_anchors[1]
				elif q == 2: # Count backwards
					sigx[0] = (1-n[-(ii+1)])*x_anchors[0]
					sigy[0] = n[-(ii+1)]*y_anchors[1]
				elif q == 3:
					sigx[0] = (1-n[ii])*x_anchors[0]
					sigy[0] = n[ii]*y_anchors[0]
				else: # Count backwards
					sigx[0] = (1-n[-(ii+1)])*x_anchors[1]
					sigy[0] = n[-(ii+1)]*y_anchors[0]
				# Scale down to start within failure envelop
				sigx[0] *= 0.5
				sigy[0] *= 0.5

				if self.verbose > 1:
					print('Start stress {}, {}'.format(sigx[0], sigy[0]))
			
				stress_dir = np.array([sigx[0], sigy[0]])/np.sqrt(sigx[0]**2+sigy[0]**2)

				# iteratively find FI=1
				i = 0
				fi = np.zeros(max_iterations)
				converged = False
				aid = np.ones(max_iterations)
				cutback = 0
				last_cutback_i = -1
				while True:
					fi[i] = failure_index_fun(sigx[i], sigy[i])
					if self.verbose > 2:
						print('Failure index: {0}'.format(fi[i]))

					# check for convergence
					if np.abs(fi[i]-1.) < tol:
						if self.verbose > 1:
							print('Converged at stress: x={0:.3f}, y={1:.3f}'.format(sigx[i], sigy[i]))
						converged = True
						break

					# Stop if not converged after max_iterations
					if i+1 >= max_iterations:
						write_nonconvergence_data()
						break

					# Check for diverging solution
					if i > 2 and i > last_cutback_i + 1:
						if np.abs(fi[i-1]-1) < np.abs(fi[i]-1):
							aid[i:] = aid[i-1]/2.							
							if cutback > 3:
								write_nonconvergence_data()
								print('---- write_nonconvergence_data()')
								break
							cutback += 1
							last_cutback_i = i
							if self.verbose > 1:
								print('Cutting back at i={}, new aid={}'.format(i,aid[i]))
							# if fi[i] == -1.:
							sigx[i] = sigx[i-2]
							sigy[i] = sigy[i-2]
							continue


					# New trial stress
					mag_last = np.sqrt(sigx[i]**2+sigy[i]**2)
					if fi[i] > 1.5:
						mag_next = 0.5*mag_last
					else:
						scale_factor = ((1-fi[i])/2.*aid[i]+1)
						mag_next = scale_factor*mag_last
						if self.verbose > 2:
							print('scale_factor: {0}'.format(scale_factor))

					sig_new = mag_next*stress_dir
					sigx[i+1] = sig_new[0]
					sigy[i+1] = sig_new[1]
					# print('sig_new: {0}'.format(sig_new))

					# increment
					i = i+1

				# Keep converged result
				if converged:
					x_values.append(sigx[i])
					y_values.append(sigy[i])

		return (x_values, y_values)


	def plot_by_criterion(self, x, y, failure_index_fun_returns_dict, base_label='', linestyle=None, color=None, marker=None, single_label=False):
		'''
		Creates a plot of the failure envelop
		x,y are the stress vectors that define the failure envelop
		'''

		grouped_ss = {}
		for i in range(0, len(x)):
			fi_dict = failure_index_fun_returns_dict(x[i], y[i])
			# print(x[i],y[i],fi_dict)
			criterion_num = fi_dict['mode']
			if criterion_num not in grouped_ss.keys():
				grouped_ss[criterion_num] = []
			grouped_ss[criterion_num].append((x[i], y[i]))

		for i, criterion_num in enumerate(grouped_ss.keys()):
			xy = np.array(grouped_ss[criterion_num])
			if i==0 and single_label:
				ll = base_label
			elif i>0 and single_label:
				ll = None
			else:
				ll = '{0} #{1}'.format(base_label,criterion_num)
			if i>0 and single_label and color==None:
				color = p[0].get_color()
			p = plt.plot(xy[:,0], xy[:,1], linestyle=linestyle, color=color, marker=marker, label=ll)


	def failure_envelop_1122(self, npts, theta=[0, 2*np.pi], criteria=[], plot=True, savetxt=None, base_label='', linestyle=None, color=None, marker=None, single_label=False):
		'''
		Implementation specific to 11 vs 22 failure envelop
		'''
		self.failure_envelop_plane('1122', npts, theta, criteria, plot, savetxt, base_label, linestyle, color, marker, single_label)


	def failure_envelop_2212(self, npts, theta=[0, np.pi], criteria=[], plot=True, savetxt=None, base_label='', linestyle=None, color=None, marker=None, single_label=False):
		'''
		Implementation specific to 11 vs 22 failure envelop
		'''
		self.failure_envelop_plane('2212', npts, theta, criteria, plot, savetxt, base_label, linestyle, color, marker, single_label)


	def failure_envelop_1112(self, npts, theta=[np.pi/2, 3*np.pi/2], criteria=[], savetxt=None, plot=True, base_label='', linestyle=None, color=None, marker=None, single_label=False):
		'''
		Implementation specific to 11 vs 12 failure envelop
		'''
		self.failure_envelop_plane('1112', npts, theta, criteria, plot, savetxt, base_label, linestyle, color, marker, single_label)


	def failure_envelop_plane(self, plane, npts, theta, criteria, plot, savetxt, base_label, linestyle, color, marker, single_label):
		'''
		Shared code among all failure envelop functions
		'''

		# Get the failure index function(s)
		fi_funs = self._get_failure_index_funs(criteria)

		# Loop through each criterion
		for fi_fun in fi_funs:
		
			# Items specific to the plane being evaluated
			if plane == '1122':
				fc = lambda x, y: fi_fun(Stress(sig11=x, sig22=y))
				x_anchors = [-self.XC, self.XT]
				y_anchors = [-self.YC, self.YT]
			elif plane == '2212':
				fc = lambda x, y: fi_fun(Stress(sig22=x, sig12=y))
				x_anchors = [-self.YC, self.YT]
				y_anchors = [-1, self.SL]
			elif plane == '1112':
				fc = lambda x, y: fi_fun(Stress(sig11=x, sig12=y))
				x_anchors = [-self.XC, self.XT]
				y_anchors = [-self.SL, self.SL]
			else:
				raise Exception('Not implemented')
			

			# Generate the failure envelop stresses and strains
			# (x, y) = self.failure_envelop_generator(x_anchors, y_anchors, theta=theta, failure_index_fun=fc, npts=npts, stress_ratios=[0.98])
			(x, y) = self.failure_envelop_generator(x_anchors, y_anchors, theta=theta, failure_index_fun=fc, npts=npts)

			# Report the results
			if plane == '1122':
				fc = lambda x, y: fi_fun(Stress(sig11=x, sig22=y), return_dict=True)
			elif plane == '2212':
				fc = lambda x, y: fi_fun(Stress(sig22=x, sig12=y), return_dict=True)
			elif plane == '1112':
				fc = lambda x, y: fi_fun(Stress(sig11=x, sig12=y), return_dict=True)
			if plot:
				self.plot_by_criterion(x, y, fc, base_label, linestyle=linestyle, color=color, marker=marker, single_label=single_label)
			if savetxt:
				with open(savetxt, 'w') as outfile:
					for i in range(0, len(x)):
						fi_dict = fc(x[i], y[i])
						outfile.write('{},{},{}\n'.format(x[i], y[i], fi_dict))


	def _get_failure_index_funs(self, criteria):
		'''
		Gets a list of function handles for the given list of criteria, where criteria are num1, num2, etc
		'''
		if criteria:
			fi_funs = []
			for c in criteria:
				criterion_fun_name = 'FI_'+c
				fi_funs.append(getattr(self, criterion_fun_name))
		else:
			fi_funs = [getattr(self, 'FI'),]
		return fi_funs