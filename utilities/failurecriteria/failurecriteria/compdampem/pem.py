"""
Functionality to use CompDam Fortran code via python extension module
"""


import matplotlib.pyplot as plt
import numpy as np
import os
import pdb
import sys

from ..stress import Stress

path_this_file = os.path.dirname(os.path.abspath(__file__))
pyextmod_path = os.path.abspath(os.path.join(path_this_file, os.pardir, os.pardir, os.pardir, os.pardir, 'pyextmod'))
sys.path.insert(0, pyextmod_path)
import CompDam_DGD

class PemLoader():

	def __init__(self, feature_flags, log_file_name='compdam_out.txt', **kwargs):
		"""
		Loads and initializes the python extension module to access CompDam Fortran code
		accepts material properties and state variables as kwargs
		"""

		# logging
		CompDam_DGD.dgd_mod.log_init(level=4, filename=log_file_name)

		# Load material properties
		m = CompDam_DGD.matprop_mod.loadmatprops('IM7-8552', 40, [
			feature_flags,
			0,
			kwargs.get('thickness', 0.2),
			0,
			0,
			0,
			0,
			0,
			kwargs.get('E1', 171420.),
			kwargs.get('E1', 9080.),
			kwargs.get('G12', 5290.),
			kwargs.get('nu12', 0.32),
			kwargs.get('nu23', 0.45),
			kwargs.get('YT', 62.3),
			kwargs.get('SL', 92.3),
			kwargs.get('GYT', 0.277),
			kwargs.get('GSL', 0.788),
			kwargs.get('eta_BK', 1.634),
			kwargs.get('YC', 199.8),
			kwargs.get('alpha0', 0.925),
			kwargs.get('E3', 8703.),
			kwargs.get('G13', 5164.),
			kwargs.get('G23', 3001.),
			kwargs.get('nu13', 0.32),
			kwargs.get('alpha11', -5.5e-06),
			kwargs.get('alpha22', 2.58e-05),
			kwargs.get('alpha_PL', 4.412e-10),
			kwargs.get('eta_PL', 5.934),
			kwargs.get('XT', 2326.2),
			kwargs.get('fXT', 0.2),
			kwargs.get('GXT', 205.),
			kwargs.get('fGXT', 0.5),
			kwargs.get('XC', 1200.),
			kwargs.get('fXC', 0.2),
			kwargs.get('GXC', 61.),
			kwargs.get('fGXC', 0.5),
			kwargs.get('cl', 10.),
			kwargs.get('w_kb', 0.05),
			0.,
			kwargs.get('mu', 0.3)])
		self.m = m

		# Load parameters from the 'CompDam.parameters' file
		p = CompDam_DGD.parameters_mod.loadparameters()
		self.p = p

		CompDam_DGD.Consistency_Mod.consistencychecks(self.m, self.p, True)
		

		# State variables
		svarray = [
			kwargs.get('CDM_d2', 0.000000000000000E+00),   # SDV1
			kwargs.get('CDM_Fb1', 0.000000000000000E+00),  # SDV2
			kwargs.get('CDM_Fb2', 0.000000000000000E+00),  # SDV3
			kwargs.get('CDM_Fb3', 0.000000000000000E+00),  # SDV4
			kwargs.get('CDM_B', 0.000000000000000E+00),  # SDV5
			kwargs.get('CDM_Lc1', 0.200000000000000E+00),  # SDV6
			kwargs.get('CDM_Lc2', 0.200000000000000E+00),  # SDV7
			kwargs.get('CDM_Lc3', kwargs.get('thickness', 0.2)),  # SDV8
			kwargs.get('CDM_FIm', 0.000000000000000E+00),  # SDV9
			kwargs.get('CDM_alpha', 0.000000000000000E+00),  # SDV10
			kwargs.get('CDM_STATUS', 1.000000000000000E+00),  # SDV11
			kwargs.get('CDM_Plas12', 0.000000000000000E+00),  # SDV12
			kwargs.get('CDM_Inel12', 0.000000000000000E+00),  # SDV13
			kwargs.get('CDM_FIfT', 0.000000000000000E+00),  # SDV14
			kwargs.get('CDM_slide1', 0.000000000000000E+00),  # SDV15
			kwargs.get('CDM_slide2', 0.000000000000000E+00),  # SDV16
			kwargs.get('CDM_FIfC', 0.000000000000000E+00),  # SDV17
			kwargs.get('CDM_d1T', 0.000000000000000E+00),  # SDV18
			kwargs.get('CDM_d1C', 0.000000000000000E+00),  # SDV19
			kwargs.get('CDM_Plas13', 0.000000000000000E+00),  # SDV20
			kwargs.get('CDM_Inel13', 0.000000000000000E+00),  # SDV21
			kwargs.get('CDM_phi0_12', 0.000000000000000E+00),  # SDV22
			kwargs.get('CDM_gamma_12', 0.000000000000000E+00),  # SDV23
			kwargs.get('CDM_phi0_13', 0.000000000000000E+00),  # SDV24
			kwargs.get('CDM_gamma_13', 0.000000000000000E+00),  # SDV25
			0.000000000000000E+00,     # SV26, CDM_reserve
		]
		sv = CompDam_DGD.statevar_mod.loadstatevars(len(svarray), svarray, m)
		self.sv = sv

		# A bit ugly, but necessary to call subroutine with in/out arguments
		position = np.zeros(3)
		phi012 = np.array([self.sv.phi0_12,])
		phi013 = np.array([self.sv.phi0_13,])
		CompDam_DGD.matprop_mod.initializephi0(self.m, self.sv.lc, position, phi012, phi013)
		self.sv.phi0_12 = phi012
		self.sv.phi0_13 = phi013

	def close(self):
		CompDam_DGD.dgd_mod.log_close()


class RambergOsgood():

	def __init__(self, compdam_obj):
		self.compdam_obj = compdam_obj

	def plot(self, shear_component=2, gamma_lims=[-0.15, 0.15], npts=100, show_constant_slope_limit=True, color='k', ax=None, show=True, legend=['ro', 'cs'], derivatives=False):
		'''
		Creates a plot of the Ramberg-Osgood curve
		'''

		if not ax:
			# ax = plt.gca()
			if derivatives:
				fig, (ax, ax2, ax3) = plt.subplots(3)
			else:
				fig, ax = plt.subplots()

		# Create the Ramberg-Osgood curve
		if show_constant_slope_limit:
			# Get gammas for RO (excluding constant slope region)
			gamma_constant_slope = self.compdam_obj.m.shear_strain_const_slope
			tau_constant_slope = CompDam_DGD.rambergosgood_mod.ramberg_osgood(self.compdam_obj.m, gamma_constant_slope, shear_component)
			gamma_ro = np.linspace(max(-1*gamma_constant_slope, gamma_lims[0]), min(gamma_constant_slope, gamma_lims[1]), npts)
			tau_ro = np.array([CompDam_DGD.rambergosgood_mod.ramberg_osgood(self.compdam_obj.m, g, shear_component) for g in gamma_ro])

			# Plot RO curve
			if 'ro' in legend:
				ll = 'Ramberg-Osgood'
			else:
				ll = None
			ax.plot(gamma_ro, tau_ro, color=color, label=ll,)

			# and for constant slope
			gamma_cs_n = np.array([])
			gamma_cs_p = np.array([])

			if gamma_lims[0] < -1*gamma_constant_slope:
				gamma_cs_n = np.linspace(gamma_lims[0], -1*gamma_constant_slope, 3)
				# Marker for the transition
				ax.plot(-1*gamma_constant_slope, -1*tau_constant_slope, 'ok')
				tau_cs_n = np.array([CompDam_DGD.rambergosgood_mod.ramberg_osgood(self.compdam_obj.m, g, shear_component) for g in gamma_cs_n])
				if 'cs' in legend:
					ll = 'Ramberg-Osgood constant slope'
				else:
					ll = None
				ax.plot(gamma_cs_n, tau_cs_n, linestyle='--', color=color, label=ll)
				
				# For output
				gamma_ro = np.concatenate((gamma_cs_n, gamma_ro))
				tau_ro = np.concatenate((tau_cs_n, tau_ro))
			
			if gamma_lims[1] > gamma_constant_slope:
				gamma_cs_p = np.linspace(gamma_constant_slope, gamma_lims[1], 3)
				# Marker for the transition
				ax.plot(gamma_constant_slope, tau_constant_slope, 'ok')
				tau_cs_p = np.array([CompDam_DGD.rambergosgood_mod.ramberg_osgood(self.compdam_obj.m, g, shear_component) for g in gamma_cs_p])
				if gamma_lims[0] < -1*gamma_constant_slope:  # Prevent having two labels in the legend
					label=None
				else:
					if 'cs' in legend:
						ll='Ramberg-Osgood constant slope'
					else:
						ll = None
				ax.plot(gamma_cs_p, tau_cs_p, linestyle='--', color=color, label=ll)

				# For output
				gamma_ro = np.concatenate((gamma_ro, gamma_cs_p))
				tau_ro = np.concatenate((tau_ro, tau_cs_p))
			
		else:
			gamma_ro = np.linspace(gamma_lims[0], gamma_lims[1], npts)
			tau_ro = np.array([CompDam_DGD.rambergosgood_mod.ramberg_osgood(self.compdam_obj.m, g, shear_component) for g in gamma_ro])
			if 'ro' in legend:
				ll = 'Ramberg-Osgood'
			else:
				ll = None
			ax.plot(gamma_ro, tau_ro, color, label=ll)

		if derivatives:
			# First derivative
			gamma = np.linspace(gamma_lims[0], gamma_lims[1], npts)
			tau = np.array([CompDam_DGD.rambergosgood_mod.ramberg_osgood(self.compdam_obj.m, g, shear_component) for g in gamma])
			dtau = np.diff(tau)
			ax2.plot(gamma[1:], dtau)
			# Second derivative
			ddtau = np.diff(dtau)
			ax3.plot(gamma[2:], ddtau)

		if show:
			plt.show()

		return (ax, np.stack((gamma_ro, tau_ro)).T)


class FiberKinking():
	def transform_stress_for_compression(self, stress, psi, phi):
		fi_kink = 0.
		if psi == None:  # Prevents the transformation from being applied twice
			res = self.stress_in_kinking_plane(stress)
			stress = res['stress']
			psi = res['psi']
		if phi == None:
			res = self.stress_in_misaligned_frame(stress)
			stress = res['stress']
			phi = res['phi']
			fi_kink = res['fi_kink']
			# splitting_state = res['splitting']
		# return (stress, psi, phi, fi_kink, splitting_state)
		return (stress, psi, phi, fi_kink)


	def stress_in_kinking_plane(self, stress, psi=None):
		# Kinking plane
		psi_orig = psi
		if psi == None:
			if (stress.sig22-stress.sig33) == 0.:
				psi = 0.
			else:
				psi = 0.5*np.arctan(2*stress.sig23/(stress.sig22-stress.sig33))

		if self.verbose > 2:
			print('psi: {:.3f} [deg]'.format(psi*180./np.pi))

		sig_psi_22 = (stress.sig22+stress.sig33)/2 + (stress.sig22-stress.sig33)/2*np.cos(2*psi) + stress.sig23*np.sin(2*psi)
		sig_psi_33 = stress.sig22 + stress.sig33 - sig_psi_22
		sig_psi_12 = stress.sig12*np.cos(psi) + stress.sig13*np.sin(psi)
		sig_psi_13 = stress.sig12*np.sin(psi) + stress.sig13*np.cos(psi)
		s = Stress.fromVec([stress.sig11, sig_psi_22, sig_psi_33, sig_psi_12, 0., sig_psi_13])
		if psi_orig == None:
			return {'stress': s, 'psi': psi}
		else:
			return s


	def stress_in_misaligned_frame(self, stress, phi=None):
		'''
		By default, computes the angle phi
		'''
		phi_orig = phi
		if phi == None:
			# Check for fiber kinking
			if CompDam_DGD.kinking_mod.is_kinking_possible(self.pem.m, stress.mat, self.pem.sv.phi0_12, shear_component=2):
				gammac, sigc = CompDam_DGD.kinking_mod.kinking_stress(self.pem.m, stress.mat, self.pem.sv.phi0_12, shear_component=2)
				fi_kink = stress.sig11/sigc
			else:
				fi_kink = -1.
			# Splitting - need to rotate stress to the misaligned frame
			(phi, found_zero) = self.current_misalignment(stress)
			# print('current_misalignment() found zero = {}'.format(found_zero))

		# Rotate stress
		sig_m_11 = (stress.sig11+stress.sig22)/2 + (stress.sig11-stress.sig22)/2*np.cos(2*phi) + stress.sig12*np.sin(2*phi)
		sig_m_22 = stress.sig11 + stress.sig22 - sig_m_11
		sig_m_12 = -(stress.sig11-stress.sig22)/2*np.sin(2*phi) + stress.sig12*np.cos(2*phi)
		sig_m_13 = stress.sig13*np.cos(phi)
		sig_m_23 = -stress.sig13*np.sin(phi)

		# TODO: check for fiber crushing failurei

		s_m = Stress.fromVec([sig_m_11, sig_m_22, stress.sig33, sig_m_12, sig_m_23, sig_m_13])
		if phi_orig == None:
			return {'stress': s_m, 'phi': phi, 'fi_kink': fi_kink}
		else:
			return s_m


	def current_misalignment(self, stress, n_max=100, tol=1e-5):
		'''
		For nonlinear
		Solves f1 = 0 to find gamma1m2m using Newton-Raphson (finds minimum if no 0 is found)

		Returns (phi, found_zero) where the bool found_zero indicates if a zero was found
		if no zero was found, failure has occurred by kinking
		'''

		# Initialize
		ff0_best = 1e5  # Large number
		gamma1m2m_best = 0.
		found_zero = False
		found_min = False

		if stress.sig11 > 0:
			print('WARNING: atttempting to find the current misalignment angle for longitudinal tension; may not converge')

		# Attempt to find f1 = 0 (intersection between driving force and shear law, linearized)
		# gamma1m2m = 0.0001  # Initial guess
		gamma1m2m = (self.pem.sv.phi0_12*(stress.sig22-stress.sig11)+stress.sig12)/(self.pem.m.g12-stress.sig22+stress.sig11)  # Initial guess
		# print('Initial guesss gamma1m2m', gamma1m2m)
		# print('current_misalignment initial guess gamma1m2m = {}'.format(gamma1m2m))
		tol = tol*self.pem.m.sl  # tolerance for convergence [shear stress]
		history_ff0 = []
		history_gamma = []
		not_coverging = False
		# print('Attempting to find current_misalignment')

		# Debug ff0 and ff1
		# fig, (ax1, ax2) = plt.subplots(2)
		# gamma = np.linspace(-0.05, 0.05, 100)
		# ff0 = np.zeros(gamma.shape)
		# ff1 = np.zeros(gamma.shape)
		# for i, g in enumerate(gamma):
		# 	resjac = CompDam_DGD.kinking_mod.kinking_stress_resjac(self.pem.m, stress.mat, self.pem.sv.phi0_12, shear_component=2, gamma=g)
		# 	ff0[i] = resjac[0,0]
		# 	ff1[i] = resjac[0,1]
		# ax1.plot(gamma, ff0)
		# ax2.plot(gamma, ff1)
		# plt.show()
		# raise ''



		for n in np.arange(0,n_max+1):
			resjac = CompDam_DGD.kinking_mod.kinking_stress_resjac(self.pem.m, stress.mat, self.pem.sv.phi0_12, shear_component=2, gamma=gamma1m2m)
			ff0 = resjac[0,0]
			ff1 = resjac[0,1]

			history_gamma.append(gamma1m2m)
			history_ff0.append(ff0)
			# print(gamma1m2m, ff0)

			if np.abs(ff0) <= tol:  # Absolute tolerance
				found_zero = True
				break
			if np.abs(ff0) < np.abs(ff0_best):
				gamma1m2m_best = gamma1m2m
				ff0_best = ff0
			elif n > 3:
				not_coverging = True
				# print('current_misalignment() is not converging')
				break
			gamma1m2m = gamma1m2m - ff0/ff1  # Newton-Raphson equation
			# gamma1m2m = gamma1m2m - ff0/ff1  # Newton-Raphson equation
		# print(history_gamma)
		# print(history_ff0)
		# print(ff0_best, gamma1m2m_best)

		# Check for valid solution
		valid_soln = True

		# Direction of shear strain
		above_left = stress.sig12 > CompDam_DGD.rambergosgood_mod.ramberg_osgood(self.pem.m, -self.pem.sv.phi0_12, shear_component=2)
		if above_left:
			# print('Above/left')
			# Above/left of RO -- gamma1m2m must be positive
			if gamma1m2m < -self.pem.sv.phi0_12:
				valid_soln = False
				# print('Converged solution is invalid: shear strain is in wrong direction')
			if gamma1m2m > self.pem.m.shear_strain_const_slope:
				valid_soln = False
				# print('Converged solution is invalid: gamma1m2m - self.pem.sv.phi0_12 > self.pem.m.shear_strain_const_slope')
		else:
			# Below/right of RO -- gamma1m2m must be negative
			# print('Below/right')
			if gamma1m2m > -self.pem.sv.phi0_12:
				valid_soln = False
				# print('Converged solution is invalid: shear strain is in wrong direction')
			if gamma1m2m < -self.pem.m.shear_strain_const_slope:
				valid_soln = False
				# print('Converged solution is invalid: gamma1m2m < -self.pem.m.shear_strain_const_slope')

			
		if n == n_max or not_coverging or not valid_soln:
			# Failed to find f1=0
			# Try to find a minimum, ie df/dg = 0
			if above_left:
				gamma1m2m = np.linspace(-self.pem.sv.phi0_12, self.pem.m.shear_strain_const_slope, n_max*2)
			else:
				gamma1m2m = np.linspace(-self.pem.m.shear_strain_const_slope, -self.pem.sv.phi0_12, n_max*2)

			ff0 = np.array([CompDam_DGD.kinking_mod.kinking_stress_resjac(self.pem.m, stress.mat, self.pem.sv.phi0_12, shear_component=2, gamma=g)[0,0] for g in gamma1m2m])
			gamma1m2m = gamma1m2m[np.argmin(np.abs(ff0))]

		# Current misalignment
		return (self.pem.sv.phi0_12 + gamma1m2m, found_zero)


	def considere(self, stresses, gamma_lims=[-0.15, 0.15], npts=100, phi0=None, shear_component=2, ax=None, show=False, plot_driving_force=True,
		legend=['ramberg_osgood', 'shear_strength', 'zero', 'kinking_critical', 'current_misalignment'], kinking_color='g', debugging_plots=False):
		'''
		Creates a plot of the considere construction
		stress is a stress object or list of stress objects
		'''

		if not isinstance(stresses, list):
			stresses = [stresses, ]

		if phi0:
			raise Exception('not implemented')
		else:
			phi0 = self.pem.sv.phi0_12
		print('Creating Considere construction for phi0 = {:.2f} deg ({:.3f} rad)'.format(phi0*180/np.pi, phi0))
		
		# Plot the Ramberg-Osgood curve
		ro = RambergOsgood(self.pem)
		(ax, ro_data) = ro.plot(shear_component=shear_component, gamma_lims=gamma_lims, npts=npts, show_constant_slope_limit=True, color='k', ax=ax, show=False, legend=[])

		# Add the shear strength to the plot
		strain_at_sl = np.interp(self.pem.m.sl, ro_data[:,1], ro_data[:,0])
		if 'shear_strength' in legend:
			ll = 'Shear strength'
		else:
			ll = None
		ax.plot(strain_at_sl, self.pem.m.sl, 'vk', label=ll)
		ax.plot(-1*strain_at_sl, -1*self.pem.m.sl, 'vk')

		# Draw curves representing the stress state
		gamma = np.linspace(gamma_lims[0], gamma_lims[1], npts)
		for i, stress in enumerate(stresses):
			# if self.verbose:
				# print('stress state:\n{}'.format(stress))
			if plot_driving_force:
				consid = np.array([CompDam_DGD.kinking_mod.kinking_driving_force(stress.mat, g, phi0, shear_component, order=0) for g in gamma])
				ax.plot(gamma, consid, '-r', label='Stress: '+stress.str_oneline)

			# Plot phi0
			if 'zero' in legend:
				ll = 'zero stress point (phi0, sig12)'
			else:
				ll = None
			ax.plot(-1*phi0, stress.sig12, 'xk', label=ll)

			# Check if kinking can occur at this (sig12, sig22)
			if CompDam_DGD.kinking_mod.is_kinking_possible(self.pem.m, stress.mat, phi0, shear_component):
				print('Kinking is possible')
				# Find kinking onset
				# stress_kk = stress.copy()
				# stress_kk.sig11 = -248.06921857627967
				# print(CompDam_DGD.kinking_mod.kinking_stress_resjac(self.pem.m, stress_kk.mat, phi0, shear_component, gamma=-0.1, debug=True))
				# raise ''
				gammac, sigc = CompDam_DGD.kinking_mod.kinking_stress(self.pem.m, stress.mat, phi0, shear_component)
				print('Kinking stress: {:.2f}, shear strain: {:.4f}'.format(sigc, gammac))
				stress_kinking = stress.copy()
				stress_kinking.sig11 = sigc
				# Plot the critical case
				consid = np.array([CompDam_DGD.kinking_mod.kinking_driving_force(stress_kinking.mat, g, phi0, shear_component, order=0) for g in gamma])
				ax.plot(gamma, consid, linestyle='-', color=kinking_color, label='Kinking: {}, phi0: {:.1f} deg'.format(stress_kinking.str_oneline, phi0*180/np.pi))
				tau12c = CompDam_DGD.kinking_mod.kinking_driving_force(stress_kinking.mat, gammac, phi0, shear_component, order=0)
				if 'kinking_critical' in legend:
					ll = 'Kinking critical point'
				else:
					ll = None
				ax.plot(gammac, tau12c, marker='s', color=kinking_color, label=ll)
			else:
				print('Kinking is not possible')

			if debugging_plots:
				# Check alternate starting points
				candidate_gammas = [
				 -1.0000000000000000E-002,
				 -2.0000000000000000E-002,
				 -3.0000000000000006E-002,
				 -4.0000000000000001E-002,
				 -5.0000000000000003E-002,
				 -6.0000000000000012E-002,
				 -7.0000000000000007E-002,
				 -8.0000000000000002E-002,
				 -8.9999999999999997E-002,
				 -0.10000000000000001
				]
				candidate_sigcs = [
				 -3429.8524775380861,
				 -1274.0454116275241,
				 -749.47666033094151,
				 -532.18023288459744,
				 -413.46970363003936,
				 -338.20051479796939,
				 -285.85567870884961,
				 -247.08830477511029,
				 -217.03313115507896,
				 -192.90664854450452    
				]

				for i in range(len(candidate_gammas)):
					gs = np.linspace(-0.2, 0.2, 100)
					stress_temp = stress.copy()
					stress_temp.sig11 = candidate_sigcs[i]
					driving_force = np.array([CompDam_DGD.kinking_mod.kinking_driving_force(stress_temp.mat, g, phi0, shear_component, order=0) for g in gamma])
					# Driving force curve for the current candidate gamma
					ax_obj = ax.plot(gs, driving_force)
					# Starting point
					color = ax_obj[-1].get_color()
					ax.plot(candidate_gammas[i], CompDam_DGD.kinking_mod.kinking_driving_force(stress_temp.mat, candidate_gammas[i], phi0, shear_component, order=0), 'x', color=color)

				# Plot error
				# eqn20 = np.zeros(gamma.shape)
				# eqn21 = np.zeros(gamma.shape)
				# residual = np.zeros(gamma.shape)
				# for i, g in enumerate(gamma):
				# 	tau = self.sn_obj.sn_fun(g)
				# 	# tau = np.sign(g)*tau
				# 	# eqn20[i] = self.kinking_driving_force(stress, g)
				# 	# eqn20[i] = tau
				# 	eqn20[i] = -0.5*(stress.sig11-stress.sig22)*np.sin(2*(self.phi0_nl+g)) + stress.sig12*np.cos(2*(self.phi0_nl+g)) - tau
				# 	eqn21[i] = -(stress.sig11-stress.sig22)*np.cos(2*(self.phi0_nl+g)) - 2*stress.sig12*np.sin(2*(self.phi0_nl+g)) - self.sn_obj.d_sn_fun(tau, order=1)
				# 	residual[i] = np.linalg.norm([eqn20[i], eqn21[i]])
				# fig1, axs = plt.subplots(3)
				# axs[0].plot(gamma, eqn20)
				# axs[1].plot(gamma, eqn21)
				# axs[2].plot(gamma, residual)
				# print('Mininum error: {}'.format(np.min(residual)))


			# Check if splitting can occur at this stress state
			(phi, found_zero) = self.current_misalignment(stress)
			gamma1m2m = phi-self.pem.sv.phi0_12
			if 'current_misalignment' in legend:
				ll = 'Current misalignment point'
			else:
				ll = None
			ax.plot(gamma1m2m, CompDam_DGD.kinking_mod.kinking_driving_force(stress.mat, gamma1m2m, phi0, shear_component, order=0), marker='^', color=kinking_color, label=ll)
			if not ll:
				consid = np.array([CompDam_DGD.kinking_mod.kinking_driving_force(stress.mat, g, phi0, shear_component, order=0) for g in gamma])
				ax.plot(gamma, consid, linestyle='--', color=kinking_color)
			print('phi: {}, gamma1m2m: {}'.format(phi, gamma1m2m))

			# stress_m = self.stress_in_misaligned_frame(stress, phi=phi)
			# print('misaligned stress:')
			# print(stress_m)
			# print(self.FI_num4(stress_m, psi=0., phi=phi, return_dict=True))
			# # print('double check failure index: ', (stress_m.sig12/(self.SL-self.etaL*stress_m.sig22)))
			# apparent_sl = self.SL - self.etaL*stress_m.sig22
			# ax.plot(self.sn_obj.sn_fun_inv(apparent_sl), apparent_sl, '+c', label='Apparent shear strength')



			# print('Running Lcf')
			# myLcf = Lcf(G12=self.G12, alpha=self.alpha, eta=self.eta, XC=self.XC, SL=self.SL, YC=self.YC, alpha0=self.alpha0)
			# x0 = [gamma1m2m, stress.sig11]
			# print('x0 ', x0)
			# print('input stress:\n', stress)
			# print(myLcf.sigc_splitting_res(x0, stress.sig22, stress.sig12))
			# if self.is_splitting_possible(stress):
			# 	print('Splitting is possible')
			# 	phi = self.current_misalignment(stress)
			# 	gamma = phi - phi0
			# 	plt.plot(gamma, self.sn_obj.sn_fun(gamma), '^k', label='Stress state for splitting')

			# if self.verbose > 1:
			# 	self.FI(stress)

		if show:
			plt.show()



