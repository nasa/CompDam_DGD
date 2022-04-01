import numpy as np
import pdb

from .compdampem.pem import RambergOsgood, FiberKinking, PemLoader
from .failurecriteria import FailureCriteria


class CompDam(FailureCriteria, FiberKinking):
	
	def __init__(self, mat, feature_flags, verbose=False, **kwargs):
		self.verbose = verbose

		self.pem = PemLoader(feature_flags, **mat, **kwargs)

		# Make strengths available as attributes for failure envelop code
		self.XT = self.pem.m.xt
		self.XC = self.pem.m.xc
		self.YT = self.pem.m.yt
		self.YC = self.pem.m.yc
		self.SL = self.pem.m.sl
		self.ST = self.pem.m.st

		# Fiber break
		self.fiber_diameter = 0.0005
		self.w_kb = 0.1
		self.fiber_rupture_strain = 0.02  # for reference: sig11 = 2571 MPa

		# if self.verbose:
		# 	print("CompDam; Calculated parameters:")
		# 	print("\tetaL: {0:.3f}".format(self.etaL))
		# 	print("\tetaT: {0:.3f}".format(self.etaT))
		# 	print("\tST: {0:.1f}".format(self.ST))
		# 	print("\tphi0 linear: {0:.2f} [deg], {1:.3f} [rad]".format(self.phi0_linear*180/np.pi, self.phi0_linear))
		# 	print("\tphi0_nl: {0:.2f} [deg], {1:.3f} [rad]".format(self.phi0_nl*180/np.pi, self.phi0_nl))


	def close(self):
		self.pem.close()

	def FI(self, stress, psi=None, alphas=np.arange(0,91,1), return_dict=False):

		calculated_failure_indices = []

		if stress.sig11 <= 0.:
			# Longitudinal compression
			(stress_m, psi, phi, fi_kink) = self.transform_stress_for_compression(stress, psi, phi=None)
			if fi_kink > 0.:
				calculated_failure_indices.append({'fi': fi_kink, 'mode': 'k', 'phi': phi, 'psi': psi})
			# if self.verbose:
			# 	print('phi: {:.3f} [deg]'.format(phi*180/np.pi))
			# 	print('gamma: {:.2f} [%]'.format((phi-self.pem.sv.phi0_12)*100))

			# Fiber break
			FIfb = self.FI_fiberbreak(stress_m, psi, phi, return_dict=True)
			calculated_failure_indices.append(FIfb)

			# Matrix cracking
			FIm = self.FI_matrix(stress_m, psi=psi, phi=phi, return_dict=True)
			calculated_failure_indices.append(FIm)
			# Re-calculate for phi0 = 0
			FIm = self.FI_matrix(stress, psi=psi, phi=0., return_dict=True)
			FIm['mode'] += '_phi=0'
			calculated_failure_indices.append(FIm)

		else:
			FIm = self.FI_matrix(stress, return_dict=True)
			calculated_failure_indices.append(FIm)
		# print(calculated_failure_indices)
		# Return the maximum failure criteria
		# print(calculated_failure_indices)
		fis = [x['fi'] for x in calculated_failure_indices]
		if return_dict:
			ii = fis.index(max(fis))
			return calculated_failure_indices[ii]
		else:
			return max(fis)


	def FI_matrix(self, stress, alphas=np.arange(0,91,1), psi=None, phi=None, return_dict=False):
		"""
		Todo: update to use PEM
		"""

		E2 = self.pem.m.e2
		YT = self.pem.m.yt
		SL = self.pem.m.sl
		ST = self.pem.m.st
		GYT = self.pem.m.gyt
		GSL = self.pem.m.gsl
		etaT = self.pem.m.etat
		etaL = self.pem.m.etal
		eta_BK = self.pem.m.eta_bk

		(stress, psi, phi, fi_kink) = self.transform_stress_for_compression(stress, psi, phi)

		# Initialize
		fis = np.zeros(alphas.shape)
		R_cr = np.zeros((3,3))
		R_cr[:,0] = np.array([1., 0., 0.])
		normal = np.zeros((3))
		normal[0] = 0.
		mode_mix_limit = 1.e-6

		# Loop through each alpha
		for i, alpha in enumerate(alphas):

			# Project stress onto crack plane
			normal[1] = np.cos(alpha)
			normal[2] = np.sin(alpha)
			R_cr[:,1] = normal
			R_cr[:,2] = np.cross(R_cr[:,0], R_cr[:,1])
			T = np.matmul(stress.mat, R_cr[:,1])

			# Get cohesive displacements
			Pen = np.zeros((3))
			penStiffMult = 1e4
			Lc = 0.2
			Pen[1] = penStiffMult*E2/Lc
			Pen[0] = Pen[1]*GYT*SL*SL/(GSL*YT*YT) # Corresponds to Turon et al (2010)
			Pen[2] = Pen[1]*GYT*ST*ST/(GSL*YT*YT)
			delta = np.matmul(np.transpose(R_cr), T) / Pen

			# Get FI
			del_s1 = delta[0]
			del_s2 = delta[2]
			del_n = delta[1]
			del_s = np.sqrt(del_s1*del_s1 + del_s2*del_s2)
			del_ = np.sqrt(np.max([0., del_n])*del_n + del_s*del_s)

			# Account for LaRC04 and calculate "mixed shear" strengths and stiffnesses
			delta_n_init = delta[1]
			if del_s > 0.:
				ds1_str = SL - etaL*Pen[1]*np.min([0., delta_n_init])  # LaRC04 longitudinal shear strength
				ds2_str = ST - etaT*Pen[1]*np.min([0., delta_n_init])  # LaRC04 transverse shear strength
				KS = np.sqrt((Pen[0]*del_s1)**2 + (Pen[2]*del_s2)**2)/del_s  # Combined shear penalty stiffness
				ds_str = KS*del_s/np.sqrt((Pen[0]*del_s1/ds1_str)**2 + (Pen[2]*del_s2/ds2_str)**2)
			else:
				ds_str = SL
				KS = Pen[0]

			# Cohesive displacements for initiation for pure mode I and mode II
			del0n = YT/Pen[1]  # mode I cohesive initiation disp.
			del0s = ds_str/KS  # mode II cohesive initiation disp.

			# Mode mixity
			beta = del_s*del_s + np.max([0., del_n])*del_n
			if beta == 0.:
				beta = 1
			else:
				beta = del_s*del_s/beta
			B = KS*beta/(KS*beta + Pen[1]*(1. - beta))

			# Mixed-mode initiation displacements
			if (B >= 1. - mode_mix_limit):
				beta = 1.
				B = 1.
				d0 = del0s
			elif (B <= mode_mix_limit):
				beta = 0.
				B = 0.
				d0 = del0n
			else:
				del0nB = np.sqrt(((1 - B**eta_BK)*del0n*del0n + KS/Pen[1]*B**eta_BK*del0s*del0s)*(1. - B))
				del0sB = np.sqrt(beta/(1. - beta))*del0nB
				d0 = np.sqrt(del0nB*del0nB + del0sB*del0sB)

			fis[i] = del_/d0
			# fis[i] = np.min([1., del_/d0])

		# Max
		fis = np.asarray(fis)
		idx = np.argmax(fis)
		fi_max = fis[idx]
		alpha_at_max = alphas[idx]

		if return_dict:
			out_dict = {'fi': fi_max, 'alpha': alpha_at_max, 'mode': 'm'}
			if fi_kink > 0.:
				out_dict['fi_kink'] = fi_kink
			return out_dict
		else:
			return fi_max


	def FI_kink(self, stress, psi=None, phi=None, return_dict=False):
		(stress, psi, phi, fi_kink) = self.transform_stress_for_compression(stress, psi, phi)
		if return_dict:
			return {'fi': fi_kink, 'mode': 'k', 'phi': phi, 'psi': psi}
		else:
			return fi_kink


	def FI_fiberbreak(self, stress, psi=None, phi=None, return_dict=False):
		(stress, psi, phi, fi_kink) = self.transform_stress_for_compression(stress, psi, phi)
		bending_strain = -1*self.fiber_diameter*np.tan(np.abs(phi/2.))/self.w_kb
		axial_strain = stress.sig11/self.pem.m.e1
		fi = (bending_strain+axial_strain)/(-1*self.fiber_rupture_strain)

		if return_dict:
			return {'fi': fi, 'mode': 'fb', 'phi': phi}
		else:
			return fi
