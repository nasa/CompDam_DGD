import os
import pdb
import numpy as np
import sys

import failurecriteria as fcrit


# def evaluate_failure_criteria(stress, fc_obj):


# 	# Get failure criterion evaluation function
# 	return getattr(fc_obj, 'FI')(stress, return_dict=True)



#-----------------------------------------
if __name__ == "__main__":

	job_name = 'OHC_delam_fkt_06'

	# Read in stresses
	files = {
		'S11_pb2': job_name + '_proc_S11_PLY_0_NOCOH_2PLYBLOCK_2.csv',
		'S11_pb5': job_name + '_proc_S11_PLY_0_NOCOH_2PLYBLOCK_5.csv',
		'S11_pb8': job_name + '_proc_S11_PLY_0_NOCOH_4PLYBLOCK_8.csv',
		'S12_pb2': job_name + '_proc_S12_PLY_0_NOCOH_2PLYBLOCK_2.csv',
		'S12_pb5': job_name + '_proc_S12_PLY_0_NOCOH_2PLYBLOCK_5.csv',
		'S12_pb8': job_name + '_proc_S12_PLY_0_NOCOH_4PLYBLOCK_8.csv',
	}

	stresses = {}
	for label, file in files.items():
		stresses[label] = np.genfromtxt(os.path.join('..', 'data', file), delimiter=',', skip_header=1)

	phi0 = 2*np.pi/180.
	mat = fcrit.Materials.materials['IM78552ACC']
	fc_obj = fcrit.CompDam(mat, feature_flags=110300, CDM_phi0_12=phi0, verbose=False)

	for pb in ('pb2', 'pb5', 'pb8'):
		s11 = stresses['S11_'+pb]
		s12 = stresses['S12_'+pb]
		(num_x, num_frames) = s11.shape
		num_frames -= 1

		fi_kink = np.zeros((num_x, num_frames))
		for j in range(num_frames):
			for i in range(num_x):
				stress = fcrit.Stress.fromVec([s11[i,j+1], 0, 0, s12[i,j+1], 0, 0])
				result =fc_obj.FI(stress, return_dict=True)
				#result = evaluate_failure_criteria(stress, phi0, material='IM78552ACC', feature_flags=110300)
				fi_kink[i,j] = result['fi']
		fi_knk = np.hstack((s11[:,0].reshape(-1,1), fi_kink))
		np.savetxt('{}_{}_phi0_2p0_fikink.csv'.format(job_name, pb), fi_knk, delimiter=',', header='x,17,23,29,33,36')
		print('Finished '+pb)

	# stress = fcrit.Stress.fromVec([-1200, 0, 0, 0, 0, 0])
	# phi0 = -1.*np.pi/180.
	# evaluate_failure_criteria(stress, phi0, material='IM78552ACC', feature_flags=110300)
