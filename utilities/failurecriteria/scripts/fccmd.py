import argparse
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.transforms
import numpy as np
import sys

import failurecriteria as fcrit


# IMPLEMENTED_FAILURE_CRITERIA = {'LaRC03': fcrit.LaRC03, 'LaRC04': fcrit.LaRC04, 'CompDam': fcrit.CompDam}
IMPLEMENTED_FAILURE_CRITERIA = {'CompDam': fcrit.CompDam}

IMPLEMENTED_FAILURE_ENVELOP_DICT = {
	'1122': {
		'x_label': r'$\sigma_{11} [MPa]$',
		'y_label': r'$\sigma_{22} [MPa]$'
	},
	'2212': {
		'x_label': r'$\sigma_{22} [MPa]$',
		'y_label': r'$\sigma_{12} [MPa]$'
	},
	'1112': {
		'x_label': r'$\sigma_{11} [MPa]$',
		'y_label': r'$\sigma_{12} [MPa]$'
	}
}

IMPLEMENTED_SHEAR_NONLINEARITY = ['linear', 'ramberg_osgood']

DEFAULT_FEATURE_FLAG = 110300

def _entry(args):
	if args.list_failurecriteria:
		for fc in IMPLEMENTED_FAILURE_CRITERIA.keys():
			print(fc)


def _set_failure_envelop_plot_limits(plane, mat):
	if plane == '1122':
		plt.xlim(left=-2.*mat['XC'], right=0.)
		plt.ylim(top=1.2*mat['YT'], bottom=-1.2*mat['YC'])
	elif plane == '2212':
		plt.xlim(left=-1.2*mat['YC'], right=1.2*mat['YT'])
		plt.ylim(top=2*mat['SL'], bottom=0.)
	elif plane == '1112':
		plt.xlim(left=-mat['G12'], right=0.)
		# plt.ylim(top=1.2*mat['SL'], bottom=0.)
		plt.ylim(top=1.2*mat['SL'], bottom=-1.2*mat['SL'])
	else:
		raise Exception('Not implemented')


def _evaluate_failure_criteria(args):
	'''
	Loop through each failure criteria and evaluate
	'''

	s = fcrit.Stress.fromVec(args.stress)

	if args.verbose:
		print('Evaluating failure criteria with the following arguments:')
		print('Failure Criteria: {}'.format(args.fc))
		print('Stress state:\n{}'.format(s))
		print('Material: {}'.format(args.material))
		if args.criteria:
			print('Criterion: {}'.format(args.criteria))

	phi0s = _initializephi0(args)
	mat = fcrit.Materials.materials[args.material]

	for phi0 in phi0s:
		for fc in args.fc:
			fc_class = IMPLEMENTED_FAILURE_CRITERIA[fc]
			# Check for shear nonlinearity
			if fc in ('LaRC04'):
				fc_class = _get_shear_nonlinearity_initializer(fc, args.shear_nonlinearity)

			if fc in ('CompDam'):
				if phi0 == None:
					phi0 = 0.
				fc_obj = fc_class(mat, feature_flags=args.feature_flags, CDM_phi0_12=phi0, verbose=args.verbose)
			else:
				fc_obj = fc_class(mat, phi0=phi0, verbose=args.verbose)

			# Get failure criterion evaluation function
			if args.criteria:
				for c in args.criteria:
					criterion_fun_name = 'FI_'+c
					fi_fun = getattr(fc_obj, criterion_fun_name)
					print(fi_fun(s, return_dict=True))
			else:
				print(getattr(fc_obj, 'FI')(s, return_dict=True))


def _plot_failure_envelops(args):
	'''
	Loop through each failure criteria and plot failure envelop
	'''

	if args.verbose:
		print('Plotting failure envelop {0} using {1} points with arguments:'.format(args.plane, args.npoints))
		print('Failure Criteria: \n{}'.format(args.fc))
		print('Material: {}'.format(args.material))
		if args.criteria:
			print('Criterion: {}'.format(args.criteria))

	phi0s = _initializephi0(args)
	mat = fcrit.Materials.materials[args.material]

	for phi0 in phi0s:
		for fc in args.fc:

			# Get the class
			fc_class = IMPLEMENTED_FAILURE_CRITERIA[fc]
			# Check for shear nonlinearity
			if fc in ('LaRC04'):
				fc_class = _get_shear_nonlinearity_initializer(fc, args.shear_nonlinearity)
			
			# Initialize the class
			if fc in ('CompDam'):
				if phi0 == None:
					phi0 = 0.
				fc_obj = fc_class(mat, feature_flags=args.feature_flags, CDM_phi0_12=phi0, verbose=args.verbose)
			else:
				fc_obj = fc_class(mat, phi0=phi0, verbose=args.verbose)
			if phi0:
				ll = r'$\phi_{0}=$' + '{:.1f}'.format(phi0*180/np.pi) + r'$^\circ$'
			else:
				ll = ''
			sl = True if len(phi0s) > 1 else False
			getattr(fc_obj, 'failure_envelop_{}'.format(args.plane))(args.npoints, criteria=args.criteria, plot=True, marker=args.marker, savetxt=args.save, base_label=ll, single_label=sl)


	plt.xlabel(IMPLEMENTED_FAILURE_ENVELOP_DICT[args.plane]['x_label'])
	plt.ylabel(IMPLEMENTED_FAILURE_ENVELOP_DICT[args.plane]['y_label'])
	_set_failure_envelop_plot_limits(args.plane, mat)
	plt.legend()
	plt.show()

def _plot_shear_curve(args):
	'''
	Creates a plot of the shear non-linearity law
	Only implemented for CompDam
	'''
	mat = fcrit.Materials.materials[args.material]
	pem = fcrit.pem.PemLoader(args.feature_flags, **mat)
	ro = fcrit.pem.RambergOsgood(pem)
	ro.plot(npts=250, gamma_lims=[0., 0.2], derivatives=False, show_constant_slope_limit=True)


def _considere_construction(args):
	'''
	Creates a considere construction
	'''

	s = fcrit.Stress.fromVec(args.stress)

	if args.verbose:
		print('Plotting considere construction with the following arguments:')
		print('Stress state:\n{}'.format(s))
		print('Material: {}'.format(args.material))

	phi0s = _initializephi0(args)
	fig = plt.figure(figsize = [14,6]) 
	ax = fig.add_subplot(111)

	mat = fcrit.Materials.materials[args.material]
	colors = list(mcolors.TABLEAU_COLORS)

	if len(phi0s) == 1:
		plot_driving_force = True
		legend = ['ramberg_osgood', 'shear_strength', 'zero', 'kinking_critical', 'current_misalignment']
	else:
		plot_driving_force = False
		legend = []

	for i, phi0 in enumerate(phi0s):
		for fc in args.fc:
			fc_class = IMPLEMENTED_FAILURE_CRITERIA[fc]
			# Check for shear nonlinearity
			if fc in ('LaRC04'):
				fc_class = _get_shear_nonlinearity_initializer(fc, args.shear_nonlinearity)
			if fc in ('CompDam'):
				if phi0 == None:
					phi0 = 0.
				fc_obj = fc_class(mat, feature_flags=args.feature_flags, CDM_phi0_12=phi0, verbose=args.verbose)
			else:
				fc_obj = fc_class(mat, phi0=phi0, verbose=args.verbose)

			if hasattr(fc_obj, 'considere'):
				fc_obj.considere(s, ax=ax, plot_driving_force=plot_driving_force, gamma_lims=[-0.2,0.2], legend=legend, kinking_color=colors[i], debugging_plots=False)
			else:
				print('SKIPPED. considere() is not implemented for {}'.format(fc))

			if fc in ('CompDam'):
				fc_obj.close()

	ax.set_ylim(bottom=-1.5*mat['SL'], top=1.5*mat['SL'])
	ax.set_xlim(left=-0.2, right=0.2)
	ax.set_xlabel('Shear strain')
	ax.set_ylabel('Shear stress')

	def legend(ax, x0=1,y0=1, direction = "v", padpoints = 3,**kwargs):
		# https://stackoverflow.com/questions/42994338/creating-figure-with-exact-size-and-no-padding-and-legend-outside-the-axes/43001737#43001737
		otrans = ax.figure.transFigure
		t = ax.legend(bbox_to_anchor=(x0,y0), loc=1, bbox_transform=otrans,**kwargs)
		plt.tight_layout(pad=0)
		ax.figure.canvas.draw()
		plt.tight_layout(pad=0)
		ppar = [0,-padpoints/72.] if direction == "v" else [-padpoints/72.,0] 
		trans2=matplotlib.transforms.ScaledTranslation(ppar[0],ppar[1],fig.dpi_scale_trans)+\
				 ax.figure.transFigure.inverted() 
		tbox = t.get_window_extent().transformed(trans2 )
		bbox = ax.get_position()
		if direction=="v":
			ax.set_position([bbox.x0, bbox.y0,bbox.width, tbox.y0-bbox.y0]) 
		else:
			ax.set_position([bbox.x0, bbox.y0,tbox.x0-bbox.x0, bbox.height])

	# Shrink current axis's height by 10% on the bottom
	legend(ax, y0=0.8, direction="h", borderaxespad=0.2)

	plt.show()


def _get_shear_nonlinearity_initializer(fc, sn):
	initializing_method = 'using_{}'.format(sn)
	return getattr(IMPLEMENTED_FAILURE_CRITERIA[fc], initializing_method)


def _initializephi0(args):
	if args.phi0_rad and args.phi0_deg:
		raise ValueError('Arguments --phi0-rad and --phi0-deg are mutually exclusive. Received both.')
	if args.phi0_deg:
		phi0s = [phi*np.pi/180. for phi in args.phi0_deg]
	elif args.phi0_rad:
		phi0s = args.phi0_rad
	else:
		phi0s = [None,]
	return phi0s


#-----------------------------------------
if __name__ == "__main__":

	'''
	Sample commands:
	python fiber_fc.py verify
	python fiber_fc.py -v eval LaRC03 -100 0 0 73 0 0 EGlassMY750
	python fiber_fc.py -v plot LaRC04 1112 EGlassMY750 --npoints 0 50 0 0 --marker +
	'''

	# mat_obj = Materials.materials['IM78552']
	# considere_dfnonlinearity(mat_obj)

	# raise ''
	# mat_obj = Materials.materials['EGlassMY750']
	# myLcf = Lcf(G12=mat_obj['G12'], alpha=mat_obj['alpha'], eta=mat_obj['eta'], XC=mat_obj['XC'], SL=mat_obj['SL'], YC=mat_obj['YC'], alpha0=mat_obj['alpha0'])
	# print(myLcf.sigc_splitting_res([0.028, -250], 0., 48))
	# myLcf.plt_1112_failure_envelope_phi0s([0.034907,], separate_curves_for_each_mode=True)
	# raise ''

	# Arguments
	parser = argparse.ArgumentParser()
	parser.add_argument('-v', '--verbose', action='store_true', default=False, help='Print information for debugging')
	parser.add_argument('-vv', '--vverbose', action='store_true', default=False, help='Print more information for debugging')
	parser.add_argument('-vvv', '--vvverbose', action='store_true', default=False, help='Print all information for debugging')
	parser.add_argument('-l', '--list_failurecriteria', action='store_true', default=False, help='List the implemented failure criteria')
	parser.set_defaults(func=_entry)
	subparsers = parser.add_subparsers()

	# # Call verifications
	# parser_verify = subparsers.add_parser('verify', help='Run verifications that show the code results overlaid on published results')
	# parser_verify.add_argument('cases', nargs='?', action='store', type=str, default=[], help='Specify specific verifications to run')
	# parser_verify.add_argument('-l', '--list_verifications', action='store_true', default=False, help='Print list of available verifications')
	# parser_verify.add_argument('--save', default=None, help='Name of file to save data used to create the plot')
	# parser_verify.set_defaults(func=Verification.entry)

	# Evaluate specific failure criteria for a given stress state
	parser_eval = subparsers.add_parser('eval', help='Evaluate a failure criteria at a given stress state')
	parser_eval.add_argument('fc', nargs='+', choices=IMPLEMENTED_FAILURE_CRITERIA.keys(), help='Specify one or more failure criteria to evaluate')
	parser_eval.add_argument('material', choices=fcrit.Materials.materials, help='Specify the material to use [run `fiber_fc.py materials --help` for info on materials]')
	parser_eval.add_argument('stress', nargs=6, type=float, action='store', help='Specify stress state [11, 22, 33, 12, 23, 13]')
	parser_eval.add_argument('--criteria', nargs='+', help='Force evaluation of a specific failure criterion')
	parser_eval.add_argument('--shear_nonlinearity', choices=IMPLEMENTED_SHEAR_NONLINEARITY, default='ramberg_osgood', help='Specify the shear nonlinear response to use')
	parser_eval.add_argument('--phi0-rad', nargs='+', type=float, help='Specify phi0[s] to use. Plots are created for each phi0. Applies only to LaRC04. [rad]')
	parser_eval.add_argument('--phi0-deg', nargs='+', type=float, help='Specify phi0[s] to use. Plots are created for each phi0. Applies only to LaRC04. [deg]')
	parser_eval.add_argument('--feature-flags', nargs=1, type=int, default=DEFAULT_FEATURE_FLAG, help='Specify the feature flags to use for CompDam')
	parser_eval.set_defaults(func=_evaluate_failure_criteria)

	# Plot failure envelop
	parser_plot = subparsers.add_parser('plot', help='Plot failure criteria envelop')
	parser_plot.add_argument('fc', nargs='+', choices=IMPLEMENTED_FAILURE_CRITERIA.keys(), help='Specify one or more failure criteria to evaluate')
	parser_plot.add_argument('plane', choices=IMPLEMENTED_FAILURE_ENVELOP_DICT.keys(), help='Specify which failure envelop to plot')
	parser_plot.add_argument('material', choices=fcrit.Materials.materials, help='Specify the material to use [run `fiber_fc.py materials --help` for info on materials]')
	parser_plot.add_argument('--criteria', nargs='+', help='Force evaluation of a specific failure criterion')
	parser_plot.add_argument('--npoints', nargs=4, type=int, default=[50]*4, help='Specify the number of points in each quadrant')
	parser_plot.add_argument('--marker', help='Specify the matplotlib marker symbol to use')
	parser_plot.add_argument('--save', default='', help='Name of file to save data used to create the plot')
	parser_plot.add_argument('--shear_nonlinearity', choices=IMPLEMENTED_SHEAR_NONLINEARITY, default='ramberg_osgood', help='Specify the shear nonlinear response to use')
	parser_plot.add_argument('--phi0-rad', nargs='+', type=float, help='Specify phi0[s] to use. Plots are created for each phi0. Applies only to LaRC04. [rad]')
	parser_plot.add_argument('--phi0-deg', nargs='+', type=float, help='Specify phi0[s] to use. Plots are created for each phi0. Applies only to LaRC04. [deg]')
	parser_plot.add_argument('--feature-flags', nargs=1, type=int, default=DEFAULT_FEATURE_FLAG, help='Specify the feature flags to use for CompDam')
	parser_plot.set_defaults(func=_plot_failure_envelops)

	parser_plot_shear = subparsers.add_parser('plot_shear', help='Plot shear curve')
	parser_plot_shear.add_argument('material', choices=fcrit.Materials.materials, help='Specify the material to use [run `fiber_fc.py materials --help` for info on materials]')
	parser_plot_shear.add_argument('--feature-flags', nargs=1, type=int, default=DEFAULT_FEATURE_FLAG, help='Specify the feature flags to use for CompDam')
	parser_plot_shear.set_defaults(func=_plot_shear_curve)


	# Evaluate considere construction
	parser_consid = subparsers.add_parser('considere', help='Create a considere construction')
	parser_consid.add_argument('fc', nargs='+', choices=IMPLEMENTED_FAILURE_CRITERIA.keys(), help='Specify one or more failure criteria to evaluate')
	parser_consid.add_argument('material', choices=fcrit.Materials.materials, help='Specify the material to use [run `fiber_fc.py materials --help` for info on materials]')
	parser_consid.add_argument('stress', nargs=6, type=float, action='store', help='Specify stress state [11, 22, 33, 12, 23, 13]')
	parser_consid.add_argument('--shear_nonlinearity', choices=IMPLEMENTED_SHEAR_NONLINEARITY, default='ramberg_osgood', help='Specify the shear nonlinear response to use')
	parser_consid.add_argument('--phi0-rad', nargs='+', type=float, help='Specify phi0[s] to use. Plots are created for each phi0. Applies only to LaRC04. [rad]')
	parser_consid.add_argument('--phi0-deg', nargs='+', type=float, help='Specify phi0[s] to use. Plots are created for each phi0. Applies only to LaRC04. [deg]')
	parser_consid.add_argument('--feature-flags', nargs=1, type=int, default=DEFAULT_FEATURE_FLAG, help='Specify the feature flags to use for CompDam')
	parser_consid.set_defaults(func=_considere_construction)

	# Parse the arguments
	if len(sys.argv)==1:
		print('\nERROR: the script expects one or more arguments. Here is the usage overview\n')
		print('Argument list: ' + str(sys.argv))
		print('')
		parser.print_help(sys.stderr)
		sys.exit(1)
	args = parser.parse_args()

	# Process arguments
	if args.vverbose:
		args.verbose = 2
	if args.vvverbose:
		args.verbose = 3

	# Call the appropriate entry point
	args.func(args)
