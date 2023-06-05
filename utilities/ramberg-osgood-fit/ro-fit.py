"""
Functions to plot and fit ramberg-osgood curves
"""

import argparse
import matplotlib.pyplot as plt
import numpy as np
import os
import scipy.optimize
import sys
import textwrap

extra_help = textwrap.dedent('''\
Sample commands:

	# Plot experimental data
	python ro-fit.py -d soden1998_eglass_MY750.csv

	# Plot experimental data and corresponding Ramberg-Osgood fit (must provide the modulus to generate a fit)
	python ro-fit.py 5830 -d soden1998_eglass_MY750.csv

	# Plot a Ramberg-Osgood curve
	python ro-fit.py 5830 4.4e-10 5.9

	# Plot experimental data, corresponding Ramberg-Osgood fit, and another Ramberg-Osgood curve
	python ro-fit.py 5830 4.4e-10 5.9 -d soden1998_eglass_MY750.csv

''')

def _main_entry(args):

	# Ramberg-Osgood equation
	ro = lambda alpha, eta, s:  (s + alpha*np.abs(s)**eta)/args.modulus

	# If experimental data argument is supplied, try to read it in and plot
	if args.expdatacsv:
		out = np.loadtxt(args.expdatacsv, delimiter=',', skiprows=args.num_header_rows)
		strain = out[:,0]
		stress = out[:,1]

		if len(stress) != len(strain):
			raise Exception('Expect the same number of data points for stress and strain.')

		if args.verbose:
			print('Successfully read in expdatacsv with {0} points'.format(len(stress)))

		if not args.noplot:
			plt.plot(strain, stress, 'kx', label='Data')

		# Calculate ramberg-osgood fit to experimental data
		if args.modulus:

			if args.verbose:
				print('Computing least-squares fit of experimental data')

			# Function to fit
			ro_lsq = lambda x, stress, strain: ro(x[0], x[1], stress) - strain
			# Initial guess for alpha, eta
			if args.eta and args.alpha:
				x0 = np.array([args.alpha, args.eta])
			else:
				x0 = np.array([4.4e-10, 5.9])
			# Solve least-squares minimization
			res_lsq = scipy.optimize.least_squares(ro_lsq, x0, args=(stress, strain))
			alpha_fit = res_lsq.x[0]
			eta_fit = res_lsq.x[1]
			
			# Print result
			print('Ramberg-Osgood fit complete:')
			print('alpha: {0}, eta: {1}'.format(alpha_fit, eta_fit))
			
			# Plot result
			if not args.noplot:
				stress_fit = np.linspace(0,max(stress),100)
				strain_fit = np.array([ro(alpha_fit, eta_fit, s) for s in stress_fit])
				plt.plot(strain_fit, stress_fit, '-b', label='Ramberg-Osgood fit')
		else:
			if args.verbose:
				print('No modulus provided, not attempting to fit Ramberg-Osgood curve')

	# If Ramberg-Osgood curve parameters are supplied, plot the curve
	if args.modulus and args.eta and args.alpha:
		stress = np.linspace(0,args.maximum_stress,100)
		strain = np.array([ro(args.alpha, args.eta, s) for s in stress])

		if not args.noplot:
			plt.plot(strain, stress, '-r', label='Ramberg-Osgood')

		if len(args.considere):
			phi0 = args.considere[0]
			sig11 = args.considere[1]
			sig22 = args.considere[2]
			sig12 = args.considere[3]

			strain_c = np.linspace(-phi0,max(strain),100)
			stress_c = np.array([(sig11-sig22)/-2*np.sin(2*(g+phi0)) + sig12*np.cos(2*(g+phi0)) for g in strain_c])
			plt.plot(strain_c, stress_c, '-k', label='Driving force')

	elif any([args.eta, args.alpha]):
		print('WARNING: must specify modulus, alpha, and eta (all three) together. Found a partial set.')

	if not args.noplot:
		plt.xlabel('Strain')
		plt.ylabel('Stress')
		plt.legend()
		plt.show()


# Main entry point
if __name__ == "__main__":

	# Arguments
	parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, epilog=extra_help)
	parser.add_argument('modulus', nargs='?', action='store', type=float, help='Ramberg-Osgood modulus.')
	parser.add_argument('alpha', nargs='?', action='store', type=float, help='Ramberg-Osgood parameter alpha. If not specified, experimental data is expected and the Ramberg-Osgood fit is calculated.')
	parser.add_argument('eta', nargs='?', action='store', type=float, help='Ramberg-Osgood parameter eta. If not specified, experimental data is expected and the Ramberg-Osgood fit is calculated.')
	parser.add_argument('-c', '--considere', nargs=4, default=[], type=float, help='Inlcude driving force curve on plot. [phi0, sig11, sig22, sig12].')
	parser.add_argument('-m', '--maximum_stress', default=100., type=float, help='Maximum stress level to use for plotting the analytical Ramberg-Osgood curve.')
	parser.add_argument('-d', '--expdatacsv', default='', type=str, help='Specify the path to a csv file containing experimental stress-strain curve as xy points. First column is strain, second column is stress. See sample file `soden1998_eglassDY063.csv`')
	parser.add_argument('-s', '--num_header_rows', default=1, type=int, help='Number of rows to skip at start of csv file containing experimental stress-strain curve data.')
	parser.add_argument('--noplot', action='store_true', default=False, help='Generate a plot of the ramberg-osgood curve.')
	parser.add_argument('-v', '--verbose', action='store_true', default=False, help='Print information for debugging.')

	parser.set_defaults(func=_main_entry)

	# Parse the args
	if len(sys.argv)<=2:
		parser.print_help(sys.stderr)
		sys.exit(1)
	args = parser.parse_args()
	args.func(args)
