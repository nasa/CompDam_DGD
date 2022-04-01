"""
TODO: add option to remove a specific package from abaqus environment
TODO: add support for 'install_requires' in local package setup files --> search for packages in conda
TODO: add option to sync packages from abaqus to conda environment
TODO: only for abaqus 2017 and 2018; add capability for earlier versions
"""

import sys
import os
import subprocess
import re
import shutil
import argparse
import stat



def _copy_contents_recursively_no_overwrite(src_dir, dest_dir):
	for root, dirs, files in os.walk(src_dir):
		# Skip top-level directories that are already in the destination
		if root == src_dir:
			existing_dest_dirs = [o for o in os.listdir(dest_dir) if os.path.isdir(os.path.join(dest_dir,o))]
			dirs[:] = [d for d in dirs if d not in existing_dest_dirs] 
		# print('--------- COUNTER ' + str(counter))
		# print('root: ' + str(root))
		# print('dirs: ' + str(dirs))
		# print('files: ' + str(files))

		for name in files:
			# Get the src, dest file paths
			if root == src_dir:
				# print('Copying: '+name)
				src = os.path.join(src_dir, name)
				dest = os.path.join(dest_dir, name)
			else:
				rs = root.split(src_dir)
				if len(rs) > 1:
					rel_path = '\\'.join(rs[1:])+'\\'+name
				else:
					rel_path = name
				src = src_dir + rel_path
				dest = dest_dir + rel_path
				# print('Copying: '+rel_path[1:])
			
			# If the file doesn't exist in the destination, add it
			if not os.path.exists(dest):
				# print('Copying')
				# print('\tsrc: ' + src)
				# print('\tdest: ' + dest)
				shutil.copy2(src, dest)

		for name in dirs:
			
			if root == src_dir:
				directory_path = os.path.join(dest_dir, name)
			else:
				rel_path = os.path.join('\\'.join(root.split(src_dir)[1:]), name)
				directory_path = dest_dir + rel_path
			if not os.path.isdir(directory_path):
				# print('Creating directory')
				# print('\t'+directory_path)
				os.makedirs(directory_path)


def _get_conda_abq_environments():
	res = subprocess.check_output('conda env list').decode('utf-8').split('\n')
	conda_abq_envs = [env.split(' ')[0] for env in res if env.startswith('abq_')]
	return conda_abq_envs


def _get_conda_packages(verbose, env):
	"""
	Uses conda list and pip list to get a complete list of python packages in the current environment
	"""

	# For return
	pkgs = []

	# Check that the conda environment exists
	ce = _get_conda_abq_environments()
	if env not in ce:
		print('\nERROR: no conda environment {0} exists, which means there are no packages to list.\n'.format(env))
		print('Maybe you want to use a different version of abaqus? Try the --abaqus-cmd option.\n')
		print('To install Abaqus packages in python try: ')
		print('> abaqus-python-addpkg.py install package_1, [package_2, ...]')
		print('\nExiting...\n')
		sys.exit()

	# Check conda
	if verbose:
		print('Generating a list of packages in the conda environment: ' + env)
		print('Running command: ' + 'conda list -n '+env)
	res = subprocess.check_output('conda list -n '+env).decode('utf-8')
	# Check for err
	if 'EnvironmentLocationNotFound' in res:
		raise ValueError('No conda environment named: ' + env)
	# Return a list of packages
	pkgs_list_conda = [pkg.split() for pkg in res.split('\n') if not pkg.startswith('#')]
	[pkgs.append({'name': pkg[0], 'version': pkg[1], 'channel': pkg[2]}) for pkg in pkgs_list_conda if len(pkg) == 3]
	conda_pkg_names = [pkg['name'] for pkg in pkgs]

	# Check pip too
	cmd = 'conda activate ' + env + '&pip list --disable-pip-version-check'
	if verbose:
		print('Checking pip')
		print('Running command: ' + cmd)
	p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
	stdout, stderr = p.communicate()
	if stdout:
		res = stdout.decode('utf-8')
		pkgs_list_pip = [pkg.split() for pkg in res.split('\n')[2:]]
		for pkg in pkgs_list_pip:
			if len(pkg) > 1:
				name = pkg[0]
				if name not in conda_pkg_names:
					if len(pkg) == 2:
						pkgs.append({'name': pkg[0], 'version': pkg[1], 'channel': ''})
					elif len(pkg) == 3:
						pkgs.append({'name': pkg[0], 'version': pkg[1], 'channel': pkg[2]})
	if stderr:
		if 'DEPRECATION' not in stderr.decode('utf-8'):
			print('\nERROR while attempting to install a pip package.')
			print(stderr)
			print('Exiting ...\n')
			sys.exit()

	return pkgs


def _get_abaqus_version(abaqus_cmd):
	ri = subprocess.check_output(abaqus_cmd + ' information=release', shell=True).decode('utf-8')
	grps = re.search(r'^Abaqus (20\d+|6.\d+)\.*(.*)$', ri, re.MULTILINE).groups()
	year = grps[0]
	hotfix = grps[1]
	return (year, hotfix)


def _get_abaqus_python_loc(abaqus_cmd):
	"""
	Get the paths to the python interpreter for CAE and the solver
	returns (path_to_cae_python, path_to_solver_python)
	"""

	rel_path_to_py = 'tools\\SMApy\\python2.7'

	# Get the abaqus version information
	(abq_year, abq_hotfix) = _get_abaqus_version(abaqus_cmd)

	# Get the abaqus python interpreter locations
	ri = subprocess.check_output(abaqus_cmd + ' information=release', shell=True).decode('utf-8')

	# CAE
	abq_cae_dirs = [l for l in ri.split('\n') if l.startswith('Abaqus is located in the directory')]
	if len(abq_cae_dirs) != 1:
		raise Exception('Parsing CAE path failed')
	else:
		try:
			abq_cae_dirs = abq_cae_dirs[0].split('Abaqus is located in the directory ')[1]
			abq_cae_dir = os.path.join(os.path.normpath(abq_cae_dirs.split()[1]), rel_path_to_py)
		except:
			raise Exception('Parsing CAE path failed')
	
	# Solver
	abq_solver_dirs = [l for l in ri.split('\n') if l.startswith('Abaqus solver stack is located in the directory')]
	if len(abq_solver_dirs) != 1:
		raise Exception('Parsing CAE path failed')
	else:
		try:
			abq_solver_dir = abq_solver_dirs[0].split('Abaqus solver stack is located in the directory ')[1]
			abq_solver_dir = os.path.join(os.path.normpath(abq_solver_dir), rel_path_to_py)
		except:
			raise Exception('Parsing solver path failed')

		return (abq_cae_dir, abq_solver_dir)


def _info(args):
	"""
	Handles 'abaqus-python-addpkg info'
	"""

	# Default behavior
	if not args.list_envs and not args.list_pkgs and not args.show_abaqus_python_directory and not args.test_pkg:
		print('Using Abaqus ' + args.abaqus_year)
		ce = _get_conda_abq_environments()
		if ce:
			if 'abq_' + args.abaqus_year in ce:
				print('')
				print('Conda environment for Abaqus ' + args.abaqus_year + ': abq_'+ args.abaqus_year)
				print('To see the packages that are installed, try:')
				print('> abaqus-python-addpkg.py info --list-pkgs')
				print('')
			else:
				print('')
				print('No conda environment for Abaqus ' + args.abaqus_year)
				print('Found the following conda environments for abaqus: ')
				[print(v) for v in ce]
				print('')
		else:
			print('')
			print('No conda environment for Abaqus')
			print('To install Abaqus packages in python try: ')
			print('> abaqus-python-addpkg.py install package_1, [package_2, ...]')
			print('')

	# List of conda environments
	if args.list_envs:
		ce = _get_conda_abq_environments()
		print('Conda environments for Abaqus:')
		[print(v) for v in ce]
		print('')

	# List packages
	if args.list_pkgs:
		pkgs = _get_conda_packages(verbose=args.verbose, env='abq_'+args.abaqus_year)

		# Format nicely
		col_name = [pkg['name'] for pkg in pkgs]
		col_name_width = max([len(u) for u in col_name]) + 4

		col_version = [pkg['version'] for pkg in pkgs]
		col_version_width = max([len(u) for u in col_version]) + 4

		print('')
		header = 'Name'.ljust(col_name_width+2)
		header += ' Version'.ljust(col_version_width+2)
		header += ' Build/Location\n'
		header += '----'.ljust(col_name_width+2)
		header += ' -------'.ljust(col_name_width)
		header += ' ----------------------'
		print(header)
		for row in pkgs:
			row_formatted = ' ' + row['name'].ljust(col_name_width+3)
			row_formatted += row['version'].ljust(col_version_width+2)
			row_formatted += row['channel']
			print(row_formatted)
		print('')

	if args.show_abaqus_python_directory:
		(cae_dir, solver_dir) = _get_abaqus_python_loc(args.abaqus_cmd)
		print('')
		print('CAE python installation: ' + cae_dir)
		print('Solver python installation: ' + solver_dir)
		print('')

	if args.test_pkg:
		(cae_dir, solver_dir) = _get_abaqus_python_loc(args.abaqus_cmd)
		with open('abaqus_python_addpkg_temp.py', 'w') as f:
			f.write('import {0}\n'.format(','.join(args.test_pkg)))
		print('Checking imports using CAE: ' + cae_dir)
		cmd = args.abaqus_cmd + ' cae nogui=abaqus_python_addpkg_temp.py'
		if args.verbose:
			print('Running command: ' + cmd)
		try:
			p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
			stdout, stderr = p.communicate()
			stdout = stdout.decode('utf-8')
			if 'Abaqus Error' in stdout:
				print('\nERROR while attempting import. Package(s) not installed properly.')
				print(stdout)
			else:
				print('Successfully imported {0}'.format(args.test_pkg))
		finally:
			os.remove('abaqus_python_addpkg_temp.py')


def _install(args):
	"""
	Handles 'abaqus-python-addpkg install'
	"""

	# Make sure at least one package was specified
	if not args.conda and not args.local and not args.pip:
		print('Must specify at least one package to install. For example, to install scipy, try:')
		print('> abaqus-python-addpkg.py install --conda scipy ')
		print('')
		print('To install a package using pip, try:')
		print('> abaqus-python-addpkg.py install --pip plotly')
		print('')
		print('To install a local package (equivalent to pip install -e), try:')
		print('> abaqus-python-addpkg.py install --local path/to/local/package')
		print('')
		return

	# Make sure the user did not specify python or numpy as packages to install
	if args.conda:
		if 'numpy' in args.conda:
			print('Ignoring request to install numpy since it is already installed by default')
			args.conda.remove('numpy')
		if 'python' in args.conda:
			print('Ignoring request to install python')
			args.conda.remove('python')
		# Make sure at least one package is still specified
		if not args.local and len(args.conda) == 0:
			print('Must specify at least one package other than numpy and python to install')
			print('')
			print('For example, to install scipy, try:')
			print('> abaqus-python-addpkg.py install --conda scipy')
			print('')
			return

	# Check for an existing conda environment
	preexisting_conda_env = False
	ce = _get_conda_abq_environments()
	conda_env_name = 'abq_'+args.abaqus_year
	if conda_env_name in ce:
		preexisting_conda_env = True
		print('Found existing conda environment: ' + 'abq_' + args.abaqus_year)

	# If not using conda, make sure conda environment is available. If not, add conda env w/ pip
	install_abq_env = False
	if (args.local and not args.conda) or (args.pip and not args.conda):
		if not preexisting_conda_env:
			args.conda = []
			install_abq_env = True
			print('Installing conda environment with python, numpy, and pip')

	# Compatibilities
	compat = {
		# '6.14': #TODO
		# '2016': #TODO
		'2017': {
			'python': '2.7.3',
			'numpy': '1.6.2',
			'matplotlib': '1.1'
		},
		'2018': {
			'python': '2.7.3',
			'numpy': '1.6.2',
			'matplotlib': '1.1'
		}
	}

	# Abaqus version check
	if args.abaqus_year not in compat.keys():
		raise ValueError('Abaqus ' + args.abaqus_year + ' is not a supported version of Abaqus')

	# Build the conda packages to install
	if args.conda or install_abq_env:

		# First, specify python and numpy versions for compatibility purposes
		pkgs_to_install = []
		pkgs_to_install.append('python='+compat[args.abaqus_year]['python'])
		pkgs_to_install.append('pip')
		pkgs_to_install.append('numpy='+compat[args.abaqus_year]['numpy'])

		# Then, install user specified pacakages
		# If specified packages are in the compatibilities list, use the known versions
		# otherwise, leave the version unspecified
		for arg in args.conda:
			if arg in compat[args.abaqus_year].keys():
				pkgs_to_install.append(arg + '=' + compat[args.abaqus_year][arg])
			else:
				pkgs_to_install.append(arg)

		# Update existing conda environment
		if preexisting_conda_env:
			cmd = 'conda install --name ' + conda_env_name + ' ' + ' '.join(pkgs_to_install)

		# Install new conda environment
		else:
			cmd = 'conda create -n '+conda_env_name+' ' + ' '.join(pkgs_to_install)

		# Create temporary conda environment
		print(cmd)
		subprocess.check_call(cmd, shell=True)

	# Add local packages
	if args.local and len(args.local) > 0:
		cmd = 'conda activate '+conda_env_name+'&cd ' + os.path.abspath(args.local[0]) + '&python setup.py develop&conda deactivate'
		p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
		stdout, stderr = p.communicate()
		if stderr:
			print('\nERROR while attempting to install a local package.')
			print(stderr)
			print('Exiting ...\n')
			sys.exit()

	# Add pip packages
	if args.pip and len(args.pip) > 0:
		for pip_package in args.pip:
			cmd = 'conda activate '+conda_env_name+'&pip install ' + pip_package + '&conda deactivate'
			if args.verbose:
				print('Executing: ' + cmd)
			p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
			stdout, stderr = p.communicate()
			if args.verbose and stdout:
				print(stdout)
			if stderr:
				print('\nERROR while attempting to install a pip package.')
				print(stderr)
				print('Exiting ...\n')
				sys.exit()


	# Copy the site-packages directory contents
	p = subprocess.Popen('conda info --base', stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
	stdout, stderr = p.communicate()
	src_dir = os.path.join(stdout.decode().rstrip()+'\\envs\\'+conda_env_name+'\\Lib\\site-packages')
	(abqpy_cae_path, abqpy_solver_path) = _get_abaqus_python_loc(args.abaqus_cmd)
	dest_dir = os.path.join(abqpy_cae_path, 'Lib\\site-packages')
	_copy_contents_recursively_no_overwrite(src_dir, dest_dir)

	print('')
	print('Python packages installed into abaqus')
	print('')
	return



def _restore(args):
	"""
	Handles 'abaqus-python-addpkg restore'
	"""

	rel_path_sp = 'Lib\\site-packages'

	# Get the paths to the python interpreters
	(abqpy_cae_path, abqpy_solver_path) = _get_abaqus_python_loc(args.abaqus_cmd)

	# print(abqpy_cae_path)
	# print(abqpy_solver_path)

	# Remove the files in the CAE site-packages
	cae_sp = os.path.join(abqpy_cae_path, rel_path_sp)
	shutil.rmtree(cae_sp)

	# Replace the CAE site-packages with the solver site-packages
	os.chmod(os.path.join(abqpy_cae_path, 'Lib'), stat.S_IRWXU| stat.S_IRWXG| stat.S_IRWXO)
	shutil.copytree(os.path.join(abqpy_solver_path, rel_path_sp), cae_sp)

	# Revert the conda environment
	conda_env_name = 'abq_'+args.abaqus_year
	ce = _get_conda_abq_environments()
	if conda_env_name in ce:
		cmd = 'conda env remove --name '+conda_env_name
		out = subprocess.check_call(cmd, shell=True)

	print('Restored abaqus ' + args.abaqus_year + ' to the original packages')
	return



# Main entry point
if __name__ == "__main__":

	# Arguments
	parser = argparse.ArgumentParser()
	parser.add_argument('-a', '--abaqus-cmd', default='abaqus', type=str, help='Specify the command to use to run abaqus')
	parser.add_argument('-v', '--verbose', action='store_true', default=False, help='Print information for debugging')
	subparsers = parser.add_subparsers()

	# info
	parser_info = subparsers.add_parser('info', help='Get information about python packages installed to the abaqus python environment')
	parser_info.add_argument('-l', '--list-envs', action='store_true', default=False, help='List the conda environments created to replicate the abaqus python environments')
	parser_info.add_argument('-L', '--list-pkgs', action='store_true', default=False, help='List the packages installed in the abaqus python environment')
	parser_info.add_argument('--show-abaqus-python-directory', action='store_true', default=False, help='Prints the directory of the abaqus python installation')
	parser_info.add_argument('--test-pkg', action='store', nargs='+', help='Attempts to import the provided package(s) to check if it is installed')
	parser_info.set_defaults(func=_info)

	# install
	parser_install = subparsers.add_parser('install', help='Install python package(s) to the abaqus python environment')
	parser_install.add_argument('-c', '--conda', action='store', nargs='+', help='Install package(s) from conda to the abaqus python environment')
	parser_install.add_argument('-l', '--local', action='store', nargs='+', help='Install package(s) stored locally to the abaqus python environment')
	parser_install.add_argument('-p', '--pip', action='store', nargs='+', help='Install package(s) using pip')
	parser_install.set_defaults(func=_install)

	# restore
	parser_install = subparsers.add_parser('restore', help='Restore the abaqus python environment to the as-installed state')
	parser_install.set_defaults(func=_restore)

	# Parse the args
	if len(sys.argv)==1:
		print('\nERROR: the script expects one or more arguments. Here is the usage overview\n')
		print('Argument list: ' + str(sys.argv))
		print('')
		parser.print_help(sys.stderr)
		sys.exit(1)
	args = parser.parse_args()

	# Get the details on the abaqus year and hot fix (for when multiple versions of abaqus are installed)
	(year, hotfix) = _get_abaqus_version(args.abaqus_cmd)
	args.abaqus_year = year
	args.abaqus_hotfix = hotfix
	if args.verbose:
		print('Using abaqus commmand: ' + args.abaqus_cmd)
		print('Using abaqus year: ' + args.abaqus_year)
		if args.abaqus_hotfix:
			print('Using abaqus hotfix: ' + args.abaqus_hotfix)

	# Call the appropriate entry point
	args.func(args)
