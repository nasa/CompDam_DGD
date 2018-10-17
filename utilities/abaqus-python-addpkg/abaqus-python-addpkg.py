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


def _get_conda_packages(env):
	res = subprocess.check_output('conda list -n '+env).decode('utf-8')
	# Check for err
	if 'EnvironmentLocationNotFound' in res:
		raise ValueError('No conda environment named: ' + env)
	# Return a list of packages
	pkgs_list = [pkg.split() for pkg in res.split('\n') if not pkg.startswith('#')]
	pkgs = []
	[pkgs.append({'name': pkg[0], 'version': pkg[1], 'channel': pkg[2]}) for pkg in pkgs_list if len(pkg) == 3]
	return pkgs


def _get_abaqus_python_loc(version):
	"""
	Get the paths to the python interpreter for CAE and the solver
	returns (path_to_cae_python, path_to_solver_python)
	"""

	rel_path_to_py = 'tools\\SMApy\\python2.7'

	# Get the abaqus python interpreter locations
	ri = subprocess.check_output('abaqus information=release', shell=True).decode('utf-8')

	# Check the abqus version
	grps = re.search(r'^Abaqus (20\d+|6.\d+)(.*)$', ri, re.MULTILINE).groups()
	abq_version_full = 'Abaqus ' + ''.join(grps)
	abq_version = grps[0]
	if abq_version != version:
		print('Fix me TODO')
		raise Exception('The specified version is not the version called with abaqus')

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
	if not args.list_envs and not args.list_pkgs:
		print('Abaqus ' + args.abaqus_version)
		ce = _get_conda_abq_environments()
		if ce:
			if 'abq_' + args.abaqus_version in ce:
				print('')
				print('Conda environment for Abaqus ' + args.abaqus_version + ': abq_'+ args.abaqus_version)
				print('To see the packages that are installed, try:')
				print('> abaqus-python-addpkg info --list-pkgs')
				print('')
			else:
				print('')
				print('No conda environment for Abaqus ' + args.abaqus_version)
				print('Found the following conda environments for abaqus: ')
				[print(v) for v in ce]
				print('')
		else:
			print('')
			print('No conda environment for Abaqus')
			print('To install Abaqus packages in python try: ')
			print('> abaqus-python-addpkg install package_1, [package_2, ...]')
			print('')

	# List of conda environments
	if args.list_envs:
		ce = _get_conda_abq_environments()
		print('Conda environments for Abaqus:')
		[print(v) for v in ce]
		print('')

	# List packages
	if args.list_pkgs:
		# TODO: format nicely
		pkgs = _get_conda_packages(env='abq_'+args.abaqus_version)
		print('')
		print('Name, Version, Build Channel')
		[print(pkg['name']+', '+pkg['version']+', '+pkg['channel']) for pkg in pkgs]
		print('')


def _install(args):
	"""
	Handles 'abaqus-python-addpkg install'
	"""

	# Make sure at least one package was specified
	if not args.conda and not args.local and not args.pip:
		print('Must specify at least one package to install. For example, to install scipy, try:')
		print('> abaqus-python-addpkg install --conda scipy ')
		print('')
		print('To install a package using pip, try:')
		print('> abaqus-python-addpkg install --pip plotly')
		print('')
		print('To install a local package (equivalent to pip install -e), try:')
		print('> abaqus-python-addpkg install --local path/to/local/package')
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
			print('> abaqus-python-addpkg install --conda scipy')
			print('')
			return

	# Check for an existing conda environment
	preexisting_conda_env = False
	ce = _get_conda_abq_environments()
	conda_env_name = 'abq_'+args.abaqus_version
	if conda_env_name in ce:
		preexisting_conda_env = True
		print('Found existing conda environment: ' + 'abq_' + args.abaqus_version)

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
	if args.abaqus_version not in compat.keys():
		raise ValueError('Abaqus ' + args.abaqus_version + ' is not a supported version of Abaqus')

	# Build the conda packages to install
	if args.conda or install_abq_env:

		# First, specify python and numpy versions for compatibility purposes
		pkgs_to_install = []
		pkgs_to_install.append('python='+compat[args.abaqus_version]['python'])
		pkgs_to_install.append('pip')
		pkgs_to_install.append('numpy='+compat[args.abaqus_version]['numpy'])

		# Then, install user specified pacakages
		# If specified packages are in the compatibilities list, use the known versions
		# otherwise, leave the version unspecified
		for arg in args.conda:
			if arg in compat[args.abaqus_version].keys():
				pkgs_to_install.append(arg + '=' + compat[args.abaqus_version][arg])
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
		cmd = 'activate '+conda_env_name+'&cd ' + os.path.abspath(args.local[0]) + '&python setup.py develop&deactivate'
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
			cmd = 'activate '+conda_env_name+'&pip install ' + pip_package + '&deactivate'
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
	src_dir = os.path.join(os.path.expanduser("~"), 'AppData\\Local\\Continuum\\anaconda3\\envs\\'+conda_env_name+'\\Lib\\site-packages')
	(abqpy_cae_path, abqpy_solver_path) = _get_abaqus_python_loc(args.abaqus_version)
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
	(abqpy_cae_path, abqpy_solver_path) = _get_abaqus_python_loc(args.abaqus_version)

	# print(abqpy_cae_path)
	# print(abqpy_solver_path)

	# Remove the files in the CAE site-packages
	cae_sp = os.path.join(abqpy_cae_path, rel_path_sp)
	shutil.rmtree(cae_sp)

	# Replace the CAE site-packages with the solver site-packages
	os.chmod(os.path.join(abqpy_cae_path, 'Lib'), stat.S_IRWXU| stat.S_IRWXG| stat.S_IRWXO)
	shutil.copytree(os.path.join(abqpy_solver_path, rel_path_sp), cae_sp)

	# Revert the conda environment
	conda_env_name = 'abq_'+args.abaqus_version
	ce = _get_conda_abq_environments()
	if conda_env_name in ce:
		cmd = 'conda env remove --name '+conda_env_name
		out = subprocess.check_call(cmd, shell=True)

	print('Restored abaqus ' + args.abaqus_version + ' to the original packages')
	return



# Main entry point
if __name__ == "__main__":

	# Arguments
	parser = argparse.ArgumentParser()
	parser.add_argument('-a', '--abaqus-version', default='2017', type=str, help='Specify the version of abaqus to use')
	parser.add_argument('-v', '--verbose', action='store_true', default=False, help='Print information for debugging')
	subparsers = parser.add_subparsers()

	# info
	parser_info = subparsers.add_parser('info', help='Get information about python packages installed to the abaqus python environment')
	parser_info.add_argument('-l', '--list-envs', action='store_true', default=False, help='List the conda environments created to replicate the abaqus python environments')
	parser_info.add_argument('-L', '--list-pkgs', action='store_true', default=False, help='List the packages installed in the abaqus python environment')
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
		parser.print_help(sys.stderr)
		sys.exit(1)
	args = parser.parse_args()
	args.func(args)
