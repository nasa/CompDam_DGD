# Remote run directory [Default is abaverify_temp]
# Command line run directory will override the value specified here if both are provided
#remote_run_directory = 'abaverify_temp'


# Add files to copy to the remote [Default is empty list]
# Paths relative to /tests directory
local_files_to_copy_to_remote = ['CompDam.parameters']


# Specify which files to copy back from remote to local after job ends
file_extensions_to_copy_to_local = ['.dat', '.inp', '.msg', '.odb', '.sta', '.py']


# Regular expression for source files to copy [Default is below]
#source_file_extensions_to_copy = r'.*\.for$'


# Copy results back to local directory? [Defualt is False]
copy_results_to_local = True


# Name of environment file to use on remote [Default is 'abaqus_v6_remote.env']
environment_file_name = 'abaqus_v6_remote.env'