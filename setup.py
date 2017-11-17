"""
Setup script
Run me after cloning the CompDam_DGD repo to finish the initial setup:
  1. Hooks for keeping for/version.for up to date
  2. ... Add other features here ...
"""

import platform
import shutil
import os
import sys
import imp
import pip

if __name__ == "__main__":

    # Check for compatible versions of Python
    if sys.version_info[0:2] != (2,7):
        raise Exception("setup.py run with Python version " + '.'.join(str(s) for s in sys.version_info[0:3]) + ". Version 2.7 is required to setup and run CompDam_DGD.")

    # Create git hooks
    shutil.copyfile(os.path.join(os.getcwd(), 'etc', 'hook.sample'), os.path.join(os.getcwd(), '.git', 'hooks', 'post-checkout'))
    shutil.copyfile(os.path.join(os.getcwd(), 'etc', 'hook.sample'), os.path.join(os.getcwd(), '.git', 'hooks', 'post-commit'))
    shutil.copyfile(os.path.join(os.getcwd(), 'etc', 'hook.sample'), os.path.join(os.getcwd(), '.git', 'hooks', 'post-update'))
    shutil.copyfile(os.path.join(os.getcwd(), 'etc', 'hook.sample'), os.path.join(os.getcwd(), '.git', 'hooks', 'post-merge'))

    # Set permissions
    if platform.system() == 'Linux':
        os.chmod(os.path.join(os.getcwd(), '.git', 'hooks', 'post-checkout'), 0775)
        os.chmod(os.path.join(os.getcwd(), '.git', 'hooks', 'post-commit'), 0775)
        os.chmod(os.path.join(os.getcwd(), '.git', 'hooks', 'post-update'), 0775)
        os.chmod(os.path.join(os.getcwd(), '.git', 'hooks', 'post-merge'), 0775)


    os.system('git checkout')

    # Check for the installation of abaverify
    try:
        imp.find_module('abaverify')
    except ImportError:
        print """
abaverify is not installed. abaverify is a Python package for running verification tests on Abaqus user subroutines.
More information can be found here: https://github.com/nasa/abaverify
        """
        install_abaverify = raw_input("Would you like to install abaverify now? (yes/no) ")
        if 'y' in install_abaverify.lower():
            pip.main(['install', '--src', os.pardir, '-e','git+https://github.com/nasa/abaverify#egg=abaverify', '--user'])


    print "Setup completed successfully"
