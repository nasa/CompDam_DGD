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
import subprocess
import json

if __name__ == "__main__":

    config = dict()
    config["git"] = False
    config["bash"] = False
    config["abaverify"] = False

    # Check for compatible versions of Python
    assert sys.version_info[0] == 3

    # Check for git
    git_available = False
    try:
        git_output = subprocess.check_output(["git", "--version"]).decode('utf-8')
        print(f"Found git version: {git_output.rstrip()}")
        git_available = True

    except subprocess.CalledProcessError:
        shutil.copy(os.path.join(os.getcwd(), 'for', 'version.for.nogit'), os.path.join(os.getcwd(), 'for', 'version.for'))
        print("Git not found")
        print("NOTICE: If you intend to make contributions to the CompDam code base, using git is required.\n\tTo setup CompDam using git, type:\n\t\tgit clone https://github.com/nasa/CompDam_DGD.git\n\tSee the readme for more details.\n")

    if git_available:
        # Create git hooks
        os.makedirs(os.path.join(os.getcwd(), '.git', 'hooks',), exist_ok=True)
        hooks = ('post-checkout', 'post-commit', 'post-update', 'post-merge')
        for hook in hooks:
            dest = os.path.join(os.getcwd(), '.git', 'hooks', hook)
            if not os.path.exists(dest):
                shutil.copyfile(os.path.join(os.getcwd(), 'etc', 'hook.sample'), dest)
        print("Copied the git hooks to update versioning on git command execution")

        # Set permissions
        if platform.system() == 'Linux':
            for hook in hooks:
                os.chmod(os.path.join(os.getcwd(), '.git', 'hooks', hook), 0o755)
            print("Changed file persions for git hooks")

        # This configures version.for
        os.system('git checkout')
        print("Checked out repository")

        # Successfully set up versioning using git
        config["git"] = True

    # Try bash (for the python extension module)
    try:
        bash_output = subprocess.check_output(["bash", "--version"])
        ver = bash_output.decode('utf-8').split('\n')[0]
        print(f"Found bash version: {ver}")
        config["bash"] = True

    except:
        print("Bash not found")
        print("NOTICE: Bash is required to run the python extension module (useful for development).\n\tIf running windows, install the windows subsystem for linux (WSL).\n\tSee the readme for more information.")



    # Check for the installation of abaverify
    try:
        # imp.find_module('abaverify')
        import abaverify
        print(f"Found abaverify version {abaverify.__version__}")
        config["abaverify"] = True
    except ImportError:
        print("abaverify is failed to import. abaverify is a Python package for running verification tests on Abaqus user subroutines.")
        print("More information can be found here: https://github.com/nasa/abaverify")

    # Write the config file
    with open('etc/config.json', 'w') as f:
        json.dump(config, f, indent=2)

    print("Setup completed successfully")
