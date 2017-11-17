"""
Generates a fortran version file. Runs on git-hook post-checkout.
"""

import subprocess


if __name__ == "__main__":
    sha = subprocess.check_output("git rev-parse HEAD", shell=True)
    t = subprocess.check_output("git show -s --format=%ci", shell=True)
    
    with open('for/version.for', 'w') as f:
        f.write('      Module version_Mod\n')
        f.write('        Character(len=40), parameter :: hash = "' + str(sha).strip() + '"\n')
        f.write('        Character(len=50), parameter :: timestamp = "' + str(t).strip() + '"\n')
        f.write('      End Module\n')
   

