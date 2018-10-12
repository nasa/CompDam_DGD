# abaqus-python-addpkg

A command line tool for managing python packages in the python installation used by abaqus/CAE

Abaqus/CAE uses a relatively old version of python (2.7.3) that is separate from other installations of python. While the abaqus python installation includes `numpy` it does not include other python packages that provide helpful functionality in python script for abaqus/CAE, e.g. `scipy`. This script utilizes [anaconda](https://www.anaconda.com/download/) to install a mirror of the abaqus python environment so that conda automatically manages the versions and dependency requirements in order to facilitate installing python packages into the abaqus/CAE python installations.


## Installation
1. Add this script to your path
2. Optionally, add '$env:PATHEXT += ";.py"' to you powershell profile (so you don't have to type the '.py' when calling the script)


## Usage
Python packages can be installed using any of three different methods: conda, pip, or local. In all cases, a conda environment named `abq_2017` is created with the same versions of python and numpy as used by abaqus.

To install a package using conda, try:
```
abaqus-python-addpkg install --conda scipy
```

To install a package using pip, try:
```
abaqus-python-addpkg install --pip plotly
```

To install a local package (equivalent to pip install -e), try:
```
abaqus-python-addpkg install --local path/to/local/package
```
Note, dependencies in local packages must already be installed to the abaqus python environment.


## Tested packages
`scipy`
`matplotlib`
`plotly`


## Caveats
The script was developed and tested using Abaqus 2017; it may not work properly with other versions of Abaqus.

The script requires conda and is currently only for use on Windows. 

While, it's helpful to have access to additional python packages, since the versions of python and numpy that Abaqus relies on are quite old, the compatible versions of other packages are also typically relatively old and may lack some functionality.

WARNING: Since this script modifies files in your abaqus installation, it is possible for the installation to be corrupted if something goes wrong. The `restore` option is provided to revert back to the original, as-installed, configuration by copying the python environment used by the Abaqus solver (which is nominally identicial to the CAE python environment). Users may wish to backup `SIMULIA\CAE\2017\win_b64\tools\SMApy\python2.7\Lib\site-packages` for an added level of redundancy.
