# abaqus-python-addpkg

A command line tool for managing python packages in the python installation used by Abaqus/CAE

Abaqus/CAE uses a relatively old version of python (2.7.3) that is separate from other installations of python. While the Abaqus python installation includes `numpy` it does not include other python packages that provide helpful functionality in python scripts for Abaqus/CAE, e.g. `scipy`. This script utilizes [anaconda](https://www.anaconda.com/download/) to install a mirror of the Abaqus python environment so that anaconda automatically manages the versions and dependency requirements in order to facilitate installing python packages into the Abaqus/CAE python installation.


## Installation
1. Add this script to your path
2. Optionally, add '$env:PATHEXT += ";.py"' to your PowerShell profile (so you don't have to type the '.py' when calling the script)


## Usage
Python packages can be installed using any of three different methods: conda, pip, or local. In all cases, a conda environment named `abq_2017` (or for Abaqus 2018, `abq_2018`) is created with the same versions of python and numpy as used by Abaqus.

To install a package using conda, try:
```
abaqus-python-addpkg install --conda scipy
```
By default, the script assumes you are using Abaqus 2017. If you use Abaqus 2018, specify this directly as:
```
abaqus-python-addpkg -a 2018 install --conda scipy
```

To install a package using pip, try:
```
abaqus-python-addpkg install --pip plotly
```

To install a local package (equivalent to pip install -e), try:
```
abaqus-python-addpkg install --local path/to/local/package
```
Note, dependencies in local packages must already be installed to the Abaqus python environment.


## Tested packages
`scipy`
`matplotlib`
`plotly`


## Caveats
The script was developed and tested using Abaqus 2017. Limited testing has been conducted using Abaqus 2018. The script may not work properly with other versions of Abaqus.

The script requires conda and is currently only for use on Windows.

While it is helpful to have access to additional python packages, the versions of python and numpy that Abaqus relies on are quite old, and the compatible versions of other packages are also typically relatively old and may lack some functionality.

WARNING: Since this script modifies files in your Abaqus installation, it is possible for the installation to be corrupted if something goes wrong. The `restore` option is provided to revert back to the original, as-installed, configuration by copying the python environment used by the Abaqus solver (which is nominally identical to the Abaqus/CAE python environment). Users may wish to backup `SIMULIA\CAE\2017\win_b64\tools\SMApy\python2.7\Lib\site-packages` for an added level of redundancy.
