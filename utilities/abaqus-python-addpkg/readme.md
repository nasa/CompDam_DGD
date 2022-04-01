# abaqus-python-addpkg

A command line tool for managing python packages in the python installation used by Abaqus/CAE

This tool is for use on Windows.

Abaqus/CAE uses a relatively old version of python (2.7.3) that is separate from other installations of python. While the Abaqus python installation includes `numpy` it does not include other python packages that provide helpful functionality in python scripts for Abaqus/CAE, e.g. `scipy`. This script utilizes [anaconda](https://www.anaconda.com/download/) to install a mirror of the Abaqus python environment so that anaconda automatically manages the versions and dependency requirements in order to facilitate installing python packages into the Abaqus/CAE python installation.

The script is implemented for `python 3`. Using the anaconda prompt or PowerShell with a python 3 environment is recommended.

## One time setup
All that is required for installation is ensuring that you can call the script (`abaqus-python-addpkg.py`) from your shell. To do so, follow these steps:
1. Add the directory containing the script (`abaqus-python-addpkg.py`) to your path. [Here's a helpful tutorial](https://www.architectryan.com/2018/03/17/add-to-the-path-on-windows-10/) on adding directories to your path.
2. Ensure that python is the default program for running `.py` files. In windows explorer, right click on the file `abaqus-python-addpkg.py` and select properties. In the properties menu, if 'Opens with' is not Python, select the `Change...` button and find the python.exe (type `where python` in the anaconda prompt to find the location of python.exe).
You can check that the script is working by typing `abaqus-python-addpkg.py --help`; if it is installed correctly, some usage information will be printed.

Optionally, most shells support recognizing specific file extensions as scripts, so you can configure your shell to allow you to omit the `.py` as well. For example, if you are using PowerShell, add '$env:PATHEXT += ";.py"' to your PowerShell profile. After completing this configuration step the command to call the script is simplified from `abaqus-python-addpkg.py` to `abaqus-python-addpkg`. The instructions below assume that you have not taken this optional step.


## Usage
Python packages can be installed using any of three different methods: conda, pip, or local. In all cases, a conda environment named `abq_2017` (or for Abaqus 2018, `abq_2018`) is created with the same versions of python and numpy as used by Abaqus.

To install a package using conda, try:
```
abaqus-python-addpkg.py install --conda scipy
```
By default, the script uses the version of abaqus called with the `abaqus` command. Alternatively, you can specify the abaqus command directly to call a specific version of abaqus:
```
abaqus-python-addpkg.py --abaqus-cmd abq2018 install --conda scipy
```

To install a package using pip, try:
```
abaqus-python-addpkg.py install --pip plotly
```

To install a local package (equivalent to pip install -e), try:
```
abaqus-python-addpkg.py install --local path/to/local/package
```
Note, dependencies in local packages must already be installed to the Abaqus python environment.

If you receive a `PackagesNotFoundError`, try `conda config --append channels conda-forge` and also `conda config --set restore_free_channel true`

For additional usage information, use the `--help` option.


## Tested packages
`scipy`
`matplotlib`
`plotly`


## Caveats
The script was developed and tested using Abaqus 2017. Limited testing has been conducted using Abaqus 2018. The script may not work properly with other versions of Abaqus.

The script requires conda and is currently only for use on Windows.

The script works with windows command prompt and PowerShell. It may not work properly with other shells.

While it is helpful to have access to additional python packages, the versions of python and numpy that Abaqus relies on are quite old, and the compatible versions of other packages are also typically relatively old and may lack some functionality.

WARNING: Since this script modifies files in your Abaqus installation, it is possible for the installation to be corrupted if something goes wrong. The `restore` option is provided to revert back to the original, as-installed, configuration by copying the python environment used by the Abaqus solver (which is nominally identical to the Abaqus/CAE python environment). Users may wish to backup the CAE python directory (which can be found using the command `python abaqus-python-addpkg.py --show-abaqus-python-directory`) for an added level of redundancy.
