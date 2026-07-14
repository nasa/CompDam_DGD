# Python extension module
CompDam can be compiled into a [Python extension module](https://docs.python.org/2/extending/extending.html), which allows many of the Fortran subroutines and functions in the `for` directory to be called from Python. The Python package [`f90wrap`](https://github.com/jameskermode/f90wrap) is used to automatically generate the Python extension modules that interface with the Fortran code. This Python extension module functionality is useful for development and debugging.

## Dependencies and setup
The python extension module requires some additional dependencies. First, the procedure only works on Linux using the bash shell. Windows users can use the [Windows Subsystem for Linux (WSL)](https://docs.microsoft.com/en-us/windows/wsl/about). In addition, `ifort` or `gfortran` 4.6+ is required. Type `gfortran --version` to check if you have this available. The remaining dependencies are python packages and can be installed as follows. The Python extension module has been tested with Python 3.10.9. `f90wrap = 0.2.16` is the most recent version found to be working properly as of April 2026.

Using [Conda](https://conda.io/projects/conda/en/latest/user-guide/getting-started.html) significantly simplifies the setup process, so it is assumed that you have a recent version of Conda available (see the [Conda installation guide](https://conda.io/projects/conda/en/latest/user-guide/install/index.html)). Further, the bash scripts described below include calls to Conda, so they will not work correctly without installing and configuring Conda as follows. Add the Conda-Forge channel by typing:
```
$ conda config --add channels conda-forge
```

Conda stores Python packages in containers called environments. Create a new environment:
```
$ conda create --name compdam python=3.10.9
```
and switch to your new environment:
```
$ source activate compdam
```
which will add `(compdam)` to the prompt. Install `numpy`, `matplotlib`, and `f90wrap` by typing:
```
(compdam) $ conda install numpy matplotlib f90wrap=0.2.16
```
After typing 'y' in response to the prompt asking if you would like to proceed, Conda will install `f90wrap` and all of its dependencies. This completes the setup process. These steps only need to be executed once.

Note, you can exit the Conda environment by typing `source deactivate`. When you open a new session, you will need to activate the Conda environment by typing `source activate compdam`.

## Example usage
This section describes how to compile CompDam into a Python extension module and run a simple example.

The relevant files are in the `pyextmod` directory, so set `pyextmod` as your current working directory.

The bash shell is required. Type `bash` to open a bash shell session if you are using a different shell. Activate the environment in which you have installed the dependencies listed above, e.g. `source activate compdam`.

Next compile the python extension module by executing `make` in the `pyextmod` directory as follows:
```
(compdam) CompDam_DGD/pyextmod> make
```

When you execute `make`, the Fortran modules in the `for` directory are compiled to object files, the shared library `_CompDam_DGD.so` is built, and the Python extension module interface `CompDam_DGD.py` is created. A large amount of output is given in the terminal. After the module is created, most of the functionality of CompDam is available in python with `import CompDam_DGD`.

The file `test_pyextmod_dgdevolve.py` shows an example of how the python extension module can be used. Just as when CompDam_DGD is used in Abaqus, `CompDam_DGD.py` expects `CompDam.parameters` and `props` files (if provided) are located in the working directory. Note the `test_pyextmod_dgdevolve.py` loads a visualization tool that shows how the element deforms as equilibrium is sought by the DGD algorithm.

It is necessary to recompile the CompDam after making changes to the Fortran code. Recompile with the `make` command. It is a good idea to run `make clean` before rerunning `make` to remove old build files.

Note that portions of CompDam that are specific to Abaqus are hidden from `f90wrap` using the preprocessor directive `#ifndef PYEXT`.

## Associated scripts
In the `tests` directory the shell scripts `pyextmod_compile.sh` and `pyextmod_run.sh` are available to help streamline execution of the python extension module. These two scripts assume that Conda environment called `compdam` is available with `abaverify` and `f90wrap`. Both must be executed with the `-i` option. The script `pyextmod_run.sh` loads a debug file and executes the specified DGD routine. The DGD routine and the debug file are specified as arguments as follows:
```
$ bash -i pyextmod_run.sh dgdevolve <job-name>
```
The `<job-name>` is the Abaqus job name where it is assumed that the debug file resides in the testOutput folder with the name `job-name-1-debug-0.py`

In the `pyextmod` directory, the `helpers.py` file includes logic to load debug.py files.
