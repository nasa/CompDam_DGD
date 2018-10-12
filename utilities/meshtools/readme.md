# meshtools
This package provides useful tools for python scripting in Abaqus CAE. The source code is located in the directory `meshtools`.

This code package is designed for use in Abaqus/CAE python scripts or in the Abaqus/CAE command line.


## Getting started

### Installation
Meshtools must be installed so that it can be found by the Abaqus CAE interpreter (which is not the same as you system-wide python installation). Each version of Abaqus installed on your machine has it's own python installation, so you will need to repeat this procedure for every version of Abaqus that you want to be able to use `meshtools`. The installation can be automated using [`abaqus-python-addpkg`](../abaqus-python-addpkg/readme.md) with the `local` option.

The following provides a short summary of how to install `meshtools` manually. Clone `meshtools` into a convient directory and cd into the directory. Then, use the abaqus CAE python executable too install the package to your abaqus python `site-packages` directory:
```
$ C:/SIMULIA/CAE/2017/win_b64/tools/SMApy/python2.7/python setup.py develop'
```
The above pertains to standard windows installations of Abaqus 2017. Adjust for other platforms; this package may not work for older versions of Abaqus.

### Usage
To use this code, import the package in parent project python scrips:
```
import meshTools as mt

...
# Other code here to create the model, mesh etc)
...

# Example usage
(e1, e2) = mt.Mesh.sortEdgesForBiasSeed(part=p, edges=e, center=mt.Point(x=0, y=0, z=0))
```
Where `p` is a reference to the part object, `e` is a reference to an `EdgeArray` that specify the edges for which a bias seed should be applied, and `center` is an instance of `Point` that specifies the location toward which the seeds should be biased. The function returns a tuple of two `EdgeArray` in the format required by the native abaqus function `seedEdgeByBias` where `e1` is passed to `end1Edges` and `e2` is passed to `end2Edges`.



## Related libraries
1. Abapy: https://github.com/lcharleux/abapy
