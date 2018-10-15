# CompDam - Deformation Gradient Decomposition (DGD)
This code is a continuum damage mechanics (CDM) material model intended for use with the Abaqus finite element code. This is a research code which aims to provide an accurate representation of mesoscale damage modes in fiber-reinforced polymer composite materials in finite element models in which each ply is discretely represented.

The CDM material model is implemented as an Abaqus/Explicit user subroutine (VUMAT) for the simulation of matrix cracks formed under tensile, compressive, and shear loading conditions and fiber fracture under tensile and compressive loading. Within CompDam, the emphasis of many recent developments has been on accurately representing the kinematics of composite damage. The kinematics of matrix cracks are represented by treating them as cohesive surfaces embedded in a deformable bulk material in accordance with the Deformation Gradient Decomposition (DGD) approach. Fiber tensile damage is modeled using conventional CDM strain-softening.

This software may be used, reproduced, and provided to others only as permitted under the terms of the agreement under which it was acquired from the U.S. Government. Neither title to, nor ownership of, the software is hereby transferred. This notice shall remain on all copies of the software.

Copyright 2016 United States Government as represented by the Administrator of the National Aeronautics and Space Administration. No copyright is claimed in the United States under Title 17, U.S. Code. All Other Rights Reserved.

Publications that describe the theories used in this code:
- Cheryl A. Rose, et al. ["Analysis Methods for Progressive Damage of Composite Structures"](https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20140001002.pdf) NASA/TM-2013-218024, July 2013.
- Frank A. Leone Jr. ["Deformation gradient tensor decomposition for representing matrix cracks in fiber-reinforced materials"](http://dx.doi.org/10.1016/j.compositesa.2015.06.014) *Composites Part A* (2015) **76**:334-341.
- Frank A. Leone Jr. ["Representing matrix cracks through decomposition of the deformation gradient tensor in continuum damage mechanics methods"](http://iccm20.org/fullpapers/file?f=Abk7n4gkWV) *Proceedings of the 20th International Conference on Composite Materials*, Copenhagen, Denmark, 19-24 July 2015.
- Andrew C. Bergan, et al., ["Development of a Mesoscale Finite Element Constitutive Model for Fiber Kinking"](https://arc.aiaa.org/doi/10.2514/6.2018-1221) *AIAA SciTech Forum*, Kissimmee, Florida, 8-12 January 2018.

Examples of this code being applied can be found in the following publications:
- Mark McElroy, et al. ["Simulation of delamination-migration and core crushing in a CFRP sandwich structure"](https://doi.org/10.1016/j.compositesa.2015.08.026) *Composites Part A* (2015) **79**:192-202.
- Frank A. Leone Jr., et al. ["Fracture-Based Mesh Size Requirements for Matrix Cracks in Continuum Damage Mechanics Models"](https://doi.org/10.2514/6.2017-0198) *AIAA SciTech Forum*, Grapevine, Texas, 9-13 January 2017.
- Imran Hyder, et al. ["Assessment of Intralaminar Progressive Damage and Failure Analysis Using an Efficient Evaluation Framework"](https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20170010326.pdf) *32nd American Society for Composites (ASC) Annual Technical Conference*, West Lafayette, Indiana, 22-25 October 2017.
- Kyongchan Song, et al. ["Continuum Damage Mechanics Models for the Analysis of Progressive Damage in Cross-Ply and Quasi-Isotropic Panels Subjected to Static Indentation"](https://arc.aiaa.org/doi/10.2514/6.2018-1466) *AIAA SciTech Forum*, Kissimmee, Florida, 8-12 January 2018.


For any questions, please contact the developers:
- Frank Leone   | [frank.a.leone@nasa.gov](mailto:frank.a.leone@nasa.gov)     | (W) 757-864-3050
- Andrew Bergan | [andrew.c.bergan@nasa.gov](mailto:andrew.c.bergan@nasa.gov) | (W) 757-864-3744
- Carlos Dávila | [carlos.g.davila@nasa.gov](mailto:carlos.g.davila@nasa.gov) | (W) 757-864-9130


## Table of contents
- [Getting started](#getting-started)
- [Model features](#model-features)
- [Elements](#elements)
- [Material properties](#material-properties)
- [State variables](#state-variables)
- [Example problems](#example-problems)
- [Advanced debugging](#advanced-debugging)
- [Python extension module](#python-extension-module)
- [Summary of tests classes](#summary-of-tests-classes)
- [Contributing](#contributing)
- [Citing CompDam](#citing-compdam)


## Getting started

### Source code
The user subroutine source code is located in the `for` directory. The main entry point is `CompDam_DGD.for`.

### Prerequisites
[Intel Fortran Compiler](https://software.intel.com/en-us/fortran-compilers) version 11.1 or newer is required to compile the code ([more information about compiler versions](usersubroutine-prerequisites.md)). It is recommended that Abaqus 2016 or newer is used with this code. Current developments and testing are conducted with Abaqus 2017. Python supporting files require Python 2.7.

### Initial setup
After cloning the CompDam_DGD git repository, it is necessary to run the setup script file `setup.py` located in the repository root directory:
```
$ python setup.py
```

### Abaqus environment file settings
The `abaqus_v6.env` file must have [`/fpp`](https://software.intel.com/en-us/node/678201), [`/Qmkl:sequential`](https://software.intel.com/en-us/node/678038), and [`/free`](https://software.intel.com/en-us/node/678227) in the `ifort` command where the format for Windows is used. The corresponding Linux format is: `-fpp`, `-free`, and `-mkl=sequential`. The `/fpp` option enables the Fortran preprocessor, which is required for the code to compile correctly. The `/free` option sets the compiler to free-formatting for the source code files. The `/Qmkl:sequential` enables the [Intel Math Kernel Library (MKL)](https://software.intel.com/en-us/intel-mkl), which provides optimized and verified functions for many mathematical operations. The MKL is used in this code for calculating eigenvalues and eigenvectors.

A sample environment file is provided in the `tests` directory for Windows and Linux systems.

### Submitting a job
This code is an Abaqus/Explicit VUMAT. Please refer to the Abaqus documentation for the general instructions on how to submit finite element analyses using user subroutines. Please see the [example input file statements](#example-input-file-statements) for details on how to interface with this particular VUMAT subroutine.

Analyses with this code **must** be run in double precision. Some of the code has double precision statements and variables hard-coded, so if Abaqus/Explicit is run in single precision, compile-time errors will arise. When submitting an Abaqus/Explicit job from the command line, double precision is specified by including the command line argument `double=both`.

For example, run the test model `test_C3D8R_elastic_fiberTension` in the `tests` directory with the following command:
```
$ abaqus job=test_C3D8R_elastic_fiberTension user=../for/CompDam_DGD.for double=both
```

### Example input file statements
Example 1, using an [external material properties file](#defining-the-material-properties-in-a-props-file):

<pre>
*Section controls, name=control_name, distortion control=YES, element deletion=YES
**
*Material, name=IM7-8552
*Density
 1.57e-09,
*Depvar, delete=11
** the above delete statement is optional
  19,
  1, CDM_d2
  2, CDM_Fb1
  3, CDM_Fb2
  4, CDM_Fb3
  5, CDM_B
  6, CDM_Lc1
  7, CDM_Lc2
  8, CDM_Lc3
  9, CDM_FIm
 10, CDM_alpha
 11, CDM_STATUS
 12, CDM_Plas12
 13, CDM_Inel12
 14, CDM_FIfT
 15, CDM_slide1
 16, CDM_slide2
 17, CDM_FIfC
 18, CDM_d1T
 19, CDM_d1C
*Characteristic Length, definition=USER, components=3
*User material, constants=3
** 1              2  3
** feature flags,  , thickness
          100001,  ,       0.1
**
*Initial Conditions, type=SOLUTION
 elset_name,  0.d0,  0.d0,  0.d0,  0.d0,  0.d0,  0.d0,  0.d0,
 0.d0,  0.d0,  -999,     1,  0.d0,  0.d0,  0.d0, 0.d0,
 0.d0,  0.d0,  0.d0,  0.d0
** In each step, NLGEOM=YES must be used. This is the default setting.
</pre>

Example 2, using an [input deck command](#defining-the-material-properties-in-the-input-deck):

<pre>
*Section controls, name=control_name, distortion control=YES, element deletion=YES
**
*Material, name=IM7-8552
*Density
 1.57e-09,
*Depvar, delete=11
** the above delete statement is optional
  19,
  1, CDM_d2
  2, CDM_Fb1
  3, CDM_Fb2
  4, CDM_Fb3
  5, CDM_B
  6, CDM_Lc1
  7, CDM_Lc2
  8, CDM_Lc3
  9, CDM_FIm
 10, CDM_alpha
 11, CDM_STATUS
 12, CDM_Plas12
 13, CDM_Inel12
 14, CDM_FIfT
 15, CDM_slide1
 16, CDM_slide2
 17, CDM_FIfC
 18, CDM_d1T
 19, CDM_d1C
*Characteristic Length, definition=USER, components=3
*User material, constants=40
** 1              2  3          4  5  6  7  8
** feature flags,  , thickness, 4, 5, 6, 7, 8
          100001,  ,       0.1,  ,  ,  ,  ,  ,
**
**  9         10        11        12        13        14        15        16
**  E1,       E2,       G12,      nu12,     nu23,     YT,       SL        GYT,
    171420.0, 9080.0,   5290.0,   0.32,     0.52,     62.3,     92.30,    0.277,
**
**  17        18        19        20        21        22        23        24
**  GSL,      eta_BK,   YC,       alpha0    E3,       G13,      G23,      nu13,
    0.788,    1.634,    199.8,    0.925,      ,          ,         ,          ,
**
**  25        26        27        28        29        30        31        32
**  alpha11,  alpha22,  alpha_PL, n_PL,     XT,       fXT,      GXT,      fGXT,
    -5.5d-6,  2.58d-5,          ,     ,     2326.2,   0.2,      133.3,    0.5,
**
**  33        34        35        36        37        38        39        40
**  XC,       fXC,      GXC,      fGXC,       cf,     w_kb,     None,     mu
    1200.1,      ,         ,          ,         ,     0.1,          ,     0.3
**
*Initial Conditions, type=SOLUTION
 elset_name,  0.d0,  0.d0,  0.d0,  0.d0,  0.d0,  0.d0,  0.d0,
 0.d0,  0.d0,  -999,     1,  0.d0,  0.d0,  0.d0, 0.d0,
 0.d0,  0.d0,  0.d0,  0.d0
** In each step, NLGEOM=YES must be used. This is the default setting.
</pre>

### Running tests
Test cases are available in the `tests` directory. The tests are useful for demonstrating the capabilities of the VUMAT as well as to verify that the code performs as intended. Try running some of the test cases to see how the code works. The test cases can be submitted as a typical Abaqus job using the Abaqus command line arguments.

### Building a shared library
CompDam_DGD can be built into a shared library file. Follow these steps:
1. Place a copy of the Abaqus environment file (with the compiler flags specified) in the `for` directory
2. In Linux, and when using Abaqus versions prior to 2017, rename `CompDam_DGD.for` to `CompDam_DGD.f`
3. From the `for` directory, run:
```
$ abaqus make library=CompDam_DGD
```
This command will create shared libraries for the operating system it is executed on (`.dll` for Windows and `.so` for Linux).

When using a pre-compiled shared library, it is only necessary to specify the location of the shared library files in the environment file (the compiler options are not required). To run an analysis using a shared library, add `usub_lib_dir = <full path to shared library file>` to the Abaqus environment file in the Abaqus working directory.

## Model features
The CompDam_DGD material model implements a variety of features that can be enabled or disabled by the user. An overview of these features is provided in this section. The material properties required for each feature are listed. References are provided to more detailed discussions of the theoretical framework for each feature.

### Fully orthotropic elasticity
The composite materials modeled with CompDam_DGD can be defined assuming either [transverse isotropy](http://www.efunda.com/formulae/solid_mechanics/mat_mechanics/hooke_iso_transverse.cfm) or [orthotropy](http://www.efunda.com/formulae/solid_mechanics/mat_mechanics/hooke_orthotropic.cfm). For a transversely isotropic material definition, the following properties must be defined: E1, E2, G12, v12, and v23. For an orthotropic material definition, the following additional properties must be defined: E2, G13, G23, and nu13.

### Matrix damage
Tensile and compressive matrix damage is modeled by embedding cohesive laws to represent cracks in the material according to the deformation gradient decomposition method of [Leone (2015)](http://doi.org/10.1016/j.compositesa.2015.06.014). The matrix crack normals can have any orientation in the 2-3 plane, defined by the angle `CDM_alpha`. The mixed-mode behavior of matrix damage initiation and evolution is defined according to the Benzeggagh-Kenane law. The initiation of compressive matrix cracks accounts for friction on the potential crack surface according to the LaRC04 failure criteria. In the notation of the [paper](http://doi.org/10.1016/j.compositesa.2015.06.014), `Q` defines the material direction that is severed by the crack. In this implementation `Q=2` except when `CDM_alpha = 90`.

The following material properties are required for the prediction of matrix damage: YT, SL, GYT, GSL, eta_BK, YC, and alpha0. The state variables related to matrix damage are `CDM_d2`, `CDM_FIm`, `CDM_B`, `CDM_alpha`, `CDM_Fb1`, `CDM_Fb2`, and `CDM_Fb3`.

### Thermal strains
The thermal strains are calculated by multiplying the 1-, 2-, and 3-direction coefficients of thermal expansion by the current &Delta;T, as provided by the Abaqus solver. The thermal strains are subtracted from the current total strain.

The required material properties are the coefficients of thermal expansion in the 1 and 2 directions. It is assumed that the 2- and 3-direction coefficients of thermal expansion are equal.

Hygroscopic strains are not accounted for. If the effects of hygroscopic expansion are to be modeled, it is recommended to smear the hygroscopic and thermal expansion coefficients to approximate the response using the solver-provided &Delta;T.

### Shear nonlinearity
Two approaches to modeling the matrix nonlinearity are available: Ramberg-Osgood plasticity and Schapery theory. These two methods are mutually exclusive and optional.

#### Ramberg-Osgood plasticity
Shear nonlinearity in the 1-2 and/or the 1-3 plane can be modeled using the [Ramberg-Osgood equation](https://en.wikipedia.org/wiki/Ramberg%E2%80%93Osgood_relationship), with its parameters selected to fit experimental data. As applied herein, the Ramberg-Ogsood equation is written in the following form for the 1-2 plane:

*&gamma;*<sub>12</sub> = [*&tau;*<sub>12</sub> + *&alpha;*<sub>PL</sub>sign(*&tau;*<sub>12</sub>)|*&tau;*<sub>12</sub>|<sup>*n*<sub>PL</sub></sup>]/*G*<sub>12</sub>

where *&gamma;*<sub>12</sub> is the shear strain and *&tau;*<sub>12</sub> is the shear stress. Likewise, the expression for the 1-3 plane is

*&gamma;*<sub>13</sub> = [*&tau;*<sub>13</sub> + *&alpha;*<sub>PL</sub>sign(*&tau;*<sub>13</sub>)|*&tau;*<sub>13</sub>|<sup>*n*<sub>PL</sub></sup>]/*G*<sub>13</sub>

Prior to the initiation of matrix damage (i.e., `CDM_d2 = 0`), the nonlinear shear response due to the above equation is plastic, and the unloading/reloading slope is unchanged. No pre-peak nonlinearity is applied to the matrix tensile or compressive responses (i.e., *&sigma;<sub>22</sub>*).

The required material inputs are the two parameters in the above equation: *&alpha;*<sub>PL</sub> and *n*<sub>PL</sub>. Note that the same constants are used for the 1-2 and 1-3 planes under the assumption of transverse isotropy (see [Seon et al. 2017](http://dpi-proceedings.com/index.php/asc32/article/view/15267)). For the 1-2 plane, the state variables `CDM_Plas12` and `CDM_Inel12` are used to track the current plastic shear strain and the total amount of inelastic plastic shear strain that has occurred through the local deformation history, respectively. For cases of monotonic loading, `CDM_Plas12` and `CDM_Inel12` should have the same magnitude. Likewise, the state variables `CDM_Plas13` and `CDM_Inel13` are utilized for the 1-3 plane. The [feature flags](#contrlling-which-features-are-enabled) can be used to enable this Ramberg-Osgood model in the 1-2 plane, 1-3 plane, or both planes.

#### Schapery micro-damage
Matrix nonlinearity in the 1-2 plane can also be modeled using Schapery theory, in which all pre-peak matrix nonlinearity is attributed to the initiation and development of micro-scale matrix damage. With this approach, the local stress/strain curves will unload to the origin, and not develop plastic strain. A simplified version of the approach of [Pineda and Waas](https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20120000914.pdf) is here applied. The micro-damage functions *e<sub>s</sub>* and *g<sub>s</sub>* are limited to third degree polynomials for ease of implementation. As such, four fitting parameters are required for each of *e<sub>s</sub>* and *g<sub>s</sub>* to define the softening of the matrix normal and shear responses to micro-damage development.

*e<sub>s</sub>*(*S<sub>r</sub>*) = *e<sub>s0</sub>* + *e<sub>s1</sub>S<sub>r</sub>* + *e<sub>s2</sub>S<sub>r</sub>*<sup>2</sup> + *e<sub>s3</sub>S<sub>r</sub>*<sup>3</sup>

*g<sub>s</sub>*(*S<sub>r</sub>*) = *g<sub>s0</sub>* + *g<sub>s1</sub>S<sub>r</sub>* + *g<sub>s2</sub>S<sub>r</sub>*<sup>2</sup> + *g<sub>s3</sub>S<sub>r</sub>*<sup>3</sup>

where *S<sub>r</sub>* is the micro-damage reduced internal state variable. These eight material properties must be defined in an [external material properties file](#defining-the-material-properties-in-a-props-file). *S<sub>r</sub>* is stored in the 12<sup>th</sup> state variable slot, replacing `CDM_Plas12`, when Schapery theory is used in a model. The 13<sup>th</sup> state variable slot is not used when Schapery micro-damage is used.

#### Schaefer
Shear nonlinearity in the 1-2 plane can be modeled using the Schaefer prepeak model. In this model, effective plastic strain is related to effective plastic stress through the power law:

*&epsilon;*<sub>plastic</sub> = *A &sigma;*<sup>*n*</sup> 

Additionally, a yield criterion function (effective stress) is defined as:

*f* = *&sigma;* = (*S*<sub>22</sub> + *a*<sub>6</sub> *S*<sub>12</sub><sup>2</sup>)<sup>1/2</sup> + *b*<sub>2</sub> *S*<sub>22</sub>

where *a*<sub>6</sub>, *b*<sub>2</sub>, *A* and *n* are material constants needed for the Schaefer prepeak material model. These four material properties must be defined in an [external material properties file](#defining-the-material-properties-in-a-props-file). Additionally, when definining  *a*<sub>6</sub>, *b*<sub>2</sub>, *A* and *n* in the external material properties file the variables are prefixed with schaefer_ (to disambiguate the otherwise nondescript material property names and symbols). 

The above two equations are used in concert to determine plastic strain through the relationship:

*&epsilon;* <sub>plastic</sub> = *n A f* <sup>*n* - 1</sup> &part;*f* / &part;*S*<sub>*i*</sub> &part;*f* / &part;*S*<sub>j</sub>

*f* (i.e., schaefer_f) and the tensorial plastic strain determined by the nonlinearity model are stored as state variables (27 through 32 for plastic strain and 33 for *f*)

### Fiber tensile damage
A continuum damage mechanics model similar to the work of [Maimí et al. (2007)](http://doi.org/10.1016/j.mechmat.2007.03.005) is used to model tensile fiber damage evolution. The model utilizes a non-interacting maximum strain failure criterion, and bilinear softening after the initiation of failure. The area under the stress-strain curve is equal to the fracture toughness divided by the element length normal to fracture, i.e., `CDM_Lc1`. The required material properties are: XT, fXT, GXT, and fGXT, where fXT and fGXT are ratios of strength and fracture toughness for bilinear softening, defined as the *n* and *m* terms in equations (25) and (26) of [Dávila et al. (2009)](http://doi.org/10.1007/s10704-009-9366-z). To model a linear softening response, both fXT and fGXT should be set equal to 0.5.

### Fiber compression damage

#### Model 1: Max strain, bilinear softening (BL)
Same model as in tension, but for compression. Assumes maximum strain failure criterion and bilinear softening. The required material properties are: XC, fXC, GXC, and fGXC.

Load reversal assumptions from [Maimí et al. (2007)](http://doi.org/10.1016/j.mechmat.2007.03.006).

#### Model 2: Placeholder

#### Model 3: Fiber kinking theory (FKT)
A model based on Budiansky's fiber kinking theory from [Budiansky 1983](https://doi.org/10.1016/0045-7949(83)90141-4), [Budiansky and Fleck 1993](https://doi.org/10.1016/0022-5096(93)90068-Q), and [Budiansky et al. 1998](https://doi.org/10.1016/S0022-5096(97)00042-2) implemented using the DGD framework. The model is described in detail in Bergan et al. 2018 [SciTech paper]. The model accounts for fiber kinking due to shear instability by considering an initial fiber misalignment, nonlinear shear stress-strain behavior via Ramberg-Osgood, and geometric nonlinearity.

The required material properties are: XC, YC, wkb, alpha0, alpha_PL, and n_PL. The [feature flag](#controlling-which-features-are-enabled) for fiber compression should be set to '3' to activate this model feature. This feature requires all 26 state variables to be defined and initialized. The relevant state variables are:
- `CDM_phi0`: initial fiber misalignment (radians).
- `CDM_gamma`: rotation of the fibers due to loading (radians). The current fiber misalignment is `CDM_phi0 + CDM_gamma`.
- `CDM_FIfC`: failure index for fiber kinking (0 to 1). The calculation for this failure index is simplistic and is accurate for unconstrained elements. It is possible for cases to arise where `CDM_FIfC = 1` but the element is not kinked. A better definition of the failure index for fiber kinking is needed. Interpret the value of `CDM_FIfC` with caution.
- `CDM_Fmi`: the components of the first column of `Fm` used for decomposing the element where `i=1,2,3`.

The initial condition for the state variable `CDM_phi0` determines the initial fiber misalignment as described in [initial conditions](#initial-conditions).

The fiber kinking theory model implemented here is preliminary and has some known shortcomings and caveats:
- The model has only been tested for C3D8R. Limited application with C3D6 demonstrated issues. No testing has been completed for other element types.
- The interaction of this model with matrix cracking has not been tested.
- No effort has been made to model progressive crushing.
- The model does not account for out-of-plane kinking.
- Other mechanisms of fiber compressive failure (e.g., shear driven fiber breaks) are not accounted for. An outcome of this is that the model predicts the material does not fail if shear deformation is fully constrained.
- No special consideration for load reversal has been included.

Relevant example problems: UNC0_C3D8R_FKT

### Friction
Friction is modeled on the damaged fraction of the cross-sectional area of DGD cracks using the approach of [Alfano and Sacco (2006)](http://doi.org/10.1002/nme.1728). The coefficient of friction *&mu;* must be defined to account for friction on the failed crack surface.

The amount of sliding which has taken place in the longitudinal and transverse directions are stored in state variables `CDM_slide1` and `CDM_slide2`, respectively.

### Strain definition
The strain is calculated using the deformation gradient tensor provided by the Abaqus solver. The default strain definition used is the Green-Lagrange strain:

*E* = (*F*<sup>T</sup>*F* - *I*)/2

Hooke's law is applied using the Green-Lagrange strain to calculate the 2<sup>nd</sup> Piola-Kirchhoff stress *S*.

### Fiber nonlinearity
Nonlinear elastic behavior in the fiber direction can be introduced with the material property c<sub>*l*</sub>. The expression used follows [Kowalski 1988](https://www.astm.org/DIGITAL_LIBRARY/STP/PAGES/STP26136S.htm):

*E<sub>1</sub>* = *E<sub>1</sub>*(1 + c<sub>*l*</sub>*&epsilon;*<sub>11</sub>)

By default, fiber nonlinearity is disabled by setting c<sub>*l*</sub> = 0.


## Elements
CompDam_DGD has been developed and tested using the Abaqus three-dimensional, reduced-integration `C3D8R` solid and `S4R` shell elements. Limited testing has been performed using the `CPS4R` plane stress element, the `SC8R` continuum shell element, and the fully-integrated `C3D8` solid element.

Because CompDam_DGD is a material model, it is expected to be compatible with structural elements generally. However, users are advised to perform tests with any previously untested element types before proceeding to use CompDam_DGD in larger structural models.


## Material properties
A set of material properties must be defined for the material of interest. This section describes how to specify the material properties.

### Defining material properties
Material properties can be defined in the input deck or in a separate `.props` file. Definition of the material properties in a `.props` file is more convenient and generally preferred since it isolates the material properties from the structural model definition.

#### Defining the material properties in a `.props` file
Using a `.props` file is a versatile means of defining material properties. The subroutine looks for a file named as `jobName_materialName` where the job name is the name of the Abaqus job (default is input deck name) and the material is name assigned as `*Material, Name=materialName` in the input deck. If no file is found, then the subroutine looks for `materialName.props`. The `.props` file must be located in the Abaqus working directory.

The `.props` should contain one property per line with the format `[NAME] = [VALUE]` where the name is symbolic name for the property and the value is assigned value for the property. Blank lines or commented text (denoted by `//`) is ignored. See the [table of material properties](#table-of-material-properties) for a complete list of material property symbolic names, acceptable values, and recommended test methods for characterizing the properties.

When a `.props` is used to define the material properties, the feature flags and thickness still must be defined in the input deck.

#### Defining the material properties in the input deck
Material properties can be defined in the input deck. Any optional material property can be left blank and the corresponding feature(s) will be disabled. The ordering of the material properties for the input deck definition is given in the first (#) column of the [table of material properties](#table-of-material-properties).

#### Table of material properties
| #  |         Symbol         |   Name   |                  Description                  |                       Units                    |            Admissible values             | Recommended Test |
|----|------------------------|----------|-----------------------------------------------|------------------------------------------------|------------------------------------------|-----------------|
|  9 | *E<sub>1</sub>*        | E1       | Longitudinal Young's modulus                  | F/L<sup>2</sup>                                | 0 < *E<sub>1</sub>* < &infin;            | ASTM D3039      |
| 10 | *E<sub>2</sub>*        | E2       | Transverse Young's modulus                    | F/L<sup>2</sup>                                | 0 < *E<sub>2</sub>* < &infin;            | ASTM D3039      |
| 11 | *G<sub>12</sub>*       | G12      | In-plane Shear modulus                        | F/L<sup>2</sup>                                | 0 < *G<sub>12</sub>* < &infin;           | ASTM D3039      |
| 12 | *&nu;<sub>12</sub>*    | nu12     | Major Poisson's ratio                         | -                                              | 0 &le; *&nu;<sub>12</sub>* &le; 1        | ASTM D3039      |
| 13 | *&nu;<sub>23</sub>*    | nu23     | Minor Poisson's ratio                         | -                                              | 0 &le; *&nu;<sub>23</sub>* &le; 1        |                 |
|    | ===                    |          |                                               |                                                |                                          |                 |
| 14 | *Y<sub>T</sub>*        | YT       | Transverse tensile strength                   | F/L<sup>2</sup>                                | 0 < *Y<sub>T</sub>* < &infin;            | ASTM D3039      |
| 15 | *S<sub>L</sub>*        | SL       | Shear strength                                | F/L<sup>2</sup>                                | 0 < *S<sub>L</sub>* < &infin;            |                 |
| 16 | *G<sub>Ic</sub>*       | GYT      | Mode I fracture toughness                     | F/L                                            | 0 < *G<sub>Ic</sub>* < &infin;           | ASTM D5528      |
| 17 | *G<sub>IIc</sub>*      | GSL      | Mode II fracture toughness                    | F/L                                            | 0 < *G<sub>IIc</sub>* < &infin;          | ASTM D7905      |
| 18 | *&eta;*                | eta_BK   | BK exponent for mode-mixity                   | -                                              | 1 &le; *&eta;* < &infin;                 |                 |
| 19 | *Y<sub>C</sub>*        | YC       | Transverse compressive strength               | F/L<sup>2</sup>                                | 0 < *Y<sub>C</sub>* < &infin;            | ASTM D3410      |
| 20 | *&alpha;<sub>0</sub>*  | alpha0   | Fracture plane angle for pure trans. comp.    | Radians                                        | 0 &le; *&alpha;<sub>0</sub>* &le; &pi;/2 |                 |
|    | ------                 |          |                                               |                                                |                                          |                 |
| 21 | *E<sub>3</sub>*        | E3       | 3-direction Young's modulus                   | F/L<sup>2</sup>                                | 0 < *E<sub>3</sub>* < &infin;            |                 |
| 22 | *G<sub>13</sub>*       | G13      | Shear modulus in 1-3 plane                    | F/L<sup>2</sup>                                | 0 < *G<sub>13</sub>* < &infin;           |                 |
| 23 | *G<sub>23</sub>*       | G23      | Shear modulus in 1-2 plane                    | F/L<sup>2</sup>                                | 0 < *G<sub>23</sub>* < &infin;           |                 |
| 24 | *&nu;<sub>13</sub>*    | nu13     | Poisson's ratio in 2-3 plane                  | -                                              | 0 &le; *&nu;<sub>13</sub>* &le; 1        |                 |
|    | ------                 |          |                                               |                                                |                                          |                 |
| 25 | *&alpha;<sub>11</sub>* | alpha11  | Coefficient of long. thermal expansion        | /&deg;                                         | -1 &le; *&alpha;<sub>11</sub>* &le; 1    |                 |
| 26 | *&alpha;<sub>22</sub>* | alpha22  | Coefficient of tran. thermal expansion        | /&deg;                                         | -1 &le; *&alpha;<sub>22</sub>* &le; 1    |                 |
|    | ------                 |          |                                               |                                                |                                          |                 |
| 27 | *&alpha;<sub>PL</sub>* | alpha_PL | Nonlinear shear parameter                     | (F/L<sup>2</sup>)<sup>1-*n*<sub>PL</sub></sup> | 0 &le; *&alpha;<sub>PL</sub>* < &infin;  |                 |
| 28 | *n<sub>PL</sub>*       | n_PL     | Nonlinear shear parameter                     | -                                              | 0 &le; *n<sub>PL</sub>* < &infin;        |                 |
|    | ------                 |          |                                               |                                                |                                          |                 |
| 29 | *X<sub>T</sub>*        | XT       | Long. tensile strength                        | F/L<sup>2</sup>                                | 0 < *X<sub>T</sub>* < &infin;            | ASTM D3039      |
| 30 | *f<sub>XT</sub>*       | fXT      | Long. tensile strength ratio                  | -                                              | 0 &le; *f<sub>XT</sub>* &le; 1           |                 |
| 31 | *G<sub>XT</sub>*       | GXT      | Long. tensile fracture toughness              | F/L                                            | 0 < *G<sub>XT</sub>* < &infin;           |                 |
| 32 | *f<sub>GXT</sub>*      | fGXT     | Long. tensile fracture toughness ratio        | -                                              | 0 &le; *f<sub>GXT</sub>* &le; 1          |                 |
|    | ------                 |          |                                               |                                                |                                          |                 |
| 33 | *X<sub>C</sub>*        | XC       | Long. compressive strength                    | F/L<sup>2</sup>                                | 0 < *X<sub>C</sub>* < &infin;            | ASTM D3410      |
| 34 | *f<sub>XC</sub>*       | fXC      | Long. compression strength ratio              | -                                              | 0 &le; *f<sub>XC</sub>* &le; 1           |                 |
| 35 | *G<sub>XC</sub>*       | GXC      | Long. compression fracture toughness          | F/L                                            | 0 < *G<sub>XC</sub>* < &infin;           |                 |
| 36 | *f<sub>GXC</sub>*      | fGXC     | Long. compression fracture toughness ratio    | -                                              | 0 &le; *f<sub>GXC</sub>* &le; 1          |                 |
|    | ------                 |          |                                               |                                                |                                          |                 |
| 37 | *c<sub>l</sub>*        | cl       | Fiber nonlinearity coefficient                | -                                              | 0 &le; *c<sub>l</sub>* 1000              |                 |
| 38 | *w<sub>kb</sub>*       | w_kb     | Width of the kink band                        | L                                              | 0 &le; *w<sub>kb</sub>* &infin;          |                 |
| 40 | *&mu;*                 | mu       | Coefficient of friction                       | -                                              | 0 &le; *&mu;* &le; 1                     |                 ||

Notes:
- The first five inputs (above the ===) are required
- Properties for each feature are grouped and separated by ------
- The number in the first column corresponds to the property order when defined in the input deck
- Properties not listed in the table above (see next section):
  1. [feature flags](#controlling-which-features-are-enabled)
  2. *reserved*
  3. [thickness](#definition-of-thickness)
  4. *reserved*
  5. *reserved*
  6. *reserved*
  7. *reserved*
  8. *reserved*
- &infin; is calculated with the Fortran intrinsic `Huge` for double precision
- In the event that both a `.props` file is found and material properties are specified in the input deck (`nprops > 8`), then the material properties from the input deck are used and a warning is used.

### Required inputs for the `*Material` data lines in the input deck
The feature flags and thickness are defined in the input deck on the material property data lines. These properties must be defined in the input deck whether the other material properties are defined via the .props file or via the input deck. While feature flags and thickness are not material properties per se, they are used in controlling the behavior of the material model.

#### Controlling which features are enabled
Model features can be enabled or disabled by two methods. The first method is specifying only the material properties required for the features you would like to enable. CompDam_DGD disables any feature for which all of the required material properties have not been assigned. If an incomplete set of material properties are defined for a feature, a warning is issued.

The second method is by specifying the status of each feature directly as a material property in the input deck. Each feature of the subroutine is controlled by a position in an integer, where 0 is disabled and 1 is enabled.
The positions correspond to the features as follows:
- Position 1: Matrix damage
- Position 2: Shear nonlinearity (1=Ramberg-Osgood 1-2 plane, 2=Schapery, 3=Ramberg-Osgood 3-D, 4=Ramberg-Osgood 1-3 plane, 5=Schaefer || more information [here](#shear-nonlinearity))
- Position 3: Fiber tensile damage
- Position 4: Fiber compression damage (1=max strain, 2=N/A, 3=FKT || more information [here](#fiber-compression-damage))
- Position 5: *reserved*
- Position 6: Friction

For example, `101000` indicates that the model will run with matrix damage and fiber tension damage enabled; `120001` indicates that the model will run with matrix damage, in-plane shear nonlinearity using Schapery theory, and friction.

#### Definition of thickness
Length along the thickness-direction associated with the current integration point.

## State variables
The table below lists all of the state variables in the model. The model requires a minimum of 18 state variables. Additional state variables are defined depending on which (if any) shear nonlinearity and fiber compression features are enabled. For fiber compression model 1: nstatev = 19 and for model 3: nstatev = 26. For shear nonlinearity models 3 or 4: nstatev = 21.

| # | Name             | Description                                                           |
|---|------------------|-----------------------------------------------------------------------|
|  1| `CDM_d2`         | d2, Matrix cohesive damage variable                                   |
|  2| `CDM_Fb1`        | Bulk material deformation gradient tensor component 12                |
|  3| `CDM_Fb2`        | Bulk material deformation gradient tensor component 22                |
|  4| `CDM_Fb3`        | Bulk material deformation gradient tensor component 32                |
|  5| `CDM_B`          | Mode Mixity (*G*<sub>II</sub> / (*G*<sub>I</sub> + *G*<sub>II</sub>)) |
|  6| `CDM_Lc1`        | Characteristic element length along 1-direction                       |
|  7| `CDM_Lc2`        | Characteristic element length along 2-direction                       |
|  8| `CDM_Lc3`        | Characteristic element length along 3-direction                       |
|  9| `CDM_FIm`        | Matrix cohesive failure criterion (del/del_0)                         |
| 10| `CDM_alpha`      | alpha, the cohesive surface normal [degrees, integer]                 |
| 11| `CDM_STATUS`     | STATUS (for element deletion)                                         |
| 12| `CDM_Plas12`     | 1-2 plastic strain                                                    |
| 13| `CDM_Inel12`     | 1-2 inelastic strain                                                  |
| 14| `CDM_FIfT`       | Failure index for fiber tension                                       |
| 15| `CDM_slide1`     | Cohesive sliding displacement, fiber direction                        |
| 16| `CDM_slide2`     | Cohesive sliding displacement, transverse direction                   |
| 17| `CDM_FIfC`       | Failure index for fiber compression                                   |
| 18| `CDM_d1T`        | Fiber tension damage variable                                         |
|---|------------------|-----------------------------------------------------------------------|
| 19| `CDM_d1C`        | Fiber compression damage variable                                     |
|---|------------------|-----------------------------------------------------------------------|
| 20| `CDM_Plas13`     | 1-3 plastic strain                                                    |
| 21| `CDM_Inel13`     | 1-3 inelastic strain                                                  |
|---|------------------|-----------------------------------------------------------------------|
| 22| `CDM_phi0`       | Initial fiber misalignment (radians)                                  |
| 23| `CDM_gamma`      | Current rotation of the fibers due to loading (radians)               |
| 24| `CDM_Fm1`        | Fm1                                                                   |
| 25| `CDM_Fm2`        | Fm2                                                                   |
| 26| `CDM_Fm3`        | Fm3                                                                   |
|---|------------------|-----------------------------------------------------------------------|
| 27| `CDM_Ep1`        | Plastic strain in 11 direction calculated using Schaefer Theory       |
| 28| `CDM_Ep2`        | Plastic strain in 22 direction calculated using Schaefer Theory       |
| 29| `CDM_Ep3`        | Plastic strain in 33 direction calculated using Schaefer Theory       |
| 30| `CDM_Ep4`        | Plastic strain in 12 direction calculated using Schaefer Theory       |
| 31| `CDM_Ep5`        | Plastic strain in 23 direction calculated using Schaefer Theory       |
| 32| `CDM_Ep6`        | Plastic strain in 31 direction calculated using Schaefer Theory       |
| 33| `CDM_fp1`        | Yield criterion (effective stress) calculated using Schaefer Theory   |

### Initial conditions
All state variables should be initialized using the `*Initial conditions` command. As a default, all state variables should be initialized as zero, except `CDM_alpha`, `CDM_STATUS`, and `CDM_phi0`.

The initial condition for `CDM_alpha` can be used to specify a predefined angle for the cohesive surface normal. To specify a predefined `CDM_alpha`, set the initial condition for `CDM_alpha` to an integer (degrees). The range of valid values for `CDM_alpha` depends on the aspect ratio of the element, but values in the range of 0 to 90 degrees are always valid. Setting `CDM_alpha` to -999 will make the subroutine evaluate cracks every 10 degrees in the 2-3 plane to find the correct crack initiation angle. Note that `CDM_alpha` is measured from the 2-axis rotating about the 1-direction. The amount by which alpha is incremented when evaluating matrix crack initiation can be changed from the default of 10 degrees by modifying `alpha_inc` in the `CompDam.parameters` file. Note that `CDM_alpha = 90` only occurs when `CDM_alpha` is initialized as 90; when `CDM_alpha` is initialized to -999, the value of 90 is ignored in the search to find the correct initiation angle since it is assumed that delaminations are handled elsewhere in the finite element model (e.g., using cohesive interface elements).

Since `CDM_STATUS` is used for element deletion, always initialize `CDM_STATUS` to 1.

The initial condition for `CDM_phi0` is used to specify the initial fiber misalignment. One of three options is used depending on the initial condition specified for `CDM_phi0` as follows:
- *&phi;<sub>0</sub>* = 0 :: The value for *&phi;<sub>0</sub>* is calculated for shear instability.
- *&phi;<sub>0</sub>* &le; 0.5 :: The value provided in the initial condition is used as the initial fiber misalignment.
- *&phi;<sub>0</sub>* = 1 :: A pseudo random uniform distribution varying spatially in the 1-direction is used. The spatial distribution algorithm relies on an uniform element size and fiber aligned mesh. The random number generator used is set to generate the same sequence of random numbers during every execution; therefore, the results are repeatable. When using the random distribution for *&phi;<sub>0</sub>*, the characteristic length must be set to include 6 components: `*Characteristic Length, definition=USER, components=6`.

Pre-existing damage can be modeled by creating an element set for the damaged region and specifying different initial conditions for this element set. For example, to create an intraply matrix crack with no out-of-plane orientation, the following initial conditions could be specified for the cracked elements:

<pre>
*Initial Conditions, type=SOLUTION
 damaged_elset,  1.d0,  0.d0,  0.d0,  0.d0,  0.d0,  0.d0,  0.d0,
          0.d0,  0.d0,     0,     1,  0.d0,  0.d0,  0.d0,  0.d0,
          0.d0,  0.d0,  0.d0,  0.d0
</pre>

## Using CompDam with Abaqus/Standard
The repository includes a developmental capability to run the CompDam VUMAT in an Abaqus/Standard analysis using a wrapper, `for/vumatWrapper.for`, that translates between the UMAT and VUMAT user subroutine interfaces. The intended usage is for Abaqus/Standard runs with little or no damage. 

### Usage
To run an analysis with CompDam in Abaqus/Standard, the following input deck template is provided. Note that 9 additional state variables are required.

<pre>
*Section controls, name=control_name, hourglass=ENHANCED
**
*Material, name=IM7-8552
*Density
 1.57e-09,
*User material, constants=40
** 1              2  3          4  5  6  7  8
** feature flags,  , thickness, 4, 5, 6, 7, 8
          100001,  ,       0.1,  ,  ,  ,  ,  ,
**
**  9         10        11        12        13        14        15        16
**  E1,       E2,       G12,      nu12,     nu23,     YT,       SL        GYT,
    171420.0, 9080.0,   5290.0,   0.32,     0.52,     62.3,     92.30,    0.277,
**
**  17        18        19        20        21        22        23        24
**  GSL,      eta_BK,   YC,       alpha0    E3,       G13,      G23,      nu13,
    0.788,    1.634,    199.8,    0.925,      ,          ,         ,          ,
**
**  25        26        27        28        29        30        31        32
**  alpha11,  alpha22,  alpha_PL, n_PL,     XT,       fXT,      GXT,      fGXT,
    -5.5d-6,  2.58d-5,          ,     ,     2326.2,   0.2,      133.3,    0.5,
**
**  33        34        35        36        37        38        39        40
**  XC,       fXC,      GXC,      fGXC,       cf,     w_kb,     None,     mu
    1200.1,      ,         ,          ,         ,     0.1,          ,     0.3
**
*Depvar, delete=11
  28,
  1, CDM_d2
  2, CDM_Fb1
  3, CDM_Fb2
  4, CDM_Fb3
  5, CDM_B
  6, CDM_Lc1
  7, CDM_Lc2
  8, CDM_Lc3
  9, CDM_FIm
 10, CDM_alpha
 11, CDM_STATUS
 12, CDM_Plas12
 13, CDM_Inel12
 14, CDM_FIfT
 15, CDM_slide1
 16, CDM_slide2
 17, CDM_FIfC
 18, CDM_d1T
 19, CDM_d1C
 20, CDM_DIRECT11
 21, CDM_DIRECT21
 22, CDM_DIRECT31
 23, CDM_DIRECT12
 24, CDM_DIRECT22
 25, CDM_DIRECT32
 26, CDM_DIRECT13
 27, CDM_DIRECT23
 28, CDM_DIRECT33
*User defined field
**
** INITIAL CONDITIONS
**
*Initial Conditions, Type=Solution
ALL_ELEMS,  0.d0,  0.d0,  0.d0,  0.d0,  0.d0,  0.d0,  0.d0,
0.d0,  0.d0,  -999,     1,  0.d0,  0.d0,  0.d0,  0.d0,
0.d0,  0.d0,  0.d0,  0.d0,  0.d0,  0.d0,  0.d0,  0.d0,
0.d0,  0.d0,  0.d0,  0.d0,  0.d0
*Initial Conditions, Type=Field, Variable=1
GLOBAL,  0.d0
** GLOBAL is an nset with all nodes attached to compdam enabled elements
** In each step, NLGEOM=YES must be used. This is NOT the default setting.
</pre>

### Current limitations 
As the `vumatWrapper` is a developmental capability, several important limitations exist at present:
1. The material Jacobian tensor is hard-coded in `for/vumatWrapper.for` for IM7/8552 elastic stiffnesses. A more general Jacobian is needed.
2. The material response can become inaccurate for large increments in rotations. If large rotations occur, small increments must be used. A cut-back scheme based on rotation increment size is needed.
3. Testing has been conducted on the C3D8R element type only.


## Example problems
The directory `examples/` includes a collection of example models that use CompDam along with corresponding files that defined the expected results (for use with Abaverify, following the same pattern as the test models in the `tests/` directory. The file `example_runner.py` can be used to automate submission of several models and/or for automatically post-processing the model results to verify that they match the expected results.


## Advanced debugging
Using an interactive debugger helps to identify issues in the Fortran code. Abaqus knowledge base article QA00000007986 describes the details involved. The following is a quick-start guide for direct application to CompDam.

Several statements for debugging need to be uncommented in the environment file. Follow these steps:
1. Copy your system environment file to your local working directory. For the example below, copy the environment file to the `tests` directory.
2. Edit the local environment file: uncomment lines that end with `# <-- Debugging`, `# <-- Debug symbols`, and `# <-- Optimization Debugging`

Run the job with the `-debug` and `-explicit` arguments. For example:
```
abaqus -j test_C3D8R_fiberTension -user ../for/CompDam_DGD.for -double both -debug -explicit
```

This command should open the [Visual Studio debugging software](https://msdn.microsoft.com/en-us/library/sc65sadd.aspx) automatically. Open the source file(s) to debug. At a minimum, open the file with the subroutine entry point `for/CompDam_DGD.for`. Set a break point by clicking in the shaded column on the left edge of the viewer. The break point will halt execution. Press <kbd>F5</kbd> to start the solver. When the break point is reached, a yellow arrow will appear and code execution will pause. Press <kbd>F5</kbd> to continue to the next break point, press <kbd>F11</kbd> to execute the next line of code following execution into function calls (Step Into), or press <kbd>F10</kbd> to execute the next line of code but not follow execution into function calls (Step Over).

To stop execution, close the Visual Studio window. Choose stop debugging and do not save your changes.

[More tips on debugging Fortran programs from Intel](https://software.intel.com/en-us/articles/tips-for-debugging-run-time-failures-in-intel-fortran-applications).

## Python extension module
CompDam can be compiled into a [Python extension module](https://docs.python.org/2/extending/extending.html), which allows many of the Fortran subroutines and functions in the `for` directory to be called from Python. The Python package [`f90wrap`](https://github.com/jameskermode/f90wrap) is used to automatically generate the Python extension modules that interface with the Fortran code. This Python extension module functionality is useful for development and debugging.

### Dependencies and setup
The python extension module requires some additional dependencies. First, the procedure only works on Linux using the bash shell. In addition, `gfortran` 4.6+ is required. Type `gfortran --version` to check if you have this available. The remaining dependencies are python packages and can be installed as follows.

Using [Conda](https://conda.io/docs/user-guide/getting-started.html) significant simplifies the setup process, so it is assumed that you have a recent version of Conda available (see the [Conda installation guide](https://conda.io/docs/user-guide/install/index.html)). Add the Conda-Forge channel by typing:
```
$ conda config --add channels conda-forge
``` 

Conda stores python packages in containers called environments. Create a new environment: 
```
$ conda create --name testing_pyextmod
```
and switch to your new environment:
```
$ source activate testing_pyextmod
```
which will add `(testing_pyextmod)` to the prompt. Install `f90wrap` by typing:
```
(testing_pyextmod) $ conda install f90wrap
```
After typing 'y' in response to the prompt asking if you would like to proceed, Conda will install `f90wrap` and all of its dependencies. This completes the setup process. These steps only need to be executed once.

Note, you can exit the Conda environment by typing `source deactivate`. When you open a new session, you will need to activate the conda environment by typing `source activate testing_pyextmod`.

### Example usage
This section describes how to compile CompDam into a python extension module and run a simple example.

The relevant files are in the `pyextmod` directory, so set `pyextmod` as your current working directory.

The bash shell is required. Type `bash` to open a bash shell session if you are using a different shell. Activate the environment in which you have installed the dependencies listed above, e.g. `source activate testing_pyextmod`.

Next compile the python extension module by executing `make` in the `pyextmod` directory as follows:
```
(testing_pyextmod) CompDam_DGD/pyextmod> make
```

When you execute `make`, the Fortran modules in the `for` directory are compiled to object files, the shared library `_CompDam_DGD.so` is built, and the Python extension module interface `CompDam_DGD.py` is created. A large amount of output is given in the terminal. After the module is created, most of the functionality of CompDam is available in python with `import CompDam_DGD`.

The file `test_pyextmod.py` shows an example of how the Python extension module can be used. Just as when CompDam_DGD is used in Abaqus, `CompDam_DGD.py` expects `CompDam.parameters` and `props` files (if provided) are located in the working directory. In the example `test_pyextmod.py`, the IM7-8552.props file is used.

It is necessary to recompile the CompDam after making changes to the Fortran code. Recompile with the `make` command. It is a good idea to run `make clean` before rerunning `make` to remove old build files.

Note that portions of CompDam that are specific to Abaqus are hidden from `f90wrap` using the preprocessor directive `#ifndef PYEXT`.


## Summary of tests classes
This section includes a brief summary of each test implemented in the `tests` folder. The input deck file names briefly describe the test. All of the input decks start with `test_<elementType>_` and end with a few words describing the test. A more detailed description for each is given below:
- *elastic_fiberTension*: Demonstrates the elastic response in the 1-direction under prescribed extension. The 1-direction stress-strain curve has the modulus E1.
- *elastic_matrixTension*: Demonstrates the elastic response in the 2-direction under prescribed extension. The 2-direction stress-strain curve has the modulus E2.
- *elementSize*: Verifies that the characteristic element lengths Lc1, Lc2, and Lc3 are being properly calculated.
- *failureEnvelope_sig11sig22*: A parametric model in which *&sigma;<sub>11</sub>* is swept from *-X<sub>C</sub>* to *X<sub>T</sub>* and *&sigma;<sub>22</sub>* is swept from *-Y<sub>C</sub>* to *Y<sub>T</sub>* in order to re-create the corresponding failure envelope.
- *failureEnvelope_sig12sig22*: A parametric model in which *&tau;<sub>12</sub>* is swept from *0* to *S<sub>L</sub>* and *&sigma;<sub>22</sub>* is swept from *-Y<sub>C</sub>* to *Y<sub>T</sub>* in order to re-create the corresponding failure envelope.
- *failureEnvelope_sig12sig23*: A parametric model in which *&tau;<sub>12</sub>* is swept from *0* to *S<sub>L</sub>* and *&tau;<sub>23</sub>* is swept from *0* to *S<sub>T</sub>* in order to re-create the corresponding failure envelope.
- *fiberCompression_CDM*: Demonstrates the constitutive response in the 1-direction under prescribed shortening. The 1-direction stress-strain curve is trilinear. A conventional CDM approach to material degradation is used.
- *fiberCompression_DGD*: Demonstrates the constitutive response in the 1-direction under prescribed shortening. The fiber kink band model is used.
- *fiberLoadReversal*: Demonstrates the constitutive response in the 1-direction under prescribed extension and shortening reversals. The 1-direction stress-strain curve shows the intended behavior under load reversal.
- *fiberTension*: Demonstrates the constitutive response in the 1-direction under prescribed extension. The 1-direction stress-strain curve is trilinear.
- *matrixCompression*: Demonstrates the constitutive response in the 2-direction under prescribed compression displacement.
- *matrixTension*: Demonstrates the constitutive response in the 2-direction under prescribed extension. The 2-direction stress-strain curve is bilinear.
- *mixedModeMatrix*: A parametric model used to stress the internal DGD convergence loop. The crack angle *&alpha;* and the direction of loading are varied. Large tensile and compressive displacements are prescribed to ensure the DGD method is able to find converged solutions under a wide variety of deformations.
- *nonlinearShear12*: Demonstrates the nonlinear Ramberg-Osgood constitutive response under prescribed simple shear deformation in the 1-2 plane with matrix damage enabled. Several cycles of loading and unloading are applied, with increasing peak displacements in each cycle.
- *nonlinearShear12_loadReversal*: Demonstrates the response of the Ramberg-Osgood model under load reversal in the 1-2 plane. Several cycles of loading and unloading are applied, with inelastic strain accumulated throughout the load history. Damage is disabled.
- *nonlinearShear13*: Demonstrates the nonlinear Ramberg-Osgood constitutive response under prescribed simple shear deformation in the 1-3 plane with matrix damage enabled. Several cycles of loading and unloading are applied, with increasing peak displacements in each cycle.
- *nonlinearShear13_loadReversal*: Demonstrates the response of the Ramberg-Osgood model under load reversal in the 1-3 plane. Several cycles of loading and unloading are applied, with inelastic strain accumulated throughout the load history. Damage is disabled.
- *schapery12*: Demonstrates the in-plane response to prescribed simple shear deformation for the Schapery micro-damage model. Several cycles of loading and unloading are applied, with increasing peak shear displacements in each cycle.
- *simpleShear12*: Demonstrates the constitutive response under prescribed simple shear. The shear stress-strain curve is bilinear.
- *simpleShear12friction*: Demonstrates the constitutive response under prescribed simple shear with friction enabled. The element is loaded under transverse compression and then sheared. Shows the friction-induced stresses.

## Contributing
We invite your contributions to CompDam_DGD! Please submit contributions (including a test case) with pull requests so that we can reproduce the behavior of interest. Commits history should be clean. Please contact the developers if you would like to make a major contribution to this repository. Here is a [checklist](contributing-checklist.md) that we use for contributions.

## Citing CompDam
If you use CompDam, please cite using the following BibTex entry:

<pre>
  @misc{CompDam,
  title={CompDam - Deformation Gradient Decomposition (DGD), v2.2.0},
  author={Leone Jr., F. A., Bergan, A. C., D\'avila, C. G. },
  note={https://github.com/nasa/CompDam_DGD},
  year={2018}
  }
</pre>
