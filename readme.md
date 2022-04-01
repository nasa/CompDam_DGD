# CompDam - Deformation Gradient Decomposition (DGD)
This code is a continuum damage mechanics (CDM) material model intended for use with the Abaqus finite element code. This is a research code which aims to provide an accurate representation of mesoscale damage modes in fiber-reinforced polymer composite materials in finite element models in which each ply is discretely represented.

The CDM material model is implemented as an Abaqus/Explicit user subroutine (VUMAT) for the simulation of matrix cracks formed under tensile, compressive, and shear loading conditions and fiber fracture under tensile and compressive loading. Within CompDam, the emphasis of many recent developments has been on accurately representing the kinematics of composite damage. The kinematics of matrix cracks are represented by treating them as cohesive surfaces embedded in a deformable bulk material in accordance with the Deformation Gradient Decomposition (DGD) approach. Fiber tensile damage is modeled using conventional CDM strain-softening.

This software may be used, reproduced, and provided to others only as permitted under the terms of the agreement under which it was acquired from the U.S. Government. Neither title to, nor ownership of, the software is hereby transferred. This notice shall remain on all copies of the software.

Copyright 2016 United States Government as represented by the Administrator of the National Aeronautics and Space Administration. No copyright is claimed in the United States under Title 17, U.S. Code. All Other Rights Reserved.

Publications that describe the theories used in this code:
- Andrew C. Bergan, ["A Three-Dimensional Mesoscale Model for In-Plane and Out-of-Plane Fiber Kinking"](https://doi.org/10.2514/6.2019-1548) *AIAA SciTech Forum*, San Diego, California, 7-11 January 2019.
- Carlos Dávila, ["From S-N to the Paris Law with a New Mixed-Mode Cohesive Fatigue Model"](https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20180004395.pdf) NASA/TP-2018-219838, June 2018.
- Andrew C. Bergan, et al., ["Development of a Mesoscale Finite Element Constitutive Model for Fiber Kinking"](https://doi.org/10.2514/6.2018-1221) *AIAA SciTech Forum*, Kissimmee, Florida, 8-12 January 2018.
- Frank A. Leone Jr. ["Deformation gradient tensor decomposition for representing matrix cracks in fiber-reinforced materials"](https://dx.doi.org/10.1016/j.compositesa.2015.06.014) *Composites Part A* (2015) **76**:334-341.
- Frank A. Leone Jr. ["Representing matrix cracks through decomposition of the deformation gradient tensor in continuum damage mechanics methods"](https://iccm20.org/fullpapers/file?f=Abk7n4gkWV) *Proceedings of the 20th International Conference on Composite Materials*, Copenhagen, Denmark, 19-24 July 2015.
- Cheryl A. Rose, et al. ["Analysis Methods for Progressive Damage of Composite Structures"](https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20140001002.pdf) NASA/TM-2013-218024, July 2013.

Examples of this code being applied can be found in the following publications:
- Andrew C. Bergan and Wade C. Jackson, ["Validation of a Mesoscale Fiber Kinking Model through Test and Analysis of Double Edge Notch Compression Specimens"](https://doi.org/10.12783/asc33/26003) *33rd American Society for Composites (ASC) Annual Technical Conference*, Seattle, Washington, 24-27 September 2018.
- Imran Hyder, et al., ["Implementation of a Matrix Crack Spacing Parameter in a Continuum Damage Mechanics Finite Element Model"](https://doi.org/10.12783/asc33/26052) *33rd American Society for Composites (ASC) Annual Technical Conference*, Seattle, Washington, 24-27 September 2018.
- Frank Leone, et al., ["Benchmarking Mixed Mode Matrix Failure in Progressive Damage and Failure Analysis Methods"](https://doi.org/10.12783/asc33/26030) *33rd American Society for Composites (ASC) Annual Technical Conference*, Seattle, Washington, 24-27 September 2018.
- Brian Justusson, et al., et al., ["Quantification of Error Associated with Using Misaligned Meshes in Continuum Damage Mechanics Material Models for Matrix Crack Growth Predictions in Composites"](https://doi.org/10.12783/asc33/26097) *33rd American Society for Composites (ASC) Annual Technical Conference*, Seattle, Washington, 24-27 September 2018.
- Kyongchan Song, et al. ["Continuum Damage Mechanics Models for the Analysis of Progressive Damage in Cross-Ply and Quasi-Isotropic Panels Subjected to Static Indentation"](https://doi.org/10.2514/6.2018-1466) *AIAA SciTech Forum*, Kissimmee, Florida, 8-12 January 2018.
- Imran Hyder, et al. ["Assessment of Intralaminar Progressive Damage and Failure Analysis Using an Efficient Evaluation Framework"](https://doi.org/10.12783/asc2017/15405) *32nd American Society for Composites (ASC) Annual Technical Conference*, West Lafayette, Indiana, 22-25 October 2017.
- Frank A. Leone, et al. ["Fracture-Based Mesh Size Requirements for Matrix Cracks in Continuum Damage Mechanics Models"](https://doi.org/10.2514/6.2017-0198) *AIAA SciTech Forum*, Grapevine, Texas, 9-13 January 2017.
- Mark McElroy, et al. ["Simulation of delamination-migration and core crushing in a CFRP sandwich structure"](https://doi.org/10.1016/j.compositesa.2015.08.026) *Composites Part A* (2015) **79**:192-202.

For any questions, please contact the developers:
- Frank Leone   | [frank.a.leone@nasa.gov](mailto:frank.a.leone@nasa.gov)     | (W) 757-864-3050
- Andrew Bergan | [andrew.c.bergan@nasa.gov](mailto:andrew.c.bergan@nasa.gov) | (W) 757-864-3744
- Carlos Dávila | [carlos.g.davila@nasa.gov](mailto:carlos.g.davila@nasa.gov) | (W) 757-864-9130


## Table of contents
- [Getting started](#getting-started)
- [Element compatibility](#element-compatibility)
- [Model features](#model-features)
- [Material properties](#material-properties)
- [State variables](#state-variables)
- [Model parameters](#model-parameters)
- [Implicit solver compatibility](#implicit-solver-compatibility)
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
[Intel Fortran Compiler](https://software.intel.com/en-us/fortran-compilers) version 11.1 or newer is required to compile the code ([more information about compiler versions](usersubroutine-prerequisites.md)). MPI must be installed and configured properly so that the MPI libraries can be linked by CompDam. It is recommended that Abaqus 2016 or newer is used with this code. Current developments and testing are conducted with Abaqus 2019. Python supporting files require Python 2.7.

### Initial setup
After cloning the CompDam_DGD git repository, it is necessary to run the setup script file `setup.py` located in the repository root directory:
```
$ python setup.py
```

The main purpose of the setup.py script is to 1) set the `for/version.for` file and 2) add git-hooks that automatically update the `for/version.for`.

In the event that you do not have access to python, rename `for/version.for.nogit` to `for/version.for` manually. The additional configuration done by `setup.py` is not strictly required.

### Abaqus environment file settings
The `abaqus_v6.env` file must have [`/fpp`](https://software.intel.com/en-us/fortran-compiler-developer-guide-and-reference-fpp), [`/Qmkl:sequential`](https://software.intel.com/en-us/fortran-compiler-developer-guide-and-reference-mkl-qmkl), and [`/free`](https://software.intel.com/en-us/fortran-compiler-developer-guide-and-reference-free) in the `ifort` command where the format for Windows is used. The corresponding Linux format is: `-fpp`, `-free`, and `-mkl=sequential`. The `/fpp` option enables the Fortran preprocessor, which is required for the code to compile correctly. The `/free` option sets the compiler to free-formatting for the source code files. The `/Qmkl:sequential` enables the [Intel Math Kernel Library (MKL)](https://software.intel.com/en-us/mkl), which provides optimized and verified functions for many mathematical operations. The MKL is used in this code for calculating eigenvalues and eigenvectors.

A sample environment file is provided in the `tests` directory for Windows and Linux systems.

### Submitting a job
This code is an Abaqus/Explicit VUMAT. Please refer to the Abaqus documentation for the general instructions on how to submit finite element analyses using user subroutines. Please see the [example input file statements](#example-input-file-statements) for details on how to interface with this particular VUMAT subroutine.

Analyses with this code **must** be run in double precision. Some of the code has double precision statements and variables hard-coded, so if Abaqus/Explicit is run in single precision, compile-time errors will arise. When submitting an Abaqus/Explicit job from the command line, double precision is specified by including the command line argument `double=both`.

Geometric nonlinearity must be used in each analysis step. Geometric nonlinearity being turned on is the default in Abaqus/Explicit analysis steps. Geometric nonlinearity can be explicitly stated in the input deck `*Step` command with the keyword option `nlgeom=YES`.

For example, run the test model `test_C3D8R_elastic_fiberTension` in the `tests` directory with the following command:
```
$ abaqus job=test_C3D8R_elastic_fiberTension user=../for/CompDam_DGD.for double=both
```

### Example input file statements
Example 1, using an [external material properties file](#defining-the-material-properties-in-a-props-file):

    *Section controls, name=control_name, distortion control=YES
    **
    *Material, name=IM7-8552
    *Density
     1.57e-09,
    *Depvar
    ** *Depvar, delete=11
    ** The delete keyword is optional. It uses a state variable as a flag to delete elements.
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
    *User material, constants=1
    ** feature flags,
              100001,
    **
    *Initial Conditions, type=SOLUTION
     elset_name,  0.d0,  0.d0,  0.d0,  0.d0,  0.d0,  0.d0,  0.d0,
     0.d0,  0.d0,  0.d0,     1,  0.d0,  0.d0,  0.d0, 0.d0,
     0.d0,  0.d0,  0.d0,  0.d0

Example 2, using an [input deck command](#defining-the-material-properties-in-the-input-deck):

    *Section controls, name=control_name, distortion control=YES
    **
    *Material, name=IM7-8552
    *Density
     1.57e-09,
    *Depvar
    ** *Depvar, delete=11
    ** The delete keyword is optional. It uses a state variable as a flag to delete elements.
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
    ** 1              2  3  4  5  6  7  8
    ** feature flags,  ,  ,  ,  ,  ,  ,  ,
              100001,  ,  ,  ,  ,  ,  ,  ,
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
    **  XC,       fXC,      GXC,      fGXC,       cl,     w_kb,     None,     mu
        1200.1,      ,         ,          ,         ,     0.1,          ,     0.3
    ** For spacing below a6=schaefer_a6, b2=schaefer_b2, n=schaefer_n and A=schaefer_A
    **  41        42        43        44
    **  a6,       b2,       n,        A
    **
    *Initial Conditions, type=SOLUTION
     elset_name,  0.d0,  0.d0,  0.d0,  0.d0,  0.d0,  0.d0,  0.d0,
     0.d0,  0.d0,  0.d0,     1,  0.d0,  0.d0,  0.d0, 0.d0,
     0.d0,  0.d0,  0.d0,  0.d0

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


## Element compatibility
CompDam_DGD has been developed and tested using the Abaqus three-dimensional (3-D), reduced-integration `C3D8R` hexahedral solid elements. Limited testing has been performed using the `CPS4R` plane stress element and the fully-integrated `C3D8` solid element. Because CompDam_DGD is a material model, it is expected to be compatible with continuum elements generally. However, users are advised to perform tests with any previously untested element types before proceeding to use CompDam_DGD in larger structural models.

CompDam is also compatible with cohesive elements for modeling delaminations, and can be used with the two-dimensional (2-D) `COH2D4` and 3-D `COH3D8` cohesive elements. Development and testing has been performed primarily with the 3-D cohesive elments.


## Model features
The CompDam_DGD material model implements a variety of features that can be enabled or disabled by the user. An overview of these features is provided in this section. The material properties required for each feature are listed. References are provided to more detailed discussions of the theoretical framework for each feature.

### Elastic response
The composite materials modeled with CompDam_DGD can be defined assuming either [transverse isotropy](https://www.efunda.com/formulae/solid_mechanics/mat_mechanics/hooke_iso_transverse.cfm) or [orthotropy](https://www.efunda.com/formulae/solid_mechanics/mat_mechanics/hooke_orthotropic.cfm). For a transversely isotropic material definition, the following properties must be defined: E1, E2, G12, v12, and v23. For an orthotropic material definition, the following additional properties must be defined: E2, G13, G23, and nu13.

### Matrix damage
Tensile and compressive matrix damage is modeled by embedding cohesive laws to represent cracks in the material according to the deformation gradient decomposition method of [Leone (2015)](https://doi.org/10.1016/j.compositesa.2015.06.014). The matrix crack normals can have any orientation in the 2-3 plane, defined by the angle `CDM_alpha`. The mixed-mode behavior of matrix damage initiation and evolution is defined according to the Benzeggagh-Kenane law. The initiation of compressive matrix cracks accounts for friction on the potential crack surface according to the LaRC04 failure criteria. In the notation of the [paper](https://doi.org/10.1016/j.compositesa.2015.06.014), `Q` defines the material direction that is severed by the crack. In this implementation `Q=2` except when `CDM_alpha = 90`.

The following material properties are required for the prediction of matrix damage: YT, SL, GYT, GSL, eta_BK, YC, and alpha0. The state variables related to matrix damage are `CDM_d2`, `CDM_FIm`, `CDM_B`, `CDM_alpha`, `CDM_Fb1`, `CDM_Fb2`, and `CDM_Fb3`.

### Thermal strains
The thermal strains are calculated by multiplying the 1-, 2-, and 3-direction coefficients of thermal expansion by the current &Delta;T, as provided by the Abaqus solver. The thermal strains are subtracted from the current total strain.

The required material properties are the coefficients of thermal expansion in the 1 and 2 directions. It is assumed that the 2- and 3-direction coefficients of thermal expansion are equal.

Hygroscopic strains are not accounted for. If the effects of hygroscopic expansion are to be modeled, it is recommended to smear the hygroscopic and thermal expansion coefficients to approximate the response using the solver-provided &Delta;T.

### Shear nonlinearity
Three approaches to modeling the matrix nonlinearity are available: Ramberg-Osgood plasticity, Schapery theory, and Schaefer plasticity. These three methods are mutually exclusive and optional.

#### Ramberg-Osgood plasticity
Shear nonlinearity in the 1-2 and/or the 1-3 plane can be modeled using the [Ramberg-Osgood equation](https://en.wikipedia.org/wiki/Ramberg%E2%80%93Osgood_relationship), with its parameters selected to fit experimental data. As applied herein, the Ramberg-Ogsood equation is written in the following form for the 1-2 plane:

*&gamma;*<sub>12</sub> = [*&tau;*<sub>12</sub> + *&alpha;*<sub>PL</sub>sign(*&tau;*<sub>12</sub>)|*&tau;*<sub>12</sub>|<sup>*n*<sub>PL</sub></sup>]/*G*<sub>12</sub>

where *&gamma;*<sub>12</sub> is the shear strain and *&tau;*<sub>12</sub> is the shear stress. Likewise, the expression for the 1-3 plane is

*&gamma;*<sub>13</sub> = [*&tau;*<sub>13</sub> + *&alpha;*<sub>PL</sub>sign(*&tau;*<sub>13</sub>)|*&tau;*<sub>13</sub>|<sup>*n*<sub>PL</sub></sup>]/*G*<sub>13</sub>

Prior to the initiation of matrix damage (i.e., `CDM_d2 = 0`), the nonlinear shear response due to the above equation is plastic, and the unloading/reloading slope is unchanged. No pre-peak nonlinearity is applied to the matrix tensile or compressive responses (i.e., *&sigma;<sub>22</sub>*).

The required material inputs are the two parameters in the above equation: *&alpha;*<sub>PL</sub> and *n*<sub>PL</sub>. Note that the same constants are used for the 1-2 and 1-3 planes under the assumption of transverse isotropy (see [Seon et al. (2017)](https://doi.org/10.12783/asc2017/15267)). For the 1-2 plane, the state variables `CDM_Plas12` and `CDM_Inel12` are used to track the current plastic shear strain and the total amount of inelastic plastic shear strain that has occurred through the local deformation history, respectively. For cases of monotonic loading, `CDM_Plas12` and `CDM_Inel12` should have the same magnitude. Likewise, the state variables `CDM_Plas13` and `CDM_Inel13` are utilized for the 1-3 plane. The [feature flags](#controlling-which-features-are-enabled) can be used to enable this Ramberg-Osgood model in the 1-2 plane, 1-3 plane, or both planes.

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

where *a*<sub>6</sub>, *b*<sub>2</sub>, *A* and *n* are material constants needed for the Schaefer prepeak material model. These four material properties must be defined in an [external material properties file](#defining-the-material-properties-in-a-props-file). Additionally, when defining  *a*<sub>6</sub>, *b*<sub>2</sub>, *A* and *n* in the external material properties file the variables are prefixed with schaefer_ (to disambiguate the otherwise nondescript material property names and symbols).

The above two equations are used in concert to determine plastic strain through the relationship:

*&epsilon;* <sub>plastic</sub> = *n A f* <sup>*n* - 1</sup> (&part;*f*/&part;*S*<sub>*i*</sub>) (&part;*f*/&part;*S*<sub>j</sub>)

*f* (i.e., schaefer_f) and the tensorial plastic strain determined by the nonlinearity model are stored as state variables (27 through 32 for plastic strain and 33 for *f*)

### Fiber tensile damage
A continuum damage mechanics model similar to the work of [Maimí et al. (2007)](https://doi.org/10.1016/j.mechmat.2007.03.005) is used to model tensile fiber damage evolution. The model utilizes a non-interacting maximum strain failure criterion, and bilinear softening after the initiation of failure. The area under the stress-strain curve is equal to the fracture toughness divided by the element length normal to fracture, i.e., `CDM_Lc1`. The required material properties are: XT, fXT, GXT, and fGXT, where fXT and fGXT are ratios of strength and fracture toughness for bilinear softening, defined as the *n* and *m* terms in equations (25) and (26) of [Dávila et al. (2009)](https://doi.org/10.1007/s10704-009-9366-z). To model a linear softening response, both fXT and fGXT should be set equal to 0.5.

### Fiber compression damage

#### Model 1: Max strain, bilinear softening (BL)
Same model as in tension, but for compression. Assumes maximum strain failure criterion and bilinear softening. The required material properties are: XC, fXC, GXC, and fGXC.

Load reversal assumptions from [Maimí et al. (2007)](https://doi.org/10.1016/j.mechmat.2007.03.006).

#### Model 2: Placeholder

#### Model 3, 4, 5: Fiber kinking theory (FKT)
A model based on Budiansky's fiber kinking theory from [Budiansky (1983)](https://doi.org/10.1016/0045-7949(83)90141-4), [Budiansky and Fleck (1993)](https://doi.org/10.1016/0022-5096(93)90068-Q), and [Budiansky et al. (1998)](https://doi.org/10.1016/S0022-5096(97)00042-2) implemented using the DGD framework to model in-plane (1-2) and/or out-of-plane (1-3) fiber kinking. The model is described in detail in [Bergan et al. (2018)](https://doi.org/10.2514/6.2018-1221), [Bergan and Jackson (2018)](https://doi.org/10.12783/asc33/26003), and [Bergan (2019)](https://doi.org/10.2514/6.2019-1548). The model accounts for fiber kinking due to shear instability by considering an initial fiber misalignment, nonlinear shear stress-strain behavior via Ramberg-Osgood, and geometric nonlinearity. Fiber failure can be introduced by specifying a critical fiber rotation angle.

The required material properties are: XC, YC, wkb, alpha0, alpha_PL, and n_PL. The [feature flag](#controlling-which-features-are-enabled) for fiber compression should be set to '3', '4', or '5' to activate this model feature. Model '3' enables in-plane (1-2 plane) fiber kinking. Model '4' enables out-of-plane (1-3 plane) fiber kinking. Model '5' enables fiber kinking in both planes (uncoupled). This feature requires 25 state variables to be defined and initialized. The relevant state variables are:
- `CDM_phi0_12`: initial fiber misalignment (radians) in the 1-2 plane.
- `CDM_phi0_13`: initial fiber misalignment (radians) in the 1-3 plane.
- `CDM_gamma_12`: rotation of the fibers due to loading (radians) in the 1-2 plane.
- `CDM_gamma_13`: rotation of the fibers due to loading (radians) in the 1-3 plane.
- `CDM_Fbi`: the components of the first column of `Fm` used for decomposing the element where `i=1,2,3`.

The current fiber misalignment is `CDM_phi0_1i + CDM_gamma_1i` where `i=2 or 3`.

The initial conditions for the state variable `CDM_phi0_12` and `CDM_phi0_13` determine the initial fiber misalignments as described in [initial conditions](#initial-conditions).

A fiber failure criterion described in [Bergan and Jackson (2018)](https://doi.org/10.12783/asc33/26003) is implemented to represent the material behavior in confined conditions under large strains (post failure). The fiber failure criterion is satisfied when

*&phi;* &ge; *&phi;*<sub>ff,c</sub>

where *&phi;* is the current fiber rotation. Once the fiber failure criterion is satisfied, the plastic shear strain is held constant. The value for *&phi;*<sub>ff,c</sub> is defined as the parameter `fkt_fiber_failure_angle` since it is not a well-defined material property. The fiber failure criterion is disabled when *&phi;*<sub>ff,c</sub> < 0. The same angle is used for in-plane and out-of-plane kinking.

The fiber kinking theory model implemented here is preliminary and has some known shortcomings and caveats:
- The model has only been tested for C3D8R. Limited application with C3D6 demonstrated issues. No testing has been completed for other element types.
- The interaction of this model with matrix cracking has not been fully tested and verified.
- No effort has been made to model progressive crushing.
- Other mechanisms of fiber compressive failure (e.g., shear driven fiber breaks) are not accounted for. An outcome of this is that the model predicts the material does not fail if shear deformation is fully constrained.
- No special consideration for load reversal has been included.

Relevant single element tests are named starting with `test_C3D8R_fiberCompression_FKT`.

### Friction
Friction is modeled on the damaged fraction of the cross-sectional area of matrix cracks and delaminations using the approach of [Alfano and Sacco (2006)](https://doi.org/10.1002/nme.1728). The coefficient of friction *&mu;* must be defined to account for friction on the failed crack surface.

The amount of sliding which has taken place in the longitudinal and transverse directions are stored in state variables `CDM_slide1` and `CDM_slide2`, respectively.

### Fatigue analyses
The cohesive fatigue constitutive model in CompDam can predict the initiation and the propagation of matrix cracks and delaminations as a function of fatigue cycles. The analyses are conducted such that the applied load (or displacement) corresponds to the maximum load of a fatigue cycle. The intended use is that the maximum load (or displacement) is held constant while fatigue damage develops with increasing step time. The constitutive model uses a specified load ratio *R*<sub>min</sub>/*R*<sub>max</sub>, the solution increment, and an automatically-calculated cycles-per-increment ratio to accumulate the damage due to fatigue loading. The cohesive fatigue model response is based on engineering approximations of the endurance limit as well as the Goodman diagram. This approach can predict the stress-life diagrams for crack initiation, the Paris law regime, as well as the transient effects of crack initiation and stable tearing.

A detailed description of the initial development of the novel cohesive fatigue law is available in a [2018 NASA technical paper by Carlos Dávila](https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20180004395.pdf). An updated fatigue damage accumulation function is derived in the 2020 NASA technical paper ["Evaluation of Fatigue Damage Accumulation Functions for Delamination Initiation and Propagation" by Dávila et al.](https://www.researchgate.net/profile/Carlos_Davila8/publication/340715727_Evaluation_of_Fatigue_Damage_Accumulation_Functions_for_Delamination_Initiation_and_Propagation/links/5e99b744a6fdcca789204fb9/Evaluation-of-Fatigue-Damage-Accumulation-Functions-for-Delamination-Initiation-and-Propagation.pdf), which is the fatigue damage accumulation function that is applied herein.

#### Usage
The fatigue capability of CompDam is disabled by default. To run a fatigue analysis, one of the analysis steps must be identified as a fatigue step. A step is identified as a fatigue step by setting the `fatigue_step` parameter to the target step number, e.g., `fatigue_step = 2` for the second analysis step to be a fatigue step. The first analysis step cannot be a fatigue step, as the model is assumed to be initially unloaded.

The load ratio *R*<sub>min</sub>/*R*<sub>max</sub> has a default value of 0.1, and can be changed using the parameter `fatigue_R_ratio`.

An example of a double cantilever beam subjected to fatigue under displacement-control is included in the `examples/` directory. The geometry and conditions of this example problem correspond to the results presented in Figure 20 of [Dávila (2018)](https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20180004395.pdf).

#### Interpreting the results of a fatigue analysis
Within a fatigue step, each solution increment represents either a number of fatigue cycles or a fractional part of a single fatigue cycle. During the solution, the number of fatigue cycles per solution increment changes based on the maximum amount of energy dissipation in any single element. If the rate of energy dissipation is too high (as defined by the parameter `fatigue_damage_max_threshold`), the increments-to-cycles ratio is decreased. If the rate of energy dissipation is too low (as defined by the parameter `fatigue_damage_min_threshold`), the increments-to-cycles ratio is increased. The parameter `cycles_per_increment_init` defines the initial ratio of fatigue cycles per solution increment. Any changes to increments-to-cycles ratio are logged in an additional output file ending in `_inc2cycles.log`, with columns for the fatigue step solution increment, the updated increments-to-cycles ratio, and the accumulated fatigue cycles.

### Finite stress and strain definitions
The strain is calculated using the deformation gradient tensor provided by the Abaqus solver. The default strain definition used is the Green-Lagrange strain:

***E*** = (***F***<sup>T</sup>***F*** - ***I***)/2

Hooke's law is applied using the Green-Lagrange strain to calculate the 2<sup>nd</sup> Piola-Kirchhoff stress ***S***.

### Fiber nonlinearity
Nonlinear elastic behavior in the fiber direction can be introduced with the material property c<sub>*l*</sub>. The expression used follows [Kowalski (1988)](https://doi.org/10.1520/STP26136S):

*E<sub>1</sub>* = *E<sub>1</sub>*(1 + c<sub>*l*</sub>*&epsilon;*<sub>11</sub>)

By default, fiber nonlinearity is disabled by setting c<sub>*l*</sub> = 0.


## Material properties
A set of material properties must be defined for the material of interest. This section describes how to specify the material properties.

### Defining material properties
Material properties can be defined in the input deck or in a separate `.props` file. Definition of the material properties in a `.props` file is more convenient and generally preferred since it isolates the material properties from the structural model definition.

#### Defining the material properties in a `.props` file
Using a `.props` file is a versatile means of defining material properties. The subroutine looks for a file named as `jobName_materialName` where the job name is the name of the Abaqus job (default is input deck name) and the material is name assigned as `*Material, name=materialName` in the input deck. If no file is found, then the subroutine looks for `materialName.props`. The `.props` file must be located in the Abaqus working directory.

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
| 18 | *&eta;*                | eta_BK   | BK exponent for mode-mixity                   | -                                              | 0 < *&eta;* < &infin;                 |                 |
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
| 37 | *c<sub>l</sub>*        | cl       | Fiber nonlinearity coefficient                | -                                              | 0 &le; *c<sub>l</sub>* 33                |                 |
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
  5. [*&gamma;*, number of fatigue cycles to reach endurance](#fatigue-properties)
  6. [*&epsilon;*, endurance limit](#fatigue-properties)
  7. [*&eta;*, brittleness](#fatigue-properties)
  8. [*p*, Paris Law curve-fitting parameter](#fatigue-properties)
- &infin; is calculated with the Fortran intrinsic `Huge` for double precision
- In the event that both a `.props` file is found and material properties are specified in the input deck (`nprops > 8`), then the material properties from the input deck are used and a warning is used.

### Inputs that can only be defined on the `*Material` data lines in the input deck

#### Controlling which features are enabled
The feature flags are defined in the input deck on the material property data lines. These properties must be defined in the input deck whether the other material properties are defined via the .props file or via the input deck. While the feature flags are not material properties per se, they are used in controlling the behavior of the material model.

Model features can be enabled or disabled by two methods. The first method is specifying only the material properties required for the features you would like to enable. CompDam_DGD disables any feature for which all of the required material properties have not been assigned. If an incomplete set of material properties are defined for a feature, a warning is issued.

The second method is by specifying the status of each feature directly as a material property in the input deck. Each feature of the subroutine is controlled by a position in an integer, where 0 is disabled and 1 is enabled. In cases where mutually exclusive options are available, numbers greater than 1 are used to specify the particular option to use.

The positions correspond to the features as follows:
- Position 1: Matrix damage (1=intra-laminar cracking in solid elements, 2=interlaminar cracking in cohesive elements)
- Position 2: Shear nonlinearity (1=Ramberg-Osgood 1-2 plane, 2=Schapery, 3=Ramberg-Osgood 3-D, 4=Ramberg-Osgood 1-3 plane, 5=Schaefer || more information [here](#shear-nonlinearity))
- Position 3: Fiber tensile damage
- Position 4: Fiber compression damage (1=max strain, 2=N/A, 3=FKT-12, 4=FKT-13, 5=FKT-3D || more information [here](#fiber-compression-damage))
- Position 5: Energy output contribution (0=all mechanisms, 1=only fracture energy, 2=only plastic energy)
- Position 6: Friction

For example, `101000` indicates that the model will run with matrix damage and fiber tension damage enabled; `120001` indicates that the model will run with matrix damage, in-plane shear nonlinearity using Schapery theory, and friction; and `200000` indicates that the material model is being applied to cohesive elements.

#### Definition of thickness
Length along the thickness-direction associated with the current integration point. This input is used only for 2-D plane stress elements and does not affect the performance of 3-D hexahedral elements or cohesive elements. The three characteristic element lengths of hexahedral elements are calculated using the VUCHARLENGTH subroutine based on the element nodal coordinates.

#### Fatigue properties
Material inputs 5 through 8 are related to the cohesive fatigue model. Each of these inputs are optional as they have the default values listed in the below table.

| # | Symbol          | Description                                                           | Default value |
|---|-----------------|-----------------------------------------------------------------------|---------------|
| 5 | *&gamma;*       | number of fatigue cycles to reach the endurance limit                 | 10,000,000    |
| 6 | *&epsilon;*     | endurance limit                                                       | 0.2           |
| 7 | *&eta;*         | brittleness                                                           | 0.95          |
| 8 | *p*             | Paris Law curve-fitting parameter                                     | 0.0           |

## State variables
The table below lists all of the state variables in the model. The model requires a minimum of 18 state variables. Additional state variables are defined depending on which (if any) shear nonlinearity and fiber compression features are enabled. For fiber compression model 1: nstatev = 19 and for model 3: nstatev = 25. For shear nonlinearity models 3 or 4: nstatev = 21.

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
| 22| `CDM_phi0_12`    | Initial fiber misalignment, 1-2 plane (radians)                       |
| 23| `CDM_gamma_12`   | Current rotation of the fibers due to loading, 1-2 plane (radians)    |
| 24| `CDM_phi0_13`    | Initial fiber misalignment, 1-3 plane (radians)                       |
| 25| `CDM_gamma_13`   | Current rotation of the fibers due to loading, 1-3 plane (radians)    |
| 26| `CDM_reserve`    | Reserved                                                              |
|---|------------------|-----------------------------------------------------------------------|
| 27| `CDM_Ep1`        | Plastic strain in 11 direction calculated using Schaefer Theory       |
| 28| `CDM_Ep2`        | Plastic strain in 22 direction calculated using Schaefer Theory       |
| 29| `CDM_Ep3`        | Plastic strain in 33 direction calculated using Schaefer Theory       |
| 30| `CDM_Ep4`        | Plastic strain in 12 direction calculated using Schaefer Theory       |
| 31| `CDM_Ep5`        | Plastic strain in 23 direction calculated using Schaefer Theory       |
| 32| `CDM_Ep6`        | Plastic strain in 31 direction calculated using Schaefer Theory       |
| 33| `CDM_fp1`        | Yield criterion (effective stress) calculated using Schaefer Theory   |

When using the material model with cohesive elements, a different set of state variables are used. Cohesive state variables with a similar continuum damage mechanics counterpart utilize the same state variable number. The cohesive element material model requires fewer state variables, however, resulting in "gaps" in the cohesive state variable numbering.

| # | Name             | Description                                                           |
|---|------------------|-----------------------------------------------------------------------|
|  1| `COH_dmg`        | Cohesive damage variable                                              |
|  2| `COH_delta_s1`   | Displacement-jump in the first shear direction                        |
|  3| `COH_delta_n`    | Displacement-jump in the normal direction                             |
|  4| `COH_delta_s2`   | Displacement-jump in the second shear direction                       |
|  5| `COH_B`          | Mode Mixity (*G*<sub>II</sub> / (*G*<sub>I</sub> + *G*<sub>II</sub>)) |
|  9| `COH_FI`         | Cohesive failure criterion (del/del_0)                                |
| 15| `COH_slide1`     | Cohesive sliding displacement, fiber direction                        |
| 16| `COH_slide2`     | Cohesive sliding displacement, transverse direction                   |

### Initial conditions
All state variables should be initialized using the `*Initial conditions` command. As a default, all state variables should be initialized as zero, except `CDM_alpha`, `CDM_STATUS`, `CDM_phi0_12`, and `CDM_phi0_13`.

The initial condition for `CDM_alpha` can be used to specify a predefined angle for the cohesive surface normal. To specify a predefined `CDM_alpha`, set the initial condition for `CDM_alpha` to an integer (degrees). The range of valid values for `CDM_alpha` depends on the aspect ratio of the element, but values in the range of 0 to 90 degrees are always valid. Setting the parameter `alpha_search` to TRUE will make the subroutine evaluate cracks every 10 degrees (by default) in the 2-3 plane to find the correct crack initiation angle. Note that `CDM_alpha` is measured from the 2-axis rotating about the 1-direction. The amount by which alpha is incremented when evaluating matrix crack initiation can be changed by modifying `alpha_inc` in the `CompDam.parameters` file. Note that `CDM_alpha = 90` only occurs when `CDM_alpha` is initialized as 90; the value of 90 is ignored in the search to find the correct initiation angle since it is assumed that delaminations are handled elsewhere in the finite element model (e.g., using cohesive interface elements).

Since `CDM_STATUS` is used for element deletion, always initialize `CDM_STATUS` to 1.

The initial condition for `CDM_phi0_12` and `CDM_phi0_13` are used to specify the initial fiber misalignment. One of the followings options is used depending on the initial condition specified for `CDM_phi0_12` and `CDM_phi0_13` as follows:
- *&phi;<sub>0</sub>* = 0 :: The value for *&phi;<sub>0</sub>* is calculated for shear instability. For 3-D kinking, *&phi;<sub>0,12</sub>* = *&phi;<sub>0,13</sub>* is required.
- *&phi;<sub>0</sub>* &le; 0.5 :: The value provided in the initial condition is used as the initial fiber misalignment.
- *&phi;<sub>0</sub>* = 1 :: A pseudo random uniform distribution varying spatially in the 1-direction is used. The spatial distribution algorithm relies on an uniform element size and fiber aligned mesh. The random number generator can be set to generate the same realizations or different realizations on multiple nominally identical analyses using the Boolean parameter `fkt_random_seed`. When using the random distribution for *&phi;<sub>0</sub>*, the characteristic length must be set to include 6 components: `*Characteristic Length, definition=USER, components=6`. For 3-D kinking, *&phi;<sub>0,12</sub>* = *&phi;<sub>0,13</sub>* is required.
- *&phi;<sub>0</sub>* = 2 :: Identical to *&phi;<sub>0</sub>* = 1, with the exception that a different realization is calculated for each ply. For 3-D kinking, *&phi;<sub>0,12</sub>* = *&phi;<sub>0,13</sub>* is required.
- *&phi;<sub>0</sub>* = 3 :: (Intended for use with 3-D FKT only) A pseudo random distribution varying spatially in the 1-direction is used with a 2-parameter lognormal distribution for the polar angle and a normal distribution for the azimuthal angle. The parameters starting with `fkt_init_misalignment` are used to control the polar and azimuthal distributions. Requires *&phi;<sub>0,12</sub>* = *&phi;<sub>0,13</sub>*.

Pre-existing damage can be modeled by creating an element set for the damaged region and specifying different initial conditions for this element set. For example, to create an intraply matrix crack with no out-of-plane orientation, the following initial conditions could be specified for the cracked elements:

    *Initial Conditions, type=SOLUTION
     damaged_elset,  1.d0,  0.d0,  0.d0,  0.d0,  0.d0,  0.d0,  0.d0,
              0.d0,  0.d0,     0,     1,  0.d0,  0.d0,  0.d0,  0.d0,
              0.d0,  0.d0,  0.d0,  0.d0


## Model parameters
A number of model parameters are used in the subroutine that are not directly related to material properties. These parameters are mostly related to internal error tolerances, iteration limits, etc. All parameters have default values and in most cases will not need to be modified.

The model parameters can be changed from their default values by including a file named `CompDam.parameters` in the working directory. An example parameters file is located in the `tests` directory. To create a job-specific parameters file, copy the example parameters file and rename it to `<job-name>.parameters`.


## Implicit solver compatibility
The repository includes a developmental capability to run the CompDam VUMAT in an Abaqus/Standard analysis using a wrapper, `for/vumatWrapper.for`, which translates between the UMAT and VUMAT user subroutine interfaces. The intended usage is for Abaqus/Standard runs with little or no damage.

### Usage
To run an analysis with CompDam in Abaqus/Standard, the following input deck template is provided. Note that 9 additional state variables are required.

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
    **  XC,       fXC,      GXC,      fGXC,       cl,     w_kb,     None,     mu
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
    0.d0,  0.d0,  0.d0,     1,  0.d0,  0.d0,  0.d0,  0.d0,
    0.d0,  0.d0,  0.d0,  0.d0,  0.d0,  0.d0,  0.d0,  0.d0,
    0.d0,  0.d0,  0.d0,  0.d0,  0.d0
    *Initial Conditions, Type=Field, Variable=1
    GLOBAL,  0.d0
    ** GLOBAL is an nset with all nodes attached to CompDam-enabled elements
    ** In each step, NLGEOM=YES must be used. This is NOT the default setting.

### Current limitations
As the `vumatWrapper` is a developmental capability, several important limitations exist at present:
1. The material Jacobian tensor is hard-coded in `for/vumatWrapper.for` for IM7/8552 elastic stiffnesses. A more general Jacobian is needed.
2. The material response can become inaccurate for large increments in rotations. If large rotations occur, small increments must be used. A cut-back scheme based on rotation increment size is needed.
3. Testing has been conducted on the C3D8R element type only.


## Example problems
The directory `examples/` includes example models that use CompDam along with corresponding files that defined the expected results (for use with Abaverify, following the same pattern as the test models in the `tests/` directory. The file `example_runner.py` can be used to automate submission of several models and/or for automatically post-processing the model results to verify that they match the expected results.


## Advanced debugging
Using an interactive debugger helps to identify issues in the Fortran code. Abaqus knowledge base article QA00000007986 describes the details involved. The following is a quick-start guide for direct application to CompDam.

Several statements for debugging need to be uncommented in the environment file. Follow these steps:
1. Copy your system environment file to your local working directory. For the example below, copy the environment file to the `tests` directory.
2. Edit the local environment file: uncomment lines that end with `# <-- Debugging`, `# <-- Debug symbols`, and `# <-- Optimization Debugging`

Run the job with the `-debug` and `-explicit` arguments. For example:
```
$ abaqus -j test_C3D8R_fiberTension -user ../for/CompDam_DGD.for -double both -debug -explicit
```

This command should open the [Visual Studio debugging software](https://msdn.microsoft.com/en-us/library/sc65sadd.aspx) automatically. Open the source file(s) to debug. At a minimum, open the file with the subroutine entry point `for/CompDam_DGD.for`. Set a break point by clicking in the shaded column on the left edge of the viewer. The break point will halt execution. Press <kbd>F5</kbd> to start the solver. When the break point is reached, a yellow arrow will appear and code execution will pause. Press <kbd>F5</kbd> to continue to the next break point, press <kbd>F11</kbd> to execute the next line of code following execution into function calls (Step Into), or press <kbd>F10</kbd> to execute the next line of code but not follow execution into function calls (Step Over).

To stop execution, close the Visual Studio window. Choose stop debugging and do not save your changes.

[More tips on debugging Fortran programs from Intel](https://software.intel.com/en-us/articles/tips-for-debugging-run-time-failures-in-intel-fortran-applications).

In case you must use remote ssh/scp access to run CompDam, a neat trick is to use the `Commands -> Keep Remote Directory up to Date...` option in [WinScp](https://winscp.net/eng/index.php). This feature can be used during development so that the CompDam files are edited locally and then automatically synced on a remote server for testing. The following mask (`Keep Remote Directory up to Date... -> Transfer Settings... -> File mask:`) can be used for syncing the source code: `*.for; *.py; *.inp; *.props; *.inc; Makefile; kind_map | tests/testOutput/; .git/; .vscode/`.

## Python extension module
CompDam can be compiled into a [Python extension module](https://docs.python.org/2/extending/extending.html), which allows many of the Fortran subroutines and functions in the `for` directory to be called from Python. The Python package [`f90wrap`](https://github.com/jameskermode/f90wrap) is used to automatically generate the Python extension modules that interface with the Fortran code. This Python extension module functionality is useful for development and debugging.

### Dependencies and setup
The python extension module requires some additional dependencies. First, the procedure only works on Linux using the bash shell. Windows users can use the [Windows Subsystem for Linux (WSL)](https://docs.microsoft.com/en-us/windows/wsl/about). In addition, `gfortran` 4.6+ is required. Type `gfortran --version` to check if you have this available. The remaining dependencies are python packages and can be installed as follows. The Python extension module works with Python 2 and 3; Python 2.7 is used for consistency with Abaqus in the following description.

Using [Conda](https://conda.io/projects/conda/en/latest/user-guide/getting-started.html) significantly simplifies the setup process, so it is assumed that you have a recent version of Conda available (see the [Conda installation guide](https://conda.io/projects/conda/en/latest/user-guide/install/index.html)). Further, the bash scripts described below include calls to Conda, so they will not work correctly without installing and configuring Conda as follows. Add the Conda-Forge channel by typing:
```
$ conda config --add channels conda-forge
```

Conda stores Python packages in containers called environments. Create a new environment:
```
$ conda create --name compdam python=2.7
```
and switch to your new environment:
```
$ source activate compdam
```
which will add `(compdam)` to the prompt. Install `numpy`, `matplotlib`, and `f90wrap` by typing:
```
(compdam) $ conda install numpy matplotlib f90wrap
```
After typing 'y' in response to the prompt asking if you would like to proceed, Conda will install `f90wrap` and all of its dependencies. This completes the setup process. These steps only need to be executed once.

Note, you can exit the Conda environment by typing `source deactivate`. When you open a new session, you will need to activate the Conda environment by typing `source activate compdam`.

### Example usage
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

### Associated scripts
In the `tests` directory the shell scripts `pyextmod_compile.sh` and `pyextmod_run.sh` are available to help streamline execution of the python extension module. These two scripts assume that Conda environment called `compdam` is available with `abaverify` and `f90wrap`. Both must be executed with the `-i` option. The script `pyextmod_run.sh` loads a debug file and executes the specified DGD routine. The DGD routine and the debug file are specified as arguments as follows:
```
$ bash -i pyextmod_run.sh dgdevolve <job-name>
```
The `<job-name>` is the Abaqus job name where it is assumed that the debug file resides in the testOutput folder with the name `job-name-1-debug-0.py`

In the `pyextmod` directory, the `helpers.py` file includes logic to load debug.py files.


## Summary of tests classes
This section includes a brief summary of each test implemented in the `tests` folder. The input deck file names briefly describe the test. All of the input decks start with `test_<elementType>_` and end with a few words describing the test. A more detailed description for each is given below:
- *elastic_fiberTension*: Demonstrates the elastic response in the 1-direction under prescribed extension. The 1-direction stress-strain curve has the modulus E1.
- *elastic_matrixTension*: Demonstrates the elastic response in the 2-direction under prescribed extension. The 2-direction stress-strain curve has the modulus E2.
- *elastic_simpeShear12*: Demonstrates the elastic response in the 1-2 plane. The 1-2 plane stress-strain curve has the module G12.
- *elementSize*: Verifies that the characteristic element lengths Lc1, Lc2, and Lc3 are being properly calculated.
- *error*: Verifies that analyses can cleanly terminate upon encountering an error within the user subroutine.
- *failureEnvelope_sig11sig22*: A parametric model in which *&sigma;<sub>11</sub>* is swept from *-X<sub>C</sub>* to *X<sub>T</sub>* and *&sigma;<sub>22</sub>* is swept from *-Y<sub>C</sub>* to *Y<sub>T</sub>* in order to re-create the corresponding failure envelope.
- *failureEnvelope_sig12sig22*: A parametric model in which *&tau;<sub>12</sub>* is swept from *0* to *S<sub>L</sub>* and *&sigma;<sub>22</sub>* is swept from *-Y<sub>C</sub>* to *Y<sub>T</sub>* in order to re-create the corresponding failure envelope.
- *failureEnvelope_sig12sig23*: A parametric model in which *&tau;<sub>12</sub>* is swept from *0* to *S<sub>L</sub>* and *&tau;<sub>23</sub>* is swept from *0* to *S<sub>T</sub>* in order to re-create the corresponding failure envelope.
- *fatigue_normal*: Demonstrates the traction-displacement curve of a cohesive law subjected to mode I fatigue loading.
- *fatigue_shear13*: Demonstrates the traction-displacement curve of a cohesive law subjected to shear fatigue loading in the 1-3 plane.
- *fiberCompression_BL*: Demonstrates the constitutive response in the 1-direction under prescribed shortening. The 1-direction stress-strain curve has a bilinear softening law. A conventional CDM approach to material degradation is used.
- *fiberCompression_FKT*: Demonstrates the constitutive response in the 1-direction under prescribed shortening. The fiber kink band model is used. `FF` indicates that the fiber failure criterion is enabled. `FN` indicates that fiber nonlinearity is enabled.
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

    @misc{CompDam,
    title={Comp{D}am - {D}eformation {G}radient {D}ecomposition ({DGD}), v2.5.0},
    author={Frank A. Leone and Andrew C. Bergan and Carlos G. D\'{a}vila},
    note={https://github.com/nasa/CompDam\_DGD},
    year={2019}
    }
