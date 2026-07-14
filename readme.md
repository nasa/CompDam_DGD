# CompDam - Deformation Gradient Decomposition (DGD)
This code is a continuum damage mechanics (CDM) material model intended for use with the Abaqus finite element code. This is a research code which aims to provide an accurate representation of mesoscale damage modes in fiber-reinforced polymer composite materials in finite element models in which each ply is discretely represented.

The CDM material model is implemented as an Abaqus/Explicit user subroutine (VUMAT) for the simulation of matrix cracks formed under tensile, compressive, and shear loading conditions and fiber fracture under tensile and compressive loading. Within CompDam, the emphasis of many recent developments has been on accurately representing the kinematics of composite damage. The kinematics of matrix cracks are represented by treating them as cohesive surfaces embedded in a deformable bulk material in accordance with the Deformation Gradient Decomposition (DGD) approach. Fiber tensile damage is modeled using conventional CDM strain-softening. This readme document provides an overview of the CompDam software. See the [List of Publications](list-of-publications.md) for references with details on the theory and use-cases.

This software may be used, reproduced, and provided to others only as permitted under the terms of the agreement under which it was acquired from the U.S. Government. Neither title to, nor ownership of, the software is hereby transferred. This notice shall remain on all copies of the software.

Copyright 2016 United States Government as represented by the Administrator of the National Aeronautics and Space Administration. No copyright is claimed in the United States under Title 17, U.S. Code. All Other Rights Reserved.

For any questions, please contact the developers:
- Andrew (Drew) Bergan | [andrew.c.bergan@nasa.gov](mailto:andrew.c.bergan@nasa.gov)
- Frank Leone   | [frank.a.leone@nasa.gov](mailto:frank.a.leone@nasa.gov)

## Table of contents
- [Getting started](#getting-started)
- [Element compatibility](#element-compatibility)
- [Model features](#model-features)
- [Material properties](#material-properties)
- [State variables](#state-variables)
- [Model parameters](#model-parameters)
- [Example problems](#example-problems)
- [Test models](#test-models)
- [Usage and best practices](#usage-and-best-practices)
- [Resources for developers](#resources-for-developers)
- [Contributing](#contributing)
- [Citing CompDam](#citing-compdam)


## Getting started

### Source code
The user subroutine source code is located in the [`for/`](for/) directory. The main entry point is [`CompDam_DGD.for`](for/CompDam_DGD.for).

### Prerequisites
Your abaqus installation must be configured to compile user subroutines, see the Abaqus documentation for details. This requires a fortran compiler and a C++ compiler, see [user subroutine prerequisites](usersubroutine-prerequisites.md) for a summary of the recommended software versions. Ryan Enos has published [a guide](https://doi.org/10.13140/RG.2.2.15382.11843) for installing and configuring Abaqus to compile user subroutines. Subroutine installation can be verified with the command `abaqus verify user_exp`. Note that CompDam requires MPI (which is not verified using `abaqus verify user_exp`). MPI must be installed and configured properly so that the MPI libraries can be linked by CompDam.

- CompDam versions <= 2.6.1 are compatible with Abaqus 2019-2023.
- CompDam versions >= 2.7.0 are compatible with Abaqus 2023+.

### Initial setup
For users that wish to run the code as-is, follow these steps:
1. Download CompDam
2. Rename the file `for/version.for.nogit` to `for/version.for`
3. [Compile the code to a shared library](#building-a-shared-library)

Note that users should compile the job on the hardware that will be used to run analyses.

For users that intend to modify the code, the following process is preferred. First `git clone` the repository.
After cloning the CompDam_DGD git repository, run the setup script file `setup.py` located in the repository root directory:
```
$ python setup.py
```

The main purpose of the setup.py script is to:
1) Set the `for/version.for` file, and
2) Add git-hooks that automatically update the `for/version.for`.

In the event that you do not have access to python, rename `for/version.for.nogit` to `for/version.for` manually. The additional configuration done by `setup.py` is not strictly required.

### Abaqus environment file settings
The `abaqus_v6.env` file must have `/fpp`, `/Qmkl:sequential`, and `/free` in the `ifort` command where the format for Windows is used.
The corresponding Linux format is: `-fpp`, `-free`, and `-mkl=sequential`.
The `/fpp` option enables the Fortran preprocessor, which is required for the code to compile correctly.
The `/free` option sets the compiler to free-formatting for the source code files.
The `/Qmkl:sequential` enables the [Intel Math Kernel Library (MKL)](https://software.intel.com/en-us/mkl), which provides optimized and verified functions for many mathematical operations. 
The MKL is used in this code for calculating eigenvalues and eigenvectors.

A sample environment file is provided: [`tests/abaqus_v6.env`](tests/abaqus_v6.env).

### Submitting a job
CompDam is an Abaqus/Explicit VUMAT. Please refer to the Abaqus documentation for the general instructions on how to submit finite element analyses using user subroutines. Please see the [example input file statements](#example-input-file-statements) for details on how to interface with this particular VUMAT subroutine.

Analyses with this code **must** be run in double precision. Some of the code has double precision statements and variables hard-coded, so if Abaqus/Explicit is run in single precision, compile-time errors will arise. When submitting an Abaqus/Explicit job from the command line, double precision is specified by including the command line argument `double=both`.

Geometric nonlinearity must be used in each analysis step. Geometric nonlinearity being turned on is the default in Abaqus/Explicit analysis steps. Geometric nonlinearity can be explicitly stated in the input deck `*Step` command with the keyword option `nlgeom=YES`.

For example, run the test model [`test_C3D8R_elastic_fiberTension`](tests/test_C3D8R_elastic_fiberTension.inp) in the [`tests/`](tests/) directory with the following command:
```
$ abaqus job=test_C3D8R_elastic_fiberTension user=../for/CompDam_DGD.for double=both
```

### Example input file statements
The snippet below includes a basic use of CompDam for continuum elements. The [`tests/`](tests/) and [`examples/`](examples/) directories include many additional cases.

    *Solid Section, elset=all, orientation=Ori-1, material=IM7-8552, controls=control_name
    *Section controls, name=control_name, distortion control=YES
    **
    *Material, name=IM7-8552
    *Density
     1.57e-09,
    *Depvar
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
    ** feature flags, 2,   thickness, 4, 5, 6, 7, 8
              100001,  ,            ,  ,  ,  ,  ,  ,
    **
    **  9         10        11        12        13        14        15        16
    **  E1,       E2,       G12,      nu12,     nu23,     YT,       SL        GYT
        171420.0, 9080.0,   5290.0,   0.32,     0.52,     62.3,     92.30,    0.277,
    **
    **  17        18        19        20        21        22        23        24
    **  GSL,      eta_BK,   YC,       alpha0    E3,       G13,      G23,      nu13,
        0.788,    1.634,    199.8,    0.925,      ,          ,         ,          ,
    **
    **  25        26        27        28        29        30        31        32
    **  alpha11,  alpha22,  alpha_PL, n_PL,     XT,       fXT,      GXT,      fGXT,
               ,         ,          ,     ,       ,          ,         ,          ,
    **
    **  33        34        35        36        37        38        39        40
    **  XC,       fXC,      GXC,      fGXC,     cl,       w_kb,     T_sf,     mu
          ,          ,         ,          ,       ,           ,         ,     0.3
    **
    *Initial Conditions, type=SOLUTION
     elset_name,  0.d0,  0.d0,  0.d0,  0.d0,  0.d0,  0.d0,  0.d0,
     0.d0,        0.d0,  0.d0,     1,  0.d0,  0.d0,  0.d0,  0.d0,
     0.d0,        0.d0,  0.d0,  0.d0

The following snippet illustrates an example of using CompDam for cohesive elements. Note that the 21st material property is the cohesive element penalty stiffess.

    *Cohesive Section, elset=COHESIVE, material=IM7-8552, response=TRACTION SEPARATION, thickness=SPECIFIED
    1.0,  ** Must be 1, or the dissipated energy will be incorrect
    *Material, name=IM7-8552
    *Density
    1.57e-09,
    *User material, constants=40, effmod
    ** feature flags, 2,   thickness,           4, 5, 6, 7, 8
              200000,  ,            ,            ,  ,  ,  ,  ,
    **
    **  9         10        11        12        13        14        15        16
    **  E1,       E2,       G12,      nu12,     nu23,     YT,       SL        GYT
          ,         ,          ,          ,         ,     62.3,     92.30,    0.277,
    **
    **  17        18        19        20        21        22        23        24
    **  GSL,      eta_BK,   YC,       alpha0    Pen,      G13,      G23,      nu13,
        0.788,    1.634,    199.8,    0.925,    1.d6,        ,         ,          ,
    **
    **  25        26        27        28        29        30        31        32
    **  alpha11,  alpha22,  alpha_PL, n_PL,     XT,       fXT,      GXT,      fGXT,
               ,         ,          ,     ,       ,          ,         ,          ,
    **
    **  33        34        35        36        37        38        39        40
    **  XC,       fXC,      GXC,      fGXC,     cl,       w_kb,     T_sf,     mu
          ,          ,         ,          ,           ,       ,         ,     0.3
    **
    *Depvar
    19,
    1, COH_dmg
    2, COH_delta_s1
    3, COH_delta_n
    4, COH_delta_s2
    5, COH_B
    9, COH_FI
    11, COH_status
    15, COH_slide1
    16, COH_slide2
    *Initial Conditions, Type=Solution
    COHESIVE,  0.d0,  0.d0,  0.d0,  0.d0,  0.d0,  0.d0,  0.d0,
    0.d0,      0.d0,  0.d0,     1,  0.d0,  0.d0,  0.d0,  0.d0,
    0.d0,      0.d0,  0.d0,  0.d0

### Running tests
Test models are available in the [`tests/`](tests/) directory. The tests are useful for demonstrating the capabilities of the VUMAT as well as to verify that the code performs as intended. Try running some of the test cases to see how the code works. The test cases can be submitted as a typical Abaqus job using the Abaqus command line arguments. Tests generally have an `_expected.py` file that defines the expected result from the analysis. If [abaverify](https://github.com/nasa/abaverify) is installed, the tests can be executed and results checked automatically using the [`test_runner.py`](tests/test_runner.py) script.

### Building a shared library
CompDam_DGD can be built into a shared library file. Follow these steps:
1. Place a copy of the Abaqus environment file in the `for/` directory (you can copy-paste [`tests/abaqus_v6.env`](tests/abaqus_v6.env))
2. In Linux, and when using Abaqus versions prior to 2017, rename `CompDam_DGD.for` to `CompDam_DGD.f`
3. From the `for/` directory, run:
```
$ abaqus make library=CompDam_DGD
```
This command will create shared libraries for the operating system it is executed on (`.dll` for Windows and `.so` for Linux).

When using a pre-compiled shared library, it is only necessary to specify the location of the shared library files in the environment file (the compiler options are not required). To run an analysis using a shared library, add `usub_lib_dir = "<full path to shared library file>"` to the Abaqus environment file in the Abaqus working directory. After that, you don't need to include `user=../for/CompDam_DGD.for` when submitting abaqus jobs.


## Element compatibility
CompDam_DGD has been developed and tested using the Abaqus three-dimensional (3-D), reduced-integration `C3D8R` hexahedral solid elements and 3-D `COH3D8` cohesive elements. Limited testing has been performed using the following element types:
- `CPS4R` plane stress element, reduced-integration
- `C3D8` 3-D hexahedral solid element
- `C3D6` 3-D wedge solid element
- `COH2D4` two-dimensional (2-D) cohesive element

Because CompDam_DGD is a material model, it is expected to be compatible with continuum elements generally. However, users are advised to perform tests with any previously untested element types before proceeding to use CompDam_DGD in larger structural models. Axisymmetric elements are not currently supported.


## Model features
The CompDam_DGD material model implements a variety of features that can be enabled or disabled by the user. An overview of these features is provided in this section. The material properties required for each feature are listed. References are provided to more detailed discussions of the theoretical framework for each feature.

### Elastic response
The composite materials modeled with CompDam_DGD can be defined assuming either [transverse isotropy](https://www.efunda.com/formulae/solid_mechanics/mat_mechanics/hooke_iso_transverse.cfm) or [orthotropy](https://www.efunda.com/formulae/solid_mechanics/mat_mechanics/hooke_orthotropic.cfm). For a transversely isotropic material definition, the following properties must be defined: E<sub>1</sub>, E<sub>2</sub>, G<sub>12</sub>, $\nu$<sub>12</sub>, and $\nu$<sub>23</sub>. For an orthotropic material definition, the following additional properties must be defined: E<sub>3</sub>, G<sub>13</sub>, G<sub>23</sub>, and $\nu$<sub>13</sub>.

### Matrix damage
Tensile and compressive matrix damage is modeled by embedding cohesive laws to represent cracks in the material according to the deformation gradient decomposition method of [Leone (2015)](https://doi.org/10.1016/j.compositesa.2015.06.014). The matrix crack normals can have any orientation in the 2-3 plane, defined by the angle `CDM_alpha`. The mixed-mode behavior of matrix damage initiation and evolution is defined according to the Benzeggagh-Kenane law. The initiation of compressive matrix cracks accounts for friction on the potential crack surface according to the LaRC04 failure criteria.

The following material properties are required for the prediction of matrix damage: YT, SL, GYT, GSL, eta_BK, YC, and alpha0. The state variables related to matrix damage are `CDM_d2`, `CDM_FIm`, `CDM_B`, `CDM_alpha`, `CDM_Fb1`, `CDM_Fb2`, and `CDM_Fb3`.

### Thermal expansion and residual stress
The thermal strain is calculated by multiplying the 1-, 2-, and 3-direction coefficients of thermal expansion by the current temperature difference, &Delta;T. The current temperature is provided by the Abaqus solver. A reference stress-free temperature `T_sf` can be defined as a material property input. The mechanical strain is found by subtracting the thermal strain from the total strain.

The required material properties are the coefficients of thermal expansion in the 1 and 2 directions. If the 3-direction coefficient of thermal expansion is not provided, it is assumed that the 2- and 3-direction coefficients of thermal expansion are equal. The stress-free temperature is assumed equal to zero unless provided as an optional input.

Hygroscopic strains are not accounted for. If the effects of hygroscopic expansion are to be modeled, it is recommended to smear the hygroscopic and thermal expansion coefficients to approximate the response using the solver-provided temperature value.

### Shear nonlinearity
Three approaches to modeling the matrix nonlinearity are available: Ramberg-Osgood plasticity, Schapery theory, and Schaefer plasticity. These three methods are mutually exclusive and optional. In most analyses, shear nonlinearity is not used.

<details>

<summary>Details on the models for shear nonlinearity</summary>

#### Ramberg-Osgood plasticity
Shear nonlinearity in the 1-2 and/or the 1-3 plane can be modeled using the [Ramberg-Osgood equation](https://en.wikipedia.org/wiki/Ramberg%E2%80%93Osgood_relationship), with its parameters selected to fit experimental data. As applied herein, the Ramberg-Ogsood equation is written in the following form for the 1-2 plane:

*&gamma;*<sub>12</sub> = [*&tau;*<sub>12</sub> + *&alpha;*<sub>PL</sub>sign(*&tau;*<sub>12</sub>)|*&tau;*<sub>12</sub>|<sup>*n*<sub>PL</sub></sup>]/*G*<sub>12</sub>

where *&gamma;*<sub>12</sub> is the shear strain and *&tau;*<sub>12</sub> is the shear stress. Likewise, the expression for the 1-3 plane is

*&gamma;*<sub>13</sub> = [*&tau;*<sub>13</sub> + *&alpha;*<sub>PL</sub>sign(*&tau;*<sub>13</sub>)|*&tau;*<sub>13</sub>|<sup>*n*<sub>PL</sub></sup>]/*G*<sub>13</sub>

Prior to the initiation of matrix damage (i.e., `CDM_d2 = 0`), the nonlinear shear response due to the above equation is plastic, and the unloading/reloading slope is unchanged. No pre-peak nonlinearity is applied to the matrix tensile or compressive responses (i.e., *&sigma;<sub>22</sub>*).

The required material inputs are the two parameters in the above equation: *&alpha;*<sub>PL</sub> and *n*<sub>PL</sub>. Note that the same constants are used for the 1-2 and 1-3 planes under the assumption of transverse isotropy (see [Seon et al. (2017)](https://doi.org/10.12783/asc2017/15267)). For the 1-2 plane, the state variables `CDM_Plas12` and `CDM_Inel12` are used to track the current plastic shear strain and the total amount of inelastic plastic shear strain that has occurred through the local deformation history, respectively. For cases of monotonic loading, `CDM_Plas12` and `CDM_Inel12` should have the same magnitude. Likewise, the state variables `CDM_Plas13` and `CDM_Inel13` are utilized for the 1-3 plane. The [feature flags](#controlling-which-features-are-enabled) can be used to enable this Ramberg-Osgood model in the 1-2 plane, 1-3 plane, or both planes.

#### Schapery micro-damage
Matrix nonlinearity in the 1-2 plane can also be modeled using Schapery theory, in which all pre-peak matrix nonlinearity is attributed to the initiation and development of micro-scale matrix damage. With this approach, the local stress/strain curves will unload to the origin, and not develop plastic strain. A simplified version of the approach of [Pineda and Waas (2011)](https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20120000914.pdf) is here applied. The micro-damage functions *e<sub>s</sub>* and *g<sub>s</sub>* are limited to third degree polynomials for ease of implementation. As such, four fitting parameters are required for each of *e<sub>s</sub>* and *g<sub>s</sub>* to define the softening of the matrix normal and shear responses to micro-damage development.

*e<sub>s</sub>*(*S<sub>r</sub>*) = *e<sub>s0</sub>* + *e<sub>s1</sub>S<sub>r</sub>* + *e<sub>s2</sub>S<sub>r</sub>*<sup>2</sup> + *e<sub>s3</sub>S<sub>r</sub>*<sup>3</sup>

*g<sub>s</sub>*(*S<sub>r</sub>*) = *g<sub>s0</sub>* + *g<sub>s1</sub>S<sub>r</sub>* + *g<sub>s2</sub>S<sub>r</sub>*<sup>2</sup> + *g<sub>s3</sub>S<sub>r</sub>*<sup>3</sup>

where *S<sub>r</sub>* is the micro-damage reduced internal state variable. *S<sub>r</sub>* is stored in the 12<sup>th</sup> state variable slot, replacing `CDM_Plas12`, when Schapery theory is used in a model. The 13<sup>th</sup> state variable slot is not used when Schapery micro-damage is used.

#### Schaefer
Shear nonlinearity in the 1-2 plane can be modeled using the Schaefer prepeak model. In this model, effective plastic strain is related to effective plastic stress through the power law:

*&epsilon;*<sub>plastic</sub> = *A &sigma;*<sup>*n*</sup>

Additionally, a yield criterion function (effective stress) is defined as:

*f* = *&sigma;* = (*S*<sub>22</sub> + *a*<sub>6</sub> *S*<sub>12</sub><sup>2</sup>)<sup>1/2</sup> + *b*<sub>2</sub> *S*<sub>22</sub>

where *a*<sub>6</sub>, *b*<sub>2</sub>, *A* and *n* are material constants needed for the Schaefer prepeak material model.

The above two equations are used in concert to determine plastic strain through the relationship:

*&epsilon;* <sub>plastic</sub> = *n A f* <sup>*n* - 1</sup> (&part;*f*/&part;*S*<sub>*i*</sub>) (&part;*f*/&part;*S*<sub>j</sub>)

*f* (i.e., schaefer_f) and the tensorial plastic strain determined by the nonlinearity model are stored as state variables (27 through 32 for plastic strain and 33 for *f*)

</details>

### Fiber tensile damage
A continuum damage mechanics model similar to the work of [Maimí et al. (2007)](https://doi.org/10.1016/j.mechmat.2007.03.005) is used to model tensile fiber damage evolution. The model utilizes a non-interacting maximum strain failure criterion, and bilinear softening after the initiation of failure. The area under the stress-strain curve is equal to the fracture toughness divided by the element length normal to fracture, i.e., `CDM_Lc1`. The required material properties are: XT, fXT, GXT, and fGXT, where fXT and fGXT are ratios of strength and fracture toughness for bilinear softening, defined as the *n* and *m* terms in equations (25) and (26) of [Dávila et al. (2009)](https://doi.org/10.1007/s10704-009-9366-z). To model a linear softening response, both fXT and fGXT should be set equal to 0.5.

### Fiber compression damage

#### Model 1: Max strain, bilinear softening (BL)
Same model as in tension, but for compression. Assumes maximum strain failure criterion and bilinear softening. The required material properties are: XC, fXC, GXC, and fGXC.

Load reversal assumptions from [Maimí et al. (2007)](https://doi.org/10.1016/j.mechmat.2007.03.006).

#### Model 2: Placeholder

#### Model 3, 4, 5: Fiber kinking theory (FKT)
A model based on Budiansky's fiber kinking theory from [Budiansky (1983)](https://doi.org/10.1016/0045-7949(83)90141-4), [Budiansky and Fleck (1993)](https://doi.org/10.1016/0022-5096(93)90068-Q), and [Budiansky et al. (1998)](https://doi.org/10.1016/S0022-5096(97)00042-2) implemented using the DGD framework to model in-plane (1-2) and/or out-of-plane (1-3) fiber kinking. The model is described in detail in [Bergan et al. (2018)](https://doi.org/10.2514/6.2018-1221), [Bergan and Jackson (2018)](https://doi.org/10.12783/asc33/26003), and [Bergan (2019)](https://doi.org/10.2514/6.2019-1548). The model accounts for fiber kinking due to shear instability by considering an initial fiber misalignment, nonlinear shear stress-strain behavior via Ramberg-Osgood, and geometric nonlinearity. Fiber failure can be introduced by specifying a critical fiber rotation angle.

<details>

<summary>Details on the FKT model</summary>

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

where *&phi;* is the current fiber rotation. Once the fiber failure criterion is satisfied, the plastic shear strain is held constant. The value for *&phi;*<sub>ff,c</sub> is defined as the [model parameter](#model-parameters) `fkt_fiber_failure_angle` since it is not a well-defined material property. The fiber failure criterion is disabled when *&phi;*<sub>ff,c</sub> < 0. The same angle is used for in-plane and out-of-plane kinking.

The fiber kinking theory model implemented here is preliminary and has some known shortcomings and caveats:
- The model has only been tested for C3D8R. Limited application with C3D6 demonstrated issues. No testing has been completed for other element types.
- The interaction of this model with matrix cracking has not been fully tested and verified.
- No effort has been made to model progressive crushing.
- Other mechanisms of fiber compressive failure (e.g., shear driven fiber breaks) are not accounted for. An outcome of this is that the model predicts the material does not fail if shear deformation is fully constrained.
- No special consideration for load reversal has been included.

Relevant single element tests are named starting with `test_C3D8R_fiberCompression_FKT`.

</details>

### Friction
Friction is modeled on the damaged fraction of the cross-sectional area of matrix cracks and delaminations using the approach of [Alfano and Sacco (2006)](https://doi.org/10.1002/nme.1728). The coefficient of friction *&mu;* must be defined to account for friction on the failed crack surface.

The amount of sliding which has taken place in the longitudinal and transverse directions are stored in state variables `CDM_slide1` and `CDM_slide2`, respectively.

### Fatigue analyses
The cohesive fatigue constitutive model in CompDam can predict the initiation and the propagation of matrix cracks and delaminations as a function of fatigue cycles. The analyses are conducted such that the applied load (or displacement) corresponds to the maximum load of a fatigue cycle. The intended use is that the maximum load (or displacement) is held constant while fatigue damage develops with increasing step time. The constitutive model uses a specified load ratio *R*<sub>min</sub>/*R*<sub>max</sub>, the solution increment, and an automatically-calculated cycles-per-increment ratio to accumulate the damage due to fatigue loading. The cohesive fatigue model response is based on engineering approximations of the endurance limit as well as the Goodman diagram. This approach can predict the stress-life diagrams for crack initiation, the Paris law regime, as well as the transient effects of crack initiation and stable tearing.

<details>

<summary>Details on the fatigue model</summary>

A detailed description of the initial development of the novel cohesive fatigue law is available in a [2018 NASA technical paper by Carlos Dávila](https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20180004395.pdf). An updated fatigue damage accumulation function is derived in the 2020 NASA technical paper ["Evaluation of Fatigue Damage Accumulation Functions for Delamination Initiation and Propagation" by Dávila et al.](https://ntrs.nasa.gov/citations/20200003113), which is the fatigue damage accumulation function that is applied herein.

#### Usage
The fatigue capability of CompDam is disabled by default. To run a fatigue analysis, one of the analysis steps must be identified as a fatigue step. A step is identified as a fatigue step by setting the [model parameter](#model-parameters) `fatigue_step` to the target step number, e.g., `fatigue_step = 2` for the second analysis step to be a fatigue step. The first analysis step cannot be a fatigue step, as the model is assumed to be initially unloaded.

The load ratio *R*<sub>min</sub>/*R*<sub>max</sub> has a default value of 0.1, and can be changed using the [model parameter](#model-parameters) `fatigue_R_ratio`.

Material properties 5-8 are specific to fatigue, see [Fatigue properties](#fatigue-properties).

An [example of a double cantilever beam subjected to fatigue](examples/interlaminar/test_DCB_3D_fatigue.inp) under displacement-control is provided. The geometry and conditions of this example problem correspond to the results presented in Figure 20 of [Dávila (2018)](https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20180004395.pdf).

#### Interpreting the results of a fatigue analysis
Within a fatigue step, each solution increment represents either a number of fatigue cycles or a fractional part of a single fatigue cycle. During the solution, the number of fatigue cycles per solution increment changes based on the maximum amount of energy dissipation in any single element. If the rate of energy dissipation is too high (as defined by the [model parameter](#model-parameters) `fatigue_damage_max_threshold`), the increments-to-cycles ratio is decreased. If the rate of energy dissipation is too low (as defined by the model parameter `fatigue_damage_min_threshold`), the increments-to-cycles ratio is increased. The model parameter `cycles_per_increment_init` defines the initial ratio of fatigue cycles per solution increment. Any changes to increments-to-cycles ratio are logged in an additional output file ending in `_inc2cycles.log`, with columns for the fatigue step solution increment, the updated increments-to-cycles ratio, and the accumulated fatigue cycles.

</details>

### Finite stress and strain definitions
The strain is calculated using the deformation gradient tensor provided by the Abaqus solver. The default strain definition used is the Green-Lagrange strain:

***E*** = (***F***<sup>T</sup>***F*** - ***I***)/2

Hooke's law is applied using the Green-Lagrange strain to calculate the 2<sup>nd</sup> Piola-Kirchhoff stress ***S***.

### Fiber nonlinearity
Nonlinear elastic behavior in the fiber direction can be introduced with the material property c<sub>*l*</sub>. The expression used follows [Kowalski (1988)](https://doi.org/10.1520/STP26136S):

*E<sub>1</sub>* = *E<sub>1</sub>*(1 + c<sub>*l*</sub>*&epsilon;*<sub>11</sub>)

By default, fiber nonlinearity is disabled by setting c<sub>*l*</sub> = 0 (or leaving it undefined).


## Material properties
A set of material properties must be defined for the material of interest. This section describes how to specify the material properties.

Material properties are defined in the input deck. Any optional material property can be left blank and the corresponding feature(s) will be disabled. The ordering of the material properties for the input deck definition is given in the first (#) column of the [table of material properties](#table-of-material-properties).

### Table of material properties
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
| 18 | *&eta;*                | eta_BK   | BK exponent for mode-mixity                   | -                                              | 0 < *&eta;* < &infin;                    |                 |
| 19 | *Y<sub>C</sub>*        | YC       | Transverse compressive strength               | F/L<sup>2</sup>                                | 0 < *Y<sub>C</sub>* < &infin;            | ASTM D3410      |
| 20 | *&alpha;<sub>0</sub>*  | alpha0   | Fracture plane angle for pure trans. comp.    | Radians                                        | 0 &le; *&alpha;<sub>0</sub>* &le; &pi;/2 |                 |
|    | ------                 |          |                                               |                                                |                                          |                 |
| 21 | *E<sub>3</sub>*        | E3       | 3-direction Young's modulus                   | F/L<sup>2</sup>                                | 0 < *E<sub>3</sub>* < &infin;            |                 |
| 22 | *G<sub>13</sub>*       | G13      | Shear modulus in 1-3 plane                    | F/L<sup>2</sup>                                | 0 < *G<sub>13</sub>* < &infin;           |                 |
| 23 | *G<sub>23</sub>*       | G23      | Shear modulus in 1-2 plane                    | F/L<sup>2</sup>                                | 0 < *G<sub>23</sub>* < &infin;           |                 |
| 24 | *&nu;<sub>13</sub>*    | nu13     | Poisson's ratio in 2-3 plane                  | -                                              | 0 &le; *&nu;<sub>13</sub>* &le; 1        |                 |
|    | ------                 |          |                                               |                                                |                                          |                 |
| 25 | *&alpha;<sub>11</sub>* | alpha11  | Coefficient of thermal expansion, fiber       | 1/&deg;                                        | -1 &le; *&alpha;<sub>11</sub>* &le; 1    |                 |
| 26 | *&alpha;<sub>22</sub>* | alpha22  | Coefficient of thermal expansion, matrix      | 1/&deg;                                        | -1 &le; *&alpha;<sub>22</sub>* &le; 1    |                 |
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
| 37 | *c<sub>l</sub>*        | cl       | Fiber nonlinearity coefficient                | -                                              | 0 &le; *c<sub>l</sub>* &le; 33           |                 |
| 38 | *w<sub>kb</sub>*       | w_kb     | Width of the kink band                        | L                                              | 0 &le; *w<sub>kb</sub>* < &infin;        |                 |
| 39 | *T<sub>sf</sub>*       | T_sf     | Stress-free temperature                       | &deg;                                          | -&infin; < *T<sub>sf</sub>* < &infin;    |                 |
| 40 | *&mu;*                 | mu       | Coefficient of friction                       | -                                              | 0 &le; *&mu;* &le; 1                     |                 |
|    | ------                 |          |                                               |                                                |                                          |                 |
| 41 | *a<sub>6</sub>*        | a6       | Schaefer a6                                   | -                                              | -&infin; < *a<sub>6</sub>* < &infin;     |                 |
| 42 | *b<sub>2</sub>*        | b2       | Schaefer b2                                   | -                                              | -&infin; < *b<sub>2</sub>* < &infin;     |                 |
| 43 | *A*                    | A        | Schaefer A                                    | -                                              | -&infin; < *A* < &infin;                 |                 |
| 44 | *n*                    | n        | Schaefer n                                    | -                                              | -&infin; < *n* < &infin;                 |                 |
|    | ------                 |          |                                               |                                                |                                          |                 |
| 45 | *e<sub>s0</sub>*       | es0      | Schapery microdamage parameter es0            | -                                              | -&infin; < *e<sub>s0</sub>*  < &infin;   |                 |
| 46 | *e<sub>s1</sub>*       | es1      | Schapery microdamage parameter es1            | -                                              | -&infin; < *e<sub>s1</sub>*  < &infin;   |                 |
| 47 | *e<sub>s2</sub>*       | es2      | Schapery microdamage parameter es2            | -                                              | -&infin; < *e<sub>s2</sub>*  < &infin;   |                 |
| 48 | *e<sub>s3</sub>*       | es3      | Schapery microdamage parameter es3            | -                                              | -&infin; < *e<sub>s3</sub>*  < &infin;   |                 |
|    | ------                 |          |                                               |                                                |                                          |                 |
| 49 | *g<sub>s0</sub>*       | gs0      | Schapery microdamage parameter gs0            | -                                              | -&infin; < *g<sub>s0</sub>*  < &infin;   |                 |
| 50 | *g<sub>s1</sub>*       | gs1      | Schapery microdamage parameter gs1            | -                                              | -&infin; < *g<sub>s1</sub>*  < &infin;   |                 |
| 51 | *g<sub>s2</sub>*       | gs2      | Schapery microdamage parameter gs2            | -                                              | -&infin; < *g<sub>s2</sub>*  < &infin;   |                 |
| 52 | *g<sub>s3</sub>*       | gs3      | Schapery microdamage parameter gs3            | -                                              | -&infin; < *g<sub>s3</sub>*  < &infin;   |                 |
|    | ------                 |          |                                               |                                                |                                          |                 |
| 53 | *&alpha;<sub>33</sub>* | alpha33  | Coefficient of thermal expansion, thickness   | 1/&deg;                                        | -1 &le; *&alpha;<sub>33</sub>* &le; 1    |                 ||

Notes:
- The inputs 9 through 13 (above the ===) are always required
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
- Providing inputs via a `.props` file is no longer supported for CompDam > 2.6.1.

Two material cards are included for the IM7/8552 material system:
 - 134 g/m^2 [`IM7-8552-C3D.inp`](tests/IM7-8552-C3D.inp), see also [`IM7-8552-COH.inp`](tests/IM7-8552-COH.inp) for how these properties are applied for cohesive elements.
 - 190 g/m^2 [`IM7-8552-C3D-190gsm.inp`](tests/IM7-8552-C3D-190gsm.inp).

### Controlling which features are enabled
The feature flags are defined in the input deck as the first material property input. While the feature flags are not material properties per se, they are used in controlling the behavior of the material model.

Model features can be enabled or disabled by two methods. The first method is specifying only the material properties required for the features you would like to enable. CompDam_DGD disables any feature for which all of the required material properties have not been assigned. If an incomplete set of material properties are defined for a feature, a warning is issued.

The second method is by specifying the status of each feature directly as a material property in the input deck. Each feature of the subroutine is controlled by a position in an integer, where 0 is disabled and >=1 is enabled. In cases where mutually exclusive options are available, numbers greater than 1 are used to specify the particular option to use.

The positions correspond to the features as follows:
- Position 1: Matrix damage (1=intralaminar cracking in solid elements, 2=interlaminar cracking in cohesive elements)
- Position 2: Shear nonlinearity (1=Ramberg-Osgood 1-2 plane, 2=Schapery, 3=Ramberg-Osgood 3-D, 4=Ramberg-Osgood 1-3 plane, 5=Schaefer || more information [here](#shear-nonlinearity))
- Position 3: Fiber tensile damage
- Position 4: Fiber compression damage (1=max strain, 2=N/A, 3=FKT-12, 4=FKT-13, 5=FKT-3D || more information [here](#fiber-compression-damage))
- Position 5: Energy output contribution (0=all mechanisms, 1=only fracture energy, 2=only plastic energy)
- Position 6: Friction

For example, `101000` indicates that the model will run with matrix damage and fiber tension damage enabled; `120001` indicates that the model will run with matrix damage, in-plane shear nonlinearity using Schapery theory, and friction; and `200000` indicates that the material model is being applied to cohesive elements.

### Finite thickness cohesive elements
The finite thickess cohesive formulation from [Sarrado et al. 2016](https://doi.org/10.1016/j.engfracmech.2016.03.020) can be activated by defining material properties `G13` and `G23` when the matrix damage feature flag = 2.
In this case, CompDam assumes the cohesive elements have some finite thickness, and a simplified DGD routine is used to decompose the element deformation into bulk and crack components.
Finite thickness cohesive elements are primarily intended for modeling cases where the elastic response of a ply interface is important, for example an adhesive layer.

### Definition of thickness
Material input #3 is the thickness, which is defined as the length along the thickness-direction associated with the current integration point.
This input is used for 2-D plane stress elements, conventional shell elements, and, in special cases, for zero-thickness cohesive elements.
Therefore, in most cases, the thickness input can be empty.

For 2-D plane stress elements (CPS4R) and shell elements (S4R), the thickness should be specified with the geometric value. See [test_CPS4R_elementSize.inp](tests/test_CPS4R_elementSize.inp) and [test_S4R_elementSize.inp](tests/test_S4R_elementSize.inp).

For cohesive elements, typically the element is defined as zero-thickness (i.e., the nodal coordinates on the top and bottom surfaces are coincident).
Additionally, the consitutive thickness is defined = 1 in `*Cohesive section`.
Therefore, CompDam internally assumes the constitutive thickness = 1 (strain=displacement) for zero-thickness cohesive elements when the material input thickness is not defined.
However, a side-effect is that the mass associated with the cohesive element is relatively large since the ply thickness is typically << 1.
In dynamic loading problems, it is possible to modify the constitutive thickness such that the mass associated with the cohesive elements is small enough to have a negligible influence on the dynamic response.
See [Dynamic loading and cohesive elements](#dynamic-loading-and-cohesive-elements).
When setting the constitutive thickness ≠ 1, the same value must be defined in both `*Cohesive section` and the material input #3 (thickness), or the dissipated energy will be incorrect.

This thickness input does not affect the performance of 3-D solid elements or finite-thickness cohesive elements.
The three characteristic element lengths of solid elements and finite-thickness cohesive elements are calculated using the VUCHARLENGTH subroutine based on the element nodal coordinates.

### Fatigue properties
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
| 26| `CDM_reserved`   | Reserved                                                              |
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

Note that the state variable names are defined using `*Depvar` and do not need to be identical to the names given here.
For example, in a model with both cohesive and continuum elements using CompDam, a user may want to name the first state variable `CD_dmg` for both element types instead of `CDM_d2` and `COH_dmg` for continuum and cohesive element, respectively.
The advantage is that damage in both element types can be visualized together in Abaqus/Viewer. 

### Initial conditions
All state variables should be initialized using the `*Initial conditions` command. As a default, all state variables should be initialized as zero, except `CDM_alpha`, `CDM_STATUS`, `CDM_phi0_12`, and `CDM_phi0_13`.

The initial condition for `CDM_alpha` can be used to specify a predefined angle for the cohesive surface normal.
`CDM_alpha` is measured from the 2-axis rotating about the 1-direction.
To specify a predefined (fixed) `CDM_alpha`, set the initial condition for `CDM_alpha` to an integer (degrees) and the [model parameter](#model-parameters) `alpha_search = FALSE`.
When the [model parameter](#model-parameters) `alpha_search = TRUE` (its default), the subroutine evaluates cracks every 10 degrees (by default) in the 2-3 plane to find the crack angle corresponding to the maximum failure criterion value.
The amount by which alpha is incremented when searching is controlled via the model parameter `alpha_inc`.
Note that `CDM_alpha = 90` only occurs when `CDM_alpha` is initialized as 90; the value of 90 is ignored in the search to find the correct initiation angle since it is assumed that delaminations are handled elsewhere in the finite element model (e.g., using cohesive interface elements).

Since `CDM_STATUS` and `COH_STATUS` are used for element deletion, always initialize to 1.

The initial condition for `CDM_phi0_12` and `CDM_phi0_13` are used to specify the initial fiber misalignment. One of the followings options is used depending on the initial condition specified for `CDM_phi0_12` and `CDM_phi0_13` as follows:
- *&phi;<sub>0</sub>* = 0 :: The value for *&phi;<sub>0</sub>* is calculated for shear instability. For 3-D kinking, *&phi;<sub>0,12</sub>* = *&phi;<sub>0,13</sub>* is required.
- *&phi;<sub>0</sub>* &le; 0.5 :: The value provided in the initial condition is used as the initial fiber misalignment.
- *&phi;<sub>0</sub>* = 1 :: A pseudo random uniform distribution varying spatially in the 1-direction is used. The spatial distribution algorithm relies on an uniform element size and fiber aligned mesh. The random number generator can be set to generate the same realizations or different realizations on multiple nominally identical analyses using the Boolean parameter `fkt_random_seed`. When using the random distribution for *&phi;<sub>0</sub>*, the characteristic length must be set to include 6 components: `*Characteristic Length, definition=USER, components=6`. For 3-D kinking, *&phi;<sub>0,12</sub>* = *&phi;<sub>0,13</sub>* is required.
- *&phi;<sub>0</sub>* = 2 :: Identical to *&phi;<sub>0</sub>* = 1, with the exception that a different realization is calculated for each ply. For 3-D kinking, *&phi;<sub>0,12</sub>* = *&phi;<sub>0,13</sub>* is required.
- *&phi;<sub>0</sub>* = 3 :: (Intended for use with 3-D FKT only) A pseudo random distribution varying spatially in the 1-direction is used with a 2-parameter lognormal distribution for the polar angle and a normal distribution for the azimuthal angle. The parameters starting with `fkt_init_misalignment` are used to control the polar and azimuthal distributions. Requires *&phi;<sub>0,12</sub>* = *&phi;<sub>0,13</sub>*.

Pre-existing damage can be modeled by creating an element set for the damaged region and specifying different initial conditions for this element set. For example, to create an intraply matrix crack with no out-of-plane inclination (`CDM_alpha = 0`), the following initial conditions could be specified for the cracked elements:

    *Initial Conditions, type=SOLUTION
     damaged_elset,  1.d0,  0.d0,  0.d0,  0.d0,  0.d0,  0.d0,  0.d0,
              0.d0,  0.d0,     0,     1,  0.d0,  0.d0,  0.d0,  0.d0,
              0.d0,  0.d0,  0.d0,  0.d0

## Model parameters
Model parameters are distinguished from material properties in that there is only one value of each model parameter in a finite element model.
In contrast, a finite element model can have many materials, each with different property values.
A number of model parameters are used in the subroutine that are not directly related to material properties.
These parameters are mostly related to internal error tolerances, iteration limits, etc.
All parameters have default values and, in general, most should not be modified.

The model parameters are defined in the input deck using `*Parameter table`. The default parameters are listed in the file [`tests/_default_parameters.inp`](tests/_default_parameters.inp). When modifying parameters, only datalines associated with non-default parameters need to be include in the input deck.

Previous versions of CompDam (<=2.6.1) used a file named `CompDam.parameters` in the working directory to define non-default parameters. The `CompDam.parameters` is not supported for ComDam versions >2.6.1.


## Example problems
The [`examples/`](examples/) directory includes example models that use CompDam along with corresponding files that defined the expected results. Currently, examples are available for a [double cantilever beam (DCB) problem](examples/interlaminar/readme.md). Additional examples are planned for future releases.


## Test models
This section includes a brief summary of each test implemented in the [`tests/`](tests/) directory. The input deck file names briefly describe the test. All of the input decks start with `test_<elementType>_` and end with a few words describing the test.

<details>

<summary>List of test models</summary>
- *elastic_fiberTension*: Demonstrates the elastic response in the 1-direction under prescribed extension. The 1-direction stress-strain curve has the modulus E1.
- *elastic_matrixTension*: Demonstrates the elastic response in the 2-direction under prescribed extension. The 2-direction stress-strain curve has the modulus E2.
- *elastic_simpeShear12*: Demonstrates the elastic response in the 1-2 plane. The 1-2 plane stress-strain curve has the module G12.
- *elementSize*: Verifies that the characteristic element lengths Lc1, Lc2, and Lc3 are being properly calculated.
- *error*: Verifies that analyses can cleanly terminate upon encountering an error within the user subroutine.
- *failureEnvelope_sig11sig22*: A parametric model in which *&sigma;<sub>11</sub>* is swept from *-X<sub>C</sub>* to *X<sub>T</sub>* and *&sigma;<sub>22</sub>* is swept from *-Y<sub>C</sub>* to *Y<sub>T</sub>* in order to re-create the corresponding failure envelope.
- *failureEnvelope_sig12sig22*: A parametric model in which *&tau;<sub>12</sub>* is swept from *0* to *S<sub>L</sub>* and *&sigma;<sub>22</sub>* is swept from *-Y<sub>C</sub>* to *Y<sub>T</sub>* in order to re-create the corresponding failure envelope.
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

</details>

## Usage and best practices
Getting the most out of CompDam requires careful model setup. Experience applying CompDam to a variety of problems has shown that choices for meshing and solution parameters are critical.
Most of the models described in the [list of publications](list-of-publications.md) utilize one layer of continuum elements per ply with fiber-aligned meshes and zero-thickness cohesive elements at the ply interfaces.
Mass scaling is usually required for modeling quasi-static loading to obtain practical solver run times.
Large mass scaling enables quick run times, which can beneficial during initial model development, however, a user should progressively reduce the mass scaling until the results are sufficiently converged, indicating dynamic effects are negligible.
The models included here in tests/ and examples utilize `mm-N` units, and that choice is recommended. Note that when using `m-N` units, the unit length is large, which will significantly influence aspects of the model (e.g., default cohesive constitutive thickness = 1).

### Element deletion
In general, element deletion is **not** recommended. It is included to provide analyst with options to avoid unwanted numerical issues when damaged elements produce nonphysical behavior. It is an advanced feature, and can introduce undesirable side effects.

Element deletion is enabled by inclduing the option `delete=11` with the `*Depvar` keyword and setting the desired behavior with the parameters as follows. One or more of the following parameters must be set to 1 to enable damage mode(s) to result in element deletion.
```
SET_STATUS_0_ON_D2
SET_STATUS_0_ON_D1T
SET_STATUS_0_ON_D1C
```

Note that element deletion can also occur when the DGD algorithm fails to converge, based on the parameter `terminate_on_no_convergence`. When `terminate_on_no_convergence = 0`, and no converged solution can be found by DGD, a warning is issued, the element is deleted, and the analysis continues.

### Cohesive element option `effmod`
As the opening becomes large in cracked elements utilizing a VUMAT, the stable time increment $\Delta t_{stable}$ reduces. This occurs because the current density at the material point is reducing (constant mass, increasing volume) and the solver assumes the initial stiffness throughout the analysis. This is an unfortunate situation in that fully damaged elements either gain mass (variable mass scaling) or reduce the stable time increment (fixed mass scaling).

This situation can be avoided in COH3D8 elements by including the `effmod` option in the `*User material` definition. With the `effmod` option, CompDam returns effective stiffness estimates to the solver that are used for computing the dilational wave speed when calculating the element's stable time increment. Since a fully damaged cohesive element has no stiffness, it can have a large stable time increment. In other words, the stiffness reduction occurs much faster than the volume change for partical input properties. Therefore, with the `effmod` option, the stable time increment of damaged cohesive elements will not govern stable time increment $\Delta t_{stable}$ for the model.

See [test_COH3D8_normal_effmod.inp](tests/test_COH3D8_normal_effmod.inp) or [test_DCB_3D_compdam.pes](examples/interlaminar/pes/test_DCB_3D_compdam.pes) for examples using `effmod`.

In cohesive elements, the damage variable completely controls the stiffness of the element, whereas in continuum elements, the effect of damage on stiffness is anisotropic. In continuum elements, hypothetically a similar calculation could be implemented so that fully damaged continuum elements do not suffer from the issue described above. However, in practice, analyses rarely proceed to the point where continuum elements are fully damaged (matrix and fiber direction) and experiencing large opening prior to catastrophic failure. Therefore, `effmod` has not been implemented for continuum elements.

### Dynamic loading and cohesive elements
For problems with dynamic loading, where it is critical to model accurately the dynamic structural response, some additional considerations are required.
The use of cohesive elements introduces non-physical mass, since cohesive elements mass defined by their density with volume calculated from the constitutive thickess (not the modeled thickness).
Therefore, zero-thickness cohesive interface layers can add substantial mass to laminate, which influences the structural response.
This issue can be mitigated with two different strategies: finite-thickness cohesive elements or prescribed constitutive thickness.

#### Finite thickness cohesive elements for dynamic problems
As first described in [Leone et al. 2023](https://ntrs.nasa.gov/citations/20230008595), finite cohesive elements can be used where the cohesive element thickness is set to a small fraction of the ply thickness (e.g., 10%) and the ply thickness is reduced accordingly so that the laminate thickness is maintained accurate.
Since the volume of the finite thickness element is calculated from its nodal coordinates, the mass is much reduced compared with zero-thickness cohesive elements.
The ply elastic stiffnesses are increased by nominal_ply_thickness/modeled_cohesive_thickness so that the laminate's elastic response is unaffected.
Finally, the density for the cohesive elements is set = material density * cohesive element thickness.
See [test_COH3D8_thick_normal.inp](tests/test_COH3D8_thick_normal.inp) or [test_DCB_3D_compdam_finite_thk.pes](examples/interlaminar/pes/test_DCB_3D_compdam_finite_thk.pes) for examples using finite thickness cohesive elements.

#### Constitutive thickness for dynamic problems
Since cohesive element density $\rho$ and constitutive thickness $t_c$ are numerical inputs with limited physical meaning, and they mainly impact the time incrementation in the explicit analysis, the user may select values that are beneficial for time incrementation.

<details>

<summary>Click to expand</summary>

The following expressions for $\rho$ and $t_c$ can be derived
$$
\rho = \frac{1}{E'} \left(\frac{m}{A_c\Delta t} \right)^2
$$
$$
t_{c} = \frac{A_c\Delta t^2 E'}{m}
$$
which allows the user to choose desired mass $m$ and time increment $\Delta t$ values. $A_c$ is the cross-sectional area and $E'$ is the effective modulus 
$$
E'/t = \frac{(K_1+K_3)}{2} + \frac{3K_2}{4}
$$
used to estimate the dilatational wave speed and $t$ is a unit thickness. $K_1, K_3, K_2$ are the penalty stiffnesses for shear in L and T and the normal direction. $K_2$ is the input material property #21 and
$$
K_1 = K_2G_{Ic}S_L^2/(G_{IIc}Y_T^2)
$$
$$
K_3 = K_2G_{Ic}S_T^2/(G_{IIc}Y_T^2)
$$
$$
S_T = Y_C*cos(\alpha_0)*(sin(\alpha_0) + cos(\alpha_0)/tan(2\alpha_0))
$$
The expression for $E'$ is only valid when `effmod` is used.
The values for $m$ and $\Delta t$ can be selected based on the values of adjacent continuum elements.
See [test_COH3D8_normal_effmod.inp](tests/test_COH3D8_normal_effmod.inp) for an example of these calculations.

In cases with high penalty stiffness (note this can be caused by low values of YT), numerical instability can occur when $t_{c} \ne 1$, and so a warning is issued by CompDam whenever $t_{c} \ne 1$ is specified.

</details>

## Resources for developers
- [Advanced debugging](advanced-debugging.md)
- [Python extension module](python-extension-module.md)


## Contributing
We invite your contributions to CompDam_DGD! Please submit contributions (including a test case) with pull requests so that we can reproduce the behavior of interest. Commits history should be clean. Please contact the developers if you would like to make a major contribution to this repository. Here is a [checklist](contributing-checklist.md) that we use for contributions.

## Citing CompDam
If you use CompDam, please cite using the following BibTex entry:

    @misc{CompDam,
    title={Comp{D}am - {D}eformation {G}radient {D}ecomposition ({DGD}), v2.7.0},
    author={Andrew C. Bergan and Frank A. Leone and Carlos G. D\'{a}vila},
    note={https://github.com/nasa/CompDam\_DGD},
    year={2026}
    }
