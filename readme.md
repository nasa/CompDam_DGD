# CompDam - Deformation Gradient Decomposition (DGD)
This code is a continuum damage mechanics (CDM) material model for use with the Abaqus finite element code. This is a research code which aims to provide an accurate representation of mesoscale damage modes in fiber-reinforced polymer materials.

The CDM material model is implemented as an Abaqus VUMAT for the simulation of matrix cracks formed under tensile, compressive, and shear loading conditions and fiber fracture under tensile and compressive loading. A smeared crack approach is used to accurately represent the kinematics of matrix cracks rotating with the deformed material.

This software may be used, reproduced, and provided to others only as permitted under the terms of the agreement under which it was acquired from the U.S. Government. Neither title to, nor ownership of, the software is hereby transferred. This notice shall remain on all copies of the software.

Copyright 2016 United States Government as represented by the Administrator of the National Aeronautics and Space Administration. No copyright is claimed in the United States under Title 17, U.S. Code. All Other Rights Reserved.

Publications that use this code include:
- Cheryl A. Rose, et al. ["Analysis Methods for Progressive Damage of Composite Structures"](http://hdl.handle.net/2060/20140001002) NASA/TM-2013-218024, July 2013.
- Frank A. Leone Jr. ["Deformation gradient tensor decomposition for representing matrix cracks in fiber-reinforced materials"](http://dx.doi.org/10.1016/j.compositesa.2015.06.014) *Composites Part A* (2015) **76**:334-341.
- Frank A. Leone Jr. ["Representing matrix cracks through decomposition of the deformation gradient tensor in continuum damage mechanics methods"](http://iccm20.org/fullpapers/file?f=Abk7n4gkWV) *ICCM20*.

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
- [Advanced debugging](#advanced-debugging)
- [Summary of tests classes](#summary-of-tests-classes)
- [Known issues](#known-issues)


## Getting started

### Source code
The user subroutine source code is located in the `for` directory. The main entry point is `CompDam_DGD.for`.

### Prerequisites
[Intel Fortran Compiler](https://software.intel.com/en-us/fortran-compilers) version 11.1 or newer is required to compile the code. It is recommended that Abaqus 6.14-1 or newer is used with this code.

### Abaqus environment file settings
The `abaqus_v6.env` file must have [`/fpp`](https://software.intel.com/en-us/node/579498), [`/Qmkl:sequential`](https://software.intel.com/en-us/node/579338), and [`/free`](https://software.intel.com/en-us/node/579524) in the `ifort` command where the format for Windows is used. The corresponding Linux format is: `-fpp`, `-free`, and `-mkl`. The `/fpp` option enables the Fortran preprocessor, which is required for the code to compile correctly. The `/free` option sets the compiler to free-formatting for the source code files. The `/Qmkl:sequential` enables the [Intel Math Kernel Library (MKL)](https://software.intel.com/en-us/intel-mkl), which provides optimized and verified functions for many mathematical operations. The MKL is used in this code for calculating eigenvalues and eigenvectors.

A sample environment file using the Abaqus 2016 environment file format is provided in the `tests` directory for Windows systems. For older versions of Abaqus, copy your system environment file (in `SIMULIA\Abaqus\6.14-1\SMA\site`) to your working directory and make sure that the `ifort` arguments are set as specified above.

### Submitting a job
This code is an Abaqus/Explicit VUMAT. Please refer to the Abaqus documentation for the general instructions on how to submit finite element analyses using user subroutines. Please see the [example input file statements](#example-input-file-statements) for details on how to interface with this particular VUMAT subroutine.

Analyses with this code **must** be run in double precision. Some of the code has double precision statements and variables hard-coded, so if Abaqus/Explicit is run in single precision, compile-time errors will arise. When submitting an Abaqus/Explicit job from the command line, double precision can be specified by including the command line argument `double=both`.

For example, run the test model `test_C3D8R_elastic_fiberTension` in the `tests` directory with the following command:
```
tests $   abaqus job=test_C3D8R_elastic_fiberTension user=../for/CompDam_DGD.for double=both
```
On Linux systems, `CompDam_DGD.for` must be replaced with `CompDam_DGD.f` and the sample environment file needs to be changed to the linux format.

### Example input file statements
Example 1, using an [external material properties file](#defining-the-material-properties-in-a-props-file):

<pre>
*Section controls, name=<name>, distortion control=YES, element deletion=YES
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
  6, CDM_Lc
  7, CDM_rfT
  8, CDM_d1
  9, CDM_FImT
 10, CDM_alpha
 11, CDM_STATUS
 12, CDM_Plas12
 13, CDM_Inel12
 14, CDM_d_comp_init
 15, CDM_slide1
 16, CDM_slide2
 17, CDM_rfC
 18, CDM_d1T
 19, CDM_d1C
*User material, constants=3
** 1              2  3
** feature flags,  , thickness
          100001,  ,       0.1
**
*Initial Conditions, type=SOLUTION
 elset_name,  0.d0,  0.d0,  0.d0,  0.d0,  0.d0,  0.d0,  0.d0,
              0.d0,  0.d0,  -999,     1,  0.d0,  0.d0,  0.d0,
              0.d0,  0.d0,  0.d0,  0.d0,  0.d0
** - alpha can be pre-set by setting SV10 to the desired angle.
**       Setting alpha = -999 will evaluate cracks every 10 degrees in
**       the 2-3 plane to find the correct crack initiation angle.
** In each step, NLGEOM=YES must be used. This is the default setting.
</pre>

Example 2, using an [input deck command](#defining-the-material-properties-in-the-input-deck):

<pre>
*Section controls, name=<name>, distortion control=YES, element deletion=YES
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
  6, CDM_Lc
  7, CDM_rfT
  8, CDM_d1
  9, CDM_FImT
 10, CDM_alpha
 11, CDM_STATUS
 12, CDM_Plas12
 13, CDM_Inel12
 14, CDM_d_comp_init
 15, CDM_slide1
 16, CDM_slide2
 17, CDM_rfC
 18, CDM_d1T
 19, CDM_d1C
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
**  XC,       fXC,      GXC,      fGXC,     rsvd,     rsvd,     rsvd,     mu
    1200.1,      ,         ,          ,         ,         ,         ,     0.3
**
*Initial Conditions, type=SOLUTION
 elset_name,  0.d0,  0.d0,  0.d0,  0.d0,  0.d0,  0.d0,  0.d0,
              0.d0,  0.d0,  -999,     1,  0.d0,  0.d0,  0.d0,
              0.d0,  0.d0,  0.d0,  0.d0,  0.d0
** - alpha can be pre-set by setting SV10 to the desired angle.
**       Setting alpha = -999 will evaluate cracks every 10 degrees in
**       the 2-3 plane to find the correct crack initiation angle.
** In each step, NLGEOM=YES must be used. This is the default setting.
</pre>

### Running tests
Test cases are available in the `tests` directory. The tests are useful for demonstrating the capabilities of the VUMAT as well as to verify that the code performs as intended. Try running some of the test cases to see how the code works. The test cases can be submitted as a typical Abaqus job using the Abaqus command line arguments.

### Building a shared library
CompDam_DGD can be built into a shared library file. Follow these steps:
1. Place a copy of the Abaqus environment file (with the compiler flags specified) in the `for` directory
2. From the `for` directory, on Windows, run:
```
for $  abaqus make library=CompDam_DGD.for
```
On Linux systems, `CompDam_DGD.for` must be replaced with `CompDam_DGD.f`.
This command will create executable files for the operating system it is executed on (`.dll` for Windows and `.so` for Linux).

When using a pre-compiled shared library, it is only necessary to specify the location of the shared library files in the environment file (the compiler options are not required). To run an analysis using a shared library, add `usub_lib_dir = <full path to shared library file>` to the Abaqus environment file in the Abaqus working directory.

## Model features
The CompDam_DGD material model implements a variety of features that can be enabled or disabled by the user. An overview of these features is provided in this section. The material properties required for each feature are listed. References are provided to more detailed discussions of the theoretical framework for each feature.

### Fully orthotropic elasticity
The composite materials modeled with CompDam_DGD can be defined assuming either [transverse isotropy](http://www.efunda.com/formulae/solid_mechanics/mat_mechanics/hooke_iso_transverse.cfm) or [orthotropy](http://www.efunda.com/formulae/solid_mechanics/mat_mechanics/hooke_orthotropic.cfm). For a transversely isotropic material definition, the following properties must be defined: E1, E2, G12, v12, and v23. For an orthotropic material definition, the following additional properties must be defined: E2, G13, G23, and v13.

### Matrix damage
Tensile and compressive matrix damage is modeled by embedding cohesive laws to represent cracks in the material according to the deformation gradient decomposition method of [Leone (2015)](http://doi.org/10.1016/j.compositesa.2015.06.014). The mixed-mode behavior of matrix cracks is defined according to the Benzeggagh-Kenane law. The initiation of compressive matrix cracks accounts for friction on the potential crack surface according to the LaRC04 failure criterion.

The following material properties are required for the prediction of matrix damage: YT, SL, GYT, GSL, eta_BK, YC, and alpha0.

### Thermal strains
The thermal strains are calculated by multiplying the 1-, 2-, and 3-direction coefficients of thermal expansion by the current &Delta;T, as provided by the Abaqus solver. The thermal strains are subtracted from the current total strain.

The required material properties are the coefficients of thermal expansion in the 1 and 2 directions. It is assumed that the 2- and 3-direction coefficients of thermal expansion are equal.

Hygroscopic strains are not accounted for. If the effects of hygroscopic expansion are to be modeled, it is recommended to smear the hygroscopic and thermal expansion coefficients to approximate the response using the solver-provided &Delta;T.

### Shear nonlinearity
Shear nonlinearity is modeled in the 1-2 plane according to the [Ramberg-Osgood equation](https://en.wikipedia.org/wiki/Ramberg%E2%80%93Osgood_relationship), with its parameters selected to fit  experimental data. As applied herein, the Ramberg-Ogsood equation is written in the following form:

*&gamma;*<sub>12</sub> = [*&tau;*<sub>12</sub> + *&alpha;*<sub>PL</sub>sign(*&tau;*<sub>12</sub>)|*&tau;*<sub>12</sub>|<sup>*n*<sub>PL</sub></sup>]/*G*<sub>12</sub>

where *&gamma;*<sub>12</sub> is the shear strain and *&tau;*<sub>12</sub> is the shear stress. Prior to the initiation of matrix damage (i.e., `d2`), the nonlinear response due to the above is equation is plastic, and the unloading/reloading slope is unchanged. The required material inputs are the two parameters in the above equation: *&alpha;*<sub>PL</sub> and *n*<sub>PL</sub>.

### Fiber tensile damage
A continuum damage mechanics model similar to the work of [Maimí et al. (2007)](http://doi.org/10.1016/j.mechmat.2007.03.005) is used. The model utilizes a non-interacting maximum strain failure criterion, and bilinear softening after the initiation of failure. The area under the stress-strain curve is equal to the fracture toughness divided by the element length. The required material properties are: XT, fXT, GXT, and fGXT, where fXT and fGXT are ratios of strength and fracture toughness for bilinear softening, analogous to the *n* and *m* terms in [Dávila et al. (2009)](http://doi.org/10.1007/s10704-009-9366-z).

### Fiber compression damage
Same model as in tension, but for compression. Assumes maximum strain failure criterion and bilinear softening. The required material properties are: XC, fXC, GXC, fGXC.

Load reversal assumptions from [Maimí et al. (2007)](http://doi.org/10.1016/j.mechmat.2007.03.006).

### Friction
Friction is modeled on the damaged fraction of the cross-sectional area of DGD cracks using the approach of [Alfano and Sacco (2006)](http://doi.org/10.1002/nme.1728). The coefficient of friction *&mu;* must be defined to account for friction on the failed crack surface.

The amount of sliding which has taken place in the longitudinal and transverse directions are stored in state variables `CDM_slide1` and `CDM_slide2`.

### Strain definition
The strain is calculated using the deformation gradient tensor provided by the Abaqus solver. The default strain definition used is the Green-Lagrange strain:

*E* = (*F*<sup>T</sup>*F* - *I*)/2

Hooke's law is applied using the Green-Lagrange strain to calculate the 2<sup>nd</sup> Piola-Kirchhoff stress *S*.


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
| #  |         Symbol         |   Name   |                  Description                  |                       Units                    |            Admissible values             | Recommeded Test |
|----|------------------------|----------|-----------------------------------------------|------------------------------------------------|------------------------------------------|-----------------|
|  9 | *E<sub>1</sub>*        | E1       | Longitudinal Young's modulus                  | F/L<sup>2</sup>                                | 0 < *E<sub>1</sub>* < &infin;            | ASTM D3039      |
| 10 | *E<sub>2</sub>*        | E2       | Transverse Young's modulus                    | F/L<sup>2</sup>                                | 0 < *E<sub>2</sub>* < &infin;            | ASTM D3039      |
| 11 | *G<sub>12</sub>*       | G12      | In-plane Shear modulus                        | F/L<sup>2</sup>                                | 0 < *G<sub>12</sub>* < &infin;           | ASTM D3039      |
| 12 | *&nu;<sub>12</sub>*    | nu12     | Major Poisson's ratio                         | -                                              | 0 &le; *&nu;<sub>12</sub>* &le; 1        | ASTM D3039      |
| 13 | *&nu;<sub>23</sub>*    | nu23     | Minor Poisson's ratio                         | -                                              | 0 &le; *&nu;<sub>23</sub>* &le; 1        |                 |
|    | ===                    |          |                                               |                                                |                                          |                 |
| 14 | *Y<sub>T</sub>*        | YT       | Transverse tensile strength, in-situ          | F/L<sup>2</sup>                                | 0 < *Y<sub>T</sub>* < &infin;            | ASTM D3039      |
| 15 | *S<sub>L</sub>*        | SL       | Shear strength, in-situ                       | F/L<sup>2</sup>                                | 0 < *S<sub>L</sub>* < &infin;            |                 |
| 16 | *G<sub>Ic</sub>*       | GYT      | Mode I fracture toughness                     | F/L                                            | 0 < *G<sub>Ic</sub>* < &infin;           | ASTM D5528      |
| 17 | *G<sub>IIc</sub>*      | GSL      | Mode II fracture toughness                    | F/L                                            | 0 < *G<sub>IIc</sub>* < &infin;          | ASTM D7905      |
| 18 | *&eta;*                | eta_BK   | BK exponent for mode-mixity                   | -                                              | 1 &le; *&eta;* < &infin;                 |                 |
| 19 | *Y<sub>C</sub>*        | YC       | Transverse compressive strength               | F/L<sup>2</sup>                                | 0 < *Y<sub>C</sub>* < &infin;            | ASTM D3410      |
| 20 | *&alpha;<sub>0</sub>*  | alpha0   | Fracture plane angle for pure trans. comp.    | Radians                                        | 0 &le; *&alpha;<sub>0</sub>* &le; &pi;/2 |                 |
|    | ------                 |          |                                               |                                                |                                          |                 |
| 21 | *E<sub>3</sub>*        | E3       | 3-Direction Young's modulus                   | F/L<sup>2</sup>                                | 0 < *E<sub>3</sub>* < &infin;            |                 |
| 22 | *G<sub>13</sub>*       | G13      | Shear modulus in 1-3 plane                    | F/L<sup>2</sup>                                | 0 < *G<sub>13</sub>* < &infin;           |                 |
| 23 | *G<sub>23</sub>*       | G23      | Shear modulus in 1-2 plane                    | F/L<sup>2</sup>                                | 0 < *G<sub>23</sub>* < &infin;           |                 |
| 24 | *&nu;<sub>13</sub>*    | nu13     | Poisson's ratio in 2-3 plane                  | -                                              | 0 &le; *&nu;<sub>13</sub>* &le; 1        |                 |
|    | ------                 |          |                                               |                                                |                                          |                 |
| 25 | *&alpha;<sub>11</sub>* | alpha11  | Coefficient of long. thermal expansion        | /&deg;                                         | -1 &le; *&alpha;<sub>11</sub>* &le; 1    |                 |
| 26 | *&alpha;<sub>22</sub>* | alpha11  | Coefficient of tran. thermal expansion        | /&deg;                                         | -1 &le; *&alpha;<sub>22</sub>* &le; 1    |                 |
|    | ------                 |          |                                               |                                                |                                          |                 |
| 27 | *&alpha;<sub>PL</sub>* | alpha_PL | Nonlinear shear parameter                     | (F/L<sup>2</sup>)<sup>1-*n*<sub>PL</sub></sup> | 0 &le; *&alpha;<sub>PL</sub>* < &infin;  |                 |
| 28 | *n<sub>PL</sub>*       | n_PL     | Nonlinear shear parameter                     | -                                              | 0 &le; *n<sub>PL</sub>* < &infin;        |                 |
|    | ------                 |          |                                               |                                                |                                          |                 |
| 29 | *X<sub>T</sub>*        | XT       | Long. tensile strength                        | F/L<sup>2</sup>                                | 0 < *X<sub>T</sub>* < &infin;            | ASTM D3039      |
| 30 | *f<sub>XT</sub>*       | fXT      | Long. tensile strength ratio                  | -                                              | 0 &le; *f<sub>XT</sub>* &le; 1           |                 |
| 31 | *G<sub>XT</sub>*       | GXT      | Long. tensile fracture toughness              | F/L                                            | 0 < *G<sub>XT</sub>* < &infin;           |                 |
| 32 | *f<sub>GXT</sub>*      | fGXT     | Long. tensile fracture toughness ratio        | -                                              | 0 &le; *f<sub>GXT</sub>* &le; 1          |                 |
|    | ------                 |          |                                               |                                                |                                          |                 |
| 33 | *X<sub>C</sub>*        | XC       | Longitudinal compressive strength             | F/L<sup>2</sup>                                | 0 < *X<sub>C</sub>* < &infin;            | ASTM D3410      |
| 34 | *f<sub>XC</sub>*       | fXC      | Long. compression strength ratio              | -                                              | 0 &le; *f<sub>XC</sub>* &le; 1           |                 |
| 35 | *G<sub>XC</sub>*       | GXC      | Long. compression fracture toughness          | F/L                                            | 0 < *G<sub>XC</sub>* < &infin;           |                 |
| 36 | *f<sub>GXC</sub>*      | fGXC     | Long. compression fracture toughness ratio    | -                                              | 0 &le; *f<sub>GXC</sub>* &le; 1          |                 |
|    | ------                 |          |                                               |                                                |                                          |                 |
| 40 | *&mu;*                 | mu       | Coefficient of friction                       | -                                              | 0 &le; *&mu;* &le; 1                     |                 ||

Notes:
- The first five inputs (above the ===) are required
- Properties for each feature are grouped and separated by ------
- The number in the first column corresponds to the property order when defined in the input deck
- Properties not listed in the table above (see next section):
  1. [feature flags](#controling-which-features-are-enabled)
  2. *reserved*
  3. [thickness](#definition-of-thickness)
  4. *reserved*
  5. *reserved*
  6. *reserved*
  7. *reserved*
  8. *reserved*
- &infin; is calculated with the Fortran intrinsic `Huge` for double precision

### Required inputs for the `*Material` data lines in the input deck
The feature flags and thickness are defined in the input deck on the material property data lines. These properties must be defined in the input deck whether the other material properties are defined via the .props file or via the input deck. While feature flags and thickness are not material properties per se, they are used in controlling the behavior of the material model.

#### Controlling which features are enabled
Model features can be enabled or disabled by two methods. The first method is specifying only the material properties required for the features you would like to enable. CompDam_DGD disables any feature for which all of the required material properties have not been assigned. If an incomplete set of material properties are defined for a feature, a warning is issued.

The second method is by specifying the status of each feature directly as a material property in the input deck. Each feature of the subroutine is controlled by a position in an integer, where 0 is disabled and 1 is enabled.
The positions correspond to the features as follows:
- Position 1: Matrix damage
- Position 2: Shear nonlinearity
- Position 3: Fiber tensile damage
- Position 4: Fiber compression damage
- Position 5: *reserved*
- Position 6: Friction

For example, `101000` indicates that the model will run with matrix damage and fiber tension damage enabled; `110001` indicates that the model will run with matrix damage, in-plane shear nonlinearity, and friction.

#### Definition of thickness
Length along the thickness-direction associated with the current integration point

## State variables
The table below lists all of the state variables in the model. The model requires a minimum of 18 state variables. Additional state variables are defined depending on which (if any) fiber compression features are enabled. When fiber compression (feature flag position 4) is used, 19 state variables are required.

| # | Name             | Description                                                           |
|---|------------------|-----------------------------------------------------------------------|
|  1| `CDM_d2`         | d2, Matrix cohesive damage variable                                   |
|  2| `CDM_Fb1`        | Cohesive Opening Displacement 1 - fiber-direction shear               |
|  3| `CDM_Fb2`        | Cohesive Opening Displacement 2 - normal                              |
|  4| `CDM_Fb3`        | Cohesive Opening Displacement 3 - transverse shear                    |
|  5| `CDM_B`          | Mode Mixity (*G*<sub>II</sub> / (*G*<sub>I</sub> + *G*<sub>II</sub>)) |
|  6| `CDM_Lc`         | Characteristic Element Length                                         |
|  7| `CDM_rfT`        | Fiber tensile threshold (maximum strain)                              |
|  8| `CDM_d1`         | d1, Fiber damage variable                                             |
|  9| `CDM_FImT`       | Matrix failure criterion (del/del_0)                                  |
| 10| `CDM_alpha`      | alpha, the cohesive surface normal [degrees, integer]                 |
| 11| `CDM_STATUS`     | STATUS (for element deletion)                                         |
| 12| `CDM_Plas12`     | 1-2 plastic strain                                                    |
| 13| `CDM_Inel12`     | 1-2 inelastic strain                                                  |
| 14| `CDM_d_comp_init`| Normal cohesive displacement at initiation                            |
| 15| `CDM_slide1`     | Cohesive sliding displacement along 1-direction                       |
| 16| `CDM_slide2`     | Cohesive sliding displacement along 2-direction                       |
| 17| `CDM_rfC`        | Fiber compression failure threshold                                   |
| 18| `CDM_d1T`        | Fiber tension damage variable                                         |
|---|------------------|-----------------------------------------------------------------------|
| 19| `CDM_d1C`        | Fiber compression damage variable                                     ||

### Initial conditions
All state variables should be initialized using the `*Initial conditions` command. All state variables should be initialized as zero, except `CDM_alpha` and `CDM_STATUS`.

The initial condition for `CDM_alpha` can be used to specify a predefined angle for the cohesive surface normal. To specify a predefined `CDM_alpha`, set the initial condition for `CDM_alpha` to an integer (degrees). The range of valid values for `CDM_alpha` depends on the aspect ratio of the element, but values in the range of 0 to 90 degrees are always valid. Setting `CDM_alpha` to -999 will make the subroutine evaluate cracks every 10 degrees in the 2-3 plane to find the correct crack initiation angle. Note that `CDM_alpha` is measured from the 2-axis rotating about the 1-direction.

Since `CDM_STATUS` is used for element deletion, always initialize `CDM_STATUS` to 1.

## Advanced debugging
Using an interactive debugger helps to identify issues in the Fortran code. Abaqus knowledge base article QA00000007986 describes the details involved. The following is a quick-start guide for direct application to CompDam.

Several statements for debugging need to be uncommented in the environment file. Follow these steps:
1. Copy your system environment file (`SIMULIA\Abaqus\6.14-1\SMA\site`) to your local working directory. For the example below, copy the environment file to the `tests` directory.
2. Edit the local environment file: uncomment lines that end with `# <-- Debugging`, `# <-- Debug symbols`, and `# <-- Optimization Debugging`

Run the job with the `-debug` and `-explicit` arguments. For example:
```
tests $ abaqus -j test_C3D8R_fiberTension -user ../for/CompDam_DGD.for -double both -debug -explicit
```

This command should open the [Visual Studio debugging software](https://msdn.microsoft.com/en-us/library/sc65sadd.aspx) automatically. Open the source file(s) to debug. At a minimum, open the file with the subroutine entry point `for/CompDam_DGD.for`. Set a break point by clicking in the shaded column on the left edge of the viewer. The break point will halt execution. Press <kbd>F5</kbd> to start the solver. When the break point is reached, a yellow arrow will appear and code execution will pause. Press <kbd>F5</kbd> to continue to the next break point, press <kbd>F11</kbd> to execute the next line of code following execution into function calls (Step Into), or press <kbd>F10</kbd> to execute the next line of code but not follow execution into function calls (Step Over).

To stop execution, close the visual studio window. Choose stop debugging and do not save your changes.

[More tips on debugging Fortran programs from Intel](https://software.intel.com/en-us/articles/tips-for-debugging-run-time-failures-in-intel-fortran-applications).

## Summary of tests classes
This section includes a brief summary of each test implemented in the `tests` folder. The input deck file names briefly describe the test. All of the input decks start with `test_<elementType>_` and end with a few words describing the test. A more detailed description for each is given below:
- *elastic_fiberTension*: Demonstrates the elastic response in the 1-direction under prescribed extension. The 1-direction stress-strain curve has the modulus E1.
- *elastic_matrixTension*: Demonstrates the elastic response in the 2-direction under prescribed extension. The 2-direction stress-strain curve has the modulus E2.
- *fiberCompression_CDM*: Demonstrates the constitutive response in the 1-direction under prescribed shortening. The 1-direction stress-strain curve is trilinear.
- *fiberLoadReversal*: Demonstrates the constitutive response in the 1-direction under prescribed extension and shortening reversals. The 1-direction stress-strain curve shows the intended behavior under load reversal.
- *fiberTension*: Demonstrates the constitutive response in the 1-direction under prescribed extension. The 1-direction stress-strain curve is trilinear.
- *matrixTension*: Demonstrates the constitutive response in the 2-direction under prescribed extension. The 2-direction stress-strain curve is bilinear.
- *nonlinearShear12*: Demonstrates the constitutive response under prescribed simple shear. Shows the response under load reversal.
- *simpleShear12*: Demonstrates the constitutive response under prescribed simple shear. The shear stress-strain curve is bilinear.
- *simpleShear12friction*: Demonstrates the constitutive response under prescribed simple shear with friction enabled. The element is loaded under transverse compression and then sheared. Shows the friction-induced stresses.

## Contributing
We invite your contributions to CompDam_DGD! Please submit contributions (including a test case) with pull requests so that we can reproduce the behavior of interest. Commits history should be clean. Please contact the developers if you would like to make a major contribution to this repository.
