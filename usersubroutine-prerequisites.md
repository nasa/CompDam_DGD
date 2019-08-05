# Abaqus Prerequisites for User Subroutines

The table below summarizes the recommended compilers for use with different versions of Abaqus. See Abaqus documentation for more details.

## Windows 7, 8, 10
| Abaqus version | C++ Compiler                | Fortran Compiler                         | MPI*            |
|----------------|-----------------------------|------------------------------------------|-----------------|
| 6.14           | Visual Studio 2010 SP1      | Intel Fortran Composer XE 2011 Update 6  | MS-MPI          |
| 2016           | Visual Studio 2012 Update 4 | Intel Fortran Composer XE 2011 Update 6  | MS-MPI          |
| 2017           | Visual Studio 2012 Update 5 | Intel Visual Fortran 16.0 update 1       | MS-MPI 4.2, 5.0 |
| 2018           | Visual Studio 2015 Update 3 | Intel Visual Fortran 16.0 update 1       | MS-MPI 5.0, 8.0 |
| 2019           | Visual Studio 2015 Update 3 | Intel Visual Fortran 16.0 update 1       | MS-MPI 9.0.1    |

* Windos only seems to allow 1 MS-MPI installation. MS-MPI 9.0.1 appears to work for Abaqus 2017 - 2019.


## Linux
| Abaqus version | C++ Compiler                | Fortran Compiler                         |
|----------------|-----------------------------|------------------------------------------|
| 6.14           | GCC 4.4.7                   | Intel Fortran Composer XE 2011 Update 6  |
| 2016           | GCC                         | Intel Fortran Composer XE 2011 Update 6  |
| 2017           | GCC                         | Intel Visual Fortran 16.0 update 1       |
| 2018           | GCC                         | Intel Visual Fortran 16.0 update 1       |
| 2019           | GCC                         | Intel Visual Fortran 16.0 update 1       |