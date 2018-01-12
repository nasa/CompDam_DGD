! This UMAT acts as a wrapper for VUMAT user material subroutines. The CompDam_DGD material subroutine was originally
! developed with the Abaqus/Explicit solver in mind, but this wrapper allows it to be used with Abaqus/Standard with
! most of its original functionality.


! Re-map the Abaqus/Explicit utility subroutines to their Abaqus/Standard equivalents via Fortran preprocessor
#define VGETOUTDIR GETOUTDIR
#define VGETJOBNAME GETJOBNAME
#define XPLB_ABQERR STDB_ABQERR
#define XPLB_EXIT XIT

! Include, via preprocessor, the source file containing the VUMAT subroutine to be wrapped. Using the preprocessor
!  include statement is necessary so that the above #define commands propagate through all relevant source files.
#include "CompDam_DGD.for"
#include "vumatWrapper.for"

Subroutine UMAT( &
  stress,statev,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT,DRPLDE,DRPLDT, &
  STRAN,DSTRAN,TIME,dtime,TEMP,DTEMP,PREDEF,DPRED,CMNAME,    &
  ndi,nshr,NTENS,nstatv,props,nprops,COORDS,DROT,PNEWDT,     &
  CELENT,DFGRD0,DFGRD1,noel,npt,layer,kspt,JSTEP,KINC)

  Use stress_Mod
  Use vumat_Wrapper_Mod

  Implicit Double Precision (a-h, o-z)
  Integer, parameter :: nprecd=2, km=1, lanneal=0, ndir=3

  ! UMAT arguments
  Dimension stress(NTENS),statev(nstatv),DDSDDE(NTENS,NTENS),  &
    DDSDDT(NTENS),DRPLDE(NTENS),STRAN(NTENS),DSTRAN(NTENS),    &
    TIME(2),PREDEF(1),DPRED(1),props(nprops),COORDS(3),        &
    DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),JSTEP(4)
  Character(len=80) :: cmname

  ! Local variables
  Double Precision :: Lc(3)
  Double Precision :: damage_old, damage_new, max_delta_damage, cutback_for_damage

  ! Parameters
  Double Precision, parameter :: zero=0.d0, one=1.d0

  ! Check and ensure that geometric nonlinearity is active for the current step
  If (JSTEP(3) /= 1) Then
    Print *, 'ERROR: CompDam_DGD requires that NLGEOM=YES for all steps.'
    Call XIT
  End If

  ! Define non-standard element lengths (a la vucharlength() calculations)
  Lc(1) = SQRT(CELENT**3/props(3))
  Lc(2) = Lc(1)
  Lc(3) = props(3)

  ! Record the maximum previous damage state
  damage_old = MAX(statev(1), statev(18), statev(19))

  Call vumat_Wrapper( &
    ! Standard UMAT arguments
    stress,statev,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT,DRPLDE,DRPLDT, &
    STRAN,DSTRAN,TIME,dtime,TEMP,DTEMP,PREDEF,DPRED,CMNAME,    &
    ndi,nshr,NTENS,nstatv,props,nprops,COORDS,DROT,PNEWDT,     &
    CELENT,DFGRD0,DFGRD1,noel,npt,layer,kspt,JSTEP,KINC,       &
    ! Optional arguments for added vumat_Wrapper functionality
    Lc,size(Lc))

  ! Define a Jacobian
  DDSDDE = StiffFunc(NTENS, 171420.d0, 9080.d0, 9080.d0, 5290.d0, 5290.d0, 2897.d0, 0.32d0, 0.32d0, 0.45d0, zero, zero, zero)

  ! Record the maximum current damage state
  damage_new = MAX(statev(1), statev(18), statev(19))

  max_delta_damage = 0.2d0
  cutback_for_damage = 0.1d0

  If (damage_old == zero) Then
    ! Cut back the time increment until the initial damage calculation is sufficiently small
    If (damage_new > max_delta_damage) Then
      PNEWDT = cutback_for_damage
    ! Else If (damage_new > zero) Then
    !   Call XIT
    End If
  End If

  ! -------------------------------------------------------------------- !
  Return
End Subroutine UMAT
