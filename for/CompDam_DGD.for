! This software may be used, reproduced, and provided to others only as permitted under the terms of the agreement
! under which it was acquired from the U.S. Government. Neither title to, nor ownership of, the software is hereby
! transferred. This notice shall remain on all copies of the software.
!
! Copyright 2016 United States Government as represented by the Administrator of the National Aeronautics and
! Space Administration. No copyright is claimed in the United States under Title 17, U.S. Code. All Other Rights
! Reserved.

#include "vumatArgs.for"
#include "version.for"
#include "forlog.for"
#include "matrixUtil.for"
#include "matProp.for"
#include "stateVar.for"
#include "parameters.for"
#include "schapery.for"
#include "stress.for"
#include "strain.for"
#include "schaefer.for"
#include "plasticity.for"
#include "fiberDamage.for"
#include "friction.for"
#include "cohesive.for"
#include "vucharlength.for"
#include "DGD.for"
#include "vexternaldb.for"


Subroutine VUMAT(  &
  ! Read only (unmodifiable) variables:
  jblock,ndir,nshr,nstatev,nfieldv,nprops,lanneal,stepTime,     &
  totalTime,dt,cmname,coordMp,charLength,props,density,         &
  strainInc,relSpinInc,tempOld,stretchOld,defgradOld,fieldOld,  &
  stressOld,stateOld,enerInternOld,enerInelasOld,tempNew,       &
  stretchNew,defgradNew,fieldNew,                               &
  ! Write only (modifiable) variables:
  stressNew,stateNew,enerInternNew,enerInelasNew)

  Implicit Double Precision (a-h, o-z)
  Integer, parameter :: j_sys_Dimension = 2, maxblk = 512

  Dimension jblock(*),props(*),density(*),coordMp(*),charLength(*),strainInc(*),relSpinInc(*),tempOld(*),tempNew(*), &
    stretchOld(*),defgradOld(*),fieldOld(*),stressOld(*),stateOld(*),enerInternOld(*),                               &
    enerInelasOld(*),stretchNew(*),defgradNew(*),fieldNew(*),stressNew(*),stateNew(*),                               &
    enerInternNew(*),enerInelasNew(*)

  Character(len=80) :: cmname

  Integer, parameter :: i_nblock=1, i_npt=2, i_layer=3, i_kspt=4, i_noel=5

  ! -------------------------------------------------------------------- !
  !    End VUMAT standard interface                                      !
  ! -------------------------------------------------------------------- !

  ! Throw error on single precision
  If (PRECISION(dt) < 15) Then
    Print *, 'ERROR: Found single precision data. CompDam must be run in double precision.'
    Call XPLB_EXIT
  Else
    Call CompDam(  &
      ! Read only (unmodifiable) variables:
      jblock(i_nblock),ndir,nshr,nstatev,nfieldv,nprops,lanneal,stepTime,     &
      totalTime,dt,cmname,coordMp,charLength,props,density,         &
      strainInc,relSpinInc,tempOld,stretchOld,defgradOld,fieldOld,  &
      stressOld,stateOld,enerInternOld,enerInelasOld,tempNew,       &
      stretchNew,defgradNew,fieldNew,                               &
      ! Write only (modifiable) variables:
      stressNew,stateNew,enerInternNew,enerInelasNew,               &
      ! Read only, extra arguments:
      jblock(i_noel), jblock(i_npt), jblock(i_layer), jblock(i_kspt))
  End If
End Subroutine VUMAT


Subroutine CompDam(  &
  ! Read only (unmodifiable) variables:
  nblock,ndir,nshr,nstatev,nfieldv,nprops,lanneal,stepTime,     &
  totalTime,dt,cmname,coordMp,charLength,props,density_abq,     &
  strainInc,relSpinInc,tempOld,stretchOld,defgradOld,fieldOld,  &
  stressOld,stateOld,enerInternOld,enerInelasOld,tempNew,       &
  stretchNew,defgradNew,fieldNew,                               &
  ! Write only (modifiable) variables:
  stressNew,stateNew,enerInternNew,enerInelasNew,               &
  ! Read only, extra arguments:
  nElement, nMatPoint, nLayer, nSecPoint)

  Use forlog_Mod
  Use matrixAlgUtil_Mod
  Use matProp_Mod
  Use stateVar_Mod
  Use parameters_Mod
  Use DGD_Mod
  Use plasticity_mod
  Use cohesive_mod
  Use friction_mod

  Implicit Double Precision (a-h, o-z)
  Integer, parameter :: j_sys_Dimension = 2, maxblk = 512

  Double Precision :: props(nprops), density_abq(nblock), coordMp(nblock,*), charLength(nblock,*),       &
    strainInc(nblock,ndir+nshr), relSpinInc(nblock,nshr), tempOld(nblock), tempNew(nblock),              &
    stretchOld(nblock,ndir+nshr), defgradOld(nblock,ndir+nshr+nshr), fieldOld(nblock,nfieldv),           &
    stressOld(nblock,ndir+nshr), stateOld(nblock,nstatev), enerInternOld(nblock), enerInelasOld(nblock), &
    stretchNew(nblock,ndir+nshr), defgradNew(nblock,ndir+nshr+nshr), fieldNew(nblock,nfieldv),           &
    stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev), enerInternNew(nblock), enerInelasNew(nblock)

  ! Extra arguments
  Integer :: nElement(nblock), nMatPoint, nLayer, nSecPoint

  Double Precision :: stepTime, totalTime, dt

  Character(len=80) :: cmname

  ! -------------------------------------------------------------------- !
  !    End VUMAT standard interface                                      !
  ! -------------------------------------------------------------------- !
  ! Damage variables, failure indices, damage thresholds
  Double Precision :: d1C_dummy, rfC_dummy

  ! Stress
  Double Precision :: Cauchy(3,3)
  Double Precision :: CauchyABQ(3,3)
  Double Precision :: T_coh(3), trial_traction_int(3), trial_traction_ele(3)  ! Cohesive tractions (shear, normal, shear)
  Double Precision :: Rot(3,3)

  ! Other
  Double Precision :: F(3,3)      ! Current deformation gradient tensor
  Double Precision :: F_old(3,3)  ! Previous deformation gradient tensor
  Double Precision :: U(3,3)      ! Current stretch tensor
  Double Precision :: density_current  ! Current density calculated from Fbulk (ignores the cohesive crack displacement)
  Logical :: Sliding
  
  ! Cohesive law
  Double Precision :: Pen(3), Pen_int(3)  ! Cohesive penalty stiffness
  Double Precision :: delta(3), delta_int(3)  ! Cohesive displacement-jump vector
  Double Precision :: dmg_penalty  ! Penalty-based damage variable
  Double Precision :: dGdGc, dGdGc_inc  ! Change in normalized energy dissipation

  ! Parameters
  Double Precision, parameter :: zero=0.d0, one=1.d0, two=2.d0

  ! Structure of all VUMAT arguments
  Type(vumatArg) :: args

  ! For access to the logger
  Type(forlog) :: logger

  ! Material
  Type(matProps) :: m
  Type(materialList), save :: user
  Type(materialList) :: user_temp
  Integer :: mats, mat_index
  Logical :: fiber_tension_snapback, fiber_compression_snapback

  ! State variables
  Type(stateVars) :: sv

  ! Solution parameters
  Type(parameters), save :: p

  ! Fatigue
  Logical :: fatigue_step
  Integer :: fatigue_parameters(2)
  Pointer(ptr_fatigue_int, fatigue_parameters)

  Interface

      Function SMAIntArrayAccess(ID)
        ! Access an existing global integer array.
        Integer(kind=8) :: SMAIntArrayAccess  ! Returns an address that can be associated with a Fortran pointer
        Integer(kind=4) :: ID       ! Array ID
      End Function SMAIntArrayAccess

  End Interface
  ! -------------------------------------------------------------------- !

  ! Initialize structure of VUMAT args
  Call args%init(nblock, ndir, nshr, nstatev, stepTime, totalTime, dt)

  ! Only write git hash on first call
  If (stepTime <= dt .AND. p%logLevel == 0) Then

    ! Initialize the logger for use in loadParameters()
    Call log%init(level=2, vumatArgStruct=args, format=1)

    ! Load the CompDam solution parameters
    p = loadParameters()

  End If

  ! Initialize the logger
  Call log%init(level=p%logLevel, vumatArgStruct=args, format=p%logFormat)

  If (stepTime <= dt) Then

    If (Allocated(user%materials)) Then
      mats = Size(user%materials)
      readMaterial: Do I = 1,mats
        If (user%materials(I)%name == trim(cmname)) Then
          m = user%materials(I)
          Exit readMaterial
        End If

        If (I == mats) Then
          Call log%info("Saving user material #"//trim(str(I+1))//": "//trim(cmname))

          ! Store all existing material data in temp array
          Allocate (user_temp%materials(mats))
          user_temp%materials(:) = user%materials(:)

          ! Increase size of user%materials by 1
          Deallocate (user%materials)
          Allocate (user%materials(mats+1))

          ! Return old material data to resized user%materials
          user%materials(1:mats) = user_temp%materials(1:mats)
          Deallocate (user_temp%materials)

          ! Add new material data to resized user%materials
          user%materials(mats+1) = loadMatProps(cmname, nprops, props)
          m = user%materials(mats+1)  ! Abbreviate

          ! Calculate alpha0
          If (m%matrixDam) Then
            m%alpha0_deg = NINT(m%alpha0*45.d0/ATAN(one))
            m%alpha0 = alpha0_DGD(m)
          End If

          ! Run consistency checks
          If (totalTime == zero) Then
            Call consistencyChecks(m, issueWarnings=.TRUE.)
          Else
            Call consistencyChecks(m, issueWarnings=.FALSE.)
          End If

          ! Make sure updated values in m are retained
          user%materials(mats+1) = m

        End If

      End Do readMaterial

    Else

      Call log%info("Saving initial user material: "//trim(cmname))
      Allocate (user%materials(1))

      ! Load the first material
      user%materials(1) = loadMatProps(cmname, nprops, props)
      m = user%materials(1)  ! Abbreviate

      ! Calculate alpha0
      If (m%matrixDam) Then
        m%alpha0_deg = NINT(m%alpha0*45.d0/ATAN(one))
        m%alpha0 = alpha0_DGD(m)
      End If

      ! Run consistency checks
      If (totalTime == zero) Then
        Call consistencyChecks(m, issueWarnings=.TRUE.)
      Else
        Call consistencyChecks(m, issueWarnings=.FALSE.)
      End If

      ! Make sure updated values in m are retained
      user%materials(1) = m

    End If

  Else  ! After user%materials(:) has been fully populated
    mat_index = getMaterialIndex(user, cmname)
    m = user%materials(mat_index)
  End If

  ! This forces the material to behave elastically during the packager and initial increment. This
  ! was motivated by issues with Intel MKL in the packager.
  If (totalTime == 0) Then
    m%matrixDam = .FALSE.
    m%shearNonlinearity12 = .FALSE.
    m%shearNonlinearity13 = .FALSE.
    m%schapery = .FALSE.
    m%schaefer = .FALSE.
    m%fiberTenDam = .FALSE.
    m%fiberCompDamBL = .FALSE.
    m%fiberCompDamFKT12 = .FALSE.
    m%fiberCompDamFKT13 = .FALSE.
    m%friction = .FALSE.
    fatigue_step = .FALSE.
  Else
    fatigue_step = .FALSE.
#ifndef PYEXT
    ! Is this a fatigue step?
    ptr_fatigue_int = SMAIntArrayAccess(1)
    If (fatigue_parameters(1) == 1) Then
      fatigue_step = .TRUE.
      Call log%debug("This is a fatigue step.")
    End If
#endif
  End If

  ! -------------------------------------------------------------------- !
  master: Do km = 1,nblock  ! Master Loop

  ! Update location information for debugging
  Call log%arg%update(nElement(km), nMatPoint, nLayer, nSecPoint)

  Call log%debug("Master loop. km = " // trim(str(km)))

  ! -------------------------------------------------------------------- !
  !    Recall previous elastic limits, plasticity, and damage states:    !
  ! -------------------------------------------------------------------- !
  sv = loadStateVars(nstatev, stateOld(km,:), m)
  
  ! -------------------------------------------------------------------- !
  !    Cohesive elements:                                                !
  ! -------------------------------------------------------------------- !
  ElementType: If (m%cohesive) Then

    ! Determine current cohesive displacement-jump. Assumes a constitutive thickness equal to one.
    If (nshr == 2) Then
      delta(1) = sv%Fb1 + two*strainInc(km,3)  ! Longitudinal shear direction (1--3)
      delta(2) = sv%Fb2 + strainInc(km,1)      ! Normal direction (3--3)
      delta(3) = sv%Fb3 + two*strainInc(km,2)  ! Transverse shear direction (2--3)
    Else
      delta(1) = sv%Fb1 + two*strainInc(km,2)  ! Shear direction (1--2)
      delta(2) = sv%Fb2 + strainInc(km,1)      ! Normal normal (2--2)
      delta(3) = zero
    End If
    
    ! Store current cohesive displacement-jump
    sv%Fb1 = delta(1)  ! Longitudinal shear displacement-jump
    sv%Fb2 = delta(2)  ! Normal displacement-jump
    sv%Fb3 = delta(3)  ! Transverse shear displacement-jump
    
    ! Define the penalty stiffnesses
    If (m%embedded_cohesive) Then
      Pen(1) = m%G13
      Pen(2) = m%E3
      Pen(3) = m%G23
    Else
      ! Define penalty stiffnesses
      ! Normal (Mode I) penalty stiffness. This is a user-input numerical constant.
      Pen(2) = m%E3
      ! Longitudinal (Pen(1)) and transverse (Pen(3)) shear penalty stiffnesses
      !  Defined according to the material property relationships in Turon (2010). Also includes the normal
      !  compression load dependence of the shear strengths from LaRC04.
      Pen(1) = Pen(2)*m%GYT*(m%SL - m%etaL*MAX(-m%YC, Pen(2)*MIN(zero, delta(2))))**2/(m%GSL*m%YT**2)
      Pen(3) = Pen(2)*m%GYT*(m%ST - m%etaT*MAX(-m%YC, Pen(2)*MIN(zero, delta(2))))**2/(m%GSL*m%YT**2)
    End If
    
    If (totalTime == 0) Then  ! Avoids calculating damage during the packager and the first solution increment
      dmg_penalty = zero
      dGdGc = zero
      T_coh = cohesive_traction(delta, Pen, dmg_penalty)
    
    Else If (m%embedded_cohesive) Then  ! Finite-thickness cohesive law with physical stiffness inputs
      ! Separate the penalty stiffness into bulk and interface penalty stiffness, a la springs in series
      ! Together, Pen(2) and Pen_int(2) have the same stiffness as the user-defined normal penalty stiffness
      Pen_int(2) = Pen(2) * p%penStiffMult
      Pen(2)     = Pen(2) * p%penStiffMult / (p%penStiffMult - one)
      
      ! Penalty-based damage variables from previous increment
      dmg_penalty = damage_area_to_penalty(sv%d2, delta_initiation(m%YT, Pen_int(2)), delta_final(m%YT, m%GYT))
      
      dGdGc = zero  ! Initialize incremental dissipated plastic energy
      
      embedded_cohesive_tolerance = (1.d-8*m%YT)**2  ! Tolerance for convergence of the EmbeddedCohesive loop
      EmbeddedCohesive: Do
      
        delta_int(2) = interface_displacement_jumps(delta(2), Pen(2), dmg_penalty, Pen_int(2), .TRUE.)
          
        Pen_int(1) = Pen_int(2)*m%GYT*(m%SL - m%etaL*MAX(-m%YC, Pen_int(2)*MIN(zero, delta_int(2))))**2/(m%GSL*m%YT**2)
        Pen_int(3) = Pen_int(2)*m%GYT*(m%ST - m%etaT*MAX(-m%YC, Pen_int(2)*MIN(zero, delta_int(2))))**2/(m%GSL*m%YT**2)
        
        delta_int(1) = interface_displacement_jumps(delta(1), Pen(1), dmg_penalty, Pen_int(1), .FALSE.)
        delta_int(3) = interface_displacement_jumps(delta(3), Pen(3), dmg_penalty, Pen_int(3), .FALSE.)
      
        Call cohesive_damage(m, p, delta_int, Pen_int, delta_int(2), sv%B, sv%FIm, fatigue_step, sv%d2, dmg_penalty, dGdGc_inc)
        dGdGc = dGdGc + dGdGc_inc

        ! Calculate the error between the interface traction and the bulk traction, in units of traction squared
        embedded_cohesive_error = zero
        trial_traction_int = cohesive_traction(delta_int, Pen_int, dmg_penalty)
        trial_traction_ele = cohesive_traction(delta - delta_int, Pen, zero)
        Do E = 1,3
          embedded_cohesive_error = embedded_cohesive_error + (trial_traction_int(E) - trial_traction_ele(E))**2
        End Do

        If (embedded_cohesive_error < embedded_cohesive_tolerance) Exit EmbeddedCohesive
      End Do EmbeddedCohesive
      
      If (m%friction .AND. delta(2) <= zero) Then  ! Closed cracks with friction
        Sliding = crack_is_sliding(delta_int, Pen_int, sv%slide, m%mu, m%mu)
        Call crack_traction_and_slip(delta_int, Pen_int, sv%slide, sv%slide, m%mu, m%mu, dmg_penalty, sv%d2, T_coh, Sliding)
      Else  ! Closed cracks without friction and open cracks
        T_coh = cohesive_traction(delta_int, Pen_int, dmg_penalty)
        sv%slide(1) = delta_int(1)
        sv%slide(2) = delta_int(3)
      End If
      
    Else  ! Zero-thickness cohesive law with numerical penalty stiffnesses
      Call cohesive_damage(m, p, delta, Pen, delta(2), sv%B, sv%FIm, fatigue_step, sv%d2, dmg_penalty, dGdGc)
      
      If (m%friction .AND. delta(2) <= zero) Then  ! Closed cracks with friction
        Sliding = crack_is_sliding(delta, Pen, sv%slide, m%mu, m%mu)
        Call crack_traction_and_slip(delta, Pen, sv%slide, sv%slide, m%mu, m%mu, dmg_penalty, sv%d2, T_coh, Sliding)
      Else  ! Closed cracks without friction and open cracks
        T_coh = cohesive_traction(delta, Pen, dmg_penalty)
        sv%slide(1) = delta(1)
        sv%slide(2) = delta(3)
      End If
    End If

    ! Update stress values
    stressNew(km,1) = T_coh(2)  ! Normal stress
    If (nshr == 2) Then  ! 3-D
      stressNew(km,3) = T_coh(1)  ! Longitudinal shear stress (1--3)
      stressNew(km,2) = T_coh(3)  ! Transverse shear stress (2--3)
    Else  ! 2-D
      stressNew(km,2) = T_coh(1)  ! Shear stress
    End If

    ! Update internal and inelastic energy terms
    enerInternNew(km) = DOT_PRODUCT(T_coh, delta)/two
    enerInternNew(km) = enerInternNew(km)/density_abq(km)
    enerInelasNew(km) = dGdGc*(m%GYT + (m%GSL - m%GYT)*sv%B**m%eta_BK)
    enerInelasNew(km) = enerInelasOld(km) + enerInelasNew(km)/density_abq(km)

  ! -------------------------------------------------------------------- !
  !    Solid elements:                                                   !
  ! -------------------------------------------------------------------- !
  Else ElementType
  
    If (totalTime > 0 .AND. sv%STATUS == 0) Cycle master

    ! -------------------------------------------------------------------- !
    !    Deformation Gradient Tensor and Right Stretch Tensor              !
    ! -------------------------------------------------------------------- !
    U = Vec2Matrix(stretchNew(km,:))
    F = Vec2Matrix(defgradNew(km,:))
    F_old = Vec2Matrix(defgradOld(km,:))

    ! -------------------------------------------------------------------- !
    ! As of Abaqus 6.16, the packager receives a defGradNew of (0.999, 0.999, 0.0, 0.001, 0.001)
    ! for S4R elements. The F(3,3) of 0.0 breaks the initial pass through the VUMAT and the model
    ! will not run. The following statement is a workaround to this problem.
    If (totalTime == 0 .AND. nshr == 1) F(3,3) = one

    ! -------------------------------------------------------------------- !
    !    Define the characteristic element lengths                         !
    ! -------------------------------------------------------------------- !
    If (sv%Lc(2) == zero) Then

      sv%Lc(1) = charLength(km, 1)
      sv%Lc(2) = charLength(km, 2)
      If (nshr == 1) Then
        sv%Lc(3) = m%thickness
      Else
        sv%Lc(3) = charLength(km, 3)
      End If
      
      If (p%check_for_snap_back) Then
        Call checkForSnapBack(m, sv%Lc, nElement(km), fiber_tension_snapback, fiber_compression_snapback)
        
        If (fiber_tension_snapback .OR. fiber_compression_snapback) Then
          mat_index = getMaterialIndex(user, cmname)
          If (fiber_tension_snapback) user%materials(mat_index)%fGXT = m%fXT
          If (fiber_compression_snapback) user%materials(mat_index)%fGXC = m%fXC
        End If
      
      End If

      ! Check for strange values of Lc2 (may indicate that vucharlength did not run)
      If (totalTime > 0) Then ! Don't run in packager
        If (sv%Lc(2) .EQ. zero) Then
          Call log%error("Found Lc2 = 0. Perhaps *Characteristic Length is missing from the input deck?")
        Else If (sv%Lc(3) .EQ. zero) Then
          Call log%error("Found Lc3 = 0. Perhaps *Characteristic Length is missing from the input deck?")
        ! Else If (sv%Lc(2) > sv%Lc(1)*10.d0 .OR. sv%Lc(2) < sv%Lc(1)/10.d0) Then
        !   Call log%warn("Found Lc = [" // trim(str(sv%Lc(1))) //','// trim(str(sv%Lc(2))) //','// trim(str(sv%Lc(3))) // "]; perhaps *Characteristic Length is missing from the input deck?")
        End If
      End If

      Call log%debug("Characteristic element lengths:")
      Call log%debug(trim(str(sv%Lc(1)))//' '//trim(str(sv%Lc(2)))//' '//trim(str(sv%Lc(3))))

    End If

    If (totalTime <= dt) Then
      sv%debugpy_count = 0  ! Initialize
    End If

    ! -------------------------------------------------------------------- !
    !    Initialize phi0                                                   !
    ! -------------------------------------------------------------------- !
    If (totalTime <= dt .AND. (m%fiberCompDamFKT12 .OR. m%fiberCompDamFKT13)) Then
      Call initializePhi0(m, sv%Lc, charLength(km, 4:6), sv%phi0_12, sv%phi0_13)
      Call log%debug("Calculated initial misalignments, phi0_12: " // trim(str(sv%phi0_12)) // " and phi0_13: " // trim(str(sv%phi0_13)))
    End If

    ! -------------------------------------------------------------------- !
    !    Initialize inelastic energy to the previous value                 !
    ! -------------------------------------------------------------------- !
    enerInelasNew(km) = enerInelasOld(km)

    ! -------------------------------------------------------------------- !
    !    Initialize fiber failure                                          !
    ! -------------------------------------------------------------------- !
    If (m%fiberCompDamFKT12 .AND. p%fkt_fiber_failure_angle > zero) Then
      sv%Inel12c = intializeFiberFailure(sv%phi0_12, p%fkt_fiber_failure_angle, m%G12, m%aPL, m%nPL)
    Else
      sv%Inel12c = Huge(zero)   ! Turn off fiber failure by setting the associate inelastic strain to a very large number
    End If
    If (m%fiberCompDamFKT13 .AND. p%fkt_fiber_failure_angle > zero) Then
      sv%Inel13c = intializeFiberFailure(sv%phi0_13, p%fkt_fiber_failure_angle, m%G13, m%aPL, m%nPL)
    Else
      sv%Inel13c = Huge(zero)   ! Turn off fiber failure by setting the associate inelastic strain to a very large number
    End If

    ! -------------------------------------------------------------------- !
    !    Damage Calculations:                                              !
    ! -------------------------------------------------------------------- !
    ! Damage initiation prediction
    If (.NOT. (m%matrixDam .AND. sv%d2 > zero) .AND. .NOT. ((m%fiberCompDamFKT12 .OR. m%fiberCompDamFKT13) .AND. sv%d1C > zero)) Then

      Call DGDInit(U,F,F_old,m,p,sv,ndir,nshr,tempNew(km)-m%T_sf,density_abq(km),Cauchy,enerInternNew(km),enerInelasNew(km),fatigue_step)

    End If

    ! Matrix crack damage evolution
    If (m%matrixDam .AND. sv%d2 > zero) Then

      Call DGDEvolve(U,F,F_old,m,p,sv,ndir,nshr,tempNew(km)-m%T_sf,density_abq(km),Cauchy,enerInternNew(km),enerInelasNew(km),fatigue_step)

    ! Fiber compression damage evolution (FKT decomposition)
    Else If ((m%fiberCompDamFKT12 .OR. m%fiberCompDamFKT13) .AND. sv%d1C > zero) Then

      Call DGDKinkband(U,F,F_old,m,p,sv,ndir,nshr,tempNew(km)-m%T_sf,density_abq(km),Cauchy,enerInternNew(km),enerInelasNew(km))

    End If

    ! -------------------------------------------------------------------- !
    !    Determine and store the stress tensor:                            !
    ! -------------------------------------------------------------------- !
    ! Rotation tensor
    Rot = MATMUL(F, MInverse(U))
    ! Cauchy stress in the current configuration
    CauchyABQ = MATMUL(TRANSPOSE(Rot), MATMUL(Cauchy, Rot))
    ! Convert to vector format
    stressNew(km,:) = Matrix2Vec(CauchyABQ, nshr)

    ! -------------------------------------------------------------------- !
    !    Excessive shear strain errors                                      !
    ! -------------------------------------------------------------------- !
    If (m%shearNonlinearity12) Then
      If (sv%Inel12 < stateOld(km,13)) Then
        Call writeDGDArgsToFile(m,p,sv,U,F,F_old,ndir,nshr,tempNew(km)-m%T_sf,density_abq(km),log%arg,'CompDam_DGD')
        Call log%terminate('Decrease in inelastic strain 12.')
      End If
    End If
    If (m%shearNonlinearity13) Then
      If (sv%Inel13 < stateOld(km,21)) Then
        Call writeDGDArgsToFile(m,p,sv,U,F,F_old,ndir,nshr,tempNew(km)-m%T_sf,density_abq(km),log%arg,'CompDam_DGD')
        Call log%terminate('Decrease in inelastic strain 13.')
      End If
    End If
    If (m%shearNonlinearity12 .OR. m%shearNonlinearity13) Then
      If (sv%Inel12 > 1.d0 .OR. sv%Inel13 > 1.d0) Then
        Call writeDGDArgsToFile(m,p,sv,U,F,F_old,ndir,nshr,tempNew(km)-m%T_sf,density_abq(km),log%arg,'CompDam_DGD')
        Call log%terminate('Excessive inelastic shear strain.')
      End If
    End If

  End If ElementType

  ! -------------------------------------------------------------------- !
  !    Store the updated state variables:                                !
  ! -------------------------------------------------------------------- !
  stateNew(km,:) = storeStateVars(sv, nstatev, m)

  ! -------------------------------------------------------------------- !
  !    Kill job for debugging purposes                                   !
  ! -------------------------------------------------------------------- !
  If (p%logLevel > 2 .AND. p%debug_kill_at_total_time > zero) Then
    If (log%arg%totalTime >= p%debug_kill_at_total_time) Then
      Call writeDGDArgsToFile(m,p,sv,U,F,F_old,ndir,nshr,tempNew(km)-m%T_sf,density_abq(km),log%arg,'CompDam_DGD')
      Call log%error('debug_kill_at_total_time condition satisfied; terminating job.')
    End If
  End If

  End Do master ! End Master Loop

  Call log%debug('End of VUMAT')

  Return
End Subroutine CompDam
