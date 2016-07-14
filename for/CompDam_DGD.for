! This software may be used, reproduced, and provided to others only as permitted under the terms of the agreement
! under which it was acquired from the U.S. Government. Neither title to, nor ownership of, the software is hereby
! transferred. This notice shall remain on all copies of the software.
!
! Copyright 2016 United States Government as represented by the Administrator of the National Aeronautics and
! Space Administration. No copyright is claimed in the United States under Title 17, U.S. Code. All Other Rights
! Reserved.

#include "VUMATArgs.for"
#include "forlog.for"
#include "matrixUtil.for"
#include "stateVar.for"
#include "clt.for"
#include "matProp.for"
#include "parameters.for"
#include "strain.for"
#include "fiberDamage.for"
#include "friction.for"
#include "cohesive.for"


Subroutine VUMAT(  &
  ! Read only (unmodifiable) variables:
  nblock,ndir,nshr,nstatev,nfieldv,nprops,lanneal,stepTime,     &
  totalTime,dt,cmname,coordMp,charLength,props,density,         &
  strainInc,relSpinInc,tempOld,stretchOld,defgradOld,fieldOld,  &
  stressOld,stateOld,enerInternOld,enerInelasOld,tempNew,       &
  stretchNew,defgradNew,fieldNew,                               &
  ! Write only (modifiable) variables:
  stressNew,stateNew,enerInternNew,enerInelasNew)

  Include 'vaba_param.inc'

  Dimension props(nprops),density(nblock),coordMp(nblock,*),charLength(nblock),strainInc(nblock,ndir+nshr),relSpinInc(nblock,nshr),tempOld(nblock),tempNew(nblock),      &
    stretchOld(nblock,ndir+nshr),defgradOld(nblock,ndir+nshr+nshr),fieldOld(nblock,nfieldv),stressOld(nblock,ndir+nshr),stateOld(nblock,nstatev),enerInternOld(nblock),  &
    enerInelasOld(nblock),stretchNew(nblock,ndir+nshr),defgradNew(nblock,ndir+nshr+nshr),fieldNew(nblock,nfieldv),stressNew(nblock,ndir+nshr),stateNew(nblock,nstatev),  &
    enerInternNew(nblock),enerInelasNew(nblock)

  Character*80 cmname

  ! -------------------------------------------------------------------- !
  !    End VUMAT standard interface                                      !
  ! -------------------------------------------------------------------- !

  ! Throw error on single precision
  If (PRECISION(dt) < 15) Then
    Print *, 'ERROR: Found single precision data. CompDam must be run in double precision.'
  Else
    Call CompDam(  &
      ! Read only (unmodifiable) variables:
      nblock,ndir,nshr,nstatev,nfieldv,nprops,lanneal,stepTime,     &
      totalTime,dt,cmname,coordMp,charLength,props,density,         &
      strainInc,relSpinInc,tempOld,stretchOld,defgradOld,fieldOld,  &
      stressOld,stateOld,enerInternOld,enerInelasOld,tempNew,       &
      stretchNew,defgradNew,fieldNew,                               &
      ! Write only (modifiable) variables:
      stressNew,stateNew,enerInternNew,enerInelasNew)
  End If
End Subroutine VUMAT


Subroutine CompDam(  &
  ! Read only (unmodifiable) variables:
  nblock,ndir,nshr,nstatev,nfieldv,nprops,lanneal,stepTime,     &
  totalTime,dt,cmname,coordMp,charLength,props,density,         &
  strainInc,relSpinInc,tempOld,stretchOld,defgradOld,fieldOld,  &
  stressOld,stateOld,enerInternOld,enerInelasOld,tempNew,       &
  stretchNew,defgradNew,fieldNew,                               &
  ! Write only (modifiable) variables:
  stressNew,stateNew,enerInternNew,enerInelasNew)

  Use forlog_Mod
  Use matrixAlgUtil_Mod
  Use matProp_Mod
  Use stateVar_Mod
  Use parameters_Mod

  Include 'vaba_param.inc'

  Double Precision :: props(nprops), density(nblock),coordMp(nblock,*), charLength(nblock), strainInc(nblock,ndir+nshr), relSpinInc(nblock,nshr), tempOld(nblock), tempNew(nblock),  &
    stretchOld(nblock,ndir+nshr), defgradOld(nblock,ndir+nshr+nshr), fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr), stateOld(nblock,nstatev), enerInternOld(nblock),         &
    enerInelasOld(nblock), stretchNew(nblock,ndir+nshr), defgradNew(nblock,ndir+nshr+nshr), fieldNew(nblock,nfieldv), stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),         &
    enerInternNew(nblock), enerInelasNew(nblock)

  Double Precision :: stepTime, totalTime, dt

  Character*80 cmname

  ! -------------------------------------------------------------------- !
  !    End VUMAT standard interface                                      !
  ! -------------------------------------------------------------------- !
  ! Element size
  Double Precision :: Lc(3)

  ! Damage variables, failure indices, damage thresholds
  Double Precision :: d1C_dummy, rfC_dummy

  ! Stress
  Double Precision :: Cauchy(3,3)
  Double Precision :: CauchyABQ(3,3)
  Double Precision :: Rot(3,3)

  ! Other
  Double Precision :: F(3,3)                                               ! Current deformation gradient tensor
  Double Precision :: F_old(3,3)                                           ! Previous deformation gradient tensor
  Double Precision :: U(3,3)

  ! Parameters
  Double Precision, parameter :: zero=0.d0, one=1.d0, two=2.d0

  ! Structure of all VUMAT argments
  Type(VUMATArg) :: args

  ! For access to the logger
  Type(forlog) :: logger

  ! Material
  Type(matProps) :: m
  Type(materialList), save :: user
  Type(materialList) :: user_temp
  Integer :: mats

  ! State variables
  Type(stateVars) :: sv

  ! Solution parameters
  Type(parameters), save :: p

  ! -------------------------------------------------------------------- !

  ! Initialize structure of VUMAT args
  Call args%init(nblock, ndir, nshr, nstatev, stepTime, totalTime, dt)

  ! Only write git hash on first call
  If ((totalTime .LE. DT) .AND. (p%logLevel .EQ. 0)) Then

    ! Initialize the logger for use in loadParameters()
    Call log%init(level=2, VUMATArgStruct=args, format=p%logFormat)

    ! Load the CompDam solution parameters
    p = loadParameters()

    Call log%info("CompDam parameters initialized")

  End If

  ! Initialize the logger
  Call log%init(level=p%logLevel, VUMATArgStruct=args, format=p%logFormat)

  If (totalTime .LE. DT) Then

    If (Allocated(user%materials)) Then
      mats = Size(user%materials)
      readMaterial: Do I = 1,mats
        If (user%materials(I)%name .EQ. trim(cmname)) Then
          m = user%materials(I)
          Exit readMaterial
        End If

        If (I .EQ. mats) Then
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
          If (totalTime .EQ. zero) Then
            user%materials(mats+1) = loadMatProps(cmname, .TRUE., nprops, props)
          Else
            user%materials(mats+1) = loadMatProps(cmname, .FALSE., nprops, props)
          End If
          m = user%materials(mats+1)
        End If

      End Do readMaterial

    Else

      Call log%info("Saving initial user material: "//trim(cmname))
      Allocate (user%materials(1))
      If (totalTime .EQ. zero) Then
        user%materials(1) = loadMatProps(cmname, .TRUE., nprops, props)
      Else
        user%materials(1) = loadMatProps(cmname, .FALSE., nprops, props)
      End If
      m = user%materials(1)

    End If

  Else  ! After user%materials(:) has been fully populated
    Do I = 1,Size(user%materials)
      If (user%materials(I)%name .EQ. trim(cmname))  m = user%materials(I)
    End Do
  End If

  ! -------------------------------------------------------------------- !
  master: Do km = 1,nblock  ! Master Loop

  Call log%debug("Master loop. km = " // trim(str(km)))

  ! -------------------------------------------------------------------- !
  !    Deformation Gradient Tensor and Right Stretch Tensor              !
  ! -------------------------------------------------------------------- !
  U = Vec2Matrix(stretchNew(km,:))
  F = Vec2Matrix(defgradNew(km,:))
  F_old = Vec2Matrix(defgradOld(km,:))

  ! -------------------------------------------------------------------- !
  ! As of Abaqus 6.16, the packager recieves a defGradNew of (0.999, 0.999, 0.0, 0.001, 0.001)
  ! for S4R elements. The F(3,3) of 0.0 breaks the initial pass through the VUMAT and the model
  ! will not run. The following statement is a workaround to this problem.
  If (totalTime .EQ. 0 .AND. nshr .EQ. 1) F(3,3) = one

  ! -------------------------------------------------------------------- !
  !    Recall previous elastic limits, plasticity, and damage states:    !
  ! -------------------------------------------------------------------- !
  sv = loadStateVars(nstatev, stateOld(km,:))

  ! The sign of the change in shear strain, used in the shear nonlinearity subroutine. This was previously a state variable.
  If (m%shearNonlinearity) sv%d_eps12 = Sign(one, (F(1,1)*F(1,2) + F(2,1)*F(2,2) + F(3,2)*F(3,1)) - (F_old(1,1)*F_old(1,2) + F_old(2,1)*F_old(2,2) + F_old(3,2)*F_old(3,1)))

  ! -------------------------------------------------------------------- !
  !    Define the characteristic element lengths                         !
  ! -------------------------------------------------------------------- !
  Lc = getCharElemLengths(m, nshr, charLength(km), sv%Lc1)
  sv%Lc1 = Lc(1)
  If (totalTime .LE. DT) Call checkForSnapBack(m, Lc)

  ! -------------------------------------------------------------------- !
  !    Damage Calculations:                                              !
  ! -------------------------------------------------------------------- !

  ! Damage initiation prediction
  If (.NOT. (m%matrixDam .AND. sv%d2 .GT. zero) .AND. .NOT. (m%fiberCompDamNew .AND. sv%d1C .GT. zero)) Then

    Call DGDInit(m,p,sv,Lc,U,F,ndir,nshr,tempNew(km),Cauchy,enerInternNew(km))

  End IF

  ! Matrix crack damage evolution
  If (m%matrixDam .AND. sv%d2 .GT. zero) Then

    Call DGDEvolve(m,p,sv,Lc,U,F,F_old,ndir,nshr,tempNew(km),Cauchy,enerInternNew(km))

  End IF

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
  !    Store the updated state variables:                                !
  ! -------------------------------------------------------------------- !
  stateNew(km,:) = storeStateVars(sv, nstatev)

  ! -------------------------------------------------------------------- !
  !    Store the internal energy                                         !
  ! -------------------------------------------------------------------- !
  enerInternNew(km) = enerInternNew(km)/density(km)

  End Do master ! End Master Loop

  Call log%debug('End of VUMAT')

  Return
End Subroutine CompDam



! <<<<<<<<<<<<<<<<<<<<<< Subroutine DGDInit >>>>>>>>>>>>>>>>>>>>>>>>>> !
!                                                                      !
!   Checks for the initiation of matrix damage, represented as a DGD   !
!    cohesive crack. If the crack orientation is a priori unknown, it  !
!    will be determined in this subroutine.                            !
!                                                                      !
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !
Subroutine DGDInit(m,p,sv,Lc,U,F,ndir,nshr,DT,Cauchy,enerIntern)

  Use forlog_Mod
  Use matrixAlgUtil_Mod
  Use matProp_Mod
  Use stateVar_Mod
  Use parameters_Mod
  Use CLT_Mod
  Use cohesive_mod
  Use strain_mod
  Use CDM_fiber_mod

  Include 'vaba_param.inc'

  ! -------------------------------------------------------------------- !
  ! Arguments
  Type(matProps), intent(IN) :: m
  Type(parameters), intent(IN) :: p
  Type(stateVars), intent(INOUT) :: sv
  Double Precision, intent(IN) :: F(3,3), U(3,3)                         ! Deformation gradient stretch tensor
  Double Precision, intent(IN) :: Lc(3)                                  ! Characteristic element size
  Integer, intent(IN) :: ndir
  Integer, intent(IN) :: nshr
  Double Precision, intent(IN) :: DT
  Double Precision, intent(OUT) :: Cauchy(ndir,ndir)                     ! Cauchy stress
  Double Precision, intent(OUT) :: enerIntern                            ! Internal energy

  ! -------------------------------------------------------------------- !
  ! Locals

  Double Precision :: Stiff(ndir+nshr,ndir+nshr)                         ! Stiffness
  Double Precision :: eps(ndir,ndir)                                     ! Strain
  Double Precision :: stress(ndir,ndir)                                  ! Stress
  Double Precision :: F_inverse_transpose(3,3)                           ! Inverse transpose of the Deformation Gradient Tensor

  ! Cohesive surface
  Double Precision :: normal(3)                                          ! Normal vector (to cohesive surface)
  Double Precision :: R_cr(3,3)                                          ! Basis coordinate system for the cohesive surface
  Double Precision :: Pen(3)                                             ! Penalty stiffnesses
  Double Precision :: T(3)                                               ! Tractions on the cohesive surface
  Double Precision :: delta(3)                                           ! Current displacement jumps in crack coordinate system
  Double Precision :: B_temp, beta                                       ! Placeholder (temp.) variables for Mode-mixity
  Double Precision :: FImT_temp

  ! Matrix crack cohesive surface normal
  Double Precision :: alpha_temp                                         ! Current alpha (used in loop through possible alphas)
  Integer :: alpha_test, alphaQ                                          ! Alpha = normal to matrix crack cohesive surface
  Integer :: Q                                                           ! Flag: Q=2 for matrix crack; Q=3 for delamination
  Integer :: A, A_min, A_max                                             ! Range through which the code searches for alpha

  ! Miscellaneous
  Double Precision :: rad_to_deg, deg_to_rad

  Double Precision, parameter :: zero=0.d0, one=1.d0, two=2.d0
  ! -------------------------------------------------------------------- !

  ! Miscellaneous constants
  rad_to_deg = 45.d0/ATAN(one)  ! Converts radians to degrees when multiplied
  deg_to_rad = one/rad_to_deg   ! Converts degrees to radians when multiplied

  ! Initialize outputs
  sv%d1    = zero
  sv%d2    = zero
  sv%FImT  = zero
  sv%B     = zero
  sv%Fb1   = zero
  sv%Fb2   = zero
  sv%Fb3   = zero

  F_inverse_transpose = MInverse(TRANSPOSE(F))

  ! Compute the strains: eps, Plas12, Inel12
  ! Note the first argument is a flag to define the strain to use
  Call Strains(m, F, U, DT, ndir, eps, sv%Plas12, sv%Inel12, sv%d_eps12)

  ! Check fiber tension or fiber compression damage
  If (eps(1,1) .GE. zero) Then  ! Fiber tension

    ! Set rfC for failure index output
    If (sv%rfC .EQ. one) Then
      sv%rfC = zero
    End If

    ! Evaluate fiber tension failure criteria and damage variable
    If (m%fiberTenDam) Then
      Call FiberTenDmg(eps, ndir, m%E1, m%XT, m%GXT, m%fXT, m%fGXT, Lc(1), sv%rfT, sv%d1T, sv%d1C, sv%STATUS)
      Call log%debug('Computed fiber damage variable, d1T ' // trim(str(sv%d1T)))

      sv%d1 = sv%d1T
    End If

    ! Build the stiffness matrix
    Stiff = StiffFunc(ndir+nshr, m%E1, m%E2, m%E3, m%G12, m%G13, m%G23, m%v12, m%v13, m%v23, sv%d1, zero, zero)

    ! Calculate stress
    stress = Hooke(Stiff, eps, nshr)
    Cauchy = convertToCauchy(stress, m%strainDef, F, U)

  Else  ! Compression in 1-dir

    ! Set rfT for failure index output
    If (sv%rfT .EQ. one) Then
      sv%rfT = zero
    End If

    If (m%fiberCompDam) Then
      Call FiberCompDmg(eps, ndir, m%E1, m%XC, m%GXC, m%fXC, m%fGXC, Lc(1), sv%rfT, sv%rfC, sv%d1T, sv%d1C, sv%STATUS)
      Call log%debug('Computed fiber damage variable, d1C ' // trim(str(sv%d1C)))

      sv%d1 = sv%d1C
    End If

    ! Build the stiffness matrix
    Stiff = StiffFunc(ndir+nshr, m%E1, m%E2, m%E3, m%G12, m%G13, m%G23, m%v12, m%v13, m%v23, sv%d1, zero, zero)

    ! Calculate stress
    stress = Hooke(Stiff, eps, nshr)
    Cauchy = convertToCauchy(stress, m%strainDef, F, U)

  End If

  ! -------------------------------------------------------------------- !
  !    Search for matrix crack initiation only when no fiber damage has occured
  ! -------------------------------------------------------------------- !
  If (m%matrixDam .AND. sv%d1T .EQ. zero .AND. sv%d1C .EQ. zero) Then
    ! Get fiber direction
    R_cr(:,1) = Norm(F(:,1))  ! fiber direction
    normal(1) = zero

    ! alphaQ is the angle between intralaminar and interlaminar oriented cracks
    alphaQ = FLOOR(ATAN(Lc(2)/Lc(3))*rad_to_deg)
    alphaQ = alphaQ - MOD(alphaQ, p%alpha_inc)

    ! Search through range of alphas to find the correct one (alpha=-999 is a flag to run this search)
    If (sv%alpha .EQ. -999) Then
      A_min = -alphaQ
      A_max = -alphaQ + 170
    ! Alpha was specified in the initial conditions, use the specified angle
    Else
      ! TODO - check that alpha is a valid angle.
      A_min = sv%alpha
      A_max = sv%alpha
    End If

    ! -------------------------------------------------------------------- !
    !    Loop through possible alphas, save the angle where the FC is max  !
    ! -------------------------------------------------------------------- !
    A = A_min
    CrackAngle: Do  ! Test various alphas

      alpha_temp = A*deg_to_rad  ! crack angle being evaluated (converted to radians)

      ! Crack normal in the reference configuration
      normal(2) = COS(alpha_temp)
      normal(3) = SIN(alpha_temp)

      ! Current crack normal direction
      R_cr(:,2) = Norm(MATMUL(F_inverse_transpose, normal))

      ! Current transverse direction
      R_cr(:,3) = CrossProduct(R_cr(:,1), R_cr(:,2))

      ! -------------------------------------------------------------------- !
      !    Determine the cohesive penalty stiffnesses                        !
      ! -------------------------------------------------------------------- !
      ! Does this DGD crack represent a crack or a delamination?
      If (A .LE. alphaQ .OR. nshr .EQ. 1) Then ! In-plane Crack or Two-dimensional
        Q = 2
      Else ! Delamination and Three-dimensional
        Q = 3
      End If
      Pen(2) = p%penStiffMult*m%E2/Lc(Q)
      Pen(1) = Pen(2)*m%GYT*m%SL*m%SL/(m%GSL*m%YT*m%YT) ! Corresponds to Turon et al (2010)
      Pen(3) = Pen(2)*m%GYT*m%ST*m%ST/(m%GSL*m%YT*m%YT)

      ! -------------------------------------------------------------------- !
      !    Determine the cohesive displacement-jump                          !
      ! -------------------------------------------------------------------- !
      T = MATMUL(Cauchy, R_cr(:,2)) ! Traction on fracture surface

      delta = MATMUL(TRANSPOSE(R_cr), T) / Pen

      ! -------------------------------------------------------------------- !
      !    Evaluate the cohesive law initiation criterion                    !
      ! -------------------------------------------------------------------- !
      Call cohesive_damage(m, delta, Pen, delta(2), B_temp, FImT_temp)

      ! -------------------------------------------------------------------- !
      !    Save the values corresponding to the maximum failure criteria     !
      ! -------------------------------------------------------------------- !
      If (FImT_temp .GT. sv%FImT) Then
        sv%FImT       = FImT_temp
        sv%B          = B_temp
        alpha_test    = A
        sv%Fb1        = F(1,Q)
        sv%Fb2        = F(2,Q)
        sv%Fb3        = F(3,Q)
        sv%mCompInit  = MIN(zero, delta(2))
      End If

      If (A .EQ. A_max) EXIT CrackAngle

      ! Advance the crack angle
      If (A .LT. m%alpha0_deg .AND. A + p%alpha_inc .GT. m%alpha0_deg) Then
        A = m%alpha0_deg
      Else If (A .EQ. m%alpha0_deg) Then
        A = A + p%alpha_inc - MOD(A + p%alpha_inc, p%alpha_inc)
      Else
        A = A + p%alpha_inc
      End If

    End Do CrackAngle

    ! -------------------------------------------------------------------- !
    !    If failure occurs, save alpha and indicate small dmg              !
    ! -------------------------------------------------------------------- !
    If (sv%FImT .GE. one) Then
      sv%d2    = 1.d-6 ! Used as a flag to call DGDEvolve
      sv%alpha = alpha_test
      Call log%info('DGDInit found FImT > one. Matrix damage initiated.')
    End If
  End If

  ! -------------------------------------------------------------------- !
  !    Update elastic energy variable.                                   !
  ! -------------------------------------------------------------------- !
  enerIntern = zero
  Do I=1,3
    Do J=1,3
      enerIntern = enerIntern + 0.5d0*stress(I,J)*eps(I,J)
    End Do
  End Do

  Return
End Subroutine DGDInit


! <<<<<<<<<<<<<<<<<<<<<< Subroutine DGDEvolve >>>>>>>>>>>>>>>>>>>>>>>> !
!                                                                      !
!   Determines the matrix damage state variable based on the current   !
!   deformation and mode mixity.                                       !
!                                                                      !
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !
Subroutine DGDEvolve(m,p,sv,Lc,U,F,F_old,ndir,nshr,DT,Cauchy,enerIntern)

  Use forlog_Mod
  Use matrixAlgUtil_Mod
  Use matProp_Mod
  Use stateVar_Mod
  Use parameters_Mod
  Use CLT_Mod
  Use cohesive_mod
  Use strain_mod
  Use CDM_fiber_mod
  Use friction_mod

  Include 'vaba_param.inc'

  ! -------------------------------------------------------------------- !
  ! Arguments
  Type(matProps), intent(IN) :: m
  Type(parameters), intent(IN) :: p
  Type(stateVars), intent(INOUT) :: sv
  Double Precision, Intent(IN) :: Lc(3)                                    ! Characteristic element size
  Double Precision, Intent(IN) :: F(3,3), U(3,3), F_old(3,3)               ! Deformation gradient, stretch tensor
  Double Precision, Intent(IN) :: DT                                       ! Delta temperature
  Integer, intent(IN) :: ndir, nshr
  Double Precision, Intent(OUT) :: Cauchy(ndir,ndir)
  Double Precision, Intent(OUT) :: enerIntern

  ! -------------------------------------------------------------------- !
  ! Locals
  Integer :: mode                                                          ! Flag for whether the crack is a standard matrix crack (0) or a fiber compression matrix crack (1)

  Double Precision :: Stiff(ndir+nshr,ndir+nshr)                           ! Stiffness
  Double Precision :: stress(ndir,ndir)                                    ! Stress (energy conjugate to strain definition)
  Double Precision :: eps(ndir,ndir)                                       ! Strain

  ! Cohesive surface
  Double Precision :: alpha_r                                              ! alpha in radians
  Double Precision :: R_cr(3,3)                                            ! Basis coordinate system for the cohesive surface
  Double Precision :: T(3)                                                 ! Traction on crack interface from bulk material stresses
  Double Precision :: T_coh(3)                                             ! Traction on crack interface from cohesive law
  Double Precision :: normal(3)                                            ! Normal to the crack plane in the reference configuration
  Double Precision :: dcr(ndir,ndir)                                       ! matrix of cohesive displacements (D in COMPA)
  Double Precision :: dMAX                                                 ! Maximum value for damage variable
  Double Precision :: Pen(3)                                               ! Penalty stiffnesses
  Double Precision :: damage_old, AdAe
  Integer :: Q, alphaQ                                                     ! Flag to specify matrix crack or delamination (Q=2 matrix crack, Q=3 delam); angle at which transition occurs (depends on element geometry)
  Integer :: MD                                                            ! MatrixDamage loop index and max iterations

  ! Equilibrium loop
  Double Precision :: F_bulk(3,3), U_bulk(3,3)
  Double Precision :: F_bulk_old(3), F_bulk_change(3), F_bulk_inverse(3,3)
  Double Precision :: Rs(3)                                                ! Residual stress vector
  Double Precision :: tol_DGD
  Double Precision :: err, err_old

  ! For jacobian
  Double Precision :: Cauchy_d(ndir,ndir,3)                                ! Derivative of the Cauchy stress tensor
  Double Precision :: stress_d(ndir,ndir,3)                                ! Derivative of the stress
  Double Precision :: eps_d(ndir,ndir,3)                                   ! Derivative of the strain
  Double Precision :: dcr_d(3,3,3), R_cr_d(3,3,3), F_bulk_d(3,3,3)
  Double Precision :: T_coh_d(3,3), T_d(3,3)
  Double Precision :: r1length, r2length                                   ! Magnitude used for computing derivative
  Double Precision :: Jac(3,3)                                             ! Jacobian
  Double Precision :: T_coh_d_den_temp

  ! Convergence tools
  Double Precision :: crack_inflection_aid, aid, crack_inflection_cutback
  Integer :: restarts, cutbacks, restarts_max
  Integer :: crack_inversions, crack_inversions_max
  Logical :: crack_inverted
  Logical :: crack_open, crack_was_open, crack_open_test
  Logical :: Restart, Cutback

  ! Misc
  Double Precision :: X(ndir,ndir)                                         ! Reference configuration
  Double Precision :: dGdGc
  Double Precision :: Plas12_temp, Inel12_temp, eps12_temp
  Double Precision :: tr, Y(ndir,ndir)
  Double Precision :: eye(ndir,ndir)
  Double Precision :: Fb_s1(3), Fb_s3(3)
  Double Precision :: L

  Double Precision :: rfT_temp, rfC_temp                                   ! Fiber damage thresholds, from previous converged solution
  Double Precision :: d1T_temp, d1C_temp                                   ! Fiber damage variables, from previous converged solution

  ! Friction
  Double Precision :: slide_old(2)
  Integer :: forced_sticking
  Logical :: Sliding

  Double Precision, parameter :: zero=0.d0, one=1.d0, two=2.d0

  ! -------------------------------------------------------------------- !
  Call log%debug('Start of DGDEvolve')

  dMAX = one ! Maximum value for damage variables

  restarts_max = 2  ! equal to the number of starting points minus 1
  crack_inflection_cutback = 0.99d0
  crack_inversions_max = 4

  ! Initialize outputs
  sv%FImT = zero

  X = zero; X(1,1) = Lc(1); X(2,2) = Lc(2); X(3,3) = Lc(3) ! Ref. Config.

  eye = zero; DO I = 1,3; eye(I,I) = one; end DO ! Identity Matrix

  tol_DGD = m%YT*p%tol_DGD_f ! Equilibrium loop tolerance [stress]

  ! crack or delamination?
  alphaQ = FLOOR(ATAN(Lc(2)/Lc(3))*45.d0/ATAN(one))
  alphaQ = alphaQ - MOD(alphaQ, p%alpha_inc)
  If (sv%alpha .LE. alphaQ .OR. nshr .EQ. 1) Then
    Q = 2 ! matrix crack
  Else
    Q = 3 ! delamination
  End If
  alpha_r = sv%alpha/45.d0*ATAN(one) ! alpha [radians]

  ! -------------------------------------------------------------------- !
  ! Penalty stiffness
  Pen(2) = p%penStiffMult*m%E2/Lc(Q)
  Pen(1) = Pen(2)*m%GYT*m%SL*m%SL/(m%GSL*m%YT*m%YT) ! Corresponds to Turon et al (2010)
  Pen(3) = Pen(2)*m%GYT*m%ST*m%ST/(m%GSL*m%YT*m%YT)

  ! -------------------------------------------------------------------- !
  !    Define a crack-based coordinate system with a basis R_cr(3,3):    !
  ! -------------------------------------------------------------------- !
  normal(1) = zero
  normal(2) = COS(alpha_r)
  normal(3) = SIN(alpha_r)

  ! Current fiber direction
  R_cr(:,1) = Norm(F(:,1))

  ! -------------------------------------------------------------------- !
  !    Initial guesses for F_bulk:                                       !
  ! -------------------------------------------------------------------- !
  F_bulk(:,:) = F(:,:)
  F_bulk(1,Q) = sv%Fb1
  F_bulk(2,Q) = sv%Fb2
  F_bulk(3,Q) = sv%Fb3

  IF (Length(F_bulk(:,Q)) .EQ. zero) F_bulk(Q,Q) = one

  ! Initialize the displ across the cohesive interface as zero
  dcr = zero

  ! -------------------------------------------------------------------- !
  !    Definition of Equilibrium alternate starting points:              !
  ! -------------------------------------------------------------------- !
  ! Starting point 1 is the Q-th column of the previous solution for F_bulk
  ! Starting point 2 is the Q-th column of F
  ! Starting point 3 is the cross product of the two non-decomposed columns of F
  Fb_s3(:) = Norm(CrossProduct(F(:,MOD(Q,3)+1), F(:,Q-1)))

  ! -------------------------------------------------------------------- !
  !    MatrixDamage Loop and solution controls definition                !
  ! -------------------------------------------------------------------- !
  MD = 0 ! Counter for MatrixDamage loop

  MatrixDamage: Do
    MD = MD + 1
    Call log%debug('MatrixDamage Start, MD: ' // trim(str(MD)))
    Fb_s1(:) = F_bulk(:,Q) ! Starting point 1
    slide_old(:) = sv%slide(:)

    AdAe = sv%d2/(sv%d2 + (one - sv%d2)*two*Pen(1)*m%GSL/(m%SL*m%SL))

    ! -------------------------------------------------------------------- !
    !    Equilibrium Loop and solution controls definition                 !
    ! -------------------------------------------------------------------- !
    Restart = .False.
    Cutback = .False.
    restarts = 0  ! Indicates the number of loop restarts with new starting points
    cutbacks = 0  ! Cut-back counter
    crack_inversions = 0
    crack_inverted = .False.

    forced_sticking = 0  ! Indicates that sliding has not been suppressed
    aid = one  ! Artificially reduces the rate of change in F_bulk
    crack_inflection_aid = one
    err = Huge(zero)

    Equilibrium: Do ! Loop to determine the current F_bulk
      Call log%debug('Equilibrium Start.')

      ! -------------------------------------------------------------------- !
      If (Restart) Then

        ! Reset the Restart and Cutback flags
        Restart = .False.
        Cutback = .False.

        ! Attempt to use the "no sliding" condition, if not converged conventionally with friction.
        If (forced_sticking .EQ. 0 .AND. m%friction) Then
          Call log%info('Attempting no sliding, MD: ' // trim(str(MD)) // ' Restart: ' // trim(str(restarts)))
          forced_sticking = 1
        Else
          ! Advance the restart counter
          restarts = restarts + 1
          Call log%info('Restarting Equilibrium loop with new start point, Restart: ' // trim(str(restarts)))
          forced_sticking = 0
        End If

        ! If all starting points have been fully used...
        If (restarts .GT. restarts_max) Then
          ! ...and if the matrix damage is already fully developed, delete the element.
          If (sv%d2 .GE. dMAX) Then
            Call log%warn('Deleting failed element for which no solution could be found.')
            sv%STATUS = 0
            Exit MatrixDamage
          End If
          ! ...raise an error and halt the subroutine.
          Call log%error('No starting points produced a valid solution.')
        End If

        cutbacks = 0
        crack_inversions = 0
        crack_inverted = .False.

        ! Restart from a starting point
        F_bulk = F
        If (restarts .EQ. 0) Then
          F_bulk(:,Q) = Fb_s1(:)  ! Use starting point 1
        Else If (restarts .EQ. 1) Then
          Continue                ! Use starting point 2
        Else If (restarts .EQ. 2) Then
          F_bulk(:,Q) = Fb_s3(:)  ! Use starting point 3
        End If

        aid = one
        crack_inflection_aid = one

        ! Reset err
        err = Huge(zero)
      ! -------------------------------------------------------------------- !
      Else If (Cutback) Then

        ! Reset the Cutback flag
        Cutback = .False.

        ! Advance the cutback counter
        cutbacks = cutbacks + 1
        Call log%info('Cutting back, Cutbacks: ' // trim(str(cutbacks)) //', Restart: ' // trim(str(restarts)))

        aid = p%cutback_amount**cutbacks

        F_bulk(:,Q) = F_bulk_old(:) - F_bulk_change(:)*aid*crack_inflection_aid

        ! Reset err to the last converged value
        err = err_old
      End If
      ! -------------------------------------------------------------------- !
      ! Initialize all temporary state variables for use in Equilibrium loop:

      ! Shear nonlinearity variables
      Plas12_temp = sv%Plas12
      Inel12_temp = sv%Inel12

      ! CDM fiber damage variables
      sv%d1 = zero
      d1T_temp = sv%d1T
      d1C_temp = sv%d1C
      rfT_temp = sv%rfT
      rfC_temp = sv%rfC

      ! Store re-used matrices
      F_bulk_inverse = MInverse(F_bulk)

      R_cr(:,2) = MATMUL(TRANSPOSE(F_bulk_inverse), normal)
      r2length = Length(R_cr(:,2))    ! Un-normalized length of R_cr(:,2), used in the Jacobian
      R_cr(:,2) = R_cr(:,2)/r2length  ! Normalized R_cr(:,2)
      R_cr(:,3) = CrossProduct(R_cr(:,1), R_cr(:,2))

      dcr(:,Q) = MATMUL(TRANSPOSE(R_cr), F(:,Q) - F_bulk(:,Q))*X(Q,Q)

      If (dcr(2,Q) .GE. zero) Then
        crack_open = .True.
      Else
        crack_open = .False.
      End If

      ! TODO for other strain component implementation: Calculate U_bulk from F_bulk
      U_bulk = zero

      ! Calculate the sign of the change in shear strain (for shear nonlinearity subroutine)
      If (m%shearNonlinearity .AND. Q .EQ. 2) Then
         eps12_temp = Sign(one, (F(1,1)*F_bulk(1,2) + F(2,1)*F_bulk(2,2) + F_bulk(3,2)*F(3,1)) - (F_old(1,1)*sv%Fb1 + F_old(2,1)*sv%Fb2 + sv%Fb3*F_old(3,1)))
      Else
         eps12_temp = sv%d_eps12
      End If

      ! Compute the strains: eps, Plas12, Inel12
      ! Note the first argument is a flag to define the strain to use
      Call Strains(m, F_bulk, U_bulk, DT, ndir, eps, Plas12_temp, Inel12_temp, eps12_temp)

      ! -------------------------------------------------------------------- !
      !    Evaluate the CDM fiber failure criteria and damage variable:      !
      ! -------------------------------------------------------------------- !
      If (eps(1,1) .GE. zero) Then
        If (m%fiberTenDam) Call FiberTenDmg(eps, ndir, m%E1, m%XT, m%GXT, m%fXT, m%fGXT, Lc(1), rfT_temp, d1T_temp, d1C_temp, sv%STATUS)
        sv%d1 = d1T_temp
        d1C_temp = sv%d1C
      Else If (m%fiberCompDam) Then
        Call FiberCompDmg(eps, ndir, m%E1, m%XC, m%GXC, m%fXC, m%fGXC, Lc(1), rfT_temp, rfC_temp, d1T_temp, d1C_temp, sv%STATUS)
        sv%d1 = d1C_temp
      End If

      ! -------------------------------------------------------------------- !
      !    Build the stiffness matrix:                                       !
      ! -------------------------------------------------------------------- !
      Stiff = StiffFunc(ndir+nshr, m%E1, m%E2, m%E3, m%G12, m%G13, m%G23, m%v12, m%v13, m%v23, sv%d1, zero, zero)

      ! -------------------------------------------------------------------- !
      !    Determine the bulk material tractions on the fracture plane:      !
      ! -------------------------------------------------------------------- !
      stress = Hooke(Stiff, eps, nshr)
      Cauchy = convertToCauchy(stress, m%strainDef, F_bulk, U_bulk)
      T = MATMUL(Cauchy, R_cr(:,2))  ! Traction on fracture surface

      ! -------------------------------------------------------------------- !
      !    Determine the cohesive tractions:                                 !
      ! -------------------------------------------------------------------- !
      If (crack_open) Then  ! Open cracks

        T_coh(:) = Pen(:)*(one - sv%d2)*dcr(:,Q)

        sv%slide(1) = dcr(1,Q)
        sv%slide(2) = dcr(3,Q)

      Else  ! Closed cracks

        If (.NOT. m%friction) Then  ! Closed cracks without friction
          T_coh(1) = Pen(1)*(one - sv%d2)*dcr(1,Q)
          T_coh(2) = Pen(2)*dcr(2,Q)
          T_coh(3) = Pen(3)*(one - sv%d2)*dcr(3,Q)

          sv%slide(1) = dcr(1,Q)
          sv%slide(2) = dcr(3,Q)

        Else  ! Closed cracks with friction
          Sliding = ALFANO_SLIP(dcr(:,Q), Pen, slide_old, m%mu, m%mu)
          If (forced_sticking .EQ. 1) Sliding = .False.
          Call ALFANO(dcr(:,Q), Pen, slide_old, sv%slide, m%mu, m%mu, sv%d2, AdAe, T_coh, Sliding)

        End If

      End If

      ! -------------------------------------------------------------------- !
      !    Define the stress residual vector, R. R is equal to the           !
      !    difference in stress between the cohesive interface and the bulk  !
      !    stress projected onto the cohesive interface                      !
      ! -------------------------------------------------------------------- !
      Rs = T_coh - MATMUL(TRANSPOSE(R_cr), T)

      ! -------------------------------------------------------------------- !
      !    Check for convergence.                                            !
      ! -------------------------------------------------------------------- !
      err_old = err
      err = Length(Rs)

      ! Check for a diverging solution
      If (err .GE. err_old) Then

        If (crack_inverted) Then
          ! A cut-back will not be performed if the crack opening state has just changed. It is necessary to
          ! allow this increase in error so the Jacobian can be calculated on "this side" of the crack state.
          Continue
        Else

          Call log%info('Solution is diverging, err: ' // trim(str(err)) // ' > ' // trim(str(err_old)))

          ! Cut-back using the current starting point
          If (cutbacks .LT. p%cutbacks_max) Then
            Cutback = .True.

          ! Restart using a new starting point, if available
          Else
            Restart = .True.
          End If

          Cycle Equilibrium

        End If

      End If

      ! Ensures that an artificial no "sliding condition" is not forced as the solution
      If (err .LT. tol_DGD .AND. forced_sticking .EQ. 1 .AND. .NOT. crack_open) Then
        If (Sliding .NE. ALFANO_SLIP(dcr(:,Q), Pen, slide_old, m%mu, m%mu)) Then
          forced_sticking = 2  ! Deactivates forced "no sliding" condition
          err = Huge(zero)  ! Resets the error. An increase in error here does not indicate divergence.
          Cycle Equilibrium
        End If
      End If

      ! Check for any inside-out deformation or an overly compressed bulk material
      If (err .LT. tol_DGD .AND. MDet(F_bulk) .LT. p%compLimit) Then

        Call log%warn('det(F_bulk) is below limit: ' // trim(str(MDet(F_bulk))) // ' Restart: ' // trim(str(restarts)))

        ! Restart using new starting point
        Restart = .True.
        Cycle Equilibrium
      End If

      ! If converged,
      If (err .LT. tol_DGD) Then

        Call log%debug('Equilibrium loop found a converged solution.')

        sv%Fb1 = F_bulk(1,Q); sv%Fb2 = F_bulk(2,Q); sv%Fb3 = F_bulk(3,Q)

        ! Update fiber damage state variables
        sv%d1T = d1T_temp
        sv%d1C = d1C_temp
        sv%rfT = rfT_temp
        sv%rfC = rfC_temp

        ! Update shear nonlinearity state variables
        sv%Plas12 = Plas12_temp
        sv%Inel12 = Inel12_temp

        ! If fully damaged
        If (sv%d2 .GE. dMAX) Then
          sv%d2 = dMAX
          sv%FImT = one
          EXIT MatrixDamage
        End If
        EXIT Equilibrium ! Check for change in sv%d2
      End If

      ! -------------------------------------------------------------------- !
      !    Find the derivative of the Cauchy stress tensor.                  !
      ! -------------------------------------------------------------------- !
      F_bulk_d = zero
      R_cr_d   = zero
      dcr_d    = zero
      eps_d    = zero
      stress_d = zero
      Cauchy_d = zero
      T_d      = zero
      T_coh_d  = zero

      DO I=1,3
        F_bulk_d(I,Q,I) = one

        ! info on derivative of x/||x||: http://blog.mmacklin.com/2012/05/
        R_cr_d(:,2,I) = -MATMUL(MATMUL(TRANSPOSE(F_bulk_inverse), MATMUL(TRANSPOSE(F_bulk_d(:,:,I)), TRANSPOSE(F_bulk_inverse))), normal)
        ! info on derivative of inverse matrix: http://planetmath.org/derivativeofinversematrix
        R_cr_d(:,2,I) = MATMUL(eye/r2length - OuterProduct(R_cr(:,2), R_cr(:,2))/r2length, R_cr_d(:,2,I))
        R_cr_d(:,3,I) = CrossProduct(R_cr(:,1), R_cr_d(:,2,I))

        dcr_d(:,Q,I) = MATMUL(TRANSPOSE(R_cr_d(:,:,I)), (F(:,Q) - F_bulk(:,Q)))*X(Q,Q) - MATMUL(TRANSPOSE(R_cr), F_bulk_d(:,Q,I))*X(Q,Q)

        If (crack_open) Then ! Open cracks
          T_coh_d(:,I) = Pen(:)*(one - sv%d2)*dcr_d(:,Q,I)

        Else If (.NOT. m%friction) Then ! Closed cracks without friction
          T_coh_d(1,I) = Pen(1)*(one - sv%d2)*dcr_d(1,Q,I)
          T_coh_d(2,I) = Pen(2)*dcr_d(2,Q,I)
          T_coh_d(3,I) = Pen(3)*(one - sv%d2)*dcr_d(3,Q,I)

        Else If (Sliding) Then ! Closed cracks with sliding friction
          T_coh_d_den_temp = (Pen(1)*(dcr(1,Q) - slide_old(1)))**2 + (Pen(3)*(dcr(3,Q) - slide_old(2)))**2

          T_coh_d(1,I) = Pen(1)*(one - sv%d2)*dcr_d(1,Q,I)
          T_coh_d(3,I) = Pen(3)*(one - sv%d2)*dcr_d(3,Q,I)

          If (T_coh_d_den_temp .NE. zero) Then
            T_coh_d(1,I) = T_coh_d(1,I) - AdAe*m%mu*Pen(2)*Pen(1)/SQRT(T_coh_d_den_temp)* &
              (dcr_d(2,Q,I)*(dcr(1,Q) - slide_old(1)) + dcr(2,Q)*dcr_d(1,Q,I) - dcr(2,Q)*(dcr(1,Q) - slide_old(1))* &
              (Pen(1)*Pen(1)*(dcr(1,Q) - slide_old(1))*dcr_d(1,Q,I) + &
               Pen(3)*Pen(3)*(dcr(3,Q) - slide_old(2))*dcr_d(3,Q,I))/T_coh_d_den_temp)

            T_coh_d(3,I) = T_coh_d(3,I) - AdAe*m%mu*Pen(2)*Pen(3)/SQRT(T_coh_d_den_temp)* &
              (dcr_d(2,Q,I)*(dcr(3,Q) - slide_old(2)) + dcr(2,Q)*dcr_d(3,Q,I) - dcr(2,Q)*(dcr(3,Q) - slide_old(2))* &
              (Pen(1)*Pen(1)*(dcr(1,Q) - slide_old(1))*dcr_d(1,Q,I) + &
               Pen(3)*Pen(3)*(dcr(3,Q) - slide_old(2))*dcr_d(3,Q,I))/T_coh_d_den_temp)
          End If
          T_coh_d(2,I) = Pen(2)*dcr_d(2,Q,I)

        Else ! Closed cracks with sticking friction
          T_coh_d(1,I) = Pen(1)*(one - sv%d2 + AdAe)*dcr_d(1,Q,I)
          T_coh_d(2,I) = Pen(2)*dcr_d(2,Q,I)
          T_coh_d(3,I) = Pen(3)*(one - sv%d2 + AdAe)*dcr_d(3,Q,I)
        End If

        eps_d(:,:,I) = (MATMUL(TRANSPOSE(F_bulk_d(:,:,I)), F_bulk) + MATMUL(TRANSPOSE(F_bulk), F_bulk_d(:,:,I)))/two
        stress_d(:,:,I) = Hooke(Stiff,eps_d(:,:,I),nshr)

        Y = MATMUL(F_bulk_inverse, F_bulk_d(:,:,I))
        tr = Y(1,1) + Y(2,2) + Y(3,3)
        Cauchy_d(:,:,I) = (MATMUL(MATMUL(F_bulk_d(:,:,I), stress) + MATMUL(F_bulk, stress_d(:,:,I)), TRANSPOSE(F_bulk)) + MATMUL(MATMUL(F_bulk, stress), TRANSPOSE(F_bulk_d(:,:,I))))/MDet(F_bulk) - Cauchy*tr

        T_d(:,I) = MATMUL(Cauchy_d(:,:,I), R_cr(:,2)) + MATMUL(Cauchy, R_cr_d(:,2,I))
      End Do

      ! -------------------------------------------------------------------- !
      !    Define the Jacobian matrix, J                                     !
      ! -------------------------------------------------------------------- !
      Jac = zero

      Jac(:,1) = T_coh_d(:,1) - MATMUL(TRANSPOSE(R_cr), T_d(:,1)) - MATMUL(TRANSPOSE(R_cr_d(:,:,1)), T)
      Jac(:,2) = T_coh_d(:,2) - MATMUL(TRANSPOSE(R_cr), T_d(:,2)) - MATMUL(TRANSPOSE(R_cr_d(:,:,2)), T)
      Jac(:,3) = T_coh_d(:,3) - MATMUL(TRANSPOSE(R_cr), T_d(:,3)) - MATMUL(TRANSPOSE(R_cr_d(:,:,3)), T)

      ! -------------------------------------------------------------------- !
      !    Calculate the new bulk deformation gradient                       !
      ! -------------------------------------------------------------------- !
      F_bulk_old(:)    = F_bulk(:,Q)
      F_bulk_change(:) = MATMUL(MInverse(Jac), Rs)
      F_bulk(:,Q)      = F_bulk_old(:) - F_bulk_change(:)*aid

      ! -------------------------------------------------------------------- !
      !    Check for a change in crack opening                               !
      ! -------------------------------------------------------------------- !

      ! Store the previous crack opening state variable
      crack_was_open = crack_open

      ! Update the crack opening state variable
      R_cr(:,2) = Norm(MATMUL(MInverse(TRANSPOSE(F_bulk)), normal))
      R_cr(:,3) = CrossProduct(R_cr(:,1), R_cr(:,2))

      dcr(:,Q) = MATMUL(TRANSPOSE(R_cr), F(:,Q) - F_bulk(:,Q))*X(Q,Q)

      If (dcr(2,Q) .GE. zero) Then
        crack_open = .True.
      Else
        crack_open = .False.
      End If

      crack_inflection_aid = one

      ! Check for a change in the crack opening state
      If (crack_open .EQ. crack_was_open) Then
        ! If there is no change, do nothing.
        crack_inverted = .False.

      Else If (crack_inversions .LT. crack_inversions_max) Then
        Call log%info('Change in crack opening detected.')

        crack_inversions = crack_inversions + 1
        crack_inverted = .True.

        ! Initialize a test variable for the crack opening state
        crack_open_test = crack_open

        ! Loop until the crack opening status changes back
        Do While (crack_open_test .EQ. crack_open)

          crack_inflection_aid = crack_inflection_aid * crack_inflection_cutback

          F_bulk(:,Q) = F_bulk_old(:) - F_bulk_change(:)*aid*crack_inflection_aid

          R_cr(:,2) = Norm(MATMUL(MInverse(TRANSPOSE(F_bulk)), normal))
          R_cr(:,3) = CrossProduct(R_cr(:,1), R_cr(:,2))

          dcr(:,Q) = MATMUL(TRANSPOSE(R_cr), F(:,Q) - F_bulk(:,Q))*X(Q,Q)

          If (dcr(2,Q) .GE. zero) Then
            crack_open_test = .True.
          Else
            crack_open_test = .False.
          End If

        End Do

        ! Return F_bulk to the state immediately before the crack state changed back
        crack_inflection_aid = crack_inflection_aid / crack_inflection_cutback
        F_bulk(:,Q) = F_bulk_old(:) - F_bulk_change(:)*aid*crack_inflection_aid

      Else
        crack_inverted = .False.
      End If
    End Do Equilibrium


    Call log%debug('Exited Equilibrium loop.')

    ! Store the old cohesive damage variable for checking for convergence
    damage_old = sv%d2

    ! Update the cohesive damage variable
    Call cohesive_damage(m, dcr(:,Q), Pen, sv%mCompInit, sv%B, sv%FImT, sv%d2, dGdGc)

    ! Check for damage advancement
    If (sv%d2 .LE. damage_old) Then  ! If there is no damage progression,
      Call log%debug('No change in matrix damage variable, d2 ' // trim(str(sv%d2)))
      EXIT MatrixDamage
    Else
      Call log%debug('Change in matrix damage variable, d2 ' // trim(str(sv%d2)))
    End If

    ! Check for convergence based on rate of energy dissipation
    If (dGdGc .LT. p%dGdGc_min) Then
      Call log%warn('Solution accepted due to small change in dmg.')
      Call log%info('MD: ' // trim(str(MD)) // ' AdAe: ' // trim(str(AdAe)))
      EXIT MatrixDamage
    End If

    ! Limit number of MatrixDamage loop iterations
    If (MD .GT. p%MD_max) Then
      Call log%warn('MatrixDamage loop limit exceeded. MD: ' // trim(str(MD)))
      Call log%info('dGdGc: ' // trim(str(dGdGc)) // ' AdAe: ' // trim(str(AdAe)))
      EXIT MatrixDamage
    End If
  End Do MatrixDamage


  Call log%info('Exited matrix damage loop, MD: ' // trim(str(MD)))

  ! -------------------------------------------------------------------- !
  !    Update elastic energy variable.                                   !
  ! -------------------------------------------------------------------- !
  enerIntern = zero
  Do I=1,3
    Do J=1,3
      enerIntern = enerIntern + 0.5d0*stress(I,J)*eps(I,J)
    End Do
  End Do

  ! -------------------------------------------------------------------- !
  Return
End Subroutine DGDEvolve
