! This software may be used, reproduced, and provided to others only as permitted under the terms of the agreement
! under which it was acquired from the U.S. Government. Neither title to, nor ownership of, the software is hereby
! transferred. This notice shall remain on all copies of the software.
!
! Copyright 2016 United States Government as represented by the Administrator of the National Aeronautics and
! Space Administration. No copyright is claimed in the United States under Title 17, U.S. Code. All Other Rights
! Reserved.

#include "vexternaldb.for"
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
  Parameter (j_sys_Dimension = 2, n_vec_Length = 136, maxblk = n_vec_Length)

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
  totalTime,dt,cmname,coordMp,charLength,props,density,         &
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

  Implicit Double Precision (a-h, o-z)
  Parameter (j_sys_Dimension = 2, n_vec_Length = 136, maxblk = n_vec_Length)

  Double Precision :: props(nprops), density(nblock), coordMp(nblock,*), charLength(nblock,*),           &
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
  Double Precision :: Rot(3,3)

  ! Other
  Double Precision :: F(3,3)                                               ! Current deformation gradient tensor
  Double Precision :: F_old(3,3)                                           ! Previous deformation gradient tensor
  Double Precision :: U(3,3)

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
  Integer :: mats

  ! State variables
  Type(stateVars) :: sv

  ! Solution parameters
  Type(parameters), save :: p

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
    Do I = 1,Size(user%materials)
      If (user%materials(I)%name == trim(cmname))  m = user%materials(I)
    End Do
  End If

  ! -------------------------------------------------------------------- !
  master: Do km = 1,nblock  ! Master Loop

  ! Update location information for debugging
  Call log%arg%update(nElement(km), nMatPoint, nLayer, nSecPoint)

  Call log%debug("Master loop. km = " // trim(str(km)))

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
  !    Recall previous elastic limits, plasticity, and damage states:    !
  ! -------------------------------------------------------------------- !
  sv = loadStateVars(nstatev, stateOld(km,:), m)

  ! The sign of the change in shear strain, used in the shear nonlinearity subroutine. Previously was a state variable.
  If (m%shearNonlinearity12) Then
    sv%d_eps12 = Sign(one, (F(1,1)*F(1,2) + F(2,1)*F(2,2) + F(3,2)*F(3,1)) - (F_old(1,1)*F_old(1,2) + F_old(2,1)*F_old(2,2) + F_old(3,2)*F_old(3,1)))
  End If
  If (m%shearNonlinearity13) Then
    sv%d_eps13 = Sign(one, (F(1,1)*F(1,3) + F(2,1)*F(2,3) + F(3,1)*F(3,3)) - (F_old(1,1)*F_old(1,3) + F_old(2,1)*F_old(2,3) + F_old(3,2)*F_old(3,3)))
  End If

  ! -------------------------------------------------------------------- !
  !    Define the characteristic element lengths                         !
  ! -------------------------------------------------------------------- !
  If (sv%Lc(1) == zero) Then

    sv%Lc(1) = charLength(km, 1)
    sv%Lc(2) = charLength(km, 2)
    If (nshr == 1) Then
      sv%Lc(3) = m%thickness
    Else
      sv%Lc(3) = charLength(km, 3)
    End If

    Call log%debug("Characteristic element lengths:")
    Call log%debug(trim(str(sv%Lc(1)))//' '//trim(str(sv%Lc(2)))//' '//trim(str(sv%Lc(3))))

  End If

  ! If (totalTime <= DT) Call checkForSnapBack(m, sv%Lc, nElement(km))

  ! -------------------------------------------------------------------- !
  !    Initialize phi0                                                   !
  ! -------------------------------------------------------------------- !
  If (totalTime <= DT .AND. m%fiberCompDamFKT) Then
    sv%phi0 = initializePhi0(sv%phi0, m%G12, m%XC, m%aPL, m%nPL, sv%Lc, charLength(km, 4:6))
  End If

  ! -------------------------------------------------------------------- !
  !    Damage Calculations:                                              !
  ! -------------------------------------------------------------------- !

  ! Damage initiation prediction
  If (.NOT. (m%matrixDam .AND. sv%d2 > zero) .AND. .NOT. (m%fiberCompDamFKT .AND. sv%d1C > zero)) Then

    Call DGDInit(U,F,m,p,sv,ndir,nshr,tempNew(km),Cauchy,enerInternNew(km), F_old)

  End IF

  ! Matrix crack damage evolution
  If (m%matrixDam .AND. sv%d2 > zero) Then

    Call DGDEvolve(U,F,F_old,m,p,sv,ndir,nshr,tempNew(km),Cauchy,enerInternNew(km))

  ! Fiber compression damage evolution (New model)
  Else If (m%fiberCompDamFKT .AND. sv%d1C > zero) Then

    Call DGDKinkband(U,F,F_old,m,p,sv,ndir,nshr,tempNew(km),Cauchy,enerInternNew(km))

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
  stateNew(km,:) = storeStateVars(sv, nstatev, m)

  ! -------------------------------------------------------------------- !
  !    Store the internal energy                                         !
  ! -------------------------------------------------------------------- !
  enerInternNew(km) = enerInternNew(km)/density(km)

  End Do master ! End Master Loop

  Call log%debug('End of VUMAT')

  Return
End Subroutine CompDam
