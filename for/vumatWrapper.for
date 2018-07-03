Module vumat_Wrapper_Mod
  ! This module contains a wrapper for VUMAT user material subroutines, allowing them to be called by Abaqus/Standard.

Contains

  Subroutine vumat_Wrapper( &
    ! Standard VUMAT arguments
    stress,statev,DDSDDE,SSE,SPD,SCD,RPL,DDSDDT,DRPLDE,DRPLDT, &
    STRAN,DSTRAN,TIME,dtime,TEMP,DTEMP,PREDEF,DPRED,cmname,    &
    ndi,nshr,NTENS,nstatv,props,nprops,COORDS,DROT,PNEWDT,     &
    CELENT,DFGRD0,DFGRD1,noel,npt,layer,kspt,JSTEP,KINC,       &
    ! Optional arguments
    vucharlength,nvucharlength)

    Use matrixAlgUtil_Mod

    Integer, parameter :: km=1, lanneal=0, ndir=3

    ! UMAT arguments
    Double Precision, intent(INOUT) :: stress(NTENS)
    Double Precision, intent(INOUT) :: statev(nstatv)
    Double Precision, intent(INOUT) :: RPL
    Double Precision, intent(OUT) :: DDSDDE(NTENS,NTENS), DDSDDT(NTENS)
    Double Precision, intent(OUT) :: DRPLDE(NTENS), DRPLDT
    Double Precision, intent(IN) :: STRAN(NTENS), DSTRAN(NTENS)
    Double Precision, intent(IN) :: TIME(2), dtime
    Double Precision, intent(IN) :: TEMP, DTEMP
    Double Precision, intent(IN) :: PREDEF(1), DPRED(1)
    Double Precision, intent(IN) :: props(nprops)
    Double Precision, intent(IN) :: COORDS(3)
    Double Precision, intent(IN) :: DROT(3,3)
    Double Precision, intent(IN) :: DFGRD0(3,3), DFGRD1(3,3)
    Double Precision, intent(INOUT) :: SSE, SPD, SCD
    Double Precision, intent(INOUT) :: PNEWDT
    Double Precision, intent(IN) :: CELENT
    Integer, intent(IN) :: JSTEP(4),ndi, nshr, NTENS, nstatv, nprops, noel, npt, layer, kspt, KINC
    Character(len=80) :: cmname

    ! Optional arguments
    Double Precision, optional :: vucharlength(*)  ! element lengths corresponding to a subroutine vucharlength() call
    Integer, intent(IN), optional :: nvucharlength ! size of vucharlength(:)

    ! Local variables
    Double Precision, allocatable :: charLength_v(:,:)
    Double Precision :: stepTime, totalTime
    Double Precision :: eig_values(3), eig_vectors(3,3)
    Double Precision :: stretch(3,3)                            ! U from polar decomposition of DFGRD1
    Double Precision :: direct(3,3)                             ! Rotation from global basis to basis at start of increment
    Double Precision :: R_13(3,3)                               ! Rotation from global basis to basis at end of increment
    Double Precision :: R(3,3)                                  ! Rotation from polar decomposition
    Double Precision :: DFGRD0rot(3,3), DFGRD1rot(3,3)          ! DFGRD in the global basis
    Double Precision :: deltaF(3,3), deltaR(3,3), deltaU(3,3)   ! Increments in F, R, U
    Double Precision :: CauchyVUMAT(3,3)                        ! Cauchy stress as returned from the vumat
    Double Precision :: CauchyRef(3,3)                          ! Cauchy stress in the global (reference) basis
    Double Precision :: CauchyABQ(3,3)                          ! Cauchy stress in the basis at the end of the increment (returned to Abaqus)
    Double Precision :: eye(3,3)
    Double Precision :: mdet_direct

    ! VUMAT arguments
    Double Precision density(km),coordMp(km,3), &
      strainInc(km,3+nshr),relSpinInc(km,nshr),tempOld(km),tempNew(km),stretchOld(km,3+nshr), &
      defgradOld(km,3+2*nshr),stressOld(km,3+nshr),stateOld(km,nstatv),charLength(km), &
      enerInternOld(km),enerInelasOld(km),stretchNew(km,3+nshr),defgradNew(km,3+2*nshr), &
      stressNew(km,3+nshr),stateNew(km,nstatv),enerInternNew(km),enerInelasNew(km), &
      fieldOld(km,1),fieldNew(km,1)
    Integer :: nfieldv, jblock(5)

    Interface
      Subroutine VUMAT( &
        jblock,ndir,nshr,nstatv,nfieldv,nprops,lanneal,stepTime,      &
        totalTime,dtime,cmname,coordMp,charLength,props,density,      &
        strainInc,relSpinInc,tempOld,stretchOld,defgradOld,fieldOld,  &
        stressOld,stateOld,enerInternOld,enerInelasOld,tempNew,       &
        stretchNew,defgradNew,fieldNew,                               &
        stressNew,stateNew,enerInternNew,enerInelasNew)

        Integer, intent(IN) :: jblock(*), ndir, nshr, nstatv, nfieldv, nprops, lanneal
        Double Precision, intent(IN) :: stepTime, totalTime, dtime
        Character(len=80) :: cmname
        Double Precision, intent(IN) :: coordMp(jblock(1),*), charLength(*), props(nprops), density(jblock(1)), &
          strainInc(jblock(1),ndir+nshr), relSpinInc(jblock(1),nshr), tempOld(jblock(1)), stretchOld(jblock(1),ndir+nshr), &
          defgradOld(jblock(1),ndir+2*nshr), fieldOld(jblock(1),nfieldv), stressOld(jblock(1),ndir+nshr), &
          stateOld(jblock(1),nstatv), enerInternOld(jblock(1)), enerInelasOld(jblock(1)), tempNew(jblock(1)), &
          stretchNew(jblock(1),ndir+nshr), defgradNew(jblock(1),ndir+2*nshr), fieldNew(jblock(1),nfieldv)
        Double Precision, intent(OUT) :: stressNew(jblock(1),ndir+nshr), stateNew(jblock(1),nstatv), enerInternNew(jblock(1)), &
          enerInelasNew(jblock(1))
      End Subroutine VUMAT
    End Interface

    ! Parameters
    Double Precision, parameter :: zero=0.d0, half=0.5d0, one=1.d0
    Integer, parameter :: i_nblock=1, i_npt=2, i_layer=3, i_kspt=4, i_noel=5

    eye = zero; DO I = 1,3; eye(I,I) = one; end DO

    ! -------------------------------------------------------------------- !
    ! Map UMAT inputs to VUMAT inputs
    ! -------------------------------------------------------------------- !
    jblock(i_nblock) = 1
    jblock(i_npt) = npt
    jblock(i_layer) = layer
    jblock(i_kspt) = kspt
    jblock(i_noel) = noel

    stepTime  = TIME(1)
    totalTime = TIME(2)

    ! Map characteristic element lengths
    If (present(nvucharlength)) Then
      Allocate(charLength_v(km,nvucharlength))
      Do i=1,nvucharlength; charLength_v(km,i) = vucharlength(i); End Do
    Else
      charLength(km) = CELENT
    End If

    coordMp(km,:) = COORDS(:)

    density(km) = one

    ! Map strains
    strainInc(km,1) = DSTRAN(1)
    strainInc(km,2) = DSTRAN(2)
    If (ndi > 2) Then
      strainInc(km,3) = DSTRAN(3)
      strainInc(km,4) = half*DSTRAN(4)
      If (nshr > 1) Then
        strainInc(km,5) = half*DSTRAN(6)
        strainInc(km,6) = half*DSTRAN(5)
      End If
    Else
      strainInc(km,3) = zero
      strainInc(km,4) = half*DSTRAN(3)
    End If

    ! TODO: relSpinInc is not currently supported in the VUMAT wrapper
    relSpinInc(km,:) = zero

    ! Map old stresses
    stressOld(km,1) = stress(1)
    stressOld(km,2) = stress(2)
    If (ndi > 2) Then
      stressOld(km,3) = stress(3)
      stressOld(km,4) = stress(4)
      If (nshr > 1) Then
        stressOld(km,5) = stress(6)
        stressOld(km,6) = stress(5)
      End If
    Else
      stressOld(km,3) = zero
      stressOld(km,4) = stress(3)
    End If

    ! Map old state variables
    stateOld(km,:) = statev(1:nstatv-9)

    ! Map energies
    enerInternOld(km) = SSE
    enerInelasOld(km) = SPD

    ! Map temperatures
    tempOld(km) = TEMP
    tempNew(km) = TEMP + DTEMP

    ! Load direct (from usdfld)
    direct = Vec2Matrix(statev(nstatv-8:nstatv))
    mdet_direct = MDet(direct)
    If (stepTime == zero .AND. ABS(mdet_direct-one)>1.d-7) Then
      direct = eye
    Else If (stepTime > zero .AND. ABS(mdet_direct-one)>1.d-7) Then
      print *, 'ERROR: det(direct) /= 1, improper rotation matrix'
      print *, 'det(direct) = ', mdet_direct
      print *, 'direct:'
      print *, transpose(direct)
      print *, 'Make sure the material card includes a "*User defined field" card'
      Call XIT
    End If

    ! Get the increment in deformation to subsequently get the increment in rotation
    deltaF = MATMUL(DFGRD1, MInverse(DFGRD0))
    ! Calculate deltaR from polar decomposition of deltaF
    Call PolarDecomp(deltaF, deltaR, deltaU)

    ! Rotation from basis 1 (global) to basis 3 (end of increment)
    R_13 = MATMUL(deltaR, direct)

    ! Rotate DFGRD to original stationary coordinate basis
    DFGRD0rot = MATMUL(MATMUL(direct, DFGRD0), TRANSPOSE(direct))
    DFGRD1rot = MATMUL(MATMUL(R_13, DFGRD1), TRANSPOSE(R_13))

    ! Map deformation gradient tensors
    defgradOld(km,:) = Matrix2Vec(DFGRD0rot, nshr, .FALSE.)
    defgradNew(km,:) = Matrix2Vec(DFGRD1rot, nshr, .FALSE.)

    ! Map stretch tensors
    Call PolarDecomp(DFGRD0rot, R, stretch)
    stretchOld(km,:) = Matrix2Vec(stretch, nshr)
    Call PolarDecomp(DFGRD1rot, R, stretch)
    stretchNew(km,:) = Matrix2Vec(stretch, nshr)

    ! TODO: Field variables are not currently supported in the VUMAT wrapper
    nfieldv = 1
    fieldOld = zero
    fieldNew = zero

    ! -------------------------------------------------------------------- !
    ! Call the VUMAT to be wrapped.
    ! -------------------------------------------------------------------- !
    If (present(nvucharlength)) Then
      ! VUMAT call with vucharlength functionality
      Call VUMAT( &
        jblock,ndir,nshr,nstatv,nfieldv,nprops,lanneal,stepTime,      &
        totalTime,dtime,cmname,coordMp,charLength_v,props,density,    &
        strainInc,relSpinInc,tempOld,stretchOld,defgradOld,fieldOld,  &
        stressOld,stateOld,enerInternOld,enerInelasOld,tempNew,       &
        stretchNew,defgradNew,fieldNew,                               &
        stressNew,stateNew,enerInternNew,enerInelasNew)
    Else
      ! VUMAT call without vucharlength functionality
      Call VUMAT( &
        jblock,ndir,nshr,nstatv,nfieldv,nprops,lanneal,stepTime,      &
        totalTime,dtime,cmname,coordMp,charLength,props,density,      &
        strainInc,relSpinInc,tempOld,stretchOld,defgradOld,fieldOld,  &
        stressOld,stateOld,enerInternOld,enerInelasOld,tempNew,       &
        stretchNew,defgradNew,fieldNew,                               &
        stressNew,stateNew,enerInternNew,enerInelasNew)
    End If

    ! CompDam specific
    ! Rotate stress to the correct reference frame
    ! Stress as a 3,3 array as returned by CompDam
    CauchyVUMAT = Vec2Matrix(stressNew(km,:))
    ! Same R used in VUMAT
    R = MATMUL(DFGRD1rot, MInverse(stretch))
    ! Stress in reference frame
    CauchyRef = MATMUL(R, MATMUL(CauchyVUMAT, TRANSPOSE(R)))
    ! Stress in abaqus current frame
    CauchyABQ = MATMUL(TRANSPOSE(R_13), MATMUL(CauchyRef, (R_13)))

    ! -------------------------------------------------------------------- !
    ! Map VUMAT outputs to UMAT inputs
    ! -------------------------------------------------------------------- !
    ! Map updated stresses
    stress(1) = CauchyABQ(1,1)
    stress(2) = CauchyABQ(2,2)
    If (ndi > 2) Then
      stress(3) = CauchyABQ(3,3)
      stress(4) = CauchyABQ(1,2)
      If (nshr > 1) Then
        stress(5) = CauchyABQ(3,1)
        stress(6) = CauchyABQ(2,3)
      End If
    Else
      stress(3) = CauchyABQ(1,2)
    End If

    ! Map new state variables
    statev(1:24) = stateNew(km,:)

    ! Map new energies
    SSE = enerInternNew(km)
    SPD = enerInelasNew(km)
    SCD = zero

    ! TODO: Calculate Jacobians via linear perturbation. A general approach requires identifying which displacement
    !  measures the VUMAT subroutine is sensitive to, i.e., strainInc, defgradNew, or stretchNew. For now, the required
    !  Jacobians must be defined in the UMAT subroutine manually.
    DDSDDE = zero
    DDSDDT = zero
    DRPLDE = zero
    DRPLDT = zero

    ! -------------------------------------------------------------------- !
    Return
  End Subroutine vumat_Wrapper

End Module vumat_Wrapper_Mod


Subroutine USDFLD(field,statev,pnewdt,direct,t,celent, &
  time,dtime,cmname,orname,nfield,nstatv,noel,npt,layer, &
  kspt,kstep,kinc,ndi,nshr,coord,jmac,jmatyp,matlayo,laccfla)

  Use matrixAlgUtil_Mod

  implicit Double Precision (a-h, o-z)

  Character*80 cmname,orname
  Character*3  flgray(15)
  Dimension field(nfield),statev(nstatv),direct(3,3),t(3,3),time(2)
  Dimension array(15),jarray(15),jmac(*),jmatyp(*),coord(*)

  ! save the direct as state variable to pass to umat
  statev(nstatv-8:nstatv) = Matrix2Vec(direct, nshr, .FALSE.)

  Return
End Subroutine USDFLD