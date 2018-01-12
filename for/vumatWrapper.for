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
    Double Precision, optional :: vucharlength(*)  ! element lengths corrresponding to a subroutine vucharlength() call
    Integer, intent(IN), optional :: nvucharlength ! size of vucharlength(:)

    ! Local variables
    Double Precision, allocatable :: charLength_v(:,:)
    Double Precision stepTime, totalTime
    Double Precision eig_values(3), eig_vectors(3,3), stretch(3,3)

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
    stateOld(km,:) = statev(:)

    ! Map energies
    enerInternOld(km) = SSE
    enerInelasOld(km) = SPD

    ! Map temperatures
    tempOld(km) = TEMP
    tempNew(km) = TEMP + DTEMP

    ! Map deformation gradient tensors
    defgradOld(km,1) = DFGRD0(1,1)
    defgradOld(km,2) = DFGRD0(2,2)
    defgradOld(km,3) = DFGRD0(3,3)
    defgradOld(km,4) = DFGRD0(1,2)

    defgradNew(km,1) = DFGRD1(1,1)
    defgradNew(km,2) = DFGRD1(2,2)
    defgradNew(km,3) = DFGRD1(3,3)
    defgradNew(km,4) = DFGRD1(1,2)
    If (nshr > 1) Then
      defgradOld(km,5) = DFGRD0(2,3)
      defgradOld(km,6) = DFGRD0(3,1)
      defgradOld(km,7) = DFGRD0(2,1)
      defgradOld(km,8) = DFGRD0(3,2)
      defgradOld(km,9) = DFGRD0(1,3)

      defgradNew(km,5) = DFGRD1(2,3)
      defgradNew(km,6) = DFGRD1(3,1)
      defgradNew(km,7) = DFGRD1(2,1)
      defgradNew(km,8) = DFGRD1(3,2)
      defgradNew(km,9) = DFGRD1(1,3)
    Else
      defgradOld(km,5) = DFGRD0(2,1)

      defgradNew(km,5) = DFGRD1(2,1)
    End If

    ! Map stretch tensors
    stretch = MATMUL(TRANSPOSE(DFGRD0), DFGRD0)
    Call SPRIND((/ stretch(1,1), stretch(2,2), stretch(3,3), stretch(1,2), stretch(1,3), stretch(2,3) /), eig_values, eig_vectors, 1, 3, 3)
    eig_values = SQRT(eig_values)
    stretch = MATMUL(TRANSPOSE(eig_vectors), MATMUL(reshape((/eig_values(1), zero, zero, zero, eig_values(2), zero, zero, zero, eig_values(3)/), (/3,3/)), eig_vectors))
    stretchOld(km,:) = Matrix2Vec(stretch, nshr)

    stretch = MATMUL(TRANSPOSE(DFGRD1), DFGRD1)
    Call SPRIND((/ stretch(1,1), stretch(2,2), stretch(3,3), stretch(1,2), stretch(1,3), stretch(2,3) /), eig_values, eig_vectors, 1, 3, 3)
    eig_values = SQRT(eig_values)
    stretch = MATMUL(TRANSPOSE(eig_vectors), MATMUL(reshape((/eig_values(1), zero, zero, zero, eig_values(2), zero, zero, zero, eig_values(3)/), (/3,3/)), eig_vectors))
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

    ! -------------------------------------------------------------------- !
    ! Map VUMAT outputs to UMAT inputs
    ! -------------------------------------------------------------------- !
    ! Map updated stresses
    stress(1) = stressNew(km,1)
    stress(2) = stressNew(km,2)
    If (ndi > 2) Then
      stress(3) = stressNew(km,3)
      stress(4) = stressNew(km,4)
      If (nshr > 1) Then
        stress(5) = stressNew(km,6)
        stress(6) = stressNew(km,5)
      End If
    Else
      stress(3) = stressNew(km,4)
    End If

    ! Map new state variables
    statev(:) = stateNew(km,:)

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
