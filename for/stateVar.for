Module stateVar_Mod
  ! Module for loading and validating internal state variables

  ! Stores the state variables
  Type stateVars

    Integer :: nstatev                ! Number of state variables

    ! Always stored
    Double Precision :: d2            ! Cohesive surface damage variable
    Double Precision :: Fb1           ! Fb = deformation gradient of the bulk material; 1,2,3 are the columns of Fb
    Double Precision :: Fb2
    Double Precision :: Fb3
    Double Precision :: Fm1
    Double Precision :: Fm2
    Double Precision :: Fm3
    Double Precision :: B             ! Mode mixity
    Double Precision :: Lc(3)         ! Characteristic element lengths
    Double Precision :: rfT           ! Fiber tension damage threshold
    Double Precision :: FIm           ! Failure index for matrix
    Integer :: alpha                  ! The cohesive surface normal [degrees, integer]. Only modified in this subroutine if matrix failure criteria is satisfied.
    Integer :: STATUS                 ! Element deletion flag (0 := delete element)
    Double Precision :: Plas12        ! Plastic shear strain
    Double Precision :: Inel12        ! Inelastic shear strain
    Double Precision :: slide(2)      ! Slip on a cohesive crack, in the fiber and transverse directions
    Double Precision :: rfC           ! Fiber compression damage threshold
    Double Precision :: d1T           ! Fiber tension damage
    Double Precision :: d1C           ! Fiber compression damage
    Double Precision :: phi0
    Double Precision :: gamma
    Double Precision :: Sr            ! Schapery micro-damage state variable, reduced
    Double Precision :: direct(9)

    ! Stored for debugging only
    Double Precision :: d_eps12

    ! Temporary values
    Double Precision :: Plas12_temp
    Double Precision :: Inel12_temp
    Double Precision :: d_eps12_temp
    Double Precision :: Sr_temp

  End Type stateVars

Contains

  Pure Function loadStateVars(nstatev, stateOld, m) result(sv)
    ! Loads state variables into named fields

    Use matProp_Mod

    ! Arguments
    Integer, intent(IN) :: nstatev
    Double Precision, intent(IN) :: stateOld(nstatev)
    Type(matProps), intent(IN) :: m

    ! Output
    Type(stateVars) :: sv

    ! Parameters
    Double Precision, parameter :: zero=0.d0, one=1.d0
    ! -------------------------------------------------------------------- !

    sv%nstatev = nstatev

    ! Global variable (not returned to abaqus)
    sv%d_eps12 = zero

    sv%d2 = MAX(zero, stateOld(1))
    sv%Fb1 = stateOld(2)
    sv%Fb2 = stateOld(3)
    sv%Fb3 = stateOld(4)
    sv%B   = stateOld(5)
    sv%Lc(1) = stateOld(6)
    sv%Lc(2) = stateOld(7)
    sv%Lc(3) = stateOld(8)
    sv%FIm = zero   ! State variable 9
    sv%alpha  = stateOld(10)
    sv%STATUS = stateOld(11)

    If (m%shearNonlinearity) Then
      sv%Plas12 = stateOld(12)
      sv%Inel12 = stateOld(13)
      sv%Sr = one
    Else If (m%schapery) Then
      sv%Plas12 = zero
      sv%Inel12 = zero
      sv%Sr = stateOld(12)
    Else
      sv%Plas12 = zero
      sv%Inel12 = zero
      sv%Sr = one
    End If

    sv%rfT = MAX(one, stateOld(14))
    sv%slide(1) = stateOld(15)
    sv%slide(2) = stateOld(16)
    sv%rfC =stateOld(17)
    IF (m%fiberCompDamBL) Then
      sv%rfC = MAX(one, sv%rfC)
    End If
    sv%d1T = MAX(zero, stateOld(18))
    sv%d1C = zero

    ! These state variables are required only for fiber compression
    If (m%fiberCompDamBL .OR. m%fiberCompDamFKT) Then
      sv%d1C = MAX(zero, stateOld(19))
      ! DGD-based fiber compression
      If (m%fiberCompDamFKT) Then
        sv%phi0 = stateOld(20)  ! State variable 20
        sv%gamma = zero  ! State variable 21
        sv%Fm1 = stateOld(22)
        sv%Fm2 = stateOld(23)
        sv%Fm3 = stateOld(24)
      End If
    End If

    Return
  End Function loadStateVars


  Pure Function storeStateVars(sv, nstatev, m) result(stateNew)
    ! Returns an array of state variables (calling code should be stateNew = store(nstatev))

    Use matProp_Mod

    ! Arguments
    Type(stateVars), intent(IN) :: sv
    Integer, intent(IN) :: nstatev
    Type(matProps), intent(IN) :: m

    ! Output
    Double Precision :: stateNew(nstatev)

    ! Parameters
    Double Precision, parameter :: zero=0.d0
    ! -------------------------------------------------------------------- !

    stateNew(1) = sv%d2
    stateNew(2) = sv%Fb1
    stateNew(3) = sv%Fb2
    stateNew(4) = sv%Fb3
    stateNew(5) = sv%B
    stateNew(6) = sv%Lc(1)
    stateNew(7) = sv%Lc(2)
    stateNew(8) = sv%Lc(3)
    stateNew(9) = sv%FIm
    stateNew(10) = sv%alpha
    stateNew(11) = sv%STATUS

    If (m%shearNonlinearity) Then
      stateNew(12) = sv%Plas12
      stateNew(13) = sv%Inel12
    Else If (m%schapery) Then
      stateNew(12) = sv%Sr
      stateNew(13) = zero
    Else
      stateNew(12) = zero
      stateNew(13) = zero
    End If

    stateNew(14) = sv%rfT
    stateNew(15) = sv%slide(1)
    stateNew(16) = sv%slide(2)
    stateNew(17) = sv%rfC
    stateNew(18) = sv%d1T

    ! Fiber compression enabled
    If (m%fiberCompDamBL .OR. m%fiberCompDamFKT) Then
      stateNew(19) = sv%d1C
      ! DGD based fiber compression enabled
      If (m%fiberCompDamFKT) Then
        stateNew(20) = sv%phi0
        stateNew(21) = sv%gamma
        stateNew(22) = sv%Fm1
        stateNew(23) = sv%Fm2
        stateNew(24) = sv%Fm3
      End If
    Else
      stateNew(19) = zero
    End If

    Return
  End Function storeStateVars


  Subroutine initializeTemp(sv, m)
    ! This function sets all temporary state variables equal to the actual state variable value

    Use matProp_Mod

    ! Arguments
    Type(stateVars), intent(INOUT) :: sv
    Type(matProps), intent(IN) :: m
    ! -------------------------------------------------------------------- !

    If (m%shearNonlinearity) Then
      sv%Plas12_temp = sv%Plas12
      sv%Inel12_temp = sv%Inel12
      sv%d_eps12_temp = sv%d_eps12
    End If

    If (m%schapery) Then
      sv%Sr_temp = sv%Sr
    End If

    Return
  End Subroutine initializeTemp

  Subroutine finalizeTemp(sv, m)
    ! This function updates actual state variables to the value of corresponding temporary state variables
    ! This function is intended to be called once iterations for which temporary state variables were needed
    ! are complete

    Use matProp_Mod

    ! Arguments
    Type(stateVars), intent(INOUT) :: sv
    Type(matProps), intent(IN) :: m
    ! -------------------------------------------------------------------- !

    If (m%shearNonlinearity) Then
      sv%Plas12 = sv%Plas12_temp
      sv%Inel12 = sv%Inel12_temp
      sv%d_eps12 = sv%d_eps12_temp
    End If

    If (m%schapery) Then
      sv%Sr = sv%Sr_temp
    End If

    Return
  End Subroutine finalizeTemp
End Module stateVar_Mod
