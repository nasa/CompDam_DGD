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
    Double Precision :: Plas13        ! Plastic shear strain
    Double Precision :: Inel13        ! Inelastic shear strain
    Double Precision :: slide(2)      ! Slip on a cohesive crack, in the fiber and transverse directions
    Double Precision :: rfC           ! Fiber compression damage threshold
    Double Precision :: d1T           ! Fiber tension damage
    Double Precision :: d1C           ! Fiber compression damage
    Double Precision :: phi0
    Double Precision :: gamma
    Double Precision :: Sr            ! Schapery micro-damage state variable, reduced
    Double Precision :: direct(9)
    Double Precision :: Ep_schaefer(6) ! Total plastic strain for schaefer theory
    Double Precision :: fp ! yield function value

    ! Stored for debugging only
    Double Precision :: d_eps12, d_eps13

    ! Temporary values
    Double Precision :: Plas12_temp
    Double Precision :: Inel12_temp
    Double Precision :: d_eps12_temp
    Double Precision :: Plas13_temp
    Double Precision :: Inel13_temp
    Double Precision :: d_eps13_temp
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

    If (m%shearNonlinearity12) Then
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

    IF (m%schaefer) THEN
      sv%Ep_schaefer(1) = stateOld(27)
      sv%Ep_schaefer(2) = stateOld(28)
      sv%Ep_schaefer(3) = stateOld(29)
      sv%Ep_schaefer(4) = stateOld(30)
      sv%Ep_schaefer(5) = stateOld(31)
      sv%Ep_schaefer(6) = stateOld(32)
      sv%fp = stateOld(33)
    ELSE
      sv%Ep_schaefer(1) = zero
      sv%Ep_schaefer(2) = zero
      sv%Ep_schaefer(3) = zero
      sv%Ep_schaefer(4) = zero
      sv%Ep_schaefer(5) = zero
      sv%Ep_schaefer(6) = zero
      sv%fp = zero
    END IF     
    sv%rfT = MAX(one, stateOld(14))
    sv%slide(1) = stateOld(15)
    sv%slide(2) = stateOld(16)
    sv%rfC =stateOld(17)
    IF (m%fiberCompDamBL) Then
      sv%rfC = MAX(one, sv%rfC)
    End If
    sv%d1T = MAX(zero, stateOld(18))
    sv%d1C = zero

    If (m%fiberCompDamBL) Then
      sv%d1C = MAX(zero, stateOld(19))
    End If

    If (m%shearNonlinearity13) Then
      sv%Plas13 = stateOld(20)
      sv%Inel13 = stateOld(21)
    Else
      sv%Plas13 = zero
      sv%Inel13 = zero
    End If

    If (m%fiberCompDamFKT) Then
      sv%d1C = MAX(zero, stateOld(19))
      sv%phi0 = stateOld(22)  ! State variable 20
      sv%gamma = zero  ! State variable 21
      sv%Fm1 = stateOld(24)
      sv%Fm2 = stateOld(25)
      sv%Fm3 = stateOld(26)
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

    If (m%shearNonlinearity12) Then
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

    If (m%fiberCompDamBL .OR. m%fiberCompDamFKT) Then
      stateNew(19) = sv%d1C
    Else If (nstatev >= 19) Then
      stateNew(19) = zero
    End If

    If (m%schaefer) Then
      stateNew(27) = sv%Ep_schaefer(1)
      stateNew(28) = sv%Ep_schaefer(2)
      stateNew(29) = sv%Ep_schaefer(3)
      stateNew(30) = sv%Ep_schaefer(4)
      stateNew(31) = sv%Ep_schaefer(5)
      stateNew(32) = sv%Ep_schaefer(6)
      stateNew(33) = sv%fp
    ELSE IF (nstatev == 33) THEN
      stateNew(27) = zero
      stateNew(28) = zero
      stateNew(29) = zero
      stateNew(30) = zero
      stateNew(31) = zero
      stateNew(32) = zero
      stateNew(33) = zero
    END IF

    If (m%shearNonlinearity13) Then
      stateNew(20) = sv%Plas13
      stateNew(21) = sv%Inel13
    Else If (nstatev > 21) Then
      stateNew(20) = zero
      stateNew(21) = zero
    End If

    If (m%fiberCompDamFKT) Then
      stateNew(22) = sv%phi0
      stateNew(23) = sv%gamma
      stateNew(24) = sv%Fm1
      stateNew(25) = sv%Fm2
      stateNew(26) = sv%Fm3
    Else If (nstatev == 26) Then
      stateNew(22) = zero
      stateNew(23) = zero
      stateNew(24) = zero
      stateNew(25) = zero
      stateNew(26) = zero
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

    If (m%shearNonlinearity12) Then
      sv%Plas12_temp = sv%Plas12
      sv%Inel12_temp = sv%Inel12
      sv%d_eps12_temp = sv%d_eps12
    End If
    If (m%shearNonlinearity13) Then
      sv%Plas13_temp = sv%Plas13
      sv%Inel13_temp = sv%Inel13
      sv%d_eps13_temp = sv%d_eps13
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

    If (m%shearNonlinearity12) Then
      sv%Plas12 = sv%Plas12_temp
      sv%Inel12 = sv%Inel12_temp
      sv%d_eps12 = sv%d_eps12_temp
    End If
    If (m%shearNonlinearity13) Then
      sv%Plas13 = sv%Plas13_temp
      sv%Inel13 = sv%Inel13_temp
      sv%d_eps13 = sv%d_eps13_temp
    End If

    If (m%schapery) Then
      sv%Sr = sv%Sr_temp
    End If

    Return
  End Subroutine finalizeTemp
End Module stateVar_Mod
