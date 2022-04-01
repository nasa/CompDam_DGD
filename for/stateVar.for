Module stateVar_Mod
  ! Module for loading and validating internal state variables

  ! Maximum possible state variables defined
  Integer, parameter :: nstatev_max = 33

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
    Double Precision :: phi0_12
    Double Precision :: gamma_12
    Double Precision :: phi0_13
    Double Precision :: gamma_13
    Double Precision :: fifbc         ! Failure index for fiber break in compression
    Double Precision :: Sr            ! Schapery micro-damage state variable, reduced
    Double Precision :: direct(9)
    Double Precision :: Ep_schaefer(6) ! Total plastic strain for schaefer theory
    Double Precision :: fp ! yield function value
    Double Precision :: enerPlasOld   ! Increment of dissipated plastic energy
    Double Precision :: enerPlasNew   ! Increment of dissipated plastic energy
    Double Precision :: enerFracOld   ! Increment of dissipated fracture energy
    Double Precision :: enerFracNew   ! Increment of dissipated fracture energy

    ! Calculated once at the start of the analysis, not passed to abaqus
    ! Since this value depends on phi0, it needs to be a state variable and not global calculated parameter
    Double Precision :: Inel12c       ! Critical plastic strain for fiber failure
    Double Precision :: Inel13c       ! Critical plastic strain for fiber failure

    ! Stored for debugging only
    Integer :: debugpy_count
    Double Precision :: old(nstatev_max)

    ! Temporary values
    Double Precision :: Plas12_temp
    Double Precision :: Inel12_temp
    Double Precision :: Plas13_temp
    Double Precision :: Inel13_temp
    Double Precision :: rfC_temp
    Double Precision :: enerPlas_temp
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
    sv%old(1:nstatev) = stateOld

    ! Global variable (not returned to abaqus)
    sv%Inel12c = Huge(zero)  ! Initialize to a large positive number so it is not reached (turns off fiber failure)
    sv%Inel13c = Huge(zero)  ! Initialize to a large positive number so it is not reached (turns off fiber failure)

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
    sv%rfT = stateOld(14)  ! MAX(one, stateOld(14))
    sv%slide(1) = stateOld(15)
    sv%slide(2) = stateOld(16)
    sv%rfC =stateOld(17)
    ! IF (m%fiberCompDamBL) Then
    !   sv%rfC = MAX(one, sv%rfC)
    ! End If
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

    If (m%fiberCompDamFKT12) Then
      sv%d1C = MAX(zero, stateOld(19))
      sv%phi0_12 = stateOld(22)
      sv%gamma_12 = stateOld(23)
      If (nstatev >= 26) Then
        sv%fifbc = stateOld(26)
      End If
    Else
      sv%phi0_12 = zero
      sv%gamma_12 = zero
    End If

    If (m%fiberCompDamFKT13) Then
      sv%d1C = MAX(zero, stateOld(19))
      sv%phi0_13 = stateOld(24)
      sv%gamma_13 = stateOld(25)
      If (nstatev >= 26) Then
        sv%fifbc = stateOld(26)
      End If
    Else
      sv%phi0_13 = zero
      sv%gamma_13 = zero
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

    If (m%fiberCompDamBL .OR. m%fiberCompDamFKT12 .OR. m%fiberCompDamFKT13) Then
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

    If (m%fiberCompDamFKT12) Then
      stateNew(22) = sv%phi0_12
      stateNew(23) = sv%gamma_12
      If (nstatev >= 26) Then
        stateNew(26) = sv%fifbc
      End If
    Else If (nstatev >= 23) Then
      stateNew(22) = zero
      stateNew(23) = zero
    End If

    If (m%fiberCompDamFKT13) Then
      stateNew(24) = sv%phi0_13
      stateNew(25) = sv%gamma_13
      If (nstatev >= 26) Then
        stateNew(26) = sv%fifbc
      End If
    Else If (nstatev >= 25) Then
      stateNew(24) = zero
      stateNew(25) = zero
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
    End If
    If (m%shearNonlinearity13) Then
      sv%Plas13_temp = sv%Plas13
      sv%Inel13_temp = sv%Inel13
    End If

    If (m%fiberCompDamFKT12 .OR. m%fiberCompDamFKT13) Then
      sv%rfC_temp = sv%rfC
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
    End If
    If (m%shearNonlinearity13) Then
      sv%Plas13 = sv%Plas13_temp
      sv%Inel13 = sv%Inel13_temp
    End If

    If (m%fiberCompDamFKT12 .OR. m%fiberCompDamFKT13) Then
      sv%rfC = sv%rfC_temp
    End If

    If (m%schapery) Then
      sv%Sr = sv%Sr_temp
    End If

    Return
  End Subroutine finalizeTemp


  Subroutine writeStateVariablesToFile(fileUnit, sv, m)
    ! Writes provided state variables to a file as a python dictionary
    ! Assumes that file opening and closing is handled elsewhere

    Use matProp_Mod

    ! Arguments
    Integer, intent(IN) :: fileUnit
    Type(stateVars), intent(IN) :: sv
    Type(matProps), intent(IN) :: m

    ! Locals
    Character(len=32) :: nameValueFmt
    Double Precision, parameter :: zero=0.d0
    ! -------------------------------------------------------------------- !

    ! Defines the format for writing the floating point numbers
    nameValueFmt = "(A,E21.15E2,A)"

    ! Write the current state variables
    write(fileUnit, "(A)") 'sv = ['
    write(fileUnit, nameValueFmt) '    ', sv%d2, ',  # d2'
    write(fileUnit, nameValueFmt) '    ', sv%Fb1, ',  # Fb1'
    write(fileUnit, nameValueFmt) '    ', sv%Fb2, ',  # Fb2'
    write(fileUnit, nameValueFmt) '    ', sv%Fb3, ',  # Fb3'
    write(fileUnit, nameValueFmt) '    ', sv%B, ',  # B'
    write(fileUnit, nameValueFmt) '    ', sv%Lc(1), ',  # Lc1'
    write(fileUnit, nameValueFmt) '    ', sv%Lc(2), ',  # Lc2'
    write(fileUnit, nameValueFmt) '    ', sv%Lc(3), ',  # Lc3'
    write(fileUnit, nameValueFmt) '    ', sv%FIm, ',  # FIm'
    write(fileUnit, "(A,I5,A)")   '    ', sv%alpha, ',  # alpha'
    write(fileUnit, "(A,I1,A)")   '    ', sv%STATUS, ',  # STATUS'
    If (m%shearNonlinearity12) Then
      write(fileUnit, nameValueFmt) '    ', sv%Plas12, ',  # Plas12'
      write(fileUnit, nameValueFmt) '    ', sv%Inel12, ',  # Inel12'
    Else If (m%schapery) Then
      write(fileUnit, nameValueFmt) '    ', sv%Sr, ',  # Sr'
      write(fileUnit, nameValueFmt) '    ', zero, ',  # SDV13'
    Else
      write(fileUnit, nameValueFmt) '    ', zero, ',  # SDV12'
      write(fileUnit, nameValueFmt) '    ', zero, ',  # SDV13'
    End If
    write(fileUnit, nameValueFmt) '    ', sv%rfT, ',  # rfT'
    write(fileUnit, nameValueFmt) '    ', sv%slide(1), ',  # slide1'
    write(fileUnit, nameValueFmt) '    ', sv%slide(2), ',  # slide2'
    write(fileUnit, nameValueFmt) '    ', sv%rfC, ',  # rfC'
    write(fileUnit, nameValueFmt) '    ', sv%d1T, ',  # d1T'
    If (m%fiberCompDamBL .OR. m%fiberCompDamFKT12 .OR. m%fiberCompDamFKT13) Then
      write(fileUnit, nameValueFmt) '    ', sv%d1C, ',  # d1C'
    Else If (sv%nstatev >= 19) Then
      write(fileUnit, nameValueFmt) '    ', zero, ',  # SDV19'
    End If
    If (m%shearNonlinearity13) Then
      write(fileUnit, nameValueFmt) '    ', sv%Plas13, ',  # Plas13'
      write(fileUnit, nameValueFmt) '    ', sv%Inel13, ',  # Inel13'
    Else If (sv%nstatev >= 21) Then
      write(fileUnit, nameValueFmt) '    ', zero, ',  # SDV20'
      write(fileUnit, nameValueFmt) '    ', zero, ',  # SDV21'
    End If
    If (m%fiberCompDamFKT12) Then
      write(fileUnit, nameValueFmt) '    ', sv%phi0_12, ',  # phi0_12'
      write(fileUnit, nameValueFmt) '    ', sv%gamma_12, ',  # gamma_12'
      If (nstatev >= 26) Then
        write(fileUnit, nameValueFmt) '    ', sv%fifbc, ',  # fifbc'
      End If
    Else If (sv%nstatev >= 23) Then
      write(fileUnit, nameValueFmt) '    ', zero, ',  # SDV22'
      write(fileUnit, nameValueFmt) '    ', zero, ',  # SDV23'
    End If
    If (m%fiberCompDamFKT13) Then
      write(fileUnit, nameValueFmt) '    ', sv%phi0_12, ',  # phi0_13'
      write(fileUnit, nameValueFmt) '    ', sv%gamma_12, ',  # gamma_13'
      If (sv%nstatev >= 26) Then
        write(fileUnit, nameValueFmt) '    ', sv%fifbc, ',  # fifbc'
      End If
    Else If (sv%nstatev >= 25) Then
      write(fileUnit, nameValueFmt) '    ', zero, ',  # SDV24'
      write(fileUnit, nameValueFmt) '    ', zero, ',  # SDV25'
    End If
    If (sv%nstatev >= 26) Then
      write(fileUnit, nameValueFmt) '    ', zero, ',  # SDV26'
    End If
    If (m%schaefer) Then
      write(fileUnit, nameValueFmt) '    ', sv%Ep_schaefer(1), ',  # Ep1'
      write(fileUnit, nameValueFmt) '    ', sv%Ep_schaefer(2), ',  # Ep2'
      write(fileUnit, nameValueFmt) '    ', sv%Ep_schaefer(3), ',  # Ep3'
      write(fileUnit, nameValueFmt) '    ', sv%Ep_schaefer(4), ',  # Ep4'
      write(fileUnit, nameValueFmt) '    ', sv%Ep_schaefer(5), ',  # Ep5'
      write(fileUnit, nameValueFmt) '    ', sv%Ep_schaefer(6), ',  # Ep6'
      write(fileUnit, nameValueFmt) '    ', sv%fp, '  # fp1'
    Else If (sv%nstatev >= 33) Then
      write(fileUnit, nameValueFmt) '    ', zero, ',  # SDV27'
      write(fileUnit, nameValueFmt) '    ', zero, ',  # SDV28'
      write(fileUnit, nameValueFmt) '    ', zero, ',  # SDV29'
      write(fileUnit, nameValueFmt) '    ', zero, ',  # SDV30'
      write(fileUnit, nameValueFmt) '    ', zero, ',  # SDV31'
      write(fileUnit, nameValueFmt) '    ', zero, ',  # SDV32'
      write(fileUnit, nameValueFmt) '    ', zero, '  #SDV33'
    End If
    write(fileUnit, "(A)") ']'

    ! Write old state variables too
#ifndef PYEXT
    write(fileUnit, "(A)") 'sv_old = ['
    write(fileUnit, nameValueFmt) '    ', sv%old(1), ',  # d2'
    write(fileUnit, nameValueFmt) '    ', sv%old(2), ',  # Fb1'
    write(fileUnit, nameValueFmt) '    ', sv%old(3), ',  # Fb2'
    write(fileUnit, nameValueFmt) '    ', sv%old(4), ',  # Fb3'
    write(fileUnit, nameValueFmt) '    ', sv%old(5), ',  # B'
    write(fileUnit, nameValueFmt) '    ', sv%old(6), ',  # Lc1'
    write(fileUnit, nameValueFmt) '    ', sv%old(7), ',  # Lc2'
    write(fileUnit, nameValueFmt) '    ', sv%old(8), ',  # Lc3'
    write(fileUnit, nameValueFmt) '    ', sv%old(9), ',  # FIm'
    write(fileUnit, "(A,I5,A)")   '    ', INT(sv%old(10)), ',  # alpha'
    write(fileUnit, "(A,I2,A)")   '    ', INT(sv%old(11)), ',  # STATUS'
    If (m%shearNonlinearity12) Then
      write(fileUnit, nameValueFmt) '    ', sv%old(12), ',  # Plas12'
      write(fileUnit, nameValueFmt) '    ', sv%old(13), ',  # Inel12'
    Else If (m%schapery) Then
      write(fileUnit, nameValueFmt) '    ', sv%old(12), ',  # Sr'
      write(fileUnit, nameValueFmt) '    ', zero, ',  # SDV13'
    Else
      write(fileUnit, nameValueFmt) '    ', zero, ',  # SDV12'
      write(fileUnit, nameValueFmt) '    ', zero, ',  # SDV13'
    End If
    write(fileUnit, nameValueFmt) '    ', sv%old(14), ',  # rfT'
    write(fileUnit, nameValueFmt) '    ', sv%old(15), ',  # slide1'
    write(fileUnit, nameValueFmt) '    ', sv%old(16), ',  # slide2'
    write(fileUnit, nameValueFmt) '    ', sv%old(17), ',  # rfC'
    write(fileUnit, nameValueFmt) '    ', sv%old(18), ',  # d1T'
    If (m%fiberCompDamBL .OR. m%fiberCompDamFKT12 .OR. m%fiberCompDamFKT13) Then
      write(fileUnit, nameValueFmt) '    ', sv%old(19), ',  # d1C'
    Else If (sv%nstatev >= 19) Then
      write(fileUnit, nameValueFmt) '    ', zero, ',  # SDV19'
    End If
    If (m%shearNonlinearity13) Then
      write(fileUnit, nameValueFmt) '    ', sv%old(20), ',  # Plas13'
      write(fileUnit, nameValueFmt) '    ', sv%old(21), ',  # Inel13'
    Else If (sv%nstatev >= 21) Then
      write(fileUnit, nameValueFmt) '    ', zero, ',  # SDV20'
      write(fileUnit, nameValueFmt) '    ', zero, ',  # SDV21'
    End If
    If (m%fiberCompDamFKT12) Then
      write(fileUnit, nameValueFmt) '    ', sv%old(22), ',  # phi0_12'
      write(fileUnit, nameValueFmt) '    ', sv%old(23), ',  # gamma_12'
    Else If (sv%nstatev >= 23) Then
      write(fileUnit, nameValueFmt) '    ', zero, ',  # SDV22'
      write(fileUnit, nameValueFmt) '    ', zero, ',  # SDV23'
    End If
    If (m%fiberCompDamFKT13) Then
      write(fileUnit, nameValueFmt) '    ', sv%old(22), ',  # phi0_13'
      write(fileUnit, nameValueFmt) '    ', sv%old(23), ',  # gamma_13'
    Else If (sv%nstatev >= 25) Then
      write(fileUnit, nameValueFmt) '    ', zero, ',  # SDV24'
      write(fileUnit, nameValueFmt) '    ', zero, ',  # SDV25'
    End If
    If (sv%nstatev >= 26) Then
      write(fileUnit, nameValueFmt) '    ', zero, ',  # SDV26'
    End If
    If (m%schaefer) Then
      write(fileUnit, nameValueFmt) '    ', sv%old(27), ',  # Ep1'
      write(fileUnit, nameValueFmt) '    ', sv%old(28), ',  # Ep2'
      write(fileUnit, nameValueFmt) '    ', sv%old(29), ',  # Ep3'
      write(fileUnit, nameValueFmt) '    ', sv%old(30), ',  # Ep4'
      write(fileUnit, nameValueFmt) '    ', sv%old(31), ',  # Ep5'
      write(fileUnit, nameValueFmt) '    ', sv%old(32), ',  # Ep6'
      write(fileUnit, nameValueFmt) '    ', sv%old(33), '  # fp1'
    Else If (sv%nstatev >= 33) Then
      write(fileUnit, nameValueFmt) '    ', zero, ',  # SDV27'
      write(fileUnit, nameValueFmt) '    ', zero, ',  # SDV28'
      write(fileUnit, nameValueFmt) '    ', zero, ',  # SDV29'
      write(fileUnit, nameValueFmt) '    ', zero, ',  # SDV30'
      write(fileUnit, nameValueFmt) '    ', zero, ',  # SDV31'
      write(fileUnit, nameValueFmt) '    ', zero, ',  # SDV32'
      write(fileUnit, nameValueFmt) '    ', zero, '  # SDV33'
    End If
    write(fileUnit, "(A)") ']'
#endif

    Return
  End Subroutine writeStateVariablesToFile
End Module stateVar_Mod
