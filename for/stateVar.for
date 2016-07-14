
Module stateVar_Mod

  ! Stores the state variables
  Type stateVars

    ! Always stored
    Double Precision :: d2            ! Cohesive surface damage variable
    Double Precision :: Fb1           ! Fb = deformation gradient of the bulk material; 1,2,3 are the columns of Fb
    Double Precision :: Fb2
    Double Precision :: Fb3
    Double Precision :: B             ! Mode mixity
    Double Precision :: Lc1
    Double Precision :: rfT           ! Fiber damage threshold
    Double Precision :: d1            ! Fiber damage variable (d1 is used for stiffness reduction; d1T and d1C are tracked for reference)
    Double Precision :: FImT          ! Failure index for matrix tension
    Integer :: alpha                  ! The cohesive surface normal [degrees, integer]. Only modified in this subroutine if matrix failure criteria is satisfied.
    Integer :: STATUS                 ! Element deletion
    Double Precision :: Plas12        ! Shear strains
    Double Precision :: Inel12
    Double Precision :: mCompInit     ! cohesive normal displacement at initiation (zero if crack is in tension)
    Double Precision :: slide(2)
    Double Precision :: rfC
    Double Precision :: d1T
    Double Precision :: d1C

    ! Stored for debugging only
    Double Precision :: d_eps12

  End Type stateVars

Contains

  Pure Function loadStateVars(nstatev, stateOld) result(sv)
    ! Loads state variables into named fields

    ! Arguments
    Integer, intent(IN) :: nstatev
    Double Precision, intent(IN) :: stateOld(nstatev)

    ! Output
    Type(stateVars) :: sv

    ! Parameters
    Double Precision, parameter :: zero=0.d0, one=1.d0
    ! -------------------------------------------------------------------- !

    ! Global variable (not returned to abaqus)
    sv%d_eps12 = zero

    sv%d2 = MAX(zero, stateOld(1))
    sv%Fb1 = stateOld(2)
    sv%Fb2 = stateOld(3)
    sv%Fb3 = stateOld(4)
    sv%B   = stateOld(5)
    sv%Lc1 = stateOld(6)
    sv%rfT = MAX(one, stateOld(7))
    sv%FImT = zero   ! State variable 9
    sv%alpha  = stateOld(10)
    sv%STATUS = stateOld(11)
    sv%Plas12 = stateOld(12)
    sv%Inel12 = stateOld(13)
    sv%mCompInit = stateOld(14)  ! the crack-normal cohesive displacement at initiation (zero if positive)
    sv%slide(1) = stateOld(15)
    sv%slide(2) = stateOld(16)
    sv%rfC = MAX(one, stateOld(17))
    sv%d1T = MAX(zero, stateOld(18))
    sv%d1C = zero
    If (nstatev .EQ. 18) Return ! Fiber compression disabled

    ! State variables required only for fiber compression
    sv%d1C = MAX(zero, stateOld(19))

    Return
  End Function loadStateVars


  Pure Function storeStateVars(sv, nstatev) result(stateNew)
    ! Returns an array of state variables (calling code should be stateNew = store(nstatev))

    ! Arguments
    Type(stateVars), intent(IN) :: sv
    Integer, intent(IN) :: nstatev

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
    stateNew(6) = sv%Lc1
    stateNew(7) = sv%rfT
    stateNew(8) = sv%d1
    stateNew(9) = sv%FImT
    stateNew(10) = sv%alpha
    stateNew(11) = sv%STATUS
    stateNew(12) = sv%Plas12
    stateNew(13) = sv%Inel12
    stateNew(14) = sv%mCompInit
    stateNew(15) = sv%slide(1)
    stateNew(16) = sv%slide(2)
    stateNew(17) = sv%rfC
    stateNew(18) = sv%d1T
    If (nstatev .EQ. 18) Return ! Fiber compression disabled
    stateNew(19) = sv%d1C

    Return
  End Function storeStateVars


End Module
