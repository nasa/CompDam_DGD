Module CDM_fiber_mod
  ! Continuum damage mechanics fiber damage calculations

Contains

  Pure Subroutine FiberTenDmg(eps,ndir,E1,XT,GXT,n,m,Lc,rfT,d1T,d1C,STATUS)
    ! The purpose of this subroutine is to evaluate fiber damage.
    ! This subroutine is for tension-only fiber damage.

    Include 'vaba_param.inc'

    ! Arguments
    Double Precision, intent(IN) :: E1                        ! Stiffness
    Double Precision, intent(IN) :: XT                        ! Strength
    Double Precision, intent(IN) :: GXT                       ! Toughness
    Double Precision, intent(IN) :: n, m                      ! Coefficients for bi-linear softening
    Double Precision, intent(IN) :: eps(ndir,ndir)            ! Strain
    Double Precision, intent(IN) :: Lc
    Double Precision, intent(IN) :: d1C                       ! CDM compressive fiber damage variable
    Integer, intent(IN) :: ndir

    Double Precision, intent(INOUT) :: rfT                    ! Fiber tension damage threshold, state variable (CDM_FIfT)
    Double Precision, intent(INOUT) :: d1T                    ! Fiber direction damage variables for tension
    Integer, intent(OUT) :: STATUS                            ! Element deletion

    ! Locals
    Double Precision :: FIfT                        ! Fiber direction failure index
    Double Precision :: d1MAX                       ! Maximum value for damage variables
    Double Precision :: eps_rfT                     ! strain corresponding to the maximum historical failure criterion
    Double Precision :: eps_0, eps_f_1, eps_f_2     ! strains for damage variable calculations
    Double Precision :: dmg_1, dmg_2                ! individual damage variables for bi-linear softening
    Double Precision :: d1T_temp                    ! Stores the old value of d1T

    Double Precision, Parameter :: zero=0.d0, one=1.d0, two=2.d0
    ! -------------------------------------------------------------------- !

    ! Initialize
    d1MAX = one
    d1T = MAX(d1T, d1C)
    d1T_temp = d1T
    STATUS = 1

    ! Fiber tension failure criterion
    FIfT = E1/XT*eps(1,1)

    ! Update damage threshold
    If (FIfT .GT. rfT) rfT = FIfT

    ! Determine d1T if damage has occurred
    If (rfT .GT. one) Then

      eps_rfT = rfT*XT/E1

      ! Calculate constants
      eps_0 = XT/E1
      eps_f_1 = two*GXT*m/(XT*n*Lc)
      eps_f_2 = two*GXT*(one - m)/(XT*(one - n)*Lc)

      ! Calculate damage in each linear softening law
      dmg_1 = MIN(one, eps_f_1*(eps_rfT - eps_0)/eps_rfT/(eps_f_1 - eps_0))
      dmg_2 = MIN(one, eps_f_2*(eps_rfT - eps_0)/eps_rfT/(eps_f_2 - eps_0))

      ! Find the combined damage variable, d1T
      d1T = n*(dmg_1 - dmg_2) + dmg_2

      ! Prevent healing; If the current damage is less than previous damage, use previous damage
      If (d1T .LT. d1T_temp) d1T = d1T_temp

      ! Delete element when it becomes fully damaged
      If (d1T .GE. d1MAX) Then
        d1T = d1MAX
        STATUS = 0
      End If

    Else
      ! This is only executed before any damage has occurred
      ! Since rfT is that state variable, it acts as the failure index before damage
      rfT = FIfT

    End If

    Return
  End Subroutine FiberTenDmg


  Pure Subroutine FiberCompDmg(eps,ndir,E1,XC,GXC,n,m,Lc,rfT,rfC,d1T,d1C,STATUS)
    ! The purpose of this subroutine is to evaluate fiber damage.
    ! This subroutine is for compression-only fiber damage.

    Include 'vaba_param.inc'

    ! Arguments
    Double Precision, intent(IN) :: E1                        ! Stiffness
    Double Precision, intent(IN) :: XC                        ! Strength
    Double Precision, intent(IN) :: GXC                       ! Toughness
    Double Precision, intent(IN) :: n, m                      ! Coefficients for bi-linear softening
    Double Precision, intent(IN) :: eps(ndir,ndir)            ! Strain
    Double Precision, intent(IN) :: Lc
    Integer, intent(IN) :: ndir

    Double Precision, intent(INOUT) :: rfC                      ! Fiber compression damage threshold, state variable
    Double Precision, intent(INOUT) :: rfT                      ! Fiber tension damage threshold, state variable
    Double Precision, intent(INOUT) :: d1C                      ! Fiber direction damage variables for compression
    Double Precision, intent(INOUT) :: d1T                      ! Fiber direction damage variables for tension
    Integer, intent(OUT) :: STATUS                              ! Element deletion

    ! Locals
    Double Precision :: FIfC                        ! Fiber direction failure index
    Double Precision :: d1MAX                       ! Maximum value for damage variables
    Double Precision :: eps_rfC                     ! strain corresponding to the maximum historical failure criterion
    Double Precision :: eps_0, eps_f_1, eps_f_2     ! strains for damage variable calculations
    Double Precision :: dmg_1, dmg_2                ! individual damage variables for bi-linear softening
    Double Precision :: d1C_temp                    ! Stores the old value of d1C

    Double Precision, Parameter :: zero=0.d0, one=1.d0, two=2.d0
    ! -------------------------------------------------------------------- !

    ! Initialize
    d1MAX = one
    d1C_temp = d1C
    STATUS = 1

    ! Fiber compression failure criterion
    FIfC = -E1/XC*eps(1,1)

    ! Update damage threshold
    If (FIfC .GT. rfC) Then
      rfC = FIfC

      ! For load reversal
      If (rfC .GT. rfT) rfT = rfC
    End If

    ! Determine d1C if damage has occurred
    If (rfC .GT. one) Then

      eps_rfC = rfC*XC/E1

      ! Calculate constants
      eps_0 = XC/E1
      eps_f_1 = two*GXC*m/(XC*n*Lc)
      eps_f_2 = two*GXC*(one - m)/(XC*(one - n)*Lc)

      ! Calculate damage in each linear softening law
      dmg_1 = MIN(one, eps_f_1*(eps_rfC - eps_0)/eps_rfC/(eps_f_1 - eps_0))
      dmg_2 = MIN(one, eps_f_2*(eps_rfC - eps_0)/eps_rfC/(eps_f_2 - eps_0))

      ! Find the combined damage variable, d1T
      d1C = n*(dmg_1 - dmg_2) + dmg_2

      ! Prevent healing; If the current damage is less than previous damage, use previous damage
      If (d1C .LT. d1C_temp) d1C = d1C_temp

      ! Delete element when it becomes fully damaged
      If (d1C .GE. d1MAX) Then
        d1C = d1MAX
        STATUS = 0
      End If

    Else
      ! This is only executed before any damage has occurred
      ! Since rfC is that state variable, it acts as the failure index before damage
      rfC = FIfC

      ! If there's no damage in tension, set rfT to zero
      If (d1T .LE. zero) rfT = zero

    End If

    Return
  End Subroutine FiberCompDmg

End Module CDM_fiber_mod
