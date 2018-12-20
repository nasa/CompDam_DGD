Module CDM_fiber_mod
  ! Continuum damage mechanics fiber damage calculations

Contains

  Subroutine FiberTenDmg(eps, ndir, E1, XT, GXT, n, m, Lc, cl, rfT, d1T, d1C, STATUS)
    ! The purpose of this subroutine is to evaluate fiber damage.
    ! This subroutine is for tension-only fiber damage.

    ! Arguments
    Double Precision, intent(IN) :: E1                        ! Stiffness
    Double Precision, intent(IN) :: XT                        ! Strength
    Double Precision, intent(IN) :: GXT                       ! Toughness
    Double Precision, intent(IN) :: n, m                      ! Coefficients for bi-linear softening
    Double Precision, intent(IN) :: eps(ndir,ndir)            ! Strain
    Double Precision, intent(IN) :: Lc
    Double Precision, intent(IN) :: d1C                       ! CDM compressive fiber damage variable
    Integer, intent(IN) :: ndir
    Double Precision, intent(IN) :: cl

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
    Double Precision :: E1_temp, E1_secant
    Double Precision :: eps_0_fn

    Double Precision, Parameter :: zero=0.d0, one=1.d0, two=2.d0, four=4.d0
    ! -------------------------------------------------------------------- !

    ! Initialize
    d1MAX = one
    d1T = MAX(d1T, d1C)
    d1T_temp = d1T
    STATUS = 1
    eps_0_fn = (-E1+SQRT(E1**two+four*E1*cl*XT))/(two*E1*cl)
    E1_secant = XT/eps_0_fn

    ! Fiber tension failure criterion
    If (cl > 0) Then  ! Account for fiber nonlinearity
      If (d1T > 0) Then   ! Use secant stiffness
        E1_temp = E1_secant
      Else  ! Use tangent stiffmess
        E1_temp = E1*(1+cl*eps(1,1))
      End If
    Else
      E1_temp = E1
    End If
    FIfT = E1_temp/XT*eps(1,1)

    ! Update damage threshold
    If (FIfT > rfT) rfT = FIfT

    ! Determine d1T if damage has occurred
    If (rfT > one) Then

      eps_rfT = rfT*XT/E1_temp

      ! Calculate constants
      If (cl > 0) Then  ! Account for fiber nonlinearity, use secant stiffness
        eps_0 = eps_0_fn
      Else
        eps_0 = XT/E1
      End If
      eps_f_1 = two*GXT*m/(XT*n*Lc)
      eps_f_2 = two*GXT*(one - m)/(XT*(one - n)*Lc)

      ! Calculate damage in each linear softening law
      dmg_1 = MIN(one, eps_f_1*(eps_rfT - eps_0)/eps_rfT/(eps_f_1 - eps_0))
      dmg_2 = MIN(one, eps_f_2*(eps_rfT - eps_0)/eps_rfT/(eps_f_2 - eps_0))

      ! Find the combined damage variable, d1T
      d1T = n*(dmg_1 - dmg_2) + dmg_2

      ! Prevent healing; If the current damage is less than previous damage, use previous damage
      If (d1T < d1T_temp) d1T = d1T_temp

      ! Delete element when it becomes fully damaged
      If (d1T >= d1MAX) Then
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


  Pure Subroutine FiberCompDmg(eps,ndir,E1,XC,GXC,n,m,Lc,cl,rfT,rfC,d1T,d1C,STATUS)
    ! The purpose of this subroutine is to evaluate fiber damage.
    ! This subroutine is for compression-only fiber damage.

    ! Arguments
    Double Precision, intent(IN) :: E1                        ! Stiffness
    Double Precision, intent(IN) :: XC                        ! Strength
    Double Precision, intent(IN) :: GXC                       ! Toughness
    Double Precision, intent(IN) :: n, m                      ! Coefficients for bi-linear softening
    Double Precision, intent(IN) :: eps(ndir,ndir)            ! Strain
    Double Precision, intent(IN) :: Lc
    Integer, intent(IN) :: ndir
    Double Precision, intent(IN) :: cl

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
    Double Precision :: E1_temp, E1_secant
    Double Precision :: eps_0_fn

    Double Precision, Parameter :: zero=0.d0, one=1.d0, two=2.d0, four=4.d0
    ! -------------------------------------------------------------------- !

    ! Initialize
    d1MAX = one
    d1C_temp = d1C
    STATUS = 1
    eps_0_fn = -(-E1+SQRT(E1**two-four*E1*cl*XC))/(two*E1*cl)
    E1_secant = XC/eps_0_fn

    ! Fiber compression failure criterion
    If (cl > 0) Then  ! Account for fiber nonlinearity
      If (d1C > 0) Then   ! Use secant stiffness
        E1_temp = E1_secant
      Else  ! Use tangent stiffmess
        E1_temp = E1*(1+cl*eps(1,1))
      End If
    Else
      E1_temp = E1
    End If
    FIfC = -E1_temp/XC*eps(1,1)

    ! Update damage threshold
    If (FIfC > rfC) Then
      rfC = FIfC

      ! For load reversal
      ! If (rfC > rfT) rfT = rfC
    End If

    ! Determine d1C if damage has occurred
    If (rfC > one) Then

      eps_rfC = rfC*XC/E1_temp

      ! Calculate constants
      If (cl > 0) Then  ! Account for fiber nonlinearity, use secant stiffness
        eps_0 = eps_0_fn
      Else
        eps_0 = XC/E1
      End If
      eps_f_1 = two*GXC*m/(XC*n*Lc)
      eps_f_2 = two*GXC*(one - m)/(XC*(one - n)*Lc)

      ! Calculate damage in each linear softening law
      dmg_1 = MIN(one, eps_f_1*(eps_rfC - eps_0)/eps_rfC/(eps_f_1 - eps_0))
      dmg_2 = MIN(one, eps_f_2*(eps_rfC - eps_0)/eps_rfC/(eps_f_2 - eps_0))

      ! Find the combined damage variable, d1T
      d1C = n*(dmg_1 - dmg_2) + dmg_2

      ! Prevent healing; If the current damage is less than previous damage, use previous damage
      If (d1C < d1C_temp) d1C = d1C_temp

      ! Delete element when it becomes fully damaged
      If (d1C >= d1MAX) Then
        d1C = d1MAX
        STATUS = 0
      End If

    Else
      ! This is only executed before any damage has occurred
      ! Since rfC is that state variable, it acts as the failure index before damage
      rfC = FIfC

      ! If there's no damage in tension, set rfT to zero
      If (d1T <= zero) rfT = zero

    End If

    Return
  End Subroutine FiberCompDmg

End Module CDM_fiber_mod
