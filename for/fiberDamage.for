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
    Integer, intent(INOUT) :: STATUS                            ! Element deletion

    ! Locals
    Double Precision :: FIfT                        ! Fiber direction failure index
    Double Precision :: d1MAX                       ! Maximum value for damage variables
    Double Precision :: eps_rfT                     ! strain corresponding to the maximum historical failure criterion
    Double Precision :: eps_0, eps_f_A, eps_f_B     ! strains for damage variable calculations
    Double Precision :: dmg_A, dmg_B                ! individual damage variables for superposed softening laws
    Double Precision :: d1T_temp                    ! Stores the old value of d1T
    Double Precision :: E1_temp, E1_secant
    Double Precision :: eps_0_fn

    Double Precision, Parameter :: zero=0.d0, one=1.d0, two=2.d0, four=4.d0
    ! -------------------------------------------------------------------- !

    ! Initialize
    d1MAX = one
    d1T = MAX(d1T, d1C)
    d1T_temp = d1T
    eps_0_fn = (-E1+SQRT(E1**two+four*E1*cl*XT))/(two*E1*cl)
    E1_secant = XT/eps_0_fn

    If (cl > 0) Then  ! Account for fiber nonlinearity
      If (d1T > 0) Then   ! Use secant stiffness
        E1_temp = E1_secant
      Else  ! Use tangent stiffmess
        E1_temp = E1*(1+cl*eps(1,1))
      End If
    Else
      E1_temp = E1
    End If

    ! Fiber tension failure criterion
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
      eps_f_A = two*GXT*m/(XT*n*Lc)
      eps_f_B = two*GXT*(one - m)/(XT*(one - n)*Lc)

      ! Calculate damage in each linear softening law
      If (eps_f_A <= eps_0) Then  ! Snap-back
        dmg_A = one
      Else
        dmg_A = MIN(one, eps_f_A*(eps_rfT - eps_0)/eps_rfT/(eps_f_A - eps_0))
      End If
      If (eps_f_B <= eps_0) Then  ! Snap-back
        dmg_B = one
      Else
        dmg_B = MIN(one, eps_f_B*(eps_rfT - eps_0)/eps_rfT/(eps_f_B - eps_0))
      End If
      
      ! Find the combined damage variable, d1T
      d1T = n*(dmg_A - dmg_B) + dmg_B

      ! Prevent healing; if the current damage is less than previous damage, use previous damage
      If (d1T < d1T_temp) d1T = d1T_temp

      ! Delete element when it becomes fully damaged
      If (d1T >= d1MAX) Then
        d1T = d1MAX
        ! STATUS = 0
      End If

    Else
      ! This is only executed before any damage has occurred
      ! If no damage has occurred, use rfT as a fiber tension failure criterion state variable
      rfT = FIfT

    End If

    Return
  End Subroutine FiberTenDmg


  Pure Subroutine FiberCompDmg(eps, ndir, E1, XC, GXC, n, m, Lc, cl, rfT, rfC, d1T, d1C, STATUS)
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
    Integer, intent(INOUT) :: STATUS                              ! Element deletion

    ! Locals
    Double Precision :: FIfC                        ! Fiber direction failure index
    Double Precision :: d1MAX                       ! Maximum value for damage variables
    Double Precision :: eps_rfC                     ! strain corresponding to the maximum historical failure criterion
    Double Precision :: eps_0, eps_f_A, eps_f_B     ! strains for damage variable calculations
    Double Precision :: dmg_A, dmg_B                ! individual damage variables for superposed softening laws
    Double Precision :: d1C_temp                    ! Stores the old value of d1C
    Double Precision :: E1_temp, E1_secant
    Double Precision :: eps_0_fn

    Double Precision, Parameter :: zero=0.d0, one=1.d0, two=2.d0, four=4.d0
    ! -------------------------------------------------------------------- !

    ! Initialize
    d1MAX = one
    d1C_temp = d1C
    eps_0_fn = -(-E1+SQRT(E1**two-four*E1*cl*XC))/(two*E1*cl)
    E1_secant = XC/eps_0_fn

    If (cl > 0) Then  ! Account for fiber nonlinearity
      If (d1C > 0) Then  ! Use secant stiffness
        E1_temp = E1_secant
      Else  ! Use tangent stiffmess
        E1_temp = E1*(1+cl*eps(1,1))
      End If
    Else
      E1_temp = E1
    End If

    ! Fiber compression failure criterion
    FIfC = -E1_temp/XC*eps(1,1)

    ! Update damage threshold
    If (FIfC > rfC) rfC = FIfC

    ! Determine d1C if damage has occurred
    If (rfC > one) Then

      eps_rfC = rfC*XC/E1_temp

      ! Calculate constants
      If (cl > 0) Then  ! Account for fiber nonlinearity, use secant stiffness
        eps_0 = eps_0_fn
      Else
        eps_0 = XC/E1
      End If
      eps_f_A = two*GXC*m/(XC*n*Lc)
      eps_f_B = two*GXC*(one - m)/(XC*(one - n)*Lc)

      ! Calculate damage in each linear softening law
      If (eps_f_A <= eps_0) Then  ! Snap-back
        dmg_A = one
      Else
        dmg_A = MIN(one, eps_f_A*(eps_rfC - eps_0)/eps_rfC/(eps_f_A - eps_0))
      End If
      If (eps_f_B <= eps_0) Then  ! Snap-back
        dmg_B = one
      Else
        dmg_B = MIN(one, eps_f_B*(eps_rfC - eps_0)/eps_rfC/(eps_f_B - eps_0))
      End If

      ! Find the combined damage variable, d1C
      d1C = n*(dmg_A - dmg_B) + dmg_B

      ! Prevent healing; If the current damage is less than previous damage, use previous damage
      If (d1C < d1C_temp) d1C = d1C_temp

      ! Delete element when it becomes fully damaged
      If (d1C >= d1MAX) Then
        d1C = d1MAX
        ! STATUS = 0
      End If

    Else
      ! This is only executed before any damage has occurred
      ! If no damage has occurred, use rfC as a fiber compression failure criterion state variable
      rfC = FIfC

      ! If there's no damage in tension, set rfT to zero
      If (d1T <= zero) rfT = zero

    End If

    Return
  End Subroutine FiberCompDmg

End Module CDM_fiber_mod
