Module cohesive_mod
  ! Module for checking for cohesive damage initiation and calculating the cohesive damage variable

Contains

  Subroutine cohesive_damage(m, p, delta, Pen, delta_n_init, B, FI, damage, dGdGc)
    ! All properties are passed in as arguments

    Use forlog_Mod
    Use matProp_Mod
    Use parameters_Mod

    ! Arguments
    Type(matProps), intent(IN) :: m
    Type(parameters), intent(IN) :: p
    Double Precision, intent(IN) :: delta(3)
    Double Precision, intent(IN) :: Pen(3), delta_n_init
    Double Precision, intent(INOUT) :: B
    Double Precision, intent(OUT) :: FI
    Double Precision, optional, intent(INOUT) :: damage
    Double Precision, optional, intent(OUT) :: dGdGc                        ! Normalized energy dissipation

    ! Locals
    Double Precision :: del                                                  ! Magnitude of current cohesive displacement
    Double Precision :: del_n, del_s1, del_s2, del_s                         ! Temporary values used to calculate del
    Double Precision :: del0n, del0s                                         ! Initiation displacements for pure Mode I and Mode II
    Double Precision :: delfn, delfs                                         ! Final displacements for pure Mode I and Mode II
    Double Precision :: d0, df                                               ! Initiation and final displacements for current mode-mixity ratio
    Double Precision :: damage_max                                           ! Maximum value for damage variable
    Double Precision :: beta                                                 ! Placeholder (temp.) variables for Mode-mixity
    Double Precision :: KS                                                   ! Shear penalty stiffness
    Double Precision :: B_old, damage_old
    Double Precision :: mode_mix_limit
    Double Precision, parameter :: zero=0.d0, one=1.d0, two=2.d0, seven=7.d0, ten=10.d0

    ! Fatigue
    Double Precision :: R_old, R_fatigue, R_fatigue_inc, R_max               ! Damage measures
    Double Precision :: damage_static                                        ! Temporary, used to compare static and fatigue damage
    Double Precision :: lambda_relative                                      ! Normalized fatigue displacement jump
    Double Precision :: fatigue_beta, fatigue_p_cf20                         ! beta and p parameters from CF20 paper
    Double Precision :: epsilon_mixed, endurance_relative                    ! relative endurance limit terms
    Integer :: fatigue_parameters(2)
    Double Precision :: cycles_per_increment(1)
    Pointer(ptr_fatigue_int, fatigue_parameters)
    Pointer(ptr_fatigue_dbl, cycles_per_increment)

#ifndef PYEXT
    Interface

      Function SMAIntArrayAccess(ID)
        ! Access an existing global integer array.
        Integer(kind=8) :: SMAIntArrayAccess  ! Returns an address that can be associated with a Fortran pointer
        Integer(kind=4) :: ID       ! Array ID
      End Function SMAIntArrayAccess

      Function SMARealArrayAccess(ID)
        ! Access an existing global real array.
        Integer(kind=8) :: SMARealArrayAccess  ! Returns an address that can be associated with a Fortran pointer
        Integer(kind=4) :: ID   ! Array ID
      End Function SMARealArrayAccess

    End Interface
#endif
    ! -------------------------------------------------------------------- !

    damage_max = one  ! Maximum value for damage variables
    mode_mix_limit = 1.d-6

    ! Cohesive displacement-jumps
    del_s1 = delta(1)
    del_s2 = delta(3)
    del_n  = delta(2)
    del_s  = SQRT(del_s1*del_s1 + del_s2*del_s2)
    del    = SQRT(MAX(zero, del_n)*del_n + del_s*del_s)

    ! Account for LaRC04 and calculate "mixed shear" strengths and stiffnesses
    If (del_s > zero) Then
      ds1_str = m%SL - m%etaL*Pen(2)*MIN(zero, delta_n_init)  ! LaRC04 longitudinal shear strength
      ds2_str = m%ST - m%etaT*Pen(2)*MIN(zero, delta_n_init)  ! LaRC04 transverse shear strength
      KS = SQRT((Pen(1)*del_s1)**2 + (Pen(3)*del_s2)**2)/del_s  ! Combined shear penalty stiffness
      ds_str = KS*del_s/SQRT((Pen(1)*del_s1/ds1_str)**2 + (Pen(3)*del_s2/ds2_str)**2)
    Else
      ds_str = m%SL
      KS = Pen(1)
    End If

    ! Cohesive displacements for initiation for pure mode I and mode II
    del0n = m%YT/Pen(2)  ! mode I cohesive initiation disp.
    del0s = ds_str/KS  ! mode II cohesive initiation disp.

    ! Mode mixity
    beta = del_s*del_s + MAX(zero, del_n)*del_n
    If (beta == zero) Then
      beta = one
    Else
      beta = del_s*del_s/beta
    End If
    B_old = B
    B = KS*beta/(KS*beta + Pen(2)*(one - beta))

    ! Mixed-mode initiation displacements
    If (B >= one - mode_mix_limit) Then
      beta = one
      B = one
      d0 = del0s
    Else If (B <= mode_mix_limit) Then
      beta = zero
      B = zero
      d0 = del0n
    Else
      d0 = SQRT( (del0n**2 + (del0s**2*KS/Pen(2) - del0n**2)*B**m%eta_BK)*(one + B*(Pen(2)/KS - one)) )
    End If

    FI = MIN(one, del/d0)

    DamageEvolution: If (present(damage)) Then
      dGdGc = zero

      ! Cohesive displacements for final failure for pure mode I and mode II
      delfn = two*m%GYT/del0n/Pen(2)  ! mode I cohesive final disp.
      delfs = two*m%GSL/del0s/KS  ! mode II cohesive final disp.

      ! Mixed-mode final failure displacements
      If (B >= one - mode_mix_limit) Then
        df = delfs
      Else If (B <= mode_mix_limit) Then
        df = delfn
      Else
        df = (del0n*delfn + (del0s*delfs*KS/Pen(2) - del0n*delfn)*B**m%eta_BK)*(one + B*(Pen(2)/KS - one))/d0
      End If

      ! Determine the new damage state
      damage_old = damage  ! Store the current damage variable
      If (del == zero) Then
        damage = zero
      Else
        damage = df*(del - d0)/del/(df - d0)  ! Calculate the new damage variable
      End If
      damage = MIN(damage, damage_max)  ! Cap damage variable
      damage = MAX(damage_old, damage)  ! Ensure no healing

#ifndef PYEXT
      ptr_fatigue_int = SMAIntArrayAccess(1)

      ! If this is a fatigue step, and static damage progression has not occurred, and the damage is not already fully developed...
      FatigueStep: If (fatigue_parameters(1) == 1 .AND. damage < damage_max) Then

        ptr_fatigue_dbl = SMARealArrayAccess(2)

        damage_static = damage

        R_old = damage_old*d0/((one - damage_old)*df + damage_old*d0)  ! Transform old damage to R form
        R_max = damage_max*d0/((one - damage_max)*df + damage_max*d0)  ! Transform damage_max to R form

        ! Relative endurance at R=-1, corrected for mode mixity B
        epsilon_mixed = m%fatigue_epsilon*(one - B*0.42d0)
        ! Relative endurance for all R ratios and mode mixities
        endurance_relative = two*epsilon_mixed/(epsilon_mixed + one + p%fatigue_R_ratio*(epsilon_mixed - one))

        ! Coefficient fatigue_beta for incremental fatigue damage calculation
        fatigue_beta = -seven*m%fatigue_eta/log10(endurance_relative)  ! Equation 23 in CF20 technical paper
        fatigue_p_cf20 = fatigue_beta + m%fatigue_p_mod

        lambda_relative = min(del/((one - R_old)*d0 + R_old*df), one)  ! Relative displacement jump
        If (lambda_relative < endurance_relative) lambda_relative = zero  ! threshold of propagation

        R_fatigue_inc = (one - R_old)**(fatigue_beta - fatigue_p_cf20) * (lambda_relative/endurance_relative)**fatigue_beta * &
                          cycles_per_increment(1)/(fatigue_p_cf20 + one)/m%fatigue_gamma  ! CF20 incremental fatigue damage

        R_fatigue_inc = MIN(R_fatigue_inc, R_max - R_old)
        R_fatigue = R_old + R_fatigue_inc  ! Updated fatigue damage variable

        ! Compare the current incremental fatigue growth to the upper and lower limit parameters
        If (R_fatigue_inc > p%fatigue_damage_max_threshold) Then
          fatigue_parameters(2) = MAX(fatigue_parameters(2), 2)
        Else If (R_fatigue_inc > p%fatigue_damage_min_threshold) Then
          fatigue_parameters(2) = MAX(fatigue_parameters(2), 1)
        End If

        damage = df*R_fatigue/(d0*(one - R_fatigue) + df*R_fatigue)  ! Transform R_fatigue to penalty stiffness form
        damage = MIN(damage, damage_max)  ! Cap the damage variable
        damage = MAX(damage_old, damage, damage_static)  ! Ensure no healing and use the greater of the static and fatigue damage

      End If FatigueStep
#endif

      ! If there is damage progression, calculate the normalized energy dissipation.
      If (damage > damage_old) dGdGc = damage/(df/d0*(one - damage) + damage) - damage_old/(df/d0*(one - damage_old) + damage_old)

    End If DamageEvolution

    Return
  End Subroutine cohesive_damage


  Pure Function cohesive_traction(delta, penalty, damage) result(traction)
    ! Calculates the cohesive tractions from displacement, penalty stiffness, and damage.

    ! Arguments
    Double Precision, intent(IN) :: delta(3)    ! cohesive displacement-jump vector
    Double Precision, intent(IN) :: penalty(3)  ! cohesive penalty stiffness
    Double Precision, intent(IN) :: damage      ! cohesive damage variable
    Double Precision :: traction(3)

    ! Locals
    Double Precision, parameter :: zero=0.d0, one=1.d0
    ! -------------------------------------------------------------------- !
    traction(:) = penalty(:)*(one - damage)*delta(:)
    If (delta(2) < zero) traction(2) = penalty(2)*delta(2)  ! normal direction

    Return
  End Function cohesive_traction


End Module cohesive_mod
