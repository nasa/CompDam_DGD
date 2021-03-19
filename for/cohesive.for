Module cohesive_mod
  ! Module for checking for cohesive damage initiation and calculating the cohesive damage variable

Contains

  Subroutine cohesive_damage(m, p, delta, Pen, delta_n_init, B, FI, evolve_fatigue, R, damage, dR)
    ! All properties are passed in as arguments

    Use forlog_Mod
    Use matProp_Mod
    Use parameters_Mod

    ! Arguments
    Type(matProps), intent(IN) :: m                                          ! Material properties
    Type(parameters), intent(IN) :: p                                        ! Solution parameters
    Double Precision, intent(IN) :: delta(3)                                 ! Cohesive displacement-jump vector
    Double Precision, intent(IN) :: Pen(3)                                   ! Penalty stiffnesses
    Double Precision, intent(IN) :: delta_n_init                             ! Normal cohesive displacement-jump at damage initiation
    Double Precision, intent(INOUT) :: B                                     ! Mode mixity
    Double Precision, intent(OUT) :: FI                                      ! Failure index
    Logical, intent(IN) :: evolve_fatigue                                    ! flag for calculating fatigue damage
    Double Precision, optional, intent(INOUT) :: R                           ! Normalized energy dissipation
    Double Precision, optional, intent(OUT) :: damage                        ! Cohesive damage state variable
    Double Precision, optional, intent(OUT) :: dR                            ! Change in normalized energy dissipation

    ! Locals
    Double Precision :: del                                                  ! Magnitude of current cohesive displacement-jump
    Double Precision :: del_n, del_s1, del_s2, del_s                         ! Temporary values used to calculate del
    Double Precision :: del0n, del0s                                         ! Initiation displacements for pure Mode I and Mode II
    Double Precision :: delfn, delfs                                         ! Final displacements for pure Mode I and Mode II
    Double Precision :: d0, df                                               ! Initiation and final displacements for current mode mixity
    Double Precision :: SL, ST, SC                                           ! Longitudinal, transverse, and combined shear strengths
    Double Precision :: R_max                                                ! Maximum value for the normalized energy dissipation
    Double Precision :: beta                                                 ! Intermediate variable for mode mixity
    Double Precision :: Pen_sh                                               ! Combined longitudinal and transverse shear penalty stiffness
    Double Precision :: B_old, R_old
    Double Precision :: R_static                                             ! Normalized energy dissipation, static
    Double Precision :: mode_mix_limit
    Double Precision, parameter :: zero=0.d0, one=1.d0, two=2.d0, seven=7.d0, ten=10.d0
    Logical :: evolve_damage

    ! Fatigue
    Double Precision :: R_fatigue, R_fatigue_inc                             ! Damage measures
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

    R_max = one  ! Maximum value for normalized energy dissipation damage variable
    mode_mix_limit = 1.d-6  ! Small variations (within this limit) from pure mode mixities are rounded to the pure cases

    ! Cohesive displacement-jumps
    del_s1 = delta(1)  ! Longitudinal shear
    del_s2 = delta(3)  ! Transverse shear
    del_n  = delta(2)  ! Normal
    del_s  = SQRT(del_s1*del_s1 + del_s2*del_s2)  ! Combined longitudinal and transverse shear
    del    = SQRT(MAX(zero, del_n)*del_n + del_s*del_s)  ! Combined normal and shear

    ! Account for LaRC04 and calculate "mixed shear" strengths and stiffnesses
    If (del_s > zero) Then
      SL = m%SL - m%etaL*MAX(-m%YC, Pen(2)*MIN(zero, delta_n_init, del_n))  ! LaRC04 longitudinal shear strength
      ST = m%ST - m%etaT*MAX(-m%YC, Pen(2)*MIN(zero, delta_n_init, del_n))  ! LaRC04 transverse shear strength
      Pen_sh = SQRT((Pen(1)*del_s1)**2 + (Pen(3)*del_s2)**2)/del_s  ! Combined shear penalty stiffness
      SC = Pen_sh*del_s/SQRT((Pen(1)*del_s1/SL)**2 + (Pen(3)*del_s2/ST)**2)  ! Combined shear strength
    Else
      SC = m%SL
      Pen_sh = Pen(1)
    End If

    ! Cohesive displacements for initiation for pure mode I and mode II
    del0n = m%YT/Pen(2)  ! mode I cohesive initiation disp.
    del0s = SC/Pen_sh  ! mode II cohesive initiation disp.

    ! Mode mixity
    beta = del_s*del_s + MAX(zero, del_n)*del_n
    If (beta == zero) Then
      beta = one
    Else
      beta = del_s*del_s/beta
    End If
    B_old = B
    B = Pen_sh*beta/(Pen_sh*beta + Pen(2)*(one - beta))

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
      d0 = SQRT( (del0n**2 + (del0s**2*Pen_sh/Pen(2) - del0n**2)*B**m%eta_BK)*(one + B*(Pen(2)/Pen_sh - one)) )
    End If

    FI = MIN(one, del/d0)  ! Failure index, static

    If (present(R)) Then
      Call log%debug("R is present. Damage can evolve.")
      evolve_damage = .TRUE.
      R_old = R  ! Store the current damage variable
    Else
      Call log%debug("R is not present.")
      evolve_damage = .FALSE.
      R_old = zero
    End If

    R_static = R_old
    R_fatigue = R_old
    R_fatigue_inc = zero

    DamageEvolution: If (evolve_damage .OR. evolve_fatigue) Then

      ! Cohesive displacements for final failure for pure mode I and mode II
      delfn = two*m%GYT/m%YT  ! mode I cohesive final disp.
      delfs = two*m%GSL/SC  ! mode II cohesive final disp.

      ! Mixed-mode final failure displacements
      If (B >= one - mode_mix_limit) Then
        df = delfs
      Else If (B <= mode_mix_limit) Then
        df = delfn
      Else
        df = (del0n*delfn + (del0s*delfs*Pen_sh/Pen(2) - del0n*delfn)*B**m%eta_BK)*(one + B*(Pen(2)/Pen_sh - one))/d0
      End If

      If (evolve_damage) Then

        ! Determine the new damage state
        R_static = (del - d0)/(df - d0)  ! Calculate the new normalized energy dissipation
        R_static = MIN(R_static, R_max)  ! Cap the normalized energy dissipation variable
        R_static = MAX(R_old, R_static)  ! Ensure no healing

      End If

#ifndef PYEXT

      ! If this is a fatigue step and the damage is not already fully developed...
      FatigueStep: If (evolve_fatigue .AND. R_old < R_max) Then

        ! Relative endurance at R=-1, corrected for mode mixity B
        epsilon_mixed = m%fatigue_epsilon*(one - B*0.42d0)
        ! Relative endurance for all R ratios and mode mixities
        endurance_relative = two*epsilon_mixed/(epsilon_mixed + one + p%fatigue_R_ratio*(epsilon_mixed - one))

        lambda_relative = del/((one - R_old)*d0 + R_old*df)  ! Relative displacement-jump
        
        FatigueEvolution: If (lambda_relative >= endurance_relative) Then  ! Fatigue damage evolution threshold

          Call log%debug("Fatigue damage is evolving.")

          ! Coefficients for incremental fatigue damage calculation
          fatigue_beta = -seven*m%fatigue_eta/log10(endurance_relative)  ! Equation 23 in CF20 technical paper
          fatigue_p_cf20 = fatigue_beta + m%fatigue_p_mod

          ptr_fatigue_dbl = SMARealArrayAccess(2)  ! Gain access to current cycles_per_increment value

          R_fatigue_inc = (one - R_old)**(fatigue_beta - fatigue_p_cf20) * (lambda_relative/endurance_relative)**fatigue_beta * &
                            cycles_per_increment(1)/(fatigue_p_cf20 + one)/m%fatigue_gamma  ! CF20 incremental fatigue damage
          R_fatigue_inc = MIN(R_fatigue_inc, R_max - R_old)  ! Cap the incremental fatigue damage so R does not exceed R_max
        
        End If FatigueEvolution

        Call log%debug("R_fatigue_inc:  " // trim(str(R_fatigue_inc)))
        
        R_fatigue = R_old + R_fatigue_inc  ! Update the fatigue normalized energy dissipation variable

        ptr_fatigue_int = SMAIntArrayAccess(1)  ! Gain access to the integer flag for controlling the cycles_per_increment
        
        ! Compare the current incremental fatigue growth to the upper and lower limit parameters
        If (R_fatigue_inc > p%fatigue_damage_max_threshold) Then
          fatigue_parameters(2) = MAX(fatigue_parameters(2), 2)
        Else If (R_fatigue_inc > p%fatigue_damage_min_threshold) Then
          fatigue_parameters(2) = MAX(fatigue_parameters(2), 1)
        End If

        R_fatigue = MIN(R_fatigue, R_max)  ! Cap the fatigue normalized energy dissipation variable

      End If FatigueStep
#endif

      If (evolve_damage) Then
        ! This corresponds to any call to cohesive_damage() in which R is an argument.

        R = MAX(R_old, R_static, R_fatigue)  ! Ensure no healing and use the greater of the static and fatigue damage
        damage = R/(R + d0/df*(one - R))  ! Penalty stiffness based cohesive damage variable
        dR = R - R_old  ! Change in normalized energy dissipation

      End If

      Call log%debug("End DamageEvolution If in cohesive_damage")

    End If DamageEvolution

    Call log%debug("End cohesive_damage")

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
