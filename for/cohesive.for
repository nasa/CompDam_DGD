Module cohesive_mod
  ! Module for checking for cohesive damage initiation and calculating the cohesive damage variable

  Interface cohesive_damage
    Module Procedure cohesive_damage_matProp, cohesive_damage_toughnessListed, cohesive_damage_propsListed
  End Interface

Contains

  Subroutine cohesive_damage_matProp(m, p, delta, Pen, delta_n_init, B, FI, damage, dGdGc)
    ! Uses the properties in a matProps struct

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
    Double Precision, optional, intent(OUT) :: dGdGc
    ! -------------------------------------------------------------------- !

    Call cohesive_damage_propsListed(m%YT,m%SL,m%ST,m%GYT,m%GSL,m%etaL,m%etaT,m%eta_BK,p,delta,Pen,delta_n_init,B,FI,damage,dGdGc)
  End Subroutine cohesive_damage_matProp


  Subroutine cohesive_damage_toughnessListed(GYT, GSL, m, p, delta, Pen, delta_n_init, B, FI, damage, dGdGc)
    ! Uses the properties in a matProps struct

    Use forlog_Mod
    Use matProp_Mod
    Use parameters_Mod

    ! Arguments
    Double Precision, intent(IN) :: GYT, GSL
    Type(matProps), intent(IN) :: m
    Type(parameters), intent(IN) :: p
    Double Precision, intent(IN) :: delta(3)
    Double Precision, intent(IN) :: Pen(3), delta_n_init
    Double Precision, intent(INOUT) :: B
    Double Precision, intent(OUT) :: FI
    Double Precision, optional, intent(INOUT) :: damage
    Double Precision, optional, intent(OUT) :: dGdGc
    ! -------------------------------------------------------------------- !

    Call cohesive_damage_propsListed(m%YT,m%SL,m%ST,GYT,GSL,m%etaL,m%etaT,m%eta_BK,p,delta,Pen,delta_n_init,B,FI,damage,dGdGc)
  End Subroutine cohesive_damage_toughnessListed


  Subroutine cohesive_damage_propsListed(YT, SL, ST, GYT, GSL, etaL, etaT, eta_BK, p, delta, Pen, delta_n_init, B, FI, damage, dGdGc)
    ! All properties are passed in as arguments

    Use forlog_Mod
    Use parameters_Mod

    ! Arguments
    Double Precision, intent(IN) :: YT,SL,ST,GYT,GSL,etaL,etaT,eta_BK        ! Material properties
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
    Double Precision :: del0nB, del0sB, delfnB, delfsB                       ! Temporary, used to calculate d0 and df
    Double Precision :: damage_max                                           ! Maximum value for damage variable
    Double Precision :: beta                                                 ! Placeholder (temp.) variables for Mode-mixity
    Double Precision :: KS                                                   ! Shear penalty stiffness
    Double Precision :: B_old, damage_old
    Double Precision :: mode_mix_limit
    Double Precision, parameter :: zero=0.d0, one=1.d0, two=2.d0, ten=10.d0

    ! Fatigue
    Double Precision :: R_old, R_fatigue, R_fatigue_inc, R_max               ! Damage measures
    Double Precision :: damage_static                                        ! Temporary, used to compare static and fatigue damage
    Double Precision :: lambda_relative                                      ! Normalized fatigue displacement jump
    Double Precision :: fatigue_beta, fatigue_gamma, brittleness, pF         ! Fatigue damage parameters
    Double Precision :: endurance_relative_B, endurance_relative             ! relative endurance limit terms
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
      ds1_str = SL - etaL*Pen(2)*MIN(zero, delta_n_init)  ! LaRC04 longitudinal shear strength
      ds2_str = ST - etaT*Pen(2)*MIN(zero, delta_n_init)  ! LaRC04 transverse shear strength
      KS = SQRT((Pen(1)*del_s1)**2 + (Pen(3)*del_s2)**2)/del_s  ! Combined shear penalty stiffness
      ds_str = KS*del_s/SQRT((Pen(1)*del_s1/ds1_str)**2 + (Pen(3)*del_s2/ds2_str)**2)
    Else
      ds_str = SL
      KS = Pen(1)
    End If

    ! Cohesive displacements for initiation for pure mode I and mode II
    del0n = YT/Pen(2)  ! mode I cohesive initiation disp.
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
      del0nB = SQRT(((one - B**eta_BK)*del0n*del0n + KS/Pen(2)*B**eta_BK*del0s*del0s)*(one - B))
      del0sB = SQRT(beta/(one - beta))*del0nB
      d0 = SQRT(del0nB*del0nB + del0sB*del0sB)
    End If

    FI = MIN(one, del/d0)

    DamageEvolution: If (present(damage)) Then
      dGdGc = zero

      ! Cohesive displacements for final failure for pure mode I and mode II
      delfn = two*GYT/del0n/Pen(2)  ! mode I cohesive final disp.
      delfs = two*GSL/del0s/KS  ! mode II cohesive final disp.

      ! Mixed-mode final failure displacements
      If (B >= one - mode_mix_limit) Then
        df = delfs
      Else If (B <= mode_mix_limit) Then
        df = delfn
      Else
        delfnB = ((one - B**eta_BK)*del0n*delfn + KS/Pen(2)*B**eta_BK*del0s*delfs)*(one - B)/del0nB
        delfsB = SQRT(beta/(one - beta))*delfnB
        df = SQRT(delfnB*delfnB + delfsB*delfsB)
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

        ! Relative endurance at R=-1, corrected for mode mixity B:
        endurance_relative_B = 0.20d0*(one - B*0.42d0)
        ! Relative endurance for all R ratios and mode mixities:
        endurance_relative = two*endurance_relative_B/(endurance_relative_B + one + p%fatigue_R_ratio*(endurance_relative_B - one))

        ! Relative displacement jump
        lambda_relative = del/((one - R_old)*d0 + R_old*df)
        If (lambda_relative < endurance_relative) Then  ! threshold of propagation
          lambda_relative = zero
        Else If (lambda_relative > one) Then
          lambda_relative = one
        End If

        brittleness = 0.95d0  ! brittleness parameter; TODO: move to either material properties or parameters

        ! Coefficients for incremental fatigue damage calculation
        fatigue_beta = -7.d0*brittleness/log10(endurance_relative)  ! Equation 23
        fatigue_gamma = 1.d7  ! number of cycles to endurance
        pF = fatigue_beta  ! TODO: move to either material properties or parameters

        R_fatigue_inc = (one - R_old)**(fatigue_beta - pF) * (lambda_relative/endurance_relative)**fatigue_beta * &
                          cycles_per_increment(1)/(pF + one)/fatigue_gamma  ! CF20 Incremental fatigue damage

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
  End Subroutine cohesive_damage_propsListed


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
