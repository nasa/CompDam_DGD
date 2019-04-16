Module cohesive_mod
  ! Module for checking for cohesive damage initiation and calculating the cohesive damage variable

  Interface cohesive_damage
    Module Procedure cohesive_damage_matProp, cohesive_damage_toughnessListed, cohesive_damage_propsListed
  End Interface

Contains

  Pure Subroutine cohesive_damage_matProp(m, delta, Pen, delta_n_init, B, FI, damage, dGdGc)
    ! Uses the properties in a matProps struct

    Use forlog_Mod
    Use matProp_Mod

    ! Arguments
    Type(matProps), intent(IN) :: m
    Double Precision, intent(IN) :: delta(3)
    Double Precision, intent(IN) :: Pen(3), delta_n_init
    Double Precision, intent(INOUT) :: B
    Double Precision, intent(OUT) :: FI
    Double Precision, optional, intent(INOUT) :: damage
    Double Precision, optional, intent(OUT) :: dGdGc
    ! -------------------------------------------------------------------- !

    Call cohesive_damage_propsListed(m%YT,m%SL,m%ST,m%GYT,m%GSL,m%etaL,m%etaT,m%eta_BK,delta,Pen,delta_n_init,B,FI,damage,dGdGc)
  End Subroutine cohesive_damage_matProp


  Pure Subroutine cohesive_damage_toughnessListed(GYT, GSL, m, delta, Pen, delta_n_init, B, FI, damage, dGdGc)
    ! Uses the properties in a matProps struct

    Use forlog_Mod
    Use matProp_Mod

    ! Arguments
    Double Precision, intent(IN) :: GYT, GSL
    Type(matProps), intent(IN) :: m
    Double Precision, intent(IN) :: delta(3)
    Double Precision, intent(IN) :: Pen(3), delta_n_init
    Double Precision, intent(INOUT) :: B
    Double Precision, intent(OUT) :: FI
    Double Precision, optional, intent(INOUT) :: damage
    Double Precision, optional, intent(OUT) :: dGdGc
    ! -------------------------------------------------------------------- !

    Call cohesive_damage_propsListed(m%YT,m%SL,m%ST,GYT,GSL,m%etaL,m%etaT,m%eta_BK,delta,Pen,delta_n_init,B,FI,damage,dGdGc)
  End Subroutine cohesive_damage_toughnessListed


  Pure Subroutine cohesive_damage_propsListed(YT, SL, ST, GYT, GSL, etaL, etaT, eta_BK, delta, Pen, delta_n_init, B, FI, damage, dGdGc)
    ! All properties are passed in as arguments

    Use forlog_Mod
    Use matProp_Mod

    ! Arguments
    Double Precision, intent(IN) :: YT,SL,ST,GYT,GSL,etaL,etaT,eta_BK        ! Material properties
    Double Precision, intent(IN) :: delta(3)
    Double Precision, intent(IN) :: Pen(3), delta_n_init
    Double Precision, intent(INOUT) :: B
    Double Precision, intent(OUT) :: FI
    Double Precision, optional, intent(INOUT) :: damage
    Double Precision, optional, intent(OUT) :: dGdGc                        ! Dissipated energy per fracture toughness

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
    Double Precision, parameter :: zero=0.d0, one=1.d0, two=2.d0
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
    dGdGc = zero

    DamageEvolution: If (present(damage) .AND. FI == one) Then

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

      ! Check for damage advancement
      If (damage <= damage_old) Then  ! If there is no damage progression...
        damage = damage_old  ! ...use the old damage variable.
      Else  ! If there is damage progression...
        ! ...calculate energy dissipated by change in damage state
        dGdGc = damage/(df/d0*(one - damage) + damage) - damage_old/(df/d0*(one - damage_old) + damage_old)
      End If

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
