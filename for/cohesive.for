Module cohesive_mod
  ! Module for calculating the cohesive damage variable

  Interface cohesive_damage
    Module Procedure cohesive_damage_matProp, cohesive_damage_toughnessListed, cohesive_damage_propsListed
  End Interface

Contains

  Pure Subroutine cohesive_damage_matProp(m,delta,Pen,mCompInit,B,FI,damage,dGdGc)
    ! Uses the properties in a matProps struct

    Use forlog_Mod
    Use matProp_Mod

    Include 'vaba_param.inc'

    ! Arguments
    Type(matProps), intent(IN) :: m
    Double Precision, intent(IN) :: delta(3)
    Double Precision, intent(IN) :: Pen(3), mCompInit
    Double Precision, intent(INOUT) :: B
    Double Precision, intent(OUT) :: FI
    Double Precision, optional, intent(INOUT) :: damage
    Double Precision, optional, intent(OUT) :: dGdGc
    ! -------------------------------------------------------------------- !

    Call cohesive_damage_propsListed(m%YT,m%SL,m%ST,m%GYT,m%GSL,m%etaL,m%etaT,m%eta_BK,delta,Pen,mCompInit,B,FI,damage,dGdGc)
  End Subroutine cohesive_damage_matProp


  Pure Subroutine cohesive_damage_toughnessListed(GYT,GSL,m,delta,Pen,mCompInit,B,FI,damage,dGdGc)
    ! Uses the properties in a matProps struct

    Use forlog_Mod
    Use matProp_Mod

    Include 'vaba_param.inc'

    ! Arguments
    Double Precision, intent(IN) :: GYT,GSL                                  ! Toughnesses
    Type(matProps), intent(IN) :: m
    Double Precision, intent(IN) :: delta(3)
    Double Precision, intent(IN) :: Pen(3), mCompInit
    Double Precision, intent(INOUT) :: B
    Double Precision, intent(OUT) :: FI
    Double Precision, optional, intent(INOUT) :: damage
    Double Precision, optional, intent(OUT) :: dGdGc
    ! -------------------------------------------------------------------- !

    Call cohesive_damage_propsListed(m%YT,m%SL,m%ST,GYT,GSL,m%etaL,m%etaT,m%eta_BK,delta,Pen,mCompInit,B,FI,damage,dGdGc)
  End Subroutine cohesive_damage_toughnessListed


  Pure Subroutine cohesive_damage_propsListed(YT,SL,ST,GYT,GSL,etaL,etaT,eta_BK,delta,Pen,mCompInit,B,FI,damage,dGdGc)
    ! All properties are passed in as arguments

    Use forlog_Mod
    Use matProp_Mod

    Include 'vaba_param.inc'

    ! Arguments
    Double Precision, intent(IN) :: YT,SL,ST,GYT,GSL,etaL,etaT,eta_BK        ! Material properties
    Double Precision, intent(IN) :: delta(3)
    Double Precision, intent(IN) :: Pen(3), mCompInit
    Double Precision, intent(INOUT) :: B
    Double Precision, intent(OUT) :: FI
    Double Precision, optional, intent(INOUT) :: damage
    Double Precision, optional, intent(OUT) :: dGdGc

    ! Locals
    Double Precision :: del                                                  ! Magnitude of current cohesive displacement
    Double Precision :: del_n, del_s1, del_s2, del_s                         ! Temporary values used to calculate del
    Double Precision :: del0n, del0s                                         ! Initation displacements for pure Mode I and Mode II
    Double Precision :: delfn, delfs                                         ! Final displacements for pure Mode I and Mode II
    Double Precision :: d0, df                                               ! Initation and final displacements for current mode-mixity ratio
    Double Precision :: del0nB, del0sB, delfnB, delfsB                       ! Temporary, used to calculate d0 and df
    Double Precision :: dMAX                                                 ! Maximum value for damage variable
    Double Precision :: beta                                                 ! Placeholder (temp.) variables for Mode-mixity
    Double Precision :: KS                                                   ! Penalty stiffnesses
    Double Precision :: B_temp, dmg_temp
    Double Precision, parameter :: zero=0.d0, one=1.d0, two=2.d0
    ! -------------------------------------------------------------------- !

    dMAX = one ! Maximum value for damage variables

    ! Cohesive displacement jumps
    del_s1 = delta(1)
    del_s2 = delta(3)
    del_n  = delta(2)
    del_s  = SQRT(del_s1*del_s1 + del_s2*del_s2)
    del    = SQRT(MAX(zero, del_n)*del_n + del_s*del_s)

    ! Account for LaRC04 and calculate "mixed shear" strengths and stiffnesses
    If (del_s .GT. zero) Then
      ds1_str = SL - etaL*Pen(2)*MIN(zero, MAX(mCompInit, del_n))  ! LaRC04 longitudinal shear strength
      ds2_str = ST - etaT*Pen(2)*MIN(zero, MAX(mCompInit, del_n))  ! LaRC04 transverse shear strength
      ds_str = SQRT((Pen(1)*del_s1/ds1_str)**2 + (Pen(3)*del_s2/ds2_str)**2)*GSL*YT*YT/(Pen(2)*del_s*GYT)
      KS = Pen(2)*GYT/GSL/YT/YT*ds_str*ds_str  ! Combined shear penalty stiffness
    Else
      ds_str = SL
      KS = Pen(1)
    End If

    ! Cohesive displacements for inititation for pure mode I and mode II
    del0n = YT/Pen(2)             ! mode I cohesive initiation disp.
    del0s = ds_str/KS             ! mode II cohesive initiation disp.

    ! Mode mixity
    beta = del_s*del_s + MAX(zero, del_n)*del_n
    If (beta .EQ. zero) Then
      beta = one
    Else
      beta = del_s*del_s/beta
    End If
    B_temp = B
    B = KS*beta/(KS*beta + Pen(2)*(one - beta))

    ! Mixed-mode initiation displacements
    If (B .GE. 0.999999d0) Then
      d0 = del0s
    Else If (B .LE. 1.d-6) Then
      d0 = del0n
    Else
      del0nB = SQRT(((one - B**eta_BK)*del0n*del0n + KS/Pen(2)*B**eta_BK*del0s*del0s)*(one - B))
      del0sB = SQRT(beta/(one - beta))*del0nB
      d0 = SQRT(del0nB*del0nB + del0sB*del0sB)
    End If

    FI = MIN(one, del/d0)

    If (present(damage)) Then

      ! Cohesive displacements for final failure for pure mode I and mode II
      delfn = two*GYT/del0n/Pen(2)  ! mode I cohesive final disp.
      delfs = two*GSL/del0s/KS  ! mode II cohesive final disp.

      ! Mixed-mode final failure displacements
      If (B .GE. 0.999999d0) Then
        df = delfs
      Else If (B .LE. 1.d-6) Then
        df = delfn
      Else
        delfnB = ((one - B**eta_BK)*del0n*delfn + KS/Pen(2)*B**eta_BK*del0s*delfs)*(one - B)/del0nB
        delfsB = SQRT(beta/(one - beta))*delfnB
        df = SQRT(delfnB*delfnB + delfsB*delfsB)
      End If

      ! Determine the new damage state
      dmg_temp = damage                     ! Store the current damage variable
      damage = df*(del - d0)/del/(df - d0)  ! Calculate the new damage variable
      damage = MIN(damage, dMAX)            ! Cap damage at dMAX

      ! Check for damage advancement

      ! If there is no damage progression...
      If (damage .LE. dmg_temp) Then

        ! ...use the old damage and mode mixity variables.
        damage = dmg_temp
        B = B_temp
      End If

      ! Energy dissipated by change in damage state
      dGdGc = damage/(df/d0*(one - damage) + damage) - dmg_temp/(df/d0*(one - dmg_temp) + dmg_temp)

    End If

    ! -------------------------------------------------------------------- !
    Return
  End Subroutine cohesive_damage_propsListed

End Module cohesive_mod
