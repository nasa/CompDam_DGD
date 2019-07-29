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
    Double Precision :: XLam                                                 ! Normalized fatigue displacement jump
    Double Precision :: fatigue_beta, fatigue_gamma                          ! Fatigue damage parameters
    Integer :: fatigue_parameters(2)
    Double Precision :: cycles_per_increment(1)
    Pointer(ptr_fatigue_int, fatigue_parameters)
    Pointer(ptr_fatigue_dbl, cycles_per_increment)

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

      ptr_fatigue_int = SMAIntArrayAccess(1)

      ! If this is a fatigue step, and static damage progression has not occurred, and the damage is not already fully developed...
      FatigueStep: If (fatigue_parameters(1) == 1 .AND. damage < damage_max) Then

        ptr_fatigue_dbl = SMARealArrayAccess(2)

        damage_static = damage

        R_old = damage_old*d0/((one - damage_old)*df + damage_old*d0)  ! Transform old damage to R form
        R_max = damage_max*d0/((one - damage_max)*df + damage_max*d0)  ! Transform damage_max to R form

        Call fatigue_coefficients(B, p%fatigue_R_ratio, fatigue_beta, fatigue_gamma)  ! Calculate fatigue damage model coefficients
        XLam = del/((one - R_old)*d0 + R_old*df)
        R_fatigue_inc = (R_old + fatigue_gamma)*XLam**fatigue_beta*cycles_per_increment(1)  ! Incremental fatigue damage
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


  Subroutine fatigue_coefficients(mode_mixity, R_Load, Bet, Gam)
    ! Returns fatigue coefficients beta and gamma as a function of mode mixity and fatigue R ratio.
    ! See Dávila, C. G., "From S-N to the Paris Law with a New Mixed-Mode Cohesive Fatigue Model" NASA/TP-2018-219838 for details.

    Use forlog_Mod

    ! Arguments
    Double Precision, intent(IN) :: mode_mixity  ! mode mixity
    Double Precision, intent(IN) :: R_Load       ! stress ratio
    Double precision, intent(OUT) :: Bet, Gam

    ! Locals
    Double Precision :: RTable(50), BetTable(50), GamTable(50)
    Double Precision :: R_effective    ! effective stress ratio
    Double Precision :: interpolation  ! lookup table interpolation coefficient
    Logical :: iFound
    Double Precision, parameter :: one=1.d0, two=2.d0
    ! -------------------------------------------------------------------- !
    Do i = 1,50; RTable(i) = -1.50d0 + 0.05d0*(i-1); End Do

    BetTable(1)=	11.916971
    BetTable(2)=	12.058467
    BetTable(3)=	12.205484
    BetTable(4)=	12.358368
    BetTable(5)=	12.517489
    BetTable(6)=	12.683265
    BetTable(7)=	12.856132
    BetTable(8)=	13.036582
    BetTable(9)=	13.225149
    BetTable(10)=	13.422418
    BetTable(11)=	13.62903
    BetTable(12)=	13.845697
    BetTable(13)=	14.073205
    BetTable(14)=	14.312421
    BetTable(15)=	14.564311
    BetTable(16)=	14.829949
    BetTable(17)=	15.110542
    BetTable(18)=	15.407437
    BetTable(19)=	15.722153
    BetTable(20)=	16.056402
    BetTable(21)=	16.41213
    BetTable(22)=	16.79155
    BetTable(23)=	17.197189
    BetTable(24)=	17.63196
    BetTable(25)=	18.099216
    BetTable(26)=	18.602858
    BetTable(27)=	19.147438
    BetTable(28)=	19.738302
    BetTable(29)=	20.381792
    BetTable(30)=	21.085445
    BetTable(31)=	21.858343
    BetTable(32)=	22.711489
    BetTable(33)=	23.658363
    BetTable(34)=	24.71566
    BetTable(35)=	25.904301
    BetTable(36)=	27.250877
    BetTable(37)=	28.789688
    BetTable(38)=	30.565781
    BetTable(39)=	32.639507
    BetTable(40)=	35.093624
    BetTable(41)=	38.044724
    BetTable(42)=	41.662594
    BetTable(43)=	46.204231
    BetTable(44)=	52.07832
    BetTable(45)=	59.976353
    BetTable(46)=	71.167618
    BetTable(47)=	88.262215
    BetTable(48)=	117.610985
    BetTable(49)=	179.723572
    BetTable(50)=	397.755554

    GamTable(1)=	0.00186604
    GamTable(2)=	0.00186981
    GamTable(3)=	0.00187372
    GamTable(4)=	0.0018778
    GamTable(5)=	0.00188205
    GamTable(6)=	0.00188649
    GamTable(7)=	0.00189112
    GamTable(8)=	0.00189596
    GamTable(9)=	0.00190103
    GamTable(10)=	0.00190634
    GamTable(11)=	0.00191191
    GamTable(12)=	0.00191777
    GamTable(13)=	0.00192392
    GamTable(14)=	0.00193041
    GamTable(15)=	0.00193726
    GamTable(16)=	0.00194449
    GamTable(17)=	0.00195216
    GamTable(18)=	0.00196028
    GamTable(19)=	0.00196892
    GamTable(20)=	0.00197812
    GamTable(21)=	0.00198794
    GamTable(22)=	0.00199845
    GamTable(23)=	0.00200971
    GamTable(24)=	0.00202183
    GamTable(25)=	0.00203491
    GamTable(26)=	0.00204907
    GamTable(27)=	0.00206444
    GamTable(28)=	0.0020812
    GamTable(29)=	0.00209955
    GamTable(30)=	0.00211972
    GamTable(31)=	0.00214202
    GamTable(32)=	0.0021668
    GamTable(33)=	0.00219451
    GamTable(34)=	0.00222572
    GamTable(35)=	0.00226112
    GamTable(36)=	0.00230166
    GamTable(37)=	0.00234853
    GamTable(38)=	0.00240337
    GamTable(39)=	0.00246841
    GamTable(40)=	0.0025468
    GamTable(41)=	0.00264314
    GamTable(42)=	0.00276437
    GamTable(43)=	0.00292155
    GamTable(44)=	0.0031333
    GamTable(45)=	0.0034336
    GamTable(46)=	0.00389128
    GamTable(47)=	0.00466899
    GamTable(48)=	0.00625619
    GamTable(49)=	0.01098453
    GamTable(50)=	0.06479965

    R_effective = (R_Load - two)/(one - mode_mixity*0.42d0) + two  ! Eq. 25 from NASA/TP-2018-219838

    iFound = .FALSE.
    Lookup: Do i=2,50
      If (R_effective > RTable(i-1) .AND. R_effective <= RTable(i)) Then
        interpolation = (R_effective - RTable(i-1))/(RTable(i) - RTable(i-1))
        Bet = BetTable(i-1) + interpolation*(BetTable(i) - BetTable(i-1))
        Gam = GamTable(i-1) + interpolation*(GamTable(i) - GamTable(i-1))
        iFound = .TRUE.
        Exit Lookup
      End If
    End Do Lookup
    If (.NOT. iFound) Then
      Call log%error("Fatigue coefficients not found. Mode mix:" // trim(str(mode_mixity)) // "R ratio:" // trim(str(R_Load)))
    End If

    Return
  End Subroutine fatigue_coefficients

End Module cohesive_mod
