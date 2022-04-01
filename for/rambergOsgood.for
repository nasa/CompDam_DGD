Module rambergOsgood_Mod
  ! Ramberg-Osgood model for nonlinear shear

Contains

	Subroutine ro_plasticity(m, shear_component, strain, strainOld, inel_max, Plas, Inel, enerPlas)
    ! The purpose of this subroutine is to calculate the plastic and inelastic strains
    ! for the given total strain state and nonlinear stress-strain curve.

    Use forlog_Mod
    Use matProp_Mod

    ! Arguments
    Type(matProps), intent(IN) :: m
    Integer, intent(IN) :: shear_component                  ! Shear component to evaluate (either 2 or 3 for in-plane or out-of-plane respectively)
    Double Precision, intent(IN) :: strain, strainOld       ! Strain
    Double Precision, intent(IN) :: inel_max                ! Maximum inelastic shear strain, for fiber failure criterion
    Double Precision, intent(INOUT) :: Plas, Inel
    Double Precision, optional, intent(OUT) :: enerPlas             ! Increment of dissipated plastic energy

    ! Locals
    Double Precision :: epsEff,tau,Inel_trial,ff0,ff1,tol
    Double Precision :: d_strain_sign
    Double Precision :: epsEffOld, tauOld
    Double Precision :: modulus
    Integer :: E, E_max
    Double Precision, parameter :: zero=0.d0, one=1.d0
    ! -------------------------------------------------------------------- !

    modulus = get_modulus(m, shear_component)

    ! Sign of the change in strain
    d_strain_sign = Sign(one, (strain - strainOld))

    ! effective strain (used to account for unloading and reloading)
    epsEffOld = ABS(strainOld - Plas) + Inel
    epsEff = ABS(strain - Plas) + Inel

    ! Compute stress following the Ramberg-Osgood curve
    tauOld = ramberg_osgood(m, epsEffOld, shear_component)
    tau = ramberg_osgood(m, epsEff, shear_component)

    ! Compute the inelastic strain (independent of unloading/reloading) and plastic strain (depends on loading direction) 
    Inel_trial = epsEff - tau/modulus
    If (Inel_trial > Inel) Then
      If (Inel_trial > inel_max) Then
        Inel_trial = inel_max
      End If
      Plas = Plas + (Inel_trial - Inel)*d_strain_sign
      Inel = Inel_trial
      If (present(enerPlas)) Then
        enerPlas = m%aPL/modulus*tau**(m%nPL + one)*m%nPL/(m%nPL + one) - m%aPL/modulus*tauOld**(m%nPL + one)*m%nPL/(m%nPL + one)
      End If
    Else
      If (present(enerPlas)) Then
        enerPlas = zero
      End If
    End If

    Return
  End Subroutine ro_plasticity


  Function ramberg_osgood(m, strain, shear_component) result(tau)
    ! Find shear stress corresponding to strain following RO law, respecting the constant slope transition
    ! strain is the "effective strain" -- load reversals are considered elsewhere

    Use forlog_Mod
    Use matProp_Mod
    Use matrixAlgUtil_Mod

    ! Input
    Type(matProps), intent(IN) :: m
    Double Precision, intent(IN) :: strain    ! Engineering shear strain
    Integer, intent(IN) :: shear_component    ! Shear component to evaluate (either 2 or 3 for in-plane or out-of-plane respectively)

    ! Output
    Double Precision :: tau

    ! Locals
    Double Precision :: epsEff, modulus
    Double Precision, parameter :: zero=0.d0, one=1.d0, two=2.d0
    ! -------------------------------------------------------------------- !

    epsEff = ABS(strain) ! effective strain
    modulus = get_modulus(m, shear_component)
    
    If (epsEff < m%shear_strain_const_slope) Then
      tau = ramberg_osgood_nocs(strain, modulus, m%aPL, m%nPL)
    Else
      tau = m%shear_const_slope_coeff(shear_component,1)*epsEff + m%shear_const_slope_coeff(shear_component,2)
    End If
    tau = SIGN(one, strain)*tau

    Return
  End Function ramberg_osgood


  Function ramberg_osgood_nocs(strain, modulus, aPL, nPL) result(tau)
    ! Find shear stress corresponding to strain following RO law
    ! Note: does not consider the constant slope transition

    Use forlog_Mod

    ! Input
    Double Precision, intent(IN) :: strain        ! engineering shear strain
    Double Precision, intent(IN) :: modulus, aPL, nPL

    ! Output
    Double Precision :: tau

    ! Locals
    Double Precision :: epsEff,ff0,ff1,tol
    Integer :: E, E_max
    Double Precision, parameter :: zero=0.d0, one=1.d0
    ! -------------------------------------------------------------------- !

    epsEff = ABS(strain) ! effective strain
    tau = Modulus*strain ! initial guess for stress

    tol = 1.d-8 ! tolerance for epsPL convergence [strain]

    E_max = 1000 ! maximum epsPL loop iterations
    epsPL: Do E = 0,E_max
      ff0 = tau + SIGN(one, tau)*aPL*ABS(tau)**nPL - Modulus*epsEff ! ff0 = 0, solve for tau
      If (ABS(ff0)/Modulus <= tol) EXIT epsPL
      If (E == E_max) Call log%error("epsPL loop in ramberg_osgood_nocs failed to converge: " // trim(str(strain)))
      ff1 = one + aPL*nPL*ABS(tau)**(nPL - one) ! derivative of ff0 w.r.t. tau
      tau = tau - ff0/ff1 ! Newton-Raphson equation
    End Do epsPL

    Return
  End Function ramberg_osgood_nocs


  Function ramberg_osgood_d(m, stress, shear_component, order, use_constant_slope) result(derivative)
    ! Find the derivative of the ramberg osgood curve (change in stress per 
    ! change in strain) at the specified stress considering the constant slope transition

    Use forlog_Mod
    Use matProp_Mod

    ! Input
    Type(matProps), intent(IN) :: m
    Double Precision, intent(IN) :: stress    ! Shear stress
    Integer, intent(IN) :: shear_component    ! Shear component to evaluate (either 2 or 3 for in-plane or out-of-plane respectively)
    Integer, intent(IN) :: order              ! Order of the derivative
    Logical, intent(IN) :: use_constant_slope

    ! Output
    Double Precision :: derivative
    Logical :: dont_use_cs

    ! Locals
    Double Precision :: shear_modulus
    Double Precision, parameter :: zero=0.d0, one=1.d0, two=2.d0, three=3.d0
    ! -------------------------------------------------------------------- !

    If (ABS(stress) < m%shear_stress_const_slope(shear_component) .OR. (.NOT. use_constant_slope)) Then
      ! Compute derivative
      shear_modulus = get_modulus(m, shear_component)
      derivative = ramberg_osgood_d_nocs(stress, shear_modulus, m%aPL, m%nPL, order)
    Else
      ! Use constant slope
      If (order == 1) Then
        derivative = m%shear_const_slope_coeff(shear_component,1)
      Else If (order == 2) Then
        derivative = zero
      Else
        Call log%terminate("ramberg_osgood_d() received invalid argument order: " // trim(str(order)))
      End If
    End If

    Return
  End Function ramberg_osgood_d


  Function ramberg_osgood_d_nocs(stress, modulus, aPL, nPL, order) result(derivative)
    ! Find the derivative of the ramberg osgood curve (change in stress per 
    ! change in strain) at the specified stress
    ! Note: does not consider the constant slope transition

    Use forlog_Mod

    ! Input
    Double Precision, intent(IN) :: stress
    Double Precision, intent(IN) :: modulus, aPL, nPL
    Integer, intent(IN) :: order

    ! Output
    Double Precision :: derivative

    ! Locals
    Double Precision :: fp,fpp
    Double Precision, parameter :: zero=0.d0, one=1.d0, two=2.d0, three=3.d0
    ! -------------------------------------------------------------------- !
    
    fp = (one + aPL*nPL*ABS(stress)**(nPL-one))/modulus
    If (order==1) Then
      derivative = 1/fp
    Else If (order==2) Then
      fpp = SIGN(one, stress)*aPL*nPL*(nPL-one)*ABS(stress)**(nPL-two)/modulus
      derivative = -fpp/(fp**three)
    Else
      Call log%terminate("ramberg_osgood_d_nocs received invalid argument order: " // trim(str(order)))
    End If

    Return
  End Function ramberg_osgood_d_nocs


  Subroutine ramberg_osgood_cs_consts(modulus, aPL, nPL, SL, strain_at_constant_slope, stress_at_constant_slope, constants)
    ! Find the constants that define the linear law once the strain_at_constant_slope
    ! threshold is reached. This is an initialization calculation that only needs to
    ! called once.

    Use forlog_Mod

    ! Input
    Double Precision, intent(IN) :: modulus, aPL, nPL, SL
    Double Precision, intent(IN) :: strain_at_constant_slope

    ! Output
    Double Precision, intent(OUT) :: stress_at_constant_slope
    Double Precision, intent(OUT) :: constants(2)  ! (slope, intercept)

    ! Locals
    Double Precision :: slope, intercept
    Double Precision, parameter :: zero=0.d0, one=1.d0, two=2.d0, three=3.d0
    ! -------------------------------------------------------------------- !

    ! Slope
    stress_at_constant_slope = ramberg_osgood_nocs(strain_at_constant_slope, modulus, aPL, nPL)
    If (stress_at_constant_slope < SL) Call log%warn('The shear curve has a constant slope before reaching SL, consider increasing the parameter `shear_strain_const_slope`')
    slope = ramberg_osgood_d_nocs(stress_at_constant_slope, modulus, aPL, nPL, 1)

    ! Intercept
    intercept = stress_at_constant_slope - slope*strain_at_constant_slope

    constants(1) = slope
    constants(2) = intercept

    Return
  End Subroutine ramberg_osgood_cs_consts


  Function get_modulus(m, shear_component) result(modulus)
    ! Get the modulus property corresponding to the shear component

    Use forlog_Mod
    Use matProp_Mod

    ! Input
    Type(matProps), intent(IN) :: m
    Integer, intent(IN) :: shear_component

    ! Output
    Double Precision :: modulus
    ! -------------------------------------------------------------------- !

    If (shear_component == 2) Then
      modulus = m%G12
    Else If (shear_component == 3) Then
      modulus = m%G13
    Else
      Call log%terminate('get_modulus() received invalid value for argument shear_component: ' // trim(str(shear_component)))
    End If

    Return
  End Function get_modulus

End Module rambergOsgood_Mod
