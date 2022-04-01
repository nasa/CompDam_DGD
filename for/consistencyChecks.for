Module consistency_Mod
  ! Module for input consistency checks

Contains

  Subroutine consistencyChecks(m, p, issueWarnings)
    ! Checks that a consistent set of properties has been defined

    Use forlog_Mod
    Use matProp_Mod
    Use parameters_Mod
    Use rambergOsgood_Mod

    ! Arguments
    Type(matProps), intent(INOUT) :: m
    Type(parameters), intent(IN) :: p
    Logical, intent(IN) :: issueWarnings

    ! Locals
    Double Precision :: eps_0tanstiff, eps_0_fn, ro_constant_slope(2)
    Double Precision, Parameter :: zero=0.d0, one=1.d0, two=2.d0, four=4.d0
    ! -------------------------------------------------------------------- !

    If (.NOT. m%cohesive) Then
      ! Check that all elastic properties have been defined for transverse isotropy
      If (.NOT. m%E1_def) Call log%error('PROPERTY ERROR: Must define a value for E1')
      If (.NOT. m%E2_def) Call log%error('PROPERTY ERROR: Must define a value for E2')
      If (.NOT. m%G12_def) Call log%error('PROPERTY ERROR: Must define a value for G12')
      If (.NOT. m%v12_def) Call log%error('PROPERTY ERROR: Must define a value for v12')
      If (.NOT. m%v23_def) Call log%error('PROPERTY ERROR: Must define a value for v23')

      ! Check if orthotropic elastic properties have been defined
      If (m%E3_def .OR. m%G13_def .OR. m%G23_def) Then
        If (.NOT. m%E3_def) Call log%error('PROPERTY ERROR: Some orthotropic elastic properties are missing. Must define a value for E3.')
        If (.NOT. m%G13_def) Call log%error('PROPERTY ERROR: Some orthotropic elastic properties are missing. Must define a value for G13.')
        If (.NOT. m%G23_def) Call log%error('PROPERTY ERROR: Some orthotropic elastic properties are missing. Must define a value for G23.')
        If (.NOT. m%v13_def) Call log%error('PROPERTY ERROR: Some orthotropic elastic properties are missing. Must define a value for v13.')
        Call log%info('PROPERTY: Orthotropic constants are defined directly')
      Else
        ! Compute these properties assuming transverse isotropy
        Call log%info('PROPERTY: Assuming transverse isotropy')
        m%E3  = m%E2
        m%G13 = m%G12
        m%v13 = m%v12
        m%G23 = m%E2/two/(one + m%v23)
      End If

    Else If (m%cohesive) Then  ! Check if cohesive element properties have been specified
      ! Check for cohesive stiffness material properties
      If (.NOT. m%E3_def) Call log%error('PROPERTY ERROR: Cohesive material laws require a definition for E3.')
      ! Check for cohesive strength material properties
      If (.NOT. m%YT_def) Call log%error('PROPERTY ERROR: Some cohesive element properties are missing. Must define a value for YT.')
      If (.NOT. m%SL_def) Call log%error('PROPERTY ERROR: Some cohesive element properties are missing. Must define a value for SL.')
      If (m%YC_def .AND. m%alpha0_def) Then
        ! Compute these properties for the transverse shear strength
        m%etaL = -m%SL*COS(two*m%alpha0)/(m%YC*COS(m%alpha0)*COS(m%alpha0))
        m%etaT = -one/TAN(two*m%alpha0)
        m%ST   = m%YC*COS(m%alpha0)*(SIN(m%alpha0) + COS(m%alpha0)/TAN(two*m%alpha0))
      Else
        Call log%info('PROPERTY: YC and/or alpha0 not defined. Assuming that S_T = S_L.')
        m%etaL = zero
        m%etaT = zero
        m%ST   = m%SL
      End If
      ! Check for cohesive fracture toughness material properties
      If (.NOT. m%GYT_def) Call log%error('PROPERTY ERROR: Some cohesive element properties are missing. Must define a value for GYT.')
      If (.NOT. m%GSL_def) Call log%error('PROPERTY ERROR: Some cohesive element properties are missing. Must define a value for GSL.')
      If (.NOT. m%eta_BK_def) Call log%error('PROPERTY ERROR: Some cohesive element properties are missing. Must define a value for eta_BK.')

      Call log%info('PROPERTY: All required cohesive element peroperties are defined.')
    End IF

    ! Check if matrix damage properties have been specified
    If (m%matrixDam) Then
      If (.NOT. m%YT_def) Call log%error('PROPERTY ERROR: Some matrix damage properties are missing. Must define a value for YT.')
      If (.NOT. m%SL_def) Call log%error('PROPERTY ERROR: Some matrix damage properties are missing. Must define a value for SL.')
      If (.NOT. m%GYT_def) Call log%error('PROPERTY ERROR: Some matrix damage properties are missing. Must define a value for GYT.')
      If (.NOT. m%GSL_def) Call log%error('PROPERTY ERROR: Some matrix damage properties are missing. Must define a value for GSL.')
      If (.NOT. m%eta_BK_def) Call log%error('PROPERTY ERROR: Some matrix damage properties are missing. Must define a value for eta_BK.')
      If (.NOT. m%YC_def) Call log%error('PROPERTY ERROR: Some matrix damage properties are missing. Must define a value for YC.')
      If (.NOT. m%alpha0_def) Call log%error('PROPERTY ERROR: Some matrix damage properties are missing. Must define a value for alpha0.')
      m%etaL = -m%SL*COS(two*m%alpha0)/(m%YC*COS(m%alpha0)*COS(m%alpha0))
      m%etaT = -one/TAN(two*m%alpha0)
      m%ST   = m%YC*COS(m%alpha0)*(SIN(m%alpha0) + COS(m%alpha0)/TAN(two*m%alpha0))

      Call log%info('PROPERTY: Matrix damage is enabled')
    Else
      Call log%info('PROPERTY: Matrix damage is disabled')
    End If

    ! Check if CTEs have been defined
    If (m%cte_def(1) .OR. m%cte_def(2)) Then
      If (.NOT. m%cte_def(1)) Call log%error('PROPERTY ERROR: Some CTE properties are missing. Must define a value for alpha11.')
      If (.NOT. m%cte_def(2)) Call log%error('PROPERTY ERROR: Some CTE properties are missing. Must define a value for alpha22.')
      Call log%info('PROPERTY: CTEs have been defined')
    Else
      If (issueWarnings) Call log%warn('PROPERTY: CTEs are being set to zero')
      m%cte(1) = zero
      m%cte(2) = zero
    End If

    ! Check if shear nonlinearity properties have been defined
    If (m%shearNonlinearity12 .OR. m%shearNonlinearity13) Then
      If (.NOT. m%aPL_def) Call log%error('PROPERTY ERROR: Some shear-nonlinearity properties are missing. Must define a value for alpha_PL.')
      If (.NOT. m%nPL_def) Call log%error('PROPERTY ERROR: Some shear-nonlinearity properties are missing. Must define a value for n_PL.')
      If (m%aPL .EQ. zero) Call log%error('PROPERTY ERROR: Shear nonlinearity is enabled, but alpha_PL is zero')
      If (m%nPL .EQ. zero) Call log%error('PROPERTY ERROR: Shear nonlinearity is enabled, but n_PL is zero')
      Call ramberg_osgood_cs_consts(m%G12, m%aPL, m%nPL, m%SL, p%shear_strain_const_slope, m%shear_stress_const_slope(2), m%shear_const_slope_coeff(2,:))
      Call ramberg_osgood_cs_consts(m%G13, m%aPL, m%nPL, m%ST, p%shear_strain_const_slope, m%shear_stress_const_slope(3), m%shear_const_slope_coeff(3,:))
      m%shear_strain_const_slope = p%shear_strain_const_slope
      Call log%info('PROPERTY: Shear-nonlinearity properties have been defined')
    Else
      Call log%info('PROPERTY: Shear-nonlinearity is disabled')
    End If

    ! Check if all Schapery micro-damage material properties have been defined
    If (m%schapery) Then
      If (.NOT. (all(m%es_def) .AND. all(m%gs_def))) Call log%error('PROPERTY ERROR: Some Schapery micro-damage are missing.')
      Call log%info('PROPERTY: Schapery micro-damage properties have been defined')
    Else
      Call log%info('PROPERTY: Schapery micro-damage is disabled')
      ! Default values that will cause no pre-peak nonlinearity
      m%es = zero
      m%es(1) = one
      m%gs = zero
      m%gs(1) = one
    End If

    ! Check is all the parameters for Schaefer theory have been defined
    If (m%schaefer) Then
      If (.NOT. m%schaefer_a6_def ) Call log%error('PROPERTY ERROR: Some schaefer theory properties are missing. Must define schaefer_a6')
      If (.NOT. m%schaefer_b2_def ) Call log%error('PROPERTY ERROR: Some schaefer theory properties are missing. Must define schaefer_b2')
      If (.NOT. m%schaefer_n_def ) Call log%error('PROPERTY ERROR: Some schaefer theory properties are missing. Must define schaefer_n')
      If (.NOT. m%schaefer_A_def ) Call log%error('PROPERTY ERROR: Some schaefer theory properties are missing. Must define schaefer_A')
      Call log%info('PROPERTY: Schaefer properties have been defined')
    Else
      Call log%info('PROPERTY: Schaefer is disabled')
      ! Default values that will cause no pre-peak nonlinearity
    End If

    ! Check if fiber tensile damage properties have been defined
    If (m%fiberTenDam) Then
      If (.NOT. m%XT_def) Call log%error('PROPERTY ERROR: Some fiber tensile damage properties are missing. Must define a value for XT.')
      If (.NOT. m%fXT_def) Call log%error('PROPERTY ERROR: Some fiber tensile damage properties are missing. Must define a value for fXT.')
      If (.NOT. m%GXT_def) Call log%error('PROPERTY ERROR: Some fiber tensile damage properties are missing. Must define a value for GXT.')
      If (.NOT. m%fGXT_def) Call log%error('PROPERTY ERROR: Some fiber tensile damage properties are missing. Must define a value for fGXT.')

      Call log%info('PROPERTY: fiber tensile damage properties have been defined')
    Else
      Call log%info('PROPERTY: fiber tensile damage is disabled')
    End If

    ! Check if the CDM fiber compression damage properties have been defined
    If (m%fiberCompDamBL) Then
      If (.NOT. m%XC_def) Call log%error('PROPERTY ERROR: Some fiber compression damage properties are missing (BL model). Must define a value for XC.')
      If (.NOT. m%fXC_def) Call log%error('PROPERTY ERROR: Some fiber compression damage properties are missing (BL model). Must define a value for fXC.')
      If (.NOT. m%GXC_def) Call log%error('PROPERTY ERROR: Some fiber compression damage properties are missing (BL model). Must define a value for GXC.')
      If (.NOT. m%fGXC_def) Call log%error('PROPERTY ERROR: Some fiber compression damage properties are missing (BL model). Must define a value for fGXC.')

      Call log%info('PROPERTY: fiber compression damage BL (model 1) properties have been defined')
    Else
      Call log%info('PROPERTY: fiber compression damage BL (model 1) is disabled')
    End If

    ! Check if the DGD fiber compression damage properties have been defined
    If (m%fiberCompDamFKT12 .OR. m%fiberCompDamFKT13) Then
      If (.NOT. m%XC_def) Call log%error('PROPERTY ERROR: Some fiber compression damage properties are missing (FKT model). Must define a value for XC.')
      If (.NOT. m%YC_def) Call log%error('PROPERTY ERROR: Some fiber compression damage properties are missing (FKT model). Must define a value for YC.')
      If (.NOT. m%w_kb_def) Call log%error('PROPERTY ERROR: Some fiber compression damage properties are missing (FKT model). Must define a value for w_kb.')
      If (.NOT. m%alpha0_def) Call log%error('PROPERTY ERROR: Some fiber compression damage properties are missing (FKT model). Must define a value for alpha0.')
      If (.NOT. m%SL_def) Call log%error('PROPERTY ERROR: Some fiber compression damage properties are missing (FKT model). Must define a value for SL.')

      ! Make sure shear nonlinearity is enabled
      If (m%fiberCompDamFKT12 .AND. (.NOT. m%shearNonlinearity12)) Call log%error('PROPERTY ERROR: Shear-nonlinearity 1-2 must be enabled with fiber compression damage FKT 1-2.')
      If (m%fiberCompDamFKT13 .AND. (.NOT. m%shearNonlinearity13)) Call log%error('PROPERTY ERROR: Shear-nonlinearity 1-3 must be enabled with fiber compression damage FKT 1-3.')

      If (m%fiberCompDamFKT12 .AND. m%fiberCompDamFKT13) Call log%info('PROPERTY: fiber compression damage FKT 3-D (model 5) properties have been defined')
      If (m%fiberCompDamFKT12) Call log%info('PROPERTY: fiber compression damage FKT in-plane (model 3) properties have been defined')
      If (m%fiberCompDamFKT13) Call log%info('PROPERTY: fiber compression damage FKT out-of-plane (model 4) properties have been defined')
    Else
      Call log%info('PROPERTY: fiber compression damage FKT is disabled')
    End If

    ! check if fiber nonlinearity has been defined
    If (m%cl_def) Then
      ! Check that the tangent stiffness does not go to zero before fiber failure
      If (m%fiberCompDamBL) Then
        eps_0tanstiff = one/(two*m%cl)
        eps_0_fn = (-m%E1+SQRT(m%E1**two-four*m%E1*m%cl*m%XC))/(two*m%E1*m%cl)  ! Strain at fiber compression failure
        If (eps_0tanstiff < eps_0_fn) Then
          Call log%error('PROPERTY ERROR: Tangent stiffness goes (eps11 =  ' // trim(str(eps_0tanstiff)) // ')to zero before fiber compression failure strain is reached(eps11 = ' // &
            trim(str(eps_0_fn)) // '). Adjust cl, E1, and/or XC.')
        End If
      End If
      Call log%info('PROPERTY: fiber nonlinearity has been defined')
    Else
      Call log%info('PROPERTY: fiber nonlinearity has not been defined. Setting cl = 0')
      m%cl = zero
    End If

    ! check if fatigue properties have been defined
    If (m%fatigue_gamma_def .AND. m%fatigue_epsilon_def .AND. m%fatigue_eta_def .AND. m%fatigue_p_mod_def) Then
      Call log%info('PROPERTY: cohesive fatigue properties have been defined')
    Else
      Call log%info('PROPERTY: cohesive fatigue properties have not been defined. Using default values.')
      If (.NOT. m%fatigue_gamma_def) m%fatigue_gamma = 1.d7
      If (.NOT. m%fatigue_epsilon_def) m%fatigue_epsilon = 0.2d0
      If (.NOT. m%fatigue_eta_def) m%fatigue_eta = 0.95d0
      If (.NOT. m%fatigue_p_mod_def) m%fatigue_p_mod = zero
    End If

    Return
  End Subroutine consistencyChecks

End Module consistency_Mod
