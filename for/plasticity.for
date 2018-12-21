Module plasticity_mod
  ! Module for all plasticity-related calculations

Contains

  Subroutine Plasticity(m, sv, p, ndir, nshr, eps, eps_old_in, use_temp)
    ! This is the main entry point for plasticity calculations

    Use matProp_Mod
    Use stateVar_Mod
    Use schaefer_Mod
    use matrixAlgUtil_Mod
    Use parameters_Mod

    ! Arguments
    Type(matProps), intent(IN) :: m
    Type(stateVars), intent(INOUT) :: sv                ! State variables related to plasticity
    Type(parameters), intent(IN) :: p                    ! parameters object (contains information like convergence tolerance which can dictate behavior of subroutines)
    Integer, intent(IN) :: ndir, nshr
    Double Precision, intent(INOUT) :: eps(ndir,ndir)    ! Strain tensor
    Double Precision, optional, intent(IN) :: eps_old_in(ndir, ndir)   ! Strain tensor
    Logical, optional :: use_temp                       ! Set to true to use temporary state variables (default=False)

    ! Locals
    Logical :: use_temporary_sv
    Double Precision:: eps_old(ndir, ndir)   ! Strain tensor

    Double Precision, parameter :: zero=0.d0, one=1.d0, two=2.d0
    ! -------------------------------------------------------------------- !

    ! Set default behavior for use_temp
    If (present(use_temp)) Then
      use_temporary_sv = use_temp
    Else
      use_temporary_sv = .FALSE.
    End If
    if (PRESENT(eps_old_in)) THEN
      eps_old = eps_old_in
    ELSE 
      eps_old = zero
    End If

    If (m%shearNonlinearity12) Then
      If (use_temporary_sv) Then
        Call ro_plasticity(two*eps(1,2), sv%d_eps12_temp, m%G12, m%aPL, m%nPL, sv%Inel12c, sv%Plas12_temp, sv%Inel12_temp)
        eps(1,2) = eps(1,2) - sv%Plas12_temp/two
        eps(2,1) = eps(1,2)
      Else
        Call ro_plasticity(two*eps(1,2), sv%d_eps12, m%G12, m%aPL, m%nPL, sv%Inel12c, sv%Plas12, sv%Inel12)
        eps(1,2) = eps(1,2) - sv%Plas12/two
        eps(2,1) = eps(1,2)
      End IF
    Else If (m%shearNonlinearity13) Then
      If (use_temporary_sv) Then
        Call ro_plasticity(two*eps(1,3), sv%d_eps13_temp, m%G13, m%aPL, m%nPL, Huge(zero), sv%Plas13_temp, sv%Inel13_temp)
        eps(1,3) = eps(1,3) - sv%Plas13_temp/two
        eps(3,1) = eps(1,3)
      Else
        Call ro_plasticity(two*eps(1,3), sv%d_eps13, m%G13, m%aPL, m%nPL, Huge(zero), sv%Plas13, sv%Inel13)
        eps(1,3) = eps(1,3) - sv%Plas13/two
        eps(3,1) = eps(1,3)
      End IF
    Else If (m%schaefer) THEN
    	! Update Ep_schaefer and f (yield function)
      CALL schaefer(m, p, m%schaefer_a6,  m%schaefer_b2,  m%schaefer_n, m%schaefer_A, eps, eps_old, ndir, nshr, sv%Ep_schaefer, sv%fp)
      !updated eps subtracting out total plastic strain eps -= sv%Ep_old
      eps = eps - Vec2Matrix(sv%Ep_schaefer)
    End If

    Return
  End Subroutine Plasticity


  Subroutine ro_plasticity(strain, d_strain_sign, modulus, aPL, nPL, inel_max, Plas, Inel)
    ! The purpose of this subroutine is to calculate the plastic and inelastic strains
    ! for the given total strain state and nonlinear stress-strain curve.

    Use forlog_Mod

    ! Arguments
    Double Precision, intent(IN) :: modulus                 ! Elastic modulus (eg G12)
    Double Precision, intent(IN) :: strain, d_strain_sign   ! Strain
    Double Precision, intent(IN) :: aPL, nPL                ! NL shear strain properties
    Double Precision, intent(IN) :: inel_max                ! Maximum inelastic shear strain, for fiber failure criterion
    Double Precision, intent(INOUT) :: Plas, Inel

    ! Locals
    Double Precision :: epsEff,tau,Inel_trial,ff0,ff1,tol
    Integer :: E, E_max
    Double Precision, parameter :: zero=0.d0, one=1.d0
    ! -------------------------------------------------------------------- !

    ! effective strain (used to account for unloading and reloading)
    epsEff = ABS(strain - Plas) + Inel

    ! Compute stress following the Ramberg-Osgood curve
    tau = ramberg_osgood(epsEff, modulus, aPL, nPL)

    ! Compute the inelastic strain (independent of unloading/reloading) and plastic strain (depends on loading direction) 
    Inel_trial = epsEff - tau/Modulus
    If (Inel_trial > Inel) Then
      If (Inel_trial > inel_max) Then
        Inel_trial = inel_max
      End If
      Plas = Plas + (Inel_trial - Inel)*d_strain_sign
      Inel = Inel_trial
    End If

    Return
  End Subroutine ro_plasticity


  Function ramberg_osgood(strain, modulus, aPL, nPL) result(tau)
    ! Find shear stress corresponding to strain following RO law

    Use forlog_Mod

    ! Input
    Double Precision, intent(IN) :: strain        ! engineering shear strain
    Double Precision, intent(IN) :: modulus, aPL, nPL

    ! Output
    Double Precision :: tau

    ! Locals
    Double Precision :: epsEff,Inel_trial,ff0,ff1,tol
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
      If (E == E_max) Call log%error("epsPL loop in ramberg_osgood failed to converge: " // trim(str(strain)))
      ff1 = one + aPL*nPL*ABS(tau)**(nPL - one) ! derivative of ff0 w.r.t. tau
      tau = tau - ff0/ff1 ! Newton-Raphson equation
    End Do epsPL

    Return
  End Function ramberg_osgood


  Function intializeFiberFailure(phi0, phiff, G12, aPL, nPL)
    ! Determines the critical plastic strain for fiber failure
    ! Since inelastic strain is used, the direction is ignored and the absolute value is returned

    Use forlog_Mod

    ! Arguments
    Double Precision, intent(IN) :: phi0                         ! Initial fiber misalignment
    Double Precision, intent(IN) :: phiff                        ! Angle for fiber failure
    Double Precision, intent(IN) :: G12, aPL, nPL                ! Shear nonlinearity

    ! Output
    Double Precision :: intializeFiberFailure

    ! Locals
    Double Precision :: phiff_radians, gamma12c, tau
    Double Precision, parameter :: one=1.d0
    ! -------------------------------------------------------------------- !

    ! Convert to radians
    phiff_radians = (ATAN(one)/45.d0)*phiff

    ! Assume the same direction as phi0
    phiff_radians = SIGN(phiff_radians, phi0)

    ! Compute the critical plastic strain for fiber failure
    gamma12c = phiff_radians - phi0
    tau = ramberg_osgood(gamma12c, G12, aPL, nPL)
    tau = SIGN(phiff_radians, tau)  ! Account for sign of tau
    intializeFiberFailure = ABS(gamma12c - tau/G12)

    Return
  End Function intializeFiberFailure

End Module plasticity_mod