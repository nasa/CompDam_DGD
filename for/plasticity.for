Module plasticity_mod
  ! Main entry point for all plasticity-related calculations

Contains

  Subroutine Plasticity(m, sv, p, ndir, nshr, eps, eps_old, enerPlas, use_temp)
    ! This is the main entry point for plasticity calculations

    Use matProp_Mod
    Use stateVar_Mod
    Use rambergOsgood_Mod
    Use schaefer_Mod
    use matrixAlgUtil_Mod
    Use parameters_Mod

    ! Arguments
    Type(matProps), intent(IN) :: m
    Type(stateVars), intent(INOUT) :: sv                ! State variables related to plasticity
    Type(parameters), intent(IN) :: p                    ! parameters object (contains information like convergence tolerance which can dictate behavior of subroutines)
    Integer, intent(IN) :: ndir, nshr
    Double Precision, intent(INOUT) :: eps(ndir,ndir)    ! Strain tensor
    Double Precision, intent(IN) :: eps_old(ndir, ndir)   ! Strain tensor, last increment
    Double Precision, intent(OUT) :: enerPlas             ! Increment of dissipated plastic energy
    Logical, optional :: use_temp                       ! Set to true to use temporary state variables (default=False)

    ! Locals
    Logical :: use_temporary_sv
    Double Precision :: enerPlas_temp                   ! Temporary dissipated plastic energy

    Double Precision, parameter :: zero=0.d0, one=1.d0, two=2.d0
    ! -------------------------------------------------------------------- !

    enerPlas = zero

    ! Set default behavior for use_temp
    If (present(use_temp)) Then
      use_temporary_sv = use_temp
    Else
      use_temporary_sv = .FALSE.
    End If

    If (m%shearNonlinearity12) Then
      If (use_temporary_sv) Then
        Call ro_plasticity(m, 2, two*eps(1,2), two*eps_old(1,2), sv%Inel12c, sv%Plas12_temp, sv%Inel12_temp, enerPlas_temp)
        eps(1,2) = eps(1,2) - sv%Plas12_temp/two
        eps(2,1) = eps(1,2)
      Else
        Call ro_plasticity(m, 2, two*eps(1,2), two*eps_old(1,2), sv%Inel12c, sv%Plas12, sv%Inel12, enerPlas_temp)
        eps(1,2) = eps(1,2) - sv%Plas12/two
        eps(2,1) = eps(1,2)
      End IF
      enerPlas = enerPlas_temp
    End If
    If (m%shearNonlinearity13) Then
      If (use_temporary_sv) Then
        Call ro_plasticity(m, 3, two*eps(1,3), two*eps_old(1,3), sv%Inel13c, sv%Plas13_temp, sv%Inel13_temp, enerPlas_temp)
        eps(1,3) = eps(1,3) - sv%Plas13_temp/two
        eps(3,1) = eps(1,3)
      Else
        Call ro_plasticity(m, 3, two*eps(1,3), two*eps_old(1,3), sv%Inel13c, sv%Plas13, sv%Inel13, enerPlas_temp)
        eps(1,3) = eps(1,3) - sv%Plas13/two
        eps(3,1) = eps(1,3)
      End If
      enerPlas = enerPlas_temp + enerPlas
    End If

    If (m%schaefer) Then
    	! Update Ep_schaefer and f (yield function)
      Call schaefer(m, p, m%schaefer_a6,  m%schaefer_b2,  m%schaefer_n, m%schaefer_A, eps, eps_old, ndir, nshr, sv%Ep_schaefer, sv%fp)
      !updated eps subtracting out total plastic strain eps -= sv%Ep_old
      eps = eps - Vec2Matrix(sv%Ep_schaefer)
    End If

    Return
  End Subroutine Plasticity


  Function intializeFiberFailure(phi0, phiff, G12, aPL, nPL)
    ! Determines the critical plastic strain for fiber failure
    ! Since inelastic strain is used, the direction is ignored and the absolute value is returned

    Use forlog_Mod
    Use rambergOsgood_Mod

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
    tau = ramberg_osgood_nocs(gamma12c, G12, aPL, nPL)
    tau = SIGN(phiff_radians, tau)  ! Account for sign of tau
    intializeFiberFailure = ABS(gamma12c - tau/G12)

    Return
  End Function intializeFiberFailure

End Module plasticity_mod