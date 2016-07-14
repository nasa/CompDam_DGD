Module strain_mod
  ! Module for all strain-related calculations

Contains

  Subroutine Strains(m,F,U,DT,ndir,eps,Plas12,Inel12,d_eps12_sign)
    ! The purpose of this subroutine is to calculate the strains (eps, plas12, inel12)
    ! for the given deformation gradient and temperature.

    Use matProp_Mod

    Include 'vaba_param.inc'

    ! Arguments
    Type(matProps), intent(IN) :: m
    Double Precision, intent(IN) :: F(3,3), U(3,3)            ! Deformation gradient, stretch tensor
    Double Precision, intent(IN) :: DT                        ! Coefficients of thermal expansion and temperature change
    Double Precision, intent(IN) :: d_eps12_sign              ! Sign of the change in shear strain from the previous increment to this one
    Integer, intent(IN) :: ndir
    Double Precision, intent(OUT) :: eps(ndir,ndir)           ! GL stain tensor, without plasticity
    Double Precision, intent(INOUT) :: Plas12, Inel12         ! State variables (CDM_Plas12, CDM_Inel12)

    ! Locals
    Parameter (zero=0.d0, two=2.d0)
    ! -------------------------------------------------------------------- !

    ! Initialize strain tensor
    eps = zero

    ! Calculate the specified definition of strain
    If (m%strainDef .EQ. 1) Then  ! Log strain
      eps = LogStrain(U,ndir)
    Else If (m%strainDef .EQ. 2) Then  ! GL strain
      eps = GLStrain(F,ndir)
    Else If (m%strainDef .EQ. 3) Then  ! Biot strain
      eps = BiotStrain(U,ndir)
    End If

    ! Account for thermal strains
    If (DT .NE. zero) Then
      eps(1,1) = eps(1,1) - m%cte(1)*DT
      eps(2,2) = eps(2,2) - m%cte(2)*DT
      eps(3,3) = eps(3,3) - m%cte(3)*DT
    End If

    ! Evaluate shear plasticity and remove plastic strain from eps
    If (m%shearNonlinearity) Then
      Call NLSHEAR(two*eps(1,2), d_eps12_sign, m%G12, m%aPL, m%nPL, Plas12, Inel12)
      eps(1,2) = eps(1,2) - Plas12/two
      eps(2,1) = eps(1,2)
    End If

    Return
  End Subroutine Strains


  Function LogStrain(U,ndir)
    ! Computes the log strain from the stretch

    Include 'vaba_param.inc'

    ! Input
    Double Precision, intent(IN) :: U(3,3)
    Integer, intent(IN) :: ndir

    ! Output
    Double Precision :: LogStrain(ndir,ndir)

    ! Locals
    Double precision :: eigVec(3,3), eigVal(3), work(1000), eigDiag(3,3)
    Integer :: ipiv(3), info, LWORK
    Double Precision, parameter :: zero=0.d0
    ! -------------------------------------------------------------------- !

    ! Get the eigenvalues and eigenvectors of stretch
    ! DSYEV is a LAPACK function; doc: http://www.netlib.org/lapack/double/dsyev.f
    eigVec = U  ! prevent overwrite
    ! Setup workspace
    LWORK = -1
    Call DSYEV( 'V', 'U', 3, eigVec, 3, eigVal, WORK, LWORK, INFO )
    LWORK = MIN( 1000, INT( WORK( 1 ) ) )
    ! Solve the eigenvalue problem
    Call DSYEV( 'V', 'U', 3, eigVec, 3, eigVal, WORK, LWORK, INFO )
    If (info .NE. 0) Then
      print *, 'WARNING'
      print *, 'Failed to compute eigenvalues of U^2. DSYEV Error.'
    End If

    ! Format log of eigenvalues as a square matrix
    eigDiag = zero
    eigDiag(1,1) = LOG(eigVal(1))
    eigDiag(2,2) = LOG(eigVal(2))
    eigDiag(3,3) = LOG(eigVal(3))

    ! Get the logarithmic strain tensor
    LogStrain = MATMUL(eigVec, MATMUL(eigDiag, TRANSPOSE(eigVec)))

    Return
  End Function LogStrain


  Pure Function GLStrain(F,ndir)
    ! Computes the green lagrange strain from the deformation gradient

    Include 'vaba_param.inc'

    ! Input
    Double Precision, intent(IN) :: F(3,3)
    Integer, intent(IN) :: ndir

    ! Output
    Double Precision :: GLStrain(ndir,ndir)

    ! Locals
    Double Precision :: eye(ndir,ndir)  ! Identity
    Double Precision, parameter :: zero=0.d0, one=1.d0, half=0.5d0
    ! -------------------------------------------------------------------- !

    ! Initialize identity matrix
    eye = zero; Do I = 1,3; eye(I,I) = one; End Do

    GLStrain = (MATMUL(TRANSPOSE(F), F) - eye)*half

    Return
  End Function GLStrain


  Pure Function BiotStrain(U,ndir)
    ! Computes the biot strain from the stetch

    Include 'vaba_param.inc'

    ! Input
    Double Precision, intent(IN) :: U(3,3)
    Integer, intent(IN) :: ndir

    ! Output
    Double Precision :: BiotStrain(ndir,ndir)

    ! Locals
    Double Precision :: eye(ndir,ndir)  ! Identity
    Double Precision, parameter :: zero=0.d0, one=1.d0
    ! -------------------------------------------------------------------- !

    ! Initialize identity matrix
    eye = zero; Do I = 1,3; eye(I,I) = one; End Do

    BiotStrain = U - eye

    Return
  End Function BiotStrain


  Subroutine NLSHEAR(strain,d_strain_sign,Modulus,aPL,nPL,Plas,Inel)
    ! The purpose of this subroutine is to calculate the plastic and inelastic strains
    ! for the given total strain state and nonlinear stress-strain curve.

    Use forlog_Mod

    Include 'vaba_param.inc'

    ! Arguments
    Double Precision, intent(IN) :: Modulus                 ! Elastic modulus (eg G12)
    Double Precision, intent(IN) :: strain, d_strain_sign   ! Strain
    Double Precision, intent(IN) :: aPL, nPL                ! NL shear strain properties
    Double Precision, intent(INOUT) :: Plas, Inel

    ! Locals
    Double Precision :: epsEff,tau,Inel_trial,ff0,ff1,tol!,EnergyDensity
    Integer :: E_max
    Double Precision, parameter :: one=1.d0
    ! -------------------------------------------------------------------- !

    epsEff = ABS(strain - Plas) + Inel ! effective strain
    tau = Modulus*(strain - Plas) ! initial guess for stress

    tol = 1.d-8 ! tolerance for epsPL convergence [strain]

    E_max = 1000 ! maximum epsPL loop iterations
    epsPL: Do E = 0,E_max
      ff0 = tau + SIGN(one, tau)*aPL*ABS(tau)**nPL - Modulus*epsEff ! ff0 = 0, solve for tau
      If (ABS(ff0)/Modulus .LE. tol) EXIT epsPL
      If (E .EQ. E_max) Call log%error("epsPL loop in NLSHEAR failed to converge: " // trim(str(Plas)) // ', ' // trim(str(Inel)) // ', ' // trim(str(strain)))
      ff1 = one + aPL*nPL*ABS(tau)**(nPL - one) ! derivative of ff0 w.r.t. tau
      tau = tau - ff0/ff1 ! Newton-Raphson equation
    End Do epsPL

    Inel_trial = epsEff - tau/Modulus
    If (Inel_trial .GT. Inel) Then
      Plas = Plas + (Inel_trial - Inel)*d_strain_sign
      Inel = Inel_trial
      !EnergyDensity = aPL/Modulus*tau**(nPL + one)*nPL/(nPL + one)
    End If

    Return
  End Subroutine NLSHEAR


End Module strain_mod
