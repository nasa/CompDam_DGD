Module DGD_Mod
  ! Module for all DGD calculations

Contains


  Subroutine DGDInit(U, F, m, p, sv, ndir, nshr, DT, Cauchy, enerIntern, F_old)
    ! Checks for the initiation of matrix damage, represented as a DGD
    ! cohesive crack. If the crack orientation is a priori unknown, it
    ! will be determined in this subroutine.

    Use forlog_Mod
    Use matrixAlgUtil_Mod
    Use matProp_Mod
    Use stateVar_Mod
    Use parameters_Mod
    Use stress_Mod
    Use cohesive_mod
    Use strain_mod
    Use plasticity_mod
    Use CDM_fiber_mod
    Use schapery_mod

    ! -------------------------------------------------------------------- !
    ! Arguments
    Type(matProps), intent(IN) :: m
    Type(parameters), intent(IN) :: p
    Type(stateVars), intent(INOUT) :: sv
    Double Precision, intent(IN) :: F(3,3), U(3,3), F_old(3,3)             ! Deformation gradient stretch tensor
    Integer, intent(IN) :: ndir
    Integer, intent(IN) :: nshr
    Double Precision, intent(IN) :: DT
    Double Precision, intent(OUT) :: Cauchy(ndir,ndir)                     ! Cauchy stress
    Double Precision, intent(OUT) :: enerIntern                            ! Internal energy

    ! -------------------------------------------------------------------- !
    ! Locals

    Double Precision :: Stiff(ndir+nshr,ndir+nshr)                         ! Stiffness
    Double Precision :: eps(ndir,ndir) !Strain
    Double Precision :: eps_old(ndir, ndir)                          ! Strain
    Double Precision :: stress(ndir,ndir)                                  ! Stress
    Double Precision :: F_inverse_transpose(3,3)                           ! Inverse transpose of the Deformation Gradient Tensor
    Double Precision :: X(3,3)                                             ! Reference configuration

    ! Cohesive surface
    Double Precision :: normal(3)                                          ! Normal vector (to cohesive surface)
    Double Precision :: R_cr(3,3)                                          ! Basis coordinate system for the cohesive surface
    Double Precision :: Pen(3)                                             ! Penalty stiffnesses
    Double Precision :: T(3)                                               ! Tractions on the cohesive surface
    Double Precision :: delta(3)                                           ! Current displacement jumps in crack coordinate system
    Double Precision :: B_temp, beta                                       ! Placeholder (temp.) variables for Mode-mixity
    Double Precision :: FIm_temp

    ! Matrix crack cohesive surface normal
    Double Precision :: alpha_temp                                         ! Current alpha (used in loop through possible alphas)
    Integer :: alpha_test, alphaQ                                          ! Alpha = normal to matrix crack cohesive surface
    Integer :: Q                                                           ! Flag: Q=2 for matrix crack; Q=3 for delamination
    Integer :: A, A_min, A_max                                             ! Range through which the code searches for alpha
    Integer :: alpha0_deg_2                                                ! negative symmetric alpha0 with normal in positive direction

    ! Fiber
    Double Precision :: FIfC                                               ! Fiber compression damage threshold
    Double Precision :: d1
    Double Precision :: normalDir(3)                                       ! Normal to the crack plane in the reference configuration
    Double Precision :: fiberDir(3)                                        ! Current fiber direction
    Double Precision :: pk2_fiberDir(3,3)                                  ! 2PK stress in the fiber direction
    Double Precision :: R_phi0(3,3)                                        ! Rotation to the misaligned frame
    Double Precision :: gamma_rphi0                                        ! Shear strain in the misaligned frame

    ! Miscellaneous
    Double Precision :: rad_to_deg, deg_to_rad
    Double Precision, parameter :: zero=0.d0, one=1.d0, two=2.d0
    Double Precision, parameter :: pert=0.0001d0                           ! Small perturbation used to compute a numerical derivative

    ! -------------------------------------------------------------------- !
    Call log%debug('Start of DGDInit')

    ! Miscellaneous constants
    rad_to_deg = 45.d0/ATAN(one)  ! Converts radians to degrees when multiplied
    deg_to_rad = one/rad_to_deg   ! Converts degrees to radians when multiplied

    ! Initialize outputs
    d1       = zero
    sv%d2    = zero
    sv%FIm   = zero
    sv%B     = zero
    sv%Fb1   = zero
    sv%Fb2   = zero
    sv%Fb3   = zero
    pk2_fiberDir = zero

    ! The sign of the change in shear strain, used in the shear nonlinearity subroutine. Previously was a state variable.
    If (m%shearNonlinearity12) Then
      sv%d_eps12 = Sign(one, (F(1,1)*F(1,2) + F(2,1)*F(2,2) + F(3,2)*F(3,1)) - (F_old(1,1)*F_old(1,2) + F_old(2,1)*F_old(2,2) + F_old(3,2)*F_old(3,1)))
    End If
    If (m%shearNonlinearity13) Then
      sv%d_eps13 = Sign(one, (F(1,1)*F(1,3) + F(2,1)*F(2,3) + F(3,1)*F(3,3)) - (F_old(1,1)*F_old(1,3) + F_old(2,1)*F_old(2,3) + F_old(3,2)*F_old(3,3)))
    End If

    ! Reference configuration
    X = zero; X(1,1) = sv%Lc(1); X(2,2) = sv%Lc(2); X(3,3) = sv%Lc(3)
    F_inverse_transpose = MInverse(TRANSPOSE(F))

    ! Compute the Green-Lagrange strain tensor: eps
    Call Strains(F, m, DT, ndir, eps)
    ! Compute the Green-Lagrange strain tensor: previous step
    Call Strains(F_old, m, DT, ndir, eps_old)

    ! Check fiber tension or fiber compression damage
    If (eps(1,1) >= -1d-6) Then    ! Fiber tension

      ! Set rfC for failure index output
      If (sv%rfC == one) Then
        sv%rfC = zero
      End If

      ! Compute the plastic strains and remove from the strain tensor
      Call Plasticity(m, sv, p, ndir, nshr, eps, eps_old, .FALSE.)

      ! Evaluate fiber tension failure criteria and damage variable
      If (m%fiberTenDam) Then
        Call FiberTenDmg(eps, ndir, m%E1, m%XT, m%GXT, m%fXT, m%fGXT, sv%Lc(1), m%cl, sv%rfT, sv%d1T, sv%d1C, sv%STATUS)
        Call log%debug('Computed fiber damage variable, d1T ' // trim(str(sv%d1T)))

        d1 = sv%d1T
      End If

      ! Build the stiffness matrix
      Call StiffFuncNL(m, ndir, nshr, d1, zero, zero, eps, Stiff, sv%Sr)

      ! Calculate stress
      stress = Hooke(Stiff, eps, nshr)
      Cauchy = convertToCauchy(stress, F)

    Else  ! Compression in 1-dir

      ! Set rfT for failure index output
      If (sv%rfT == one) Then
        sv%rfT = zero
      End If

      ! Check for fiber compression damage initiation
      If (m%fiberCompDamFKT) Then

        ! -------------------------------------------------------------------- !
        !    Compute stress in the material (considering phi0)                 !
        ! -------------------------------------------------------------------- !

        ! Rotation from reference frame to fiber misaligned frame
        fiberDir = (/cos(sv%phi0), sin(sv%phi0), zero/)
        normalDir = (/-sin(sv%phi0), cos(sv%phi0), zero/)
        R_phi0(:,1) = fiberDir
        R_phi0(:,2) = normalDir
        R_phi0(:,3) = (/zero, zero, one/)

        ! Transform strain to the fiber frame
        eps = MATMUL(TRANSPOSE(R_phi0), MATMUL(eps, R_phi0))

        ! Compute the plastic strains and remove from the strain tensor
        Call Plasticity(m, sv, p, ndir, nshr, eps, eps_old)

        ! Get total 1,2 strain component
        gamma_rphi0 = two*(eps(1,2) + sv%Plas12/two)

        ! Only decompose element if the plastic strain is nonnegligible and the kink band is smaller than the element size
        If (sv%Inel12 > 0.00001d0) Then
          If (m%w_kb/sv%Lc(1) < p%kb_decompose_thres) Then
            Call log%debug('DGDInit triggering DGDKinkband.')
            sv%d1C   = 1.d-6    ! Used as a flag to call DGDEvolve
            sv%alpha = zero     ! Assume an in-plane kink band
            sv%Fb1 = F(1,1)
            sv%Fb2 = F(2,1)
            sv%Fb3 = F(3,1)
          End If
        End If

        ! Compute the undamaged stiffness matrix
        Call StiffFuncNL(m, ndir, nshr, d1, zero, zero, eps, Stiff, sv%Sr)

        ! Calculate stress
        pk2_fiberDir = Hooke(Stiff, eps, nshr)  ! 2PK in the fiber direction
        stress = MATMUL(R_phi0, MATMUL(pk2_fiberDir, TRANSPOSE(R_phi0)))  ! 2PK rotated back to the reference direction
        Cauchy = convertToCauchy(stress, F)  ! Cauchy stress in the reference frame

        ! Failure index for fiber kinking
        sv%rfC = abs((-1*Cauchy(1,1)*cos(two*(sv%phi0+gamma_rphi0)))/((ramberg_osgood(gamma_rphi0 + pert, m%G12, m%aPL, m%nPL) - ramberg_osgood(gamma_rphi0 - pert, m%G12, m%aPL, m%nPL)) / (two * pert)))

        ! Rotation to the crack frame / fiber-aligned frame
        R_cr(:,1) = Norm(MATMUL(F, fiberDir))
        R_cr(:,2) = Norm(MATMUL(F_inverse_transpose, normalDir))
        R_cr(:,3) = CrossProduct(R_cr(:,1), R_cr(:,2))

        ! Calculate the angle that define the rotation of the fibers
        sv%gamma = ATAN(R_cr(2,1)/R_cr(1,1)) - sv%phi0

        d1 = zero

      Else

        ! Compute the plastic strains and remove from the strain tensor
        Call Plasticity(m, sv, p, ndir, nshr, eps, eps_old)

        If (m%fiberCompDamBL) Then
          Call FiberCompDmg(eps, ndir, m%E1, m%XC, m%GXC, m%fXC, m%fGXC, sv%Lc(1), m%cl, sv%rfT, sv%rfC, sv%d1T, sv%d1C, sv%STATUS)
          Call log%debug('Computed fiber damage variable, d1C ' // trim(str(sv%d1C)))

          d1 = sv%d1C

          ! Load reversal
          If (d1 > sv%d1T) Then
            sv%d1T = d1
          End If
        End If

        ! Build the stiffness matrix
        Call StiffFuncNL(m, ndir, nshr, d1, zero, zero, eps, Stiff, sv%Sr)

        ! Calculate stress
        stress = Hooke(Stiff, eps, nshr)
        Cauchy = convertToCauchy(stress, F)

      End If

    End If

    ! -------------------------------------------------------------------- !
    !    Search for matrix crack initiation only when no fiber damage has occurred
    ! -------------------------------------------------------------------- !
    If (m%matrixDam .AND. sv%d1T == zero .AND. sv%d1C == zero) Then
      ! Get fiber direction
      R_cr(:,1) = Norm(F(:,1)) ! fiber direction
      normal(1) = zero

      ! alphaQ is the angle between intralaminar and interlaminar oriented cracks
      alphaQ = FLOOR(ATAN(sv%Lc(2)/sv%Lc(3))*rad_to_deg)
      alphaQ = alphaQ - MOD(alphaQ, p%alpha_inc)

      ! Search through range of alphas to find the correct one (alpha=-999 is a flag to run this search)
      If (sv%alpha == -999) Then
        A_min = -alphaQ
        A_max = alphaQ

        If (-m%alpha0_deg < A_min) Then
          alpha0_deg_2 = 180 - m%alpha0_deg
        Else
          alpha0_deg_2 = -m%alpha0_deg
        End If

      ! Alpha was specified in the initial conditions, use the specified angle
      Else
        ! TODO - check that alpha is a valid angle.
        A_min = sv%alpha
        A_max = sv%alpha
      End If

      ! -------------------------------------------------------------------- !
      !    Loop through possible alphas, save the angle where the FC is max  !
      ! -------------------------------------------------------------------- !
      A = A_min
      CrackAngle: Do  ! Test various alphas

        alpha_temp = A*deg_to_rad  ! crack angle being evaluated (converted to radians)

        ! Crack normal in the reference configuration
        normal(2) = COS(alpha_temp)
        normal(3) = SIN(alpha_temp)

        ! Current crack normal direction
        R_cr(:,2) = Norm(MATMUL(F_inverse_transpose, normal))

        ! Current transverse direction
        R_cr(:,3) = CrossProduct(R_cr(:,1), R_cr(:,2))

        ! -------------------------------------------------------------------- !
        !    Determine the cohesive penalty stiffnesses                        !
        ! -------------------------------------------------------------------- !
        ! Does this DGD crack represent a crack or a delamination?
        If (A == 90) Then
          Q = 3
        Else
          Q = 2
        End If
        Pen(2) = p%penStiffMult*m%E2/sv%Lc(Q)
        Pen(1) = Pen(2)*m%GYT*m%SL*m%SL/(m%GSL*m%YT*m%YT) ! Corresponds to Turon et al (2010)
        Pen(3) = Pen(2)*m%GYT*m%ST*m%ST/(m%GSL*m%YT*m%YT)

        ! -------------------------------------------------------------------- !
        !    Determine the cohesive displacement-jump                          !
        ! -------------------------------------------------------------------- !
        T = MATMUL(Cauchy, R_cr(:,2)) ! Traction on fracture surface

        delta = MATMUL(TRANSPOSE(R_cr), T) / Pen

        ! -------------------------------------------------------------------- !
        !    Evaluate the cohesive law initiation criterion                    !
        ! -------------------------------------------------------------------- !
        Call cohesive_damage(m, delta, Pen, delta(2), B_temp, FIm_temp)

        ! -------------------------------------------------------------------- !
        !    Save the values corresponding to the maximum failure criteria     !
        ! -------------------------------------------------------------------- !
        If (FIm_temp > sv%FIm) Then
          sv%FIm        = FIm_temp
          sv%B          = B_temp
          alpha_test    = A
          sv%Fb1        = F(1,Q)
          sv%Fb2        = F(2,Q)
          sv%Fb3        = F(3,Q)
        End If

        If (A == A_max) EXIT CrackAngle

        ! Advance the crack angle
        NextAngle: Do
          ! Check to see if incrementing alpha would pass over +alpha0
          If (A < m%alpha0_deg .AND. A + p%alpha_inc > m%alpha0_deg) Then
            A = m%alpha0_deg
          ! If already at +alpha0, increment to next nearest multiple of alpha_inc
          Else If (A == m%alpha0_deg) Then
            A = A + p%alpha_inc - MOD(A + p%alpha_inc, p%alpha_inc)
          ! Check to see if incrementing alpha would pass over -alpha0
          Else If (A < alpha0_deg_2 .AND. A + p%alpha_inc > alpha0_deg_2) Then
            A = alpha0_deg_2
          ! If already at -alpha0, increment to next nearest multiple of alpha_inc
          Else If (A == alpha0_deg_2) Then
            A = A + p%alpha_inc - MOD(A + p%alpha_inc, p%alpha_inc)
          Else
            A = A + p%alpha_inc
          End If

          If (A /= 90) Exit NextAngle  ! Only evaluate 90 if is set via initial condition
        End Do NextAngle

      End Do CrackAngle

      ! -------------------------------------------------------------------- !
      !    If failure occurs, save alpha and indicate small dmg              !
      ! -------------------------------------------------------------------- !
      If (sv%FIm >= one) Then
        sv%d2    = 1.d-8 ! Used as a flag to call DGDEvolve
        sv%alpha = alpha_test
        Call log%info('DGDInit found FIm > one. Matrix damage initiated.')
        Call log%debug('alpha = ' // trim(str(sv%alpha)))
      End If
    End If

    ! -------------------------------------------------------------------- !
    !    Update elastic energy variable.                                   !
    ! -------------------------------------------------------------------- !
    enerIntern = zero
    Do I=1,3
      Do J=1,3
        enerIntern = enerIntern + 0.5d0*stress(I,J)*eps(I,J)
      End Do
    End Do

    Return
  End Subroutine DGDInit


  Subroutine DGDEvolve(U, F, F_old, m, p, sv, ndir, nshr, DT, Cauchy, enerIntern, enerInelas)
    ! Determines the matrix damage state variable based on the current   !
    ! deformation and mode mixity.                                       !

    Use forlog_Mod
    Use matrixAlgUtil_Mod
    Use matProp_Mod
    Use stateVar_Mod
    Use parameters_Mod
    Use stress_Mod
    Use cohesive_mod
    Use strain_mod
    Use plasticity_mod
    Use CDM_fiber_mod
    Use friction_mod
    Use schapery_mod

    ! -------------------------------------------------------------------- !
    ! Arguments
    Type(matProps), intent(IN) :: m
    Type(parameters), intent(IN) :: p
    Type(stateVars), intent(INOUT) :: sv
    Double Precision, Intent(IN) :: F(3,3), U(3,3), F_old(3,3)               ! Deformation gradient, stretch tensor
    Double Precision, Intent(IN) :: DT                                       ! Delta temperature
    Integer, intent(IN) :: ndir, nshr
    Double Precision, Intent(OUT) :: Cauchy(ndir,ndir)
    Double Precision, Intent(OUT) :: enerIntern, enerInelas

    ! -------------------------------------------------------------------- !
    ! Locals
    Integer :: mode                                                          ! Flag for whether the crack is a standard matrix crack (0) or a fiber compression matrix crack (1)

    Double Precision :: Stiff(ndir+nshr,ndir+nshr)                           ! Stiffness
    Double Precision :: stress(ndir,ndir)                                    ! Stress (energy conjugate to strain definition)
    Double Precision :: eps(ndir,ndir)                                       ! Strain
    Double Precision :: eps_old(ndir,ndir)                                   ! Strain

    ! Cohesive surface
    Double Precision :: alpha_rad                                            ! alpha in radians
    Double Precision :: R_cr(3,3)                                            ! Basis coordinate system for the cohesive surface
    Double Precision :: T(3), T_bulk(3)                                      ! Traction on crack interface from bulk material stresses
    Double Precision :: T_coh(3)                                             ! Traction on crack interface from cohesive law
    Double Precision :: normal(3)                                            ! Normal to the crack plane in the reference configuration
    Double Precision :: delta_coh(ndir,ndir)                                 ! matrix of cohesive displacements
    Double Precision :: damage_max                                           ! Maximum value for damage variable
    Double Precision :: Pen(3)                                               ! Penalty stiffnesses
    Double Precision :: damage_old, AdAe
    Double Precision :: delta_n_init
    Integer :: Q, alphaQ                                                     ! Flag to specify matrix crack or delamination (Q=2 matrix crack, Q=3 delam); angle at which transition occurs (depends on element geometry)
    Integer :: MD, EQk                                                       ! MatrixDamage and equilibrium loop indices

    ! Equilibrium loop
    Double Precision :: F_bulk(3,3)
    Double Precision :: F_bulk_old(3), F_bulk_change(3), F_bulk_inverse(3,3)
    Double Precision :: Residual(3)                                          ! Residual stress vector
    Double Precision :: tol_DGD
    Double Precision :: err, err_old

    ! For jacobian
    Double Precision :: Cauchy_d(ndir,ndir,3)                                ! Derivative of the Cauchy stress tensor
    Double Precision :: stress_d(ndir,ndir,3)                                ! Derivative of the stress
    Double Precision :: eps_d(ndir,ndir,3)                                   ! Derivative of the strain
    Double Precision :: delta_coh_d(3,3,3), R_cr_d(3,3,3), F_bulk_d(3,3,3)
    Double Precision :: T_coh_d(3,3), T_d(3,3)
    Double Precision :: r1length, r2length                                   ! Magnitude used for computing derivative
    Double Precision :: Jac(3,3)                                             ! Jacobian
    Double Precision :: T_coh_d_den_temp

    ! Convergence tools
    Double Precision :: crack_inflection_aid, aid, crack_inflection_cutback
    Integer :: restarts, cutbacks, restarts_max
    Integer :: crack_inversions, crack_inversions_max
    Logical :: crack_inverted
    Logical :: crack_open, crack_was_open, crack_open_test
    Logical :: Restart, Cutback

    ! Misc
    Double Precision :: X(ndir,ndir)                                         ! Reference configuration
    Double Precision :: dGdGc
    Double Precision :: tr, Y(ndir,ndir)
    Double Precision :: eye(ndir,ndir)
    Double Precision :: Fb_s1(3), Fb_s3(3)
    Double Precision :: L

    ! Added by Drew for fiber compression
    !Double Precision, Intent(IN) :: XC,w_kb,phi_ff
    Double Precision :: fiberDir(3)                                          ! Misaligned fiber direction in the reference configuration
    Double Precision :: rfT_temp, rfC_temp                                   ! Fiber damage thresholds, from previous converged solution
    Double Precision :: d1T_temp, d1C_temp                                   ! Fiber damage variables, from previous converged solution
    Double Precision :: phi0
    Double Precision :: d1

    ! Friction
    Double Precision :: slide_old(2)
    Integer :: forced_sticking
    Logical :: Sliding

    Double Precision, parameter :: zero=0.d0, one=1.d0, two=2.d0

    ! -------------------------------------------------------------------- !
    Call log%debug('Start of DGDEvolve')

    damage_max = one ! Maximum value for damage variables

    restarts_max = 2  ! equal to the number of starting points minus 1
    crack_inflection_cutback = 0.99d0
    crack_inversions_max = 4

    ! Initialize outputs
    sv%FIm = zero
    enerInelas = zero

    X = zero; X(1,1) = sv%Lc(1); X(2,2) = sv%Lc(2); X(3,3) = sv%Lc(3) ! Ref. Config.

    eye = zero; Do I = 1,3; eye(I,I) = one; End Do ! Identity Matrix

    tol_DGD = m%YT*p%tol_DGD_f ! Equilibrium loop tolerance [stress]

    ! crack or delamination?
    alphaQ = FLOOR(ATAN(sv%Lc(2)/sv%Lc(3))*45.d0/ATAN(one))
    alphaQ = alphaQ - MOD(alphaQ, p%alpha_inc)
    If (sv%alpha /= 90) Then
      Q = 2 ! matrix crack
    Else
      Q = 3 ! delamination
    End If
    alpha_rad = sv%alpha/45.d0*ATAN(one) ! alpha [radians]

    ! -------------------------------------------------------------------- !
    ! Penalty stiffness
    Pen(2) = p%penStiffMult*m%E2/sv%Lc(Q)
    Pen(1) = Pen(2)*m%GYT*m%SL*m%SL/(m%GSL*m%YT*m%YT) ! Corresponds to Turon et al (2010)
    Pen(3) = Pen(2)*m%GYT*m%ST*m%ST/(m%GSL*m%YT*m%YT)

    ! -------------------------------------------------------------------- !
    !    Define a crack-based coordinate system with a basis R_cr(3,3):    !
    ! -------------------------------------------------------------------- !
    normal(1) = zero
    normal(2) = COS(alpha_rad)
    normal(3) = SIN(alpha_rad)

    ! Current fiber direction
    R_cr(:,1) = Norm(F(:,1))

    ! -------------------------------------------------------------------- !
    !    Initial guesses for F_bulk:                                       !
    ! -------------------------------------------------------------------- !
    F_bulk(:,:) = F(:,:)
    F_bulk(1,Q) = sv%Fb1
    F_bulk(2,Q) = sv%Fb2
    F_bulk(3,Q) = sv%Fb3

    If (Length(F_bulk(:,Q)) == zero) F_bulk(Q,Q) = one

    ! Initialize the displ across the cohesive interface as zero
    delta_coh = zero

    ! -------------------------------------------------------------------- !
    !    Definition of Equilibrium alternate starting points:              !
    ! -------------------------------------------------------------------- !
    ! Starting point 1 is the Q-th column of the previous solution for F_bulk
    ! Starting point 2 is the Q-th column of F
    ! Starting point 3 is the cross product of the two non-decomposed columns of F
    Fb_s3(:) = Norm(CrossProduct(F(:,MOD(Q,3)+1), F(:,Q-1)))

    ! -------------------------------------------------------------------- !
    !    MatrixDamage Loop and solution controls definition                !
    ! -------------------------------------------------------------------- !
    MD = 0 ! Counter for MatrixDamage loop

    MatrixDamage: Do
      MD = MD + 1
      Call log%debug('MD', MD)
      Fb_s1(:) = F_bulk(:,Q) ! Starting point 1
      slide_old(:) = sv%slide(:)

      AdAe = sv%d2/(sv%d2 + (one - sv%d2)*two*Pen(1)*m%GSL/(m%SL*m%SL))

      ! -------------------------------------------------------------------- !
      !    Equilibrium Loop and solution controls definition                 !
      ! -------------------------------------------------------------------- !
      Restart = .False.
      Cutback = .False.
      restarts = 0  ! Indicates the number of loop restarts with new starting points
      cutbacks = 0  ! Cut-back counter
      crack_inversions = 0
      crack_inverted = .False.

      forced_sticking = 0  ! Indicates that sliding has not been suppressed
      aid = one  ! Artificially reduces the rate of change in F_bulk
      crack_inflection_aid = one
      err = Huge(zero)
      EQk = 0

      Equilibrium: Do ! Loop to determine the current F_bulk
        EQk = EQk + 1
        Call log%debug('EQk', EQk)

        ! -------------------------------------------------------------------- !
        If (Restart) Then

          ! Reset the Restart and Cutback flags
          Restart = .False.
          Cutback = .False.

          ! Attempt to use the "no sliding" condition, if not converged conventionally with friction.
          If (forced_sticking == 0 .AND. m%friction) Then
            Call log%info('Attempting no sliding, MD: ' // trim(str(MD)) // ' Restart: ' // trim(str(restarts)))
            forced_sticking = 1
          Else
            ! Advance the restart counter
            restarts = restarts + 1
            Call log%info('Restarting Equilibrium loop with new start point, Restart: ' // trim(str(restarts)))
            forced_sticking = 0
          End If

          ! If all starting points have been fully used...
          If (restarts > restarts_max) Then
            ! ...and if the matrix damage is already fully developed, delete the element.
            If (sv%d2 >= damage_max) Then
              If (sv%alpha == -999) Then
                Call writeDGDArgsToFile(m,p,sv,U,F,F_old,ndir,nshr,DT,log%arg,'DGDEvolve')
                Call log%terminate('Invalid alpha. Check value for alpha in the initial conditions.')
              End If
              If (p%terminate_on_no_convergence) Then
                Call writeDGDArgsToFile(m,p,sv,U,F,F_old,ndir,nshr,DT)
                Call log%terminate('DGDEvolve nonconvergence, terminating analysis.')
              Else
                Call log%warn('DGDEvolve nonconvergence, deleting failed element.')
                sv%STATUS = 0
              End If
              Exit MatrixDamage
            End If
            ! ...raise an error and halt the subroutine.
            Call writeDGDArgsToFile(m,p,sv,U,F,F_old,ndir,nshr,DT,log%arg,'DGDEvolve')
            Call log%terminate('No starting points produced a valid solution (DGDEvolve).')
            Exit MatrixDamage
          End If

          cutbacks = 0
          crack_inversions = 0
          crack_inverted = .False.

          ! Restart from a starting point
          F_bulk = F
          If (restarts == 0) Then
            F_bulk(:,Q) = Fb_s1(:)  ! Use starting point 1
          Else If (restarts == 1) Then
            Continue                ! Use starting point 2
          Else If (restarts == 2) Then
            F_bulk(:,Q) = Fb_s3(:)  ! Use starting point 3
          End If

          aid = one
          crack_inflection_aid = one

          ! Reset err
          err = Huge(zero)
        ! -------------------------------------------------------------------- !
        Else If (Cutback) Then

          ! Reset the Cutback flag
          Cutback = .False.

          ! Advance the cutback counter
          cutbacks = cutbacks + 1
          Call log%info('Cutting back, Cutbacks: ' // trim(str(cutbacks)) //', Restart: ' // trim(str(restarts)))

          aid = p%cutback_amount**cutbacks

          F_bulk(:,Q) = F_bulk_old(:) - F_bulk_change(:)*aid*crack_inflection_aid

          ! Reset err to the last converged value
          err = err_old
        End If
        ! -------------------------------------------------------------------- !
        ! Initialize all temporary state variables for use in Equilibrium loop:

        ! Shear nonlinearity temporary variables
        Call initializeTemp(sv, m)

        ! CDM fiber damage variables
        d1 = zero
        d1T_temp = sv%d1T
        d1C_temp = sv%d1C
        rfT_temp = sv%rfT
        rfC_temp = sv%rfC

        ! Store re-used matrices
        F_bulk_inverse = MInverse(F_bulk)

        R_cr(:,2) = MATMUL(TRANSPOSE(F_bulk_inverse), normal)
        r2length = Length(R_cr(:,2))    ! Un-normalized length of R_cr(:,2), used in the Jacobian
        R_cr(:,2) = R_cr(:,2)/r2length  ! Normalized R_cr(:,2)
        R_cr(:,3) = CrossProduct(R_cr(:,1), R_cr(:,2))

        delta_coh(:,Q) = MATMUL(TRANSPOSE(R_cr), F(:,Q) - F_bulk(:,Q))*X(Q,Q)

        If (delta_coh(2,Q) >= zero) Then
          crack_open = .True.
        Else
          crack_open = .False.
        End If

        ! Calculate the sign of the change in shear strain (for shear nonlinearity subroutine)
        If (m%shearNonlinearity12 .AND. Q == 2) Then
           sv%d_eps12_temp = Sign(one, (F(1,1)*F_bulk(1,2) + F(2,1)*F_bulk(2,2) + F_bulk(3,2)*F(3,1)) - (F_old(1,1)*sv%Fb1 + F_old(2,1)*sv%Fb2 + sv%Fb3*F_old(3,1)))
        End If
        If (m%shearNonlinearity13 .AND. Q == 3) Then
           sv%d_eps13_temp = Sign(one, (F(1,1)*F_bulk(1,3) + F(2,1)*F_bulk(2,3) + F_bulk(3,1)*F(3,3)) - (F_old(1,1)*sv%Fb1 + F_old(2,1)*sv%Fb2 + F_old(3,2)*sv%Fb3))
        End If

        ! Compute the Green-Lagrange strain tensor for the bulk material: eps
        Call Strains(F_bulk, m, DT, ndir, eps)

        ! Compute the plastic strains and remove from the strain tensor
        Call Plasticity(m, sv, p, ndir, nshr, eps, use_temp=.TRUE.)

        ! -------------------------------------------------------------------- !
        !    Evaluate the CDM fiber failure criteria and damage variable:      !
        ! -------------------------------------------------------------------- !
        If (eps(1,1) >= zero) Then
          If (m%fiberTenDam) Call FiberTenDmg(eps, ndir, m%E1, m%XT, m%GXT, m%fXT, m%fGXT, sv%Lc(1), m%cl, rfT_temp, d1T_temp, d1C_temp, sv%STATUS)
          d1 = d1T_temp
          d1C_temp = sv%d1C
        Else If (m%fiberCompDamBL) Then
          Call FiberCompDmg(eps, ndir, m%E1, m%XC, m%GXC, m%fXC, m%fGXC, sv%Lc(1), m%cl, rfT_temp, rfC_temp, d1T_temp, d1C_temp, sv%STATUS)
          d1 = d1C_temp
        End If

        ! -------------------------------------------------------------------- !
        !    Build the stiffness matrix:                                       !
        ! -------------------------------------------------------------------- !
        Call StiffFuncNL(m, ndir, nshr, d1, zero, zero, eps, Stiff, sv%Sr_temp)

        ! -------------------------------------------------------------------- !
        !    Determine the bulk material tractions on the fracture plane:      !
        ! -------------------------------------------------------------------- !
        stress = Hooke(Stiff, eps, nshr)
        Cauchy = convertToCauchy(stress, F_bulk)
        T = MATMUL(Cauchy, R_cr(:,2))  ! Traction on fracture surface
        T_bulk = MATMUL(TRANSPOSE(R_cr), T)

        ! -------------------------------------------------------------------- !
        !    Determine the cohesive tractions:                                 !
        ! -------------------------------------------------------------------- !
        If (crack_open) Then  ! Open cracks

          T_coh(:) = cohesive_traction(delta_coh(:,Q), Pen(:), sv%d2)
          sv%slide(1) = delta_coh(1,Q)
          sv%slide(2) = delta_coh(3,Q)

        Else  ! Closed cracks

          If (.NOT. m%friction) Then  ! Closed cracks without friction
            T_coh(:) = cohesive_traction(delta_coh(:,Q), Pen(:), sv%d2)
            sv%slide(1) = delta_coh(1,Q)
            sv%slide(2) = delta_coh(3,Q)

          Else  ! Closed cracks with friction
            Sliding = crack_is_sliding(delta_coh(:,Q), Pen, slide_old, m%mu, m%mu)
            If (forced_sticking == 1) Sliding = .False.
            Call crack_traction_and_slip(delta_coh(:,Q), Pen, slide_old, sv%slide, m%mu, m%mu, sv%d2, AdAe, T_coh, Sliding)

          End If

        End If

        ! -------------------------------------------------------------------- !
        !    Define the stress residual vector, R. R is equal to the           !
        !    difference in stress between the cohesive interface and the bulk  !
        !    stress projected onto the cohesive interface                      !
        ! -------------------------------------------------------------------- !
        Residual = T_coh - T_bulk

        ! Error in terms of the stress residual on the crack surface
        err_old = err
        err = Length(Residual)

        ! Output for visualizations
        Call log%debug('err', err)
        Call log%debug('Fb', F_bulk)

        ! -------------------------------------------------------------------- !
        !    Check for convergence.                                            !
        ! -------------------------------------------------------------------- !

        ! Check for a diverging solution
        If (err >= err_old) Then

          If (crack_inverted) Then
            ! A cut-back will not be performed if the crack opening state has just changed. It is necessary to
            ! allow this increase in error so the Jacobian can be calculated on "this side" of the crack state.
            Continue
          Else

            Call log%info('Solution is diverging, err: ' // trim(str(err)) // ' > ' // trim(str(err_old)))

            ! Cut-back using the current starting point
            If (cutbacks < p%cutbacks_max) Then
              Cutback = .True.

            ! Restart using a new starting point, if available
            Else
              Restart = .True.
            End If

            Cycle Equilibrium

          End If

        End If

        ! Ensures that an artificial no "sliding condition" is not forced as the solution
        If (err < tol_DGD .AND. forced_sticking == 1 .AND. .NOT. crack_open) Then
          If (Sliding .NEQV. crack_is_sliding(delta_coh(:,Q), Pen, slide_old, m%mu, m%mu)) Then
            forced_sticking = 2  ! Deactivates forced "no sliding" condition
            err = Huge(zero)  ! Resets the error. An increase in error here does not indicate divergence.
            Cycle Equilibrium
          End If
        End If

        ! Check for any inside-out deformation or an overly compressed bulk material
        If (err < tol_DGD .AND. MDet(F_bulk) < p%compLimit) Then

          Call log%warn('det(F_bulk) is below limit: ' // trim(str(MDet(F_bulk))) // ' Restart: ' // trim(str(restarts)))

          ! Restart using new starting point
          Restart = .True.
          Cycle Equilibrium
        End If

        ! If converged,
        If (err < tol_DGD) Then

          Call log%debug('Equilibrium loop found a converged solution.')

          sv%Fb1 = F_bulk(1,Q); sv%Fb2 = F_bulk(2,Q); sv%Fb3 = F_bulk(3,Q)

          ! Update fiber damage state variables
          sv%d1T = d1T_temp
          sv%d1C = d1C_temp
          sv%rfT = rfT_temp
          sv%rfC = rfC_temp

          ! Update shear nonlinearity state variables
          Call finalizeTemp(sv, m)

          ! If fully damaged
          If (sv%d2 >= damage_max) Then
            sv%d2 = damage_max
            sv%FIm = one
            EXIT MatrixDamage
          End If
          EXIT Equilibrium ! Check for change in sv%d2
        End If

        ! -------------------------------------------------------------------- !
        !    Find the derivative of the Cauchy stress tensor.                  !
        ! -------------------------------------------------------------------- !
        F_bulk_d    = zero
        R_cr_d      = zero
        delta_coh_d = zero
        eps_d       = zero
        stress_d    = zero
        Cauchy_d    = zero
        T_d         = zero
        T_coh_d     = zero

        Do I=1,3
          F_bulk_d(I,Q,I) = one

          ! info on derivative of x/||x||: http://blog.mmacklin.com/2012/05/
          R_cr_d(:,2,I) = -MATMUL(MATMUL(TRANSPOSE(F_bulk_inverse), MATMUL(TRANSPOSE(F_bulk_d(:,:,I)), TRANSPOSE(F_bulk_inverse))), normal)
          ! info on derivative of inverse matrix: http://planetmath.org/derivativeofinversematrix
          R_cr_d(:,2,I) = MATMUL(eye/r2length - OuterProduct(R_cr(:,2), R_cr(:,2))/r2length, R_cr_d(:,2,I))
          R_cr_d(:,3,I) = CrossProduct(R_cr(:,1), R_cr_d(:,2,I))

          delta_coh_d(:,Q,I) = MATMUL(TRANSPOSE(R_cr_d(:,:,I)), (F(:,Q) - F_bulk(:,Q)))*X(Q,Q) - MATMUL(TRANSPOSE(R_cr), F_bulk_d(:,Q,I))*X(Q,Q)

          If (crack_open) Then ! Open cracks
            T_coh_d(:,I) = Pen(:)*(one - sv%d2)*delta_coh_d(:,Q,I)

          Else If (.NOT. m%friction) Then ! Closed cracks without friction
            T_coh_d(1,I) = Pen(1)*(one - sv%d2)*delta_coh_d(1,Q,I)
            T_coh_d(2,I) = Pen(2)*delta_coh_d(2,Q,I)
            T_coh_d(3,I) = Pen(3)*(one - sv%d2)*delta_coh_d(3,Q,I)

          Else If (Sliding) Then ! Closed cracks with sliding friction
            T_coh_d_den_temp = (Pen(1)*(delta_coh(1,Q) - slide_old(1)))**2 + (Pen(3)*(delta_coh(3,Q) - slide_old(2)))**2

            T_coh_d(1,I) = Pen(1)*(one - sv%d2)*delta_coh_d(1,Q,I)
            T_coh_d(3,I) = Pen(3)*(one - sv%d2)*delta_coh_d(3,Q,I)

            If (T_coh_d_den_temp /= zero) Then
              T_coh_d(1,I) = T_coh_d(1,I) - AdAe*m%mu*Pen(2)*Pen(1)/SQRT(T_coh_d_den_temp)* &
                (delta_coh_d(2,Q,I)*(delta_coh(1,Q) - slide_old(1)) + delta_coh(2,Q)*delta_coh_d(1,Q,I) - delta_coh(2,Q)*(delta_coh(1,Q) - slide_old(1))* &
                (Pen(1)*Pen(1)*(delta_coh(1,Q) - slide_old(1))*delta_coh_d(1,Q,I) + &
                 Pen(3)*Pen(3)*(delta_coh(3,Q) - slide_old(2))*delta_coh_d(3,Q,I))/T_coh_d_den_temp)

              T_coh_d(3,I) = T_coh_d(3,I) - AdAe*m%mu*Pen(2)*Pen(3)/SQRT(T_coh_d_den_temp)* &
                (delta_coh_d(2,Q,I)*(delta_coh(3,Q) - slide_old(2)) + delta_coh(2,Q)*delta_coh_d(3,Q,I) - delta_coh(2,Q)*(delta_coh(3,Q) - slide_old(2))* &
                (Pen(1)*Pen(1)*(delta_coh(1,Q) - slide_old(1))*delta_coh_d(1,Q,I) + &
                 Pen(3)*Pen(3)*(delta_coh(3,Q) - slide_old(2))*delta_coh_d(3,Q,I))/T_coh_d_den_temp)
            End If
            T_coh_d(2,I) = Pen(2)*delta_coh_d(2,Q,I)

          Else ! Closed cracks with sticking friction
            T_coh_d(1,I) = Pen(1)*(one - sv%d2 + AdAe)*delta_coh_d(1,Q,I)
            T_coh_d(2,I) = Pen(2)*delta_coh_d(2,Q,I)
            T_coh_d(3,I) = Pen(3)*(one - sv%d2 + AdAe)*delta_coh_d(3,Q,I)
          End If

          eps_d(:,:,I) = (MATMUL(TRANSPOSE(F_bulk_d(:,:,I)), F_bulk) + MATMUL(TRANSPOSE(F_bulk), F_bulk_d(:,:,I)))/two
          stress_d(:,:,I) = Hooke(Stiff,eps_d(:,:,I),nshr)

          Y = MATMUL(F_bulk_inverse, F_bulk_d(:,:,I))
          tr = Y(1,1) + Y(2,2) + Y(3,3)
          Cauchy_d(:,:,I) = (MATMUL(MATMUL(F_bulk_d(:,:,I), stress) + MATMUL(F_bulk, stress_d(:,:,I)), TRANSPOSE(F_bulk)) + MATMUL(MATMUL(F_bulk, stress), TRANSPOSE(F_bulk_d(:,:,I))))/MDet(F_bulk) - Cauchy*tr

          T_d(:,I) = MATMUL(Cauchy_d(:,:,I), R_cr(:,2)) + MATMUL(Cauchy, R_cr_d(:,2,I))
        End Do

        ! -------------------------------------------------------------------- !
        !    Define the Jacobian matrix, J                                     !
        ! -------------------------------------------------------------------- !
        Jac = zero

        Jac(:,1) = T_coh_d(:,1) - MATMUL(TRANSPOSE(R_cr), T_d(:,1)) - MATMUL(TRANSPOSE(R_cr_d(:,:,1)), T)
        Jac(:,2) = T_coh_d(:,2) - MATMUL(TRANSPOSE(R_cr), T_d(:,2)) - MATMUL(TRANSPOSE(R_cr_d(:,:,2)), T)
        Jac(:,3) = T_coh_d(:,3) - MATMUL(TRANSPOSE(R_cr), T_d(:,3)) - MATMUL(TRANSPOSE(R_cr_d(:,:,3)), T)

        ! -------------------------------------------------------------------- !
        !    Calculate the new bulk deformation gradient                       !
        ! -------------------------------------------------------------------- !
        F_bulk_old(:)    = F_bulk(:,Q)
        F_bulk_change(:) = MATMUL(MInverse(Jac), Residual)
        F_bulk(:,Q)      = F_bulk_old(:) - F_bulk_change(:)*aid

        ! -------------------------------------------------------------------- !
        !    Check for a change in crack opening                               !
        ! -------------------------------------------------------------------- !

        ! Store the previous crack opening state variable
        crack_was_open = crack_open

        ! Update the crack opening state variable
        R_cr(:,2) = Norm(MATMUL(MInverse(TRANSPOSE(F_bulk)), normal))
        R_cr(:,3) = CrossProduct(R_cr(:,1), R_cr(:,2))

        delta_coh(:,Q) = MATMUL(TRANSPOSE(R_cr), F(:,Q) - F_bulk(:,Q))*X(Q,Q)

        If (delta_coh(2,Q) >= zero) Then
          crack_open = .True.
        Else
          crack_open = .False.
        End If

        crack_inflection_aid = one

        ! Check for a change in the crack opening state
        If (crack_open .EQV. crack_was_open) Then
          ! If there is no change, do nothing.
          crack_inverted = .False.

        Else If (crack_inversions < crack_inversions_max) Then
          If (crack_open) Then
            Call log%info('Change in crack opening. Crack now open.')
          Else
            Call log%info('Change in crack opening. Crack now closed.')
          End If
          crack_inversions = crack_inversions + 1
          crack_inverted = .True.

          ! Initialize a test variable for the crack opening state
          crack_open_test = crack_open

          ! Loop until the crack opening status changes back
          Do While (crack_open_test .EQV. crack_open)

            crack_inflection_aid = crack_inflection_aid * crack_inflection_cutback

            F_bulk(:,Q) = F_bulk_old(:) - F_bulk_change(:)*aid*crack_inflection_aid

            R_cr(:,2) = Norm(MATMUL(MInverse(TRANSPOSE(F_bulk)), normal))
            R_cr(:,3) = CrossProduct(R_cr(:,1), R_cr(:,2))

            delta_coh(:,Q) = MATMUL(TRANSPOSE(R_cr), F(:,Q) - F_bulk(:,Q))*X(Q,Q)

            If (delta_coh(2,Q) >= zero) Then
              crack_open_test = .True.
            Else
              crack_open_test = .False.
            End If

          End Do

          ! Return F_bulk to the state immediately before the crack state changed back
          crack_inflection_aid = crack_inflection_aid / crack_inflection_cutback
          F_bulk(:,Q) = F_bulk_old(:) - F_bulk_change(:)*aid*crack_inflection_aid

        Else
          crack_inverted = .False.
        End If
      End Do Equilibrium


      Call log%debug('Exited Equilibrium loop.')

      ! Store the old cohesive damage variable for checking for convergence
      damage_old = sv%d2

      If (MD == 1) delta_n_init = MIN(zero, delta_coh(2,Q))

      ! Update the cohesive damage variable
      Call cohesive_damage(m, delta_coh(:,Q), Pen, delta_n_init, sv%B, sv%FIm, sv%d2, dGdGc)

      ! Check for damage advancement
      If (sv%d2 <= damage_old) Then  ! If there is no damage progression,
        Call log%debug('No change in matrix damage variable, d2 ' // trim(str(sv%d2)))
        EXIT MatrixDamage
      Else
        Call log%debug('Change in matrix damage variable, d2 ' // trim(str(sv%d2)))
        enerInelas = enerInelas + dGdGc*(m%GYT + (m%GSL - m%GYT)*sv%B**m%eta_BK)
      End If

      ! Check for convergence based on rate of energy dissipation
      If (dGdGc < p%dGdGc_min) Then
        Call log%info('Solution accepted due to small change in dmg.')
        Call log%info('MD: ' // trim(str(MD)) // ' AdAe: ' // trim(str(AdAe)))
        EXIT MatrixDamage
      End If

      ! Limit number of MatrixDamage loop iterations
      If (MD > p%MD_max) Then
        Call log%info('MatrixDamage loop limit exceeded. MD: ' // trim(str(MD)))
        Call log%info('dGdGc: ' // trim(str(dGdGc)) // ' AdAe: ' // trim(str(AdAe)))
        EXIT MatrixDamage
      End If

    End Do MatrixDamage


    Call log%debug('Exited matrix damage loop, MD: ' // trim(str(MD)))

    ! -------------------------------------------------------------------- !
    !    Update elastic energy variable.                                   !
    ! -------------------------------------------------------------------- !
    enerIntern = zero
    Do I=1,3
      Do J=1,3
        enerIntern = enerIntern + 0.5d0*stress(I,J)*eps(I,J)
      End Do
    End Do

    ! -------------------------------------------------------------------- !
    Return
  End Subroutine DGDEvolve


  Subroutine DGDKinkband(U, F, F_old, m, p, sv, ndir, nshr, DT, Cauchy, enerIntern)

    Use forlog_Mod
    Use matrixAlgUtil_Mod
    Use matProp_Mod
    Use stateVar_Mod
    Use parameters_Mod
    Use stress_Mod
    Use strain_mod
    Use plasticity_mod

    ! -------------------------------------------------------------------- !
    ! Arguments
    Type(matProps), intent(IN) :: m
    Type(parameters), intent(IN) :: p
    Type(stateVars), intent(INOUT) :: sv
    Double Precision, Intent(IN) :: F(3,3), U(3,3), F_old(3,3)               ! Deformation gradient, stretch tensor
    Double Precision, Intent(IN) :: DT                                       ! Delta temperature
    Integer, intent(IN) :: ndir, nshr
    Double Precision, Intent(OUT) :: Cauchy(ndir,ndir)
    Double Precision, Intent(OUT) :: enerIntern

    ! -------------------------------------------------------------------- !
    ! Locals
    Double Precision :: Stiff(ndir+nshr,ndir+nshr)                           ! Stiffness
    Double Precision :: R_cr(3,3)                                            ! Basis coordinate system for the cohesive surface
    Double Precision :: normalDir(3)                                         ! Normal to the crack plane in the reference configuration
    Double Precision :: fiberDir(3)                                          ! Current fiber direction
    Double Precision :: gamma_rphi0                                          ! Shear strain in the misaligned frame
    Double Precision :: X(ndir,ndir)                                         ! Reference configuration
    Double Precision :: eye(ndir,ndir)
    Double Precision :: R_phi0(3,3)                                          ! Transformation to the misaligned frame from the reference frame
    Double Precision :: d1

    ! Kink band region
    Double Precision :: Fkb(3,3)
    Double Precision :: Fkb_inverse(3,3)
    Double Precision :: Fkb_old(3,3)
    Double Precision :: epskb(ndir,ndir)
    Double Precision :: pk2_fiberDirkb(3,3)
    Double Precision :: stresskb(ndir,ndir)
    Double Precision :: Cauchykb(ndir,ndir)
    Double Precision :: Tkb(3)
    Double Precision :: eps12kb_dir

    ! Material region
    Double Precision :: Fm(3,3)
    Double Precision :: epsm(ndir,ndir)
    Double Precision :: pk2_fiberDirm(3,3)
    Double Precision :: stressm(ndir,ndir)
    Double Precision :: Cauchym(ndir,ndir)
    Double Precision :: Tm(3)

    ! Constants
    Double Precision, parameter :: zero=0.d0, one=1.d0, two=2.d0
    Double Precision, parameter :: pert=0.0001d0                             ! Small perturbation used to compute a numerical derivative

    ! Internal convergence
    Integer :: EQk                                                           ! Equilibrium loop index
    Double Precision :: Residual(3)                                          ! Residual stress vector
    Double Precision :: x_change(3)
    Double Precision :: tol_DGD
    Double Precision :: err, err_old
    Double Precision :: aid
    Logical :: Cutback, Restart
    Double Precision :: cutbacks, restarts

    ! Jacobian
    Double Precision :: Jac(3,3)
    Double Precision :: dFkb_dFkb1(3,3,3)
    Double Precision :: dR_cr_dFkb1(3,3,3)
    Double Precision :: dFm_dFkb1(3,3,3)
    Double Precision :: dEm_dFkb1(3,3,3)
    Double Precision :: dSm_dFkb1(3,3,3)
    Double Precision :: dCauchym_dFkb1(3,3,3)
    Double Precision :: dTm_dFkb1(3,3)
    Double Precision :: dEkb_dFkb1(3,3,3)
    Double Precision :: dSkb_dFkb1(3,3,3)
    Double Precision :: dCauchykb_dFkb1(3,3,3)
    Double Precision :: dTkb_dFkb1(3,3)
    Double Precision :: Ym(3,3), Ykb(3,3)
    Double Precision :: dplas12, dinel12

    ! -------------------------------------------------------------------- !

    ! Initialize
    aid = one
    d1 = zero
    sv%gamma = zero
    X = zero; X(1,1) = sv%Lc(1); X(2,2) = sv%Lc(2); X(3,3) = sv%Lc(3) ! Ref. Config.
    eye = zero; DO I = 1,3; eye(I,I) = one; end DO ! Identity Matrix
    tol_DGD = m%YT*p%tol_DGD_f ! Equilibrium loop tolerance [stress]

    ! -------------------------------------------------------------------- !
    !    Starting point                                                    !
    ! -------------------------------------------------------------------- !
    Fm(:,:) = F(:,:)
    Fkb(:,:) = F(:,:); Fkb(1,1) = sv%Fb1; Fkb(2,1) = sv%Fb2; Fkb(3,1) = sv%Fb3

    ! Make sure that we have a valid starting point
    If (Length(Fkb(:,1)) == zero) Fkb(1,1) = one

    ! Normal to misaligned fibers (reference config)
    normalDir = (/-sin(sv%phi0), cos(sv%phi0), zero/)

    ! Misaligned fiber direction (reference config)
    fiberDir = (/cos(sv%phi0), sin(sv%phi0), zero/)

    ! Build R_phi0 matrix
    R_phi0(:,1) = fiberDir
    R_phi0(:,2) = normalDir
    R_phi0(:,3) = (/zero, zero, one/)

    ! Initialize
    Fkb_old = F_old
    Fkb_old(:,1) = (/sv%Fb1, sv%Fb2, sv%Fb3/)

    ! -------------------------------------------------------------------- !
    !    Equilibrium Loop and solution controls definition                 !
    ! -------------------------------------------------------------------- !
    err = Huge(zero)
    EQk = 0
    Cutback = .False.
    cutbacks = 0  ! Cut-back counter
    Restart = .False.
    restarts = 0  ! Restart counter

    Equilibrium: Do ! Loop to determine the current Fkb(:,1)
      EQk = EQk + 1

      ! Reset counters if restart
      If (Restart) Then
        EQk = 0
        cutbacks = 0
        Restart = .False.
      End If

      ! Reduce aid if cutting back
      If (Cutback) Then

        ! Reset the Cutback flag
          Cutback = .False.

          ! Advance the cutback counter
          cutbacks = cutbacks + 1
          Call log%info('Cutting back, Cutbacks: ' // trim(str(cutbacks)))

          aid = p%cutback_amount**cutbacks
      End If

      Call log%debug('Equilibrium Start EQk: ' // trim(str(EQk)) // ' cutbacks: ' // trim(str(cutbacks)) // ' restarts: ' // trim(str(restarts)))


      ! -------------------------------------------------------------------- !
      ! Initialize all temporary state variables for use in Equilibrium loop:
      Call initializeTemp(sv, m)



      ! -------------------------------------------------------------------- !
      ! Compatibility constraints
      Fm(:,1) = (F(:,1)*sv%Lc(1) - Fkb(:,1)*m%w_kb)/(sv%Lc(1)-m%w_kb)

      ! Store re-used matrices
      Fkb_inverse = MInverse(Fkb)

      ! -------------------------------------------------------------------- !
      ! Define the current crack coordinate system
      R_cr(:,1) = MATMUL(Fkb, fiberDir)
      r1length = Length(R_cr(:,1))    ! Un-normalized length of R_cr(:,1), used in the Jacobian
      R_cr(:,1) = R_cr(:,1)/r1length  ! Normalized R_cr(:,2)
      R_cr(:,2) = MATMUL(TRANSPOSE(Fkb_inverse), normalDir)
      r2length = Length(R_cr(:,2))    ! Un-normalized length of R_cr(:,2), used in the Jacobian
      R_cr(:,2) = R_cr(:,2)/r2length  ! Normalized R_cr(:,2)
      R_cr(:,3) = CrossProduct(R_cr(:,1), R_cr(:,2))

      ! Calculate the angle that define the rotation of the fibers
      sv%gamma = ATAN(R_cr(2,1)/R_cr(1,1)) - sv%phi0

      ! Calculate the sign of the change in shear strain (for shear nonlinearity subroutine)
      sv%d_eps12_temp = Sign(one, (Fkb(1,1)*Fkb(1,2) + Fkb(2,1)*Fkb(2,2) + Fkb(3,1)*Fkb(3,2)) - (sv%Fb1*Fkb_old(1,2) + sv%Fb2*Fkb_old(2,2) + sv%Fb3*Fkb_old(3,2)))

      ! -------------------------------------------------------------------- !
      !    Calculate the stress in the bulk material region:                 !
      ! -------------------------------------------------------------------- !
      ! epsm = GLStrain(Fm,ndir)
      Call Strains(Fm, m, DT, ndir, epsm)
      epsm = MATMUL(TRANSPOSE(R_phi0), MATMUL(epsm, R_phi0))
      Call StiffFuncNL(m, ndir, nshr, d1, zero, zero, epsm, Stiff, sv%Sr)
      pk2_fiberDirm = Hooke(Stiff, epsm, nshr) ! 2PK
      stressm = MATMUL(R_phi0, MATMUL(pk2_fiberDirm, TRANSPOSE(R_phi0)))  ! 2PK rotated back to the reference direction
      Cauchym = convertToCauchy(stressm, Fm)
      Tm = MATMUL(Cauchym, R_cr(:,1))  ! Traction on surface with normal in fiber direction

      ! -------------------------------------------------------------------- !
      !    Calculate the stress in the kink band bulk region:                 !
      ! -------------------------------------------------------------------- !
      Call Strains(Fkb, m, DT, ndir, epskb)
      epskb = MATMUL(TRANSPOSE(R_phi0), MATMUL(epskb, R_phi0))
      Call Plasticity(m, sv, p, ndir, nshr, epskb, use_temp=.TRUE.)
      gamma_rphi0 = two*(epskb(1,2) + sv%Plas12_temp/two)
      Call StiffFuncNL(m, ndir, nshr, d1, zero, zero, epskb, Stiff, sv%Sr)
      pk2_fiberDirkb = Hooke(Stiff, epskb, nshr) ! 2PK in the fiber direction
      stresskb = MATMUL(R_phi0, MATMUL(pk2_fiberDirkb, TRANSPOSE(R_phi0)))  ! 2PK rotated back to the reference direction
      Cauchykb = convertToCauchy(stresskb, Fkb)
      Cauchy = Cauchykb
      Tkb = MATMUL(Cauchykb, R_cr(:,1))  ! Traction on surface with normal in fiber direction


      ! Failure index for fiber kinking
      sv%rfC = abs((-1*Cauchy(1,1)*cos(two*(sv%phi0+gamma_rphi0)))/((ramberg_osgood(gamma_rphi0 + pert, m%G12, m%aPL, m%nPL) - ramberg_osgood(gamma_rphi0 - pert, m%G12, m%aPL, m%nPL)) / (two * pert)))

      ! -------------------------------------------------------------------- !
      !    Define the stress residual vector, R. R is equal to the           !
      !    difference in stress between the cohesive interface and the bulk  !
      !    stress projected onto the cohesive interface                      !
      ! -------------------------------------------------------------------- !
      Residual(1:3) = MATMUL(TRANSPOSE(R_cr), Tm) - MATMUL(TRANSPOSE(R_cr), Tkb)

      ! -------------------------------------------------------------------- !
      !    Check for convergence.                                            !
      ! -------------------------------------------------------------------- !
      err_old = err
      err = Length(Residual)
      percentChangeInErr = (err - err_old)/err_old
      Call log%debug('err: ' // trim(str(err)) // ', errchange: ' // trim(str(percentChangeInErr)) // ', aid: ' // trim(str(aid)))

      ! Do not bother to attempt to find equilibrium if the element was deleted in Strains
      If (sv%STATUS == 0) Then
        EXIT Equilibrium
      End If

      ! Check for any inside-out deformation or an overly compressed bulk material
      If (err < tol_DGD .AND. MDet(Fkb) < p%compLimit) Then
        Call writeDGDArgsToFile(m,p,sv,U,F,F_old,ndir,nshr,DT,log%arg,'DGDKinkband')
        If (p%terminate_on_no_convergence) Then
          Call log%terminate('Highly distorted element (DGDKinkband).')
        Else
          Call log%warn('Deleting highly distorted element (DGDKinkband).')
          sv%STATUS = 0
        End If
        EXIT Equilibrium
      End If

      ! If converged,
      If (err < tol_DGD) Then

        ! Save starting point
        sv%Fb1 = Fkb(1,1); sv%Fb2 = Fkb(2,1); sv%Fb3 = Fkb(3,1)

        ! Update shear nonlinearity state variables
        Call finalizeTemp(sv, m)

        Call log%debug('Successfully converged. Iterations: ' // trim(str(EQk)))

        EXIT Equilibrium
      End If

      ! Check for maximum number of iterations
      If (EQk > p%EQ_max) Then
        If (restarts < 1) Then
          restarts = restarts + 1
          Restart = .True.  ! One restart is allowed; modified jacobian is used
          Call log%info('Restarting with modified jacobian.')
          Cycle Equilibrium
        Else
          Call writeDGDArgsToFile(m,p,sv,U,F,F_old,ndir,nshr,DT,log%arg,'DGDKinkband')
          If (p%terminate_on_no_convergence) Then
            Call log%terminate('Equilibrium loop reached maximum number of iterations (DGDKinkband).')
          Else
            Call log%warn('Deleting element for which equilibrium loop reached maximum number of iterations (DGDKinkband).')
            sv%STATUS = 0
          End If
          Exit Equilibrium
        End If
      End If

      ! Check for a diverging solution
      If (percentChangeInErr > 0.1d0) Then
        Call log%info('Solution is diverging, err: ' // trim(str(err)) // ' > ' // trim(str(err_old)))

        ! Cut-back using the current starting point
        If (cutbacks < p%cutbacks_max) Then
          Cutback = .True.

        ! Quit
        Else If (restarts < 1) Then
          restarts = restarts + 1
          Restart = .True.  ! One restart is allowed; modified jacobian is used
          Call log%info('Restarting with modified jacobian.')
        Else
          Call writeDGDArgsToFile(m,p,sv,U,F,F_old,ndir,nshr,DT,log%arg,'DGDKinkband')
          If (p%terminate_on_no_convergence) Then
            Call log%terminate('Reached max cutback limit (DGDKinkband).')
          Else
            Call log%warn('Reached max cutback limit, deleting element (DGDKinkband).')
            sv%STATUS = 0
          End If
          EXIT Equilibrium
        End If

        Cycle Equilibrium
      End If


      ! -------------------------------------------------------------------- !
      !    Find the derivative of the Cauchy stress tensor.                  !
      ! -------------------------------------------------------------------- !
      ! Initialize
      dFkb_dFkb1 = zero
      dR_cr_dFkb1 = zero
      dFm_dFkb1 = zero
      dEm_dFkb1 = zero
      dSm_dFkb1 = zero
      dCauchym_dFkb1 = zero
      dTm_dFkb1 = zero
      dEkb_dFkb1 = zero
      dplas12 = zero
      dinel12 = zero
      dSkb_dFkb1 = zero
      dCauchykb_dFkb1 = zero
      dTkb_dFkb1 = zero
      Jac = zero

      DO I=1,3

        dFkb_dFkb1(I,1,I) = one

        ! Coordinate system
        ! dR_cr_dFkb1
        dR_cr_dFkb1(:,1,I) = MATMUL(eye/r1length - OuterProduct(R_cr(:,1), R_cr(:,1))/r1length, MATMUL(dFkb_dFkb1(:,:,I), fiberDir))
        dR_cr_dFkb1(:,2,I) = -MATMUL(MATMUL(TRANSPOSE(Fkb_inverse), MATMUL(TRANSPOSE(dFkb_dFkb1(:,:,I)), TRANSPOSE(Fkb_inverse))), normalDir)
        dR_cr_dFkb1(:,2,I) = MATMUL(eye/r2length - OuterProduct(R_cr(:,2), R_cr(:,2))/r2length, dR_cr_dFkb1(:,2,I))
        ! cross product
        dR_cr_dFkb1(:,3,I) =  CrossProduct(dR_cr_dFkb1(:,1,I), R_cr(:,2)) +  CrossProduct(R_cr(:,1), dR_cr_dFkb1(:,2,I))

        ! --- dsigm/dFkb1
        ! Bulk material
        dFm_dFkb1(I,1,I) = -m%w_kb/(X(1,1)-m%w_kb);
        dEm_dFkb1(:,:,I) = (MATMUL(TRANSPOSE(dFm_dFkb1(:,:,I)), Fm) + MATMUL(TRANSPOSE(Fm), dFm_dFkb1(:,:,I)))/two
        dEm_dFkb1(:,:,I) = MATMUL(TRANSPOSE(R_phi0), MATMUL(dEm_dFkb1(:,:,I), R_phi0))
        dSm_dFkb1(:,:,I) = Hooke(Stiff, dEm_dFkb1(:,:,I), nshr)
        dSm_dFkb1(:,:,I) = MATMUL(R_phi0, MATMUL(dSm_dFkb1(:,:,I), TRANSPOSE(R_phi0)))
        Ym = MATMUL(MInverse(Fm), dFm_dFkb1(:,:,I))
        tr = Ym(1,1) + Ym(2,2) + Ym(3,3)
        dCauchym_dFkb1(:,:,I) = (MATMUL(MATMUL(dFm_dFkb1(:,:,I), stressm) + MATMUL(Fm, dSm_dFkb1(:,:,I)), TRANSPOSE(Fm)) + MATMUL(MATMUL(Fm, stressm), TRANSPOSE(dFm_dFkb1(:,:,I))))/MDet(Fm) - Cauchym*tr
        dTm_dFkb1(:,I) = MATMUL(dCauchym_dFkb1(:,:,I), R_cr(:,1)) + MATMUL(Cauchym, dR_cr_dFkb1(:,1,I))

        ! Kink band bulk
        dEkb_dFkb1(:,:,I) = (MATMUL(TRANSPOSE(dFkb_dFkb1(:,:,I)), Fkb) + MATMUL(TRANSPOSE(Fkb), dFkb_dFkb1(:,:,I)))/two
        dEkb_dFkb1(:,:,I) = MATMUL(TRANSPOSE(R_phi0), MATMUL(dEkb_dFkb1(:,:,I), R_phi0))
        If (restarts .EQ. 1) Then
          dplas12 = sv%Plas12_temp
          dinel12 = sv%Inel12_temp
          Call ro_plasticity(two*dEkb_dFkb1(1,2,I), sv%d_eps12_temp, m%G12, m%aPL, m%nPL, sv%Inel12c, dplas12, dinel12)
          dEkb_dFkb1(1,2,I) = dEkb_dFkb1(1,2,I) - dplas12/two
          dEkb_dFkb1(2,1,I) = dEkb_dFkb1(1,2,I)
        End If
        dSkb_dFkb1(:,:,I) = Hooke(Stiff, dEkb_dFkb1(:,:,I), nshr)
        dSkb_dFkb1(:,:,I) = MATMUL(R_phi0, MATMUL(dSkb_dFkb1(:,:,I), TRANSPOSE(R_phi0)))
        Ykb = MATMUL(MInverse(Fkb), dFkb_dFkb1(:,:,I))
        tr = Ykb(1,1) + Ykb(2,2) + Ykb(3,3)
        dCauchykb_dFkb1(:,:,I) = (MATMUL(MATMUL(dFkb_dFkb1(:,:,I), stresskb) + MATMUL(Fkb, dSkb_dFkb1(:,:,I)), TRANSPOSE(Fkb)) + MATMUL(MATMUL(Fkb, stresskb), TRANSPOSE(dFkb_dFkb1(:,:,I))))/MDet(Fkb) - Cauchykb*tr
        dTkb_dFkb1(:,I) = MATMUL(dCauchykb_dFkb1(:,:,I), R_cr(:,1)) + MATMUL(Cauchykb, dR_cr_dFkb1(:,1,I))

        ! Contribution to jacobian
        Jac(:,i) = MATMUL(TRANSPOSE(dR_cr_dFkb1(:,:,i)), Tm) + MATMUL(TRANSPOSE(R_cr), dTm_dFkb1(:,i)) - (MATMUL(TRANSPOSE(dR_cr_dFkb1(:,:,i)), Tkb) + MATMUL(TRANSPOSE(R_cr), dTkb_dFkb1(:,i)));

      End Do

      ! -------------------------------------------------------------------- !
      ! Calculate the new Fkb
      x_change(:) = MATMUL(MInverse(Jac), Residual)
      Fkb(:,1) = Fkb(:,1) - x_change*aid

    End Do Equilibrium

    ! -------------------------------------------------------------------- !
    !    Update elastic energy variable.                                   !
    ! -------------------------------------------------------------------- !
    enerIntern = zero
    Do I=1,3
      Do J=1,3
        enerIntern = enerIntern + 0.5d0*stressm(I,J)*epsm(I,J)
      End Do
    End Do

    ! -------------------------------------------------------------------- !
    Return
  End Subroutine DGDKinkband


  Function alpha0_DGD(m)
    ! Determines the orientation of the angle alpha0 when subject to sigma22 = -Yc

    Use forlog_Mod
    Use matrixAlgUtil_Mod
    Use stress_Mod
    Use matProp_Mod

    ! Arguments
    Type(matProps), intent(IN) :: m
    ! Double Precision, Intent(IN) :: alpha0, E1, E2, E3, G12, G13, G23, v12, v13, v23, Yc

    ! Locals
    Double Precision :: E3, G13, v13, G23    ! For transverse isotropy assumption
    Double Precision :: C(6,6)               ! 3-D Stiffness
    Double Precision :: F(3)                 ! Represents diagonal of deformation gradient tensor
    Double Precision :: Residual(3)          ! Residual vector
    Double Precision :: tolerance
    Double Precision :: err
    Double Precision :: Jac(3,3)             ! Jacobian
    Integer :: counter
    Double Precision, parameter :: zero=0.d0, one=1.d0, two=2.d0
    ! -------------------------------------------------------------------- !

    Call log%debug('Start of Function alpha0_DGD')

    tolerance = 1.d-4  ! alphaLoop tolerance

    ! Assume transverse isotropy if needed
    If (m%E3 < one) Then
      E3 = m%E2
    Else
      E3 = m%E3
    End If
    If (m%G13 < one) Then
      G13 = m%G12
    Else
      G13 = m%G13
    End If
    If (m%v13 == zero) Then
      v13 = m%v12
    Else
      v13 = m%v13
    End If
    If (m%G23 < one) Then
      G23 = m%E2/two/(one + m%v23)
    Else
      G23 = m%G23
    End If

    ! Build the stiffness matrix
    C = StiffFunc(6, m%E1, m%E2, E3, m%G12, G13, G23, m%v12, v13, m%v23, zero, zero, zero)

    ! Make an initial guess
    F = (/ one, one, one /)

    ! -------------------------------------------------------------------- !
    !    alphaLoop Loop and solution controls definition                   !
    ! -------------------------------------------------------------------- !
    counter = 0  ! Counter for alphaLoop
    counter_max = 100

    alphaLoop: Do  ! Loop to determine the F which corresponds to sigma22 = -Yc
      counter = counter + 1
      ! -------------------------------------------------------------------- !
      !    Define the stress residual vector, R. R is equal to the           !
      !    difference in stress between the cohesive interface and the bulk  !
      !    stress projected onto the cohesive interface                      !
      ! -------------------------------------------------------------------- !
      Residual(1) = (C(1,1)*(F(1)*F(1) - one) + C(1,2)*(F(2)*F(2) - one) + C(1,3)*(F(3)*F(3) - one))/two
      Residual(2) = (C(2,1)*(F(1)*F(1) - one) + C(2,2)*(F(2)*F(2) - one) + C(2,3)*(F(3)*F(3) - one))/two + m%Yc*F(1)*F(3)/F(2)
      Residual(3) = (C(3,1)*(F(1)*F(1) - one) + C(3,2)*(F(2)*F(2) - one) + C(3,3)*(F(3)*F(3) - one))/two

      ! Check for convergence
      err = Length(Residual)

      ! If converged,
      If (err < tolerance) Then
        alpha0_DGD = ATAN(F(2)/F(3)*TAN(m%alpha0))
        EXIT alphaLoop
      End If
      IF (counter == counter_max) Call log%error('Function alpha0_DGD failed to converge')

      ! Define the Jacobian matrix, J
      Jac = zero

      Jac(1,1) = C(1,1)*F(1)
      Jac(1,2) = C(1,2)*F(2)
      Jac(1,3) = C(1,3)*F(3)

      Jac(2,1) = C(2,1)*F(1) + two*m%Yc*F(3)/F(2)
      Jac(2,2) = C(2,2)*F(1) - two*m%Yc*F(3)*F(1)/(F(2)*F(2))
      Jac(2,3) = C(2,3)*F(1) + two*m%Yc*F(1)/F(2)

      Jac(3,1) = C(3,1)*F(1)
      Jac(3,2) = C(3,2)*F(2)
      Jac(3,3) = C(3,3)*F(3)

      ! Calculate the new diagonal deformation gradient
      F = F - MATMUL(MInverse(Jac), Residual)

    End Do alphaLoop

    Return
  End Function alpha0_DGD


  Subroutine writeDGDArgsToFile(m, p, sv, U, F, F_old, ndir, nshr, DT, args, called_from)
    ! Print DGDEvolve args at error

    Use matProp_Mod
    Use stateVar_Mod
    Use parameters_Mod
    Use vumatArg_Mod

    ! Arguments
    Type(matProps), intent(IN) :: m
    Type(parameters), intent(IN) :: p
    Type(stateVars), intent(IN) :: sv
    Double Precision, Intent(IN) :: F(3,3), U(3,3), F_old(3,3)               ! Deformation gradient, stretch tensor
    Double Precision, Intent(IN) :: DT                                       ! Delta temperature
    Type(vumatArg) :: args
    Character(*), intent(IN) :: called_from
    Integer, intent(IN) :: ndir, nshr

    ! Locals
    Integer :: lenOutputDir, lenJobName, debugpy_count_local
    Character(len=256) :: outputDir, fileName, jobName
    Character(len=32) :: debugpy_count_str, nElement_str
    Integer, parameter :: file_unit = 101
    ! -------------------------------------------------------------------- !

#ifndef PYEXT
    Call VGETJOBNAME(jobName, lenJobName)
    Call VGETOUTDIR(outputDir, lenOutputDir)

    ! Sometimes the debug file counter is going to -1; unclear why
    debugpy_count_local = sv%debugpy_count
    If (debugpy_count_local < 0) debugpy_count_local = 0

    ! Build up the filename
    write (nElement_str, *) args%nElement
    nElement_str = adjustl(nElement_str)
    write (debugpy_count_str, *) debugpy_count_local
    debugpy_count_str = adjustl(debugpy_count_str)
    fileName = trim(outputDir) // '/' // trim(jobName) // '-' // trim(nElement_str) // '-debug-'  ! Name of output file
    open(unit = file_unit, file = trim(fileName)//trim(debugpy_count_str)// '.py')
#else
    fileName = 'pyextmod_debug.py'
    open(unit = file_unit, file = fileName, status='replace', recl=1000)
#endif
    Call writeMaterialPropertiesToFile(file_unit, m)
    Call writeParametersToFile(file_unit, p)
    Call writeStateVariablesToFile(file_unit, sv, m)
    Call write3x3ToFile(file_unit, U, 'U')
    Call write3x3ToFile(file_unit, F, 'F')
    Call write3x3ToFile(file_unit, F_old, 'F_old')
    write(file_unit, "(A,E22.15E2)") 'thickness = ', sv%Lc(3)
    write(file_unit, "(A,E22.15E2)") 'temp_change = ', DT
    write(file_unit, "(A,I1,A)") 'ndir = ', ndir
    write(file_unit, "(A,I1,A)") 'nshr = ', nshr
    write(file_unit, "(A,I10,A)") 'element = ', args%nElement
    write(file_unit, "(A,E22.15E2)") 'total_time = ', args%totalTime
    write(file_unit, "(A,A,A)") 'called_from = "', called_from, '"'

    close(file_unit)

    Return
  End Subroutine writeDGDArgsToFile


  Subroutine write3x3ToFile(fileUnit, array, label)
    ! Write a 3x3 tensor to a file as a python dictionary
    ! Assumes that file opening and closing is handled elsewhere

    ! Arguments
    Integer, intent(IN) :: fileUnit
    Double Precision, intent(IN) :: array(3,3)
    Character(*) :: label

    ! Locals
    Character(len=32), parameter :: nameValueFmt = "(A,E22.15E2,A)"
    ! -------------------------------------------------------------------- !

    write(101, "(A)") trim(label) // ' = ['
    write(101, nameValueFmt) '    ', array(1,1), ','
    write(101, nameValueFmt) '    ', array(2,2), ','
    write(101, nameValueFmt) '    ', array(3,3), ','
    write(101, nameValueFmt) '    ', array(1,2), ','
    write(101, nameValueFmt) '    ', array(2,3), ','
    write(101, nameValueFmt) '    ', array(3,1), ','
    write(101, nameValueFmt) '    ', array(2,1), ','
    write(101, nameValueFmt) '    ', array(3,2), ','
    write(101, nameValueFmt) '    ', array(1,3)
    write(101, "(A)") ']'

  End Subroutine write3x3ToFile


#ifdef PYEXT
  Subroutine log_init(level, fileName, totalTime)

    Use forlog_Mod

    ! Arguments
    Integer, intent(IN) :: level
    Character(*), intent(IN) :: fileName
    Double Precision, optional, intent(IN) :: totalTime
    ! -------------------------------------------------------------------- !

    log%fileUnit = 107
    log%level = level

    open(log%fileUnit, file=trim(fileName), status='replace', recl=1000)

    ! VUMATArgs
    If (present(totalTime)) Then
      log%arg%totalTime = totalTime
    End If

  End Subroutine log_init

  Subroutine log_close()
    Use forlog_Mod
    close(log%fileUnit)
  End Subroutine log_close

#endif

End Module DGD_Mod
