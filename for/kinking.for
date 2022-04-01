Module kinking_Mod
  ! Module for all fiber kinking calculations

Contains

  Function is_kinking_possible(m, stress, phi0, shear_component) result(kp)
    ! Determines if kinking is possible
    ! returns boolean

    Use forlog_Mod
    Use matProp_Mod
    Use parameters_Mod
    Use rambergOsgood_Mod

    ! Arguments
    Type(matProps), intent(IN) :: m
    Double Precision, intent(IN) :: stress(3,3)  ! Current stress state
    Double Precision, intent(IN) :: phi0
    Integer, intent(IN) :: shear_component       ! Shear component to evaluate (either 2 or 3 for in-plane or out-of-plane respectively)

    ! Output
    Logical :: kp

    ! Locals
    Double Precision :: shear_modulus
    Double Precision :: sig_shear_const_slope    ! Shear stress @ strain where constant slope starts
    Double Precision :: sig_const_slope(3,3)     ! Stress tensor for the driving force = shear law @ strain where constant slope starts
    Double Precision :: driving_force_slope      ! Slope of driving force curve @ strain where driving force = shear law (positive)
    Double Precision :: driving_force_slope2     ! Slope of driving force curve @ strain where driving force = shear law (negative)
    Double Precision, parameter :: zero=0.d0, one=1.d0, two=2.d0
    ! -------------------------------------------------------------------- !

    ! Initialize
    sig_const_slope = stress

    ! Calculate the sig11 stress value that makes the driving force curve = shear curve at the point where the shear stiffness becomes constant
    ! Positive
    sig_const_slope(1,1) = stress(shear_component,shear_component) - &
      2*(m%shear_stress_const_slope(shear_component) - stress(1,shear_component)*COS(2*(phi0 + m%shear_strain_const_slope))) &
      /(SIN(2*(phi0 + m%shear_strain_const_slope)))

    ! Find the slope of the driving force curve using the stress calculated above
    driving_force_slope = kinking_driving_force(sig_const_slope, m%shear_strain_const_slope, phi0, shear_component, 1)

    ! Negative
    sig_const_slope(1,1) = stress(shear_component,shear_component) - &
      2*(-m%shear_stress_const_slope(shear_component) - stress(1,shear_component)*COS(2*(phi0 - m%shear_strain_const_slope))) &
      /(SIN(2*(phi0 - m%shear_strain_const_slope)))
    driving_force_slope2 = kinking_driving_force(sig_const_slope, m%shear_strain_const_slope, phi0, shear_component, 1)

    ! If the driving force slope is less than the shear law minimum slope, kinking cannot occur
    If (driving_force_slope < m%shear_const_slope_coeff(shear_component,1) .OR. driving_force_slope2 < m%shear_const_slope_coeff(shear_component,1)) Then
      kp = .FALSE.
    Else
      kp = .TRUE.
    End If

    Return
  End Function is_kinking_possible


  Function kinking_driving_force(stress, gamma, phi0, shear_component, order) result(df)
    ! Calculates the kinking driving force for the given stress state and shear strain

    ! The kinking driving force is the curve in the Considere construction representing
    ! the applied stress state. The driving force is units of stress

    ! The argument order is as follows:
    ! order=0: comput the drive force
    ! order=1: first derivative
    ! order=2: second derivative 
    ! returns boolean

    Use forlog_Mod

    ! Arguments
    Double Precision, intent(IN) :: stress(3,3)
    Double Precision, intent(IN) :: gamma
    Double Precision, intent(IN) :: phi0
    Integer, intent(IN) :: shear_component       ! Shear component to evaluate (either 2 or 3 for in-plane or out-of-plane respectively)
    Integer, intent(IN) :: order

    ! Output
    Double Precision :: df

    ! Locals
    Double Precision, parameter :: zero=0.d0, one=1.d0, two=2.d0, four=4.d0
    ! -------------------------------------------------------------------- !

    If (order == 0) Then
      df = (stress(1,1)-stress(shear_component,shear_component))/(-two)*SIN(two*(gamma+phi0)) + stress(1,shear_component)*COS(two*(gamma+phi0))

    Else If (order == 1) Then
      df = -one*(stress(1,1)-stress(shear_component,shear_component))*COS(two*(gamma+phi0)) - two*stress(1,shear_component)*SIN(two*(gamma+phi0))

    Else If (order == 2) Then
      df = two*(stress(1,1)-stress(shear_component,shear_component))*SIN(two*(gamma+phi0)) - four*stress(1,shear_component)*COS(two*(gamma+phi0))
    
    Else
      Call log%terminate('kinking_driving_force() received invalid value for argument input: ' // trim(str(order)))
    End If

    Return
  End Function kinking_driving_force


  Function kinking_stress(m, stress, phi0, shear_component) result(kinking_soln)
    ! Finds the compressive stress component at which kinking will occur. Uses Newton-Raphson to 
    ! solve the Considere construction for the current stress state. Assumes that the calling
    ! code has already checked that kinking is possible - this code triggers an error if 
    ! no solution is found.
    ! 
    ! returns the shear strain and compressive stress sig11 at kinking
    !
    ! TODO: tol, max_increments should be code parameters

    Use forlog_Mod
    Use matrixAlgUtil_Mod
    Use matProp_Mod
    Use rambergOsgood_Mod

    ! Arguments
    Type(matProps), intent(IN) :: m
    Double Precision, intent(IN) :: stress(3,3)  ! Current stress state
    Double Precision, intent(IN) :: phi0
    Integer, intent(IN) :: shear_component       ! Shear component to evaluate (either 2 or 3 for in-plane or out-of-plane respectively)

    ! Output
    Double Precision :: kinking_soln(2)

    ! Locals
    Double Precision :: shear_modulus            ! Shear modulus in the current plane (governed by shear_component)
    Double Precision :: tau                      ! Shear stress (temporary value)
    ! Locals - for Newton Raphson
    Double Precision :: x(2)                     ! Vector of unknowns (gamma, sig11)
    Double Precision :: gamma, sig               ! Current trial solution values
    Double Precision :: residual(2)              ! Vector result of evaluating the considere construction at the current iteration
    Double Precision :: stress_temp(3,3)         ! Temporary stress tensor
    Double Precision :: err
    Double Precision :: jac(2,2)                 ! Jacobian for Newton-Raphson
    Double Precision :: resjac(2,3)              ! Horizontal concatenation of residual and jacobian
    Integer :: counter                           ! Number of Newton-Raphson iterations
    Logical :: admissible_soln                   ! For checking validity of solution
    ! Locals - altenate starting point
    Double Precision :: sigc_min                 ! Minimum possible value for sig11 for kinking to be possible
    Double Precision :: candidate_gammas(10)
    Double Precision :: candidate_sigcs(10)
    Double Precision :: candidate_sp(2)
    ! Double Precision :: candidate_sp(10,2)       ! Potential (gamma, sigc) to consider when finding alt starting point
    Double Precision :: candidate_sp_errs(10)    ! Scalar error calculated for each potential starting point
    Double Precision :: temp(11)
    Integer :: i
    Double Precision :: gamma_guess

    Double Precision, parameter :: zero=0.d0, one=1.d0, two=2.d0, tol=1.d-4
    Integer, parameter :: max_increments=10
    ! -------------------------------------------------------------------- !

    ! print *, 'Starting kinking_stress()'

    ! Depending on whether we are working in plane 12 or plane 13
    shear_modulus = get_modulus(m, shear_component)

    ! Initialization
    stress_temp = stress

    ! Find starting point by sweeping through a range of tentative starting points
    If (stress(1,2) > -shear_modulus*phi0) Then
      Call linspace(MAX(zero, stress(1,2)/shear_modulus), m%shear_strain_const_slope, temp)
    Else
      Call linspace(ABS(MIN(zero, stress(1,2)/shear_modulus)), m%shear_strain_const_slope, temp)
      temp = -1*temp
    End If
    candidate_gammas = temp(2:11)
    
    StartPoint: Do i=1,10

      candidate_sigcs(i) = (ramberg_osgood_d(m, ramberg_osgood(m, candidate_gammas(i), shear_component), shear_component, 1, .TRUE.) &
        + two*stress(1,2)*SIN(two*(candidate_gammas(i)+phi0)))/(-COS(two*(candidate_gammas(i)+phi0))) + stress(2,2)
      
      ! Find the error at each candidate solution
      stress_temp(1,1) = candidate_sigcs(i)
      resjac = kinking_stress_resjac(m, stress_temp, phi0, shear_component, candidate_gammas(i))
      candidate_sp_errs(i) = Length(resjac(:,1))
    End Do StartPoint

    ! Take the minimum error as the new starting point
    i = minloc(candidate_sp_errs, 1)
    ! Find solution using Newton-Raphson
    candidate_sp = (/candidate_gammas(i), candidate_sigcs(i)/)
    x = kinking_stress_nr(m, candidate_sp, stress, phi0, shear_component)

    ! Final solution
    kinking_soln = x

    Return
  End Function kinking_stress

  Function kinking_stress_resjac(m, stress, phi0, shear_component, gamma) result(resjac)
    ! Local function to help find the residual and jacobian for kinking_stress 
    ! 
    ! returns 2x3 where first col is residual and the 2x2 is the jacobian

    Use forlog_Mod
    Use matrixAlgUtil_Mod
    Use matProp_Mod
    Use rambergOsgood_Mod

    ! Arguments
    Type(matProps), intent(IN) :: m
    Double Precision, intent(IN) :: stress(3,3)  ! Current stress state
    Double Precision, intent(IN) :: phi0
    Integer, intent(IN) :: shear_component       ! Shear component to evaluate (either 2 or 3 for in-plane or out-of-plane respectively)
    Double Precision, intent(IN) :: gamma

    ! Output
    Double Precision :: resjac(2,3)

    ! Locals
    Double Precision :: tau                      ! Shear stress (temporary value)

    Double Precision, parameter :: zero=0.d0, one=1.d0, two=2.d0
    ! -------------------------------------------------------------------- !


    ! If (debug) Then
    !   print *, 'Calling kinking_stress_resjac()'
    !   print *, 'stress', Transpose(stress)
    !   print *, 'phi0', phi0
    !   print *, 'shear_component', shear_component
    !   print *, 'gamma', gamma
    ! End If

    tau = ramberg_osgood(m, gamma, shear_component)
    resjac(1,1) = kinking_driving_force(stress, gamma, phi0, shear_component, 0) - tau                                            ! Driving force and shear law have equal value
    resjac(2,1) = kinking_driving_force(stress, gamma, phi0, shear_component, 1) - ramberg_osgood_d(m, tau, shear_component, 1, .FALSE.)   ! Driving force and shear law have equal slope

    resjac(1,2) = resjac(2,1)
    resjac(1,3) = -SIN(two*(phi0+gamma))/two
    resjac(2,2) = kinking_driving_force(stress, gamma, phi0, shear_component, 2) - ramberg_osgood_d(m, tau, shear_component, 2, .FALSE.)
    ! print *, 'term1', kinking_driving_force(stress, gamma, phi0, shear_component, 2)
    ! print *, 'term2', ramberg_osgood_d(m, tau, shear_component, 2, .FALSE.)
    resjac(2,3) = -COS(two*(phi0+gamma))

    Return
  End Function kinking_stress_resjac

  Function kinking_stress_nr(m, x0, stress, phi0, shear_component) result(x)
    ! Local function: use Newton-Raphson to find the kinking stress and shear strain

    Use forlog_Mod
    Use matrixAlgUtil_Mod
    Use matProp_Mod
    Use rambergOsgood_Mod

    ! Arguments
    Type(matProps), intent(IN) :: m
    Double Precision, intent(IN) :: x0(2)        ! Starting point
    Double Precision, intent(IN) :: stress(3,3)  ! Current stress state
    Double Precision, intent(IN) :: phi0
    Integer, intent(IN) :: shear_component       ! Shear component to evaluate (either 2 or 3 for in-plane or out-of-plane respectively)

    ! Output
    Double Precision :: x(2)

    ! Locals
    Double Precision :: gamma, sig               ! Current trial solution values
    Double Precision :: stress_temp(3,3)
    Double Precision :: resjac(2,3)              ! Horizontal concatenation of residual and jacobian
    Double Precision :: residual(2)
    Double Precision :: err
    Double Precision :: jac(2,2)                 ! Jacobian for Newton-Raphson
    Integer :: counter

    Double Precision, parameter :: zero=0.d0, one=1.d0, two=2.d0, tol=1.d-6
    Integer, parameter :: max_increments=100
    ! -------------------------------------------------------------------- !

    ! Initialize
    x = x0
    stress_temp = stress
    counter = 0

    ! If (.TRUE.) Then
    !   print *, 'Calling kinking_stress_nr()'
    !   print *, 'x0', x0
    !   print *, 'stress', Transpose(stress)
    !   print *, 'phi0', phi0
    !   print *, 'shear_component', shear_component
    ! End If

    ! Newton-Raphson loop
    NR: Do

      ! If (.TRUE.) Then
      !   print *, counter
      !   print *, x
      ! End If
      gamma = x(1)
      sig = x(2)
      stress_temp(1,1) = sig
      
      ! Compute the residual
      resjac = kinking_stress_resjac(m, stress_temp, phi0, shear_component, gamma)
      residual = resjac(:,1)
      ! print *, 'residual'
      ! print *, residual
      err = Length(residual)
      ! print *, 'err', err

      ! Check for convergence or other limits
      If (err < tol) Then
        Call log%debug('Number kinking stress iterations ' // trim(str(counter)))
        Exit NR
      
      Else If (counter > max_increments) Then
        Call log%debug('kinking_stress() reached maximum number of iterations')
        Exit NR
        ! print *, 'x0', x0
        ! print *, 'stress', stress
        ! print *, 'phi0', phi0
        ! print *, 'shear_component', shear_component
        ! Call log%error('kinking_stress() reached maximum number of iterations')

      End If
     
      ! Compute the jacobian
      jac = resjac(:,2:3)
      ! print *, 'jac'
      ! print *, TRANSPOSE(jac)

      ! update the guess
      x = x - MATMUL(MInverse2x2(jac), residual)
      
      counter = counter + 1
    End Do NR

    Return
  End Function kinking_stress_nr

  Function kinking_stress_nrslope(m, gamma0, stress, phi0, shear_component) result(gamma)
    ! Local function: use Newton-Raphson to find the gamma that solves the slope eqn

    Use forlog_Mod
    Use matrixAlgUtil_Mod
    Use matProp_Mod
    Use rambergOsgood_Mod

    ! Arguments
    Type(matProps), intent(IN) :: m
    Double Precision, intent(IN) :: gamma0
    Double Precision, intent(IN) :: stress(3,3)  ! Current stress state
    Double Precision, intent(IN) :: phi0
    Integer, intent(IN) :: shear_component       ! Shear component to evaluate (either 2 or 3 for in-plane or out-of-plane respectively)

    ! Output
    Double Precision :: gamma

    ! Locals
    Double Precision :: err
    Integer :: counter
    Double Precision :: tau                      ! Shear stress (temporary value)
    Double Precision :: residual
    Double Precision :: ff1

    Double Precision, parameter :: zero=0.d0, one=1.d0, two=2.d0, tol=1.d-6
    Integer, parameter :: max_increments=200
    ! -------------------------------------------------------------------- !

    gamma = gamma0
    counter = 0

    NR: Do counter=1,max_increments
      tau = ramberg_osgood(m, gamma, shear_component)
      residual = kinking_driving_force(stress, gamma, phi0, shear_component, 1) - ramberg_osgood_d(m, tau, shear_component, 1, .TRUE.)
      err = ABS(residual)
      If (err < tol) Then
        Exit NR
      Else If (counter > max_increments) Then
        Call log%error('kinking_stress_nrslope() reached maximum number of iterations')
        Exit NR
      End If
      ff1 = kinking_driving_force(stress, gamma, phi0, shear_component, 2) - ramberg_osgood_d(m, tau, shear_component, 2, .TRUE.)
      gamma = gamma - residual/ff1 ! Newton-Raphson equation
    End Do NR

    Return
  End Function kinking_stress_nrslope


  Function kinking_failure_criterion(m, stress, shear_strain, phi0, shear_component) result(fc)
    ! Computes a simple failure index for fiber kinking based on the shear strain

    Use forlog_Mod
    Use matrixAlgUtil_Mod
    Use matProp_Mod
    Use rambergOsgood_Mod

    ! Arguments
    Type(matProps), intent(IN) :: m
    Double Precision, intent(IN) :: stress(3,3)  ! Current stress state
    Double Precision, intent(IN) :: shear_strain  ! Current engineering shear strain (in-misaligned frame)
    Double Precision, intent(IN) :: phi0
    Integer, intent(IN) :: shear_component       ! Shear component to evaluate (either 2 or 3 for in-plane or out-of-plane respectively)

    ! Output
    Double Precision :: fc

    ! Locals
    Double Precision :: kink_stress_res(2)       ! (sig11, gamma) at which kinking will occur
    Double Precision :: sig_kink, gamma_kink     ! stress and shear strain when kinking will occur
    Double Precision :: test_stress(3,3)

    Double Precision, parameter :: zero=0.d0
    ! -------------------------------------------------------------------- !

    If (is_kinking_possible(m, stress, phi0, shear_component)) Then
      kink_stress_res = kinking_stress(m, stress, phi0, shear_component)
      ! print *, 'kink_stress_res', kink_stress_res
      gamma_kink = kink_stress_res(1)
      fc = shear_strain/gamma_kink
      ! print *, 'fc', fc
      ! Call log%error('stop')
    Else
      fc = zero
    End If

    Return
  End Function kinking_failure_criterion


  Function fiber_break_fi(m, stress11, phi) result(fc)
    ! Calculates the failure index for fiber breaks under compression
    ! Strictly informational at this point

    Use matProp_Mod
    Use rambergOsgood_Mod

    ! Arguments
    Type(matProps), intent(IN) :: m
    Double Precision, intent(IN) :: stress11
    Double Precision, intent(IN) :: phi          ! Current rotation

    ! Output
    Double Precision :: fc

    ! Locals
    Double Precision :: bending_strain
    Double Precision :: axial_strain

    Double Precision, parameter :: zero=0.d0, two=2.d0
    Double Precision, parameter :: fiber_diameter=0.0005d0
    Double Precision, parameter :: fiber_rupture_strain=-0.015d0
    ! -------------------------------------------------------------------- !

    bending_strain = -fiber_diameter*TAN(ABS(phi/two))/m%w_kb
    axial_strain = stress11/m%E1
    fc = (bending_strain+axial_strain)/fiber_rupture_strain
    
    Return
  End Function fiber_break_fi

End Module kinking_Mod