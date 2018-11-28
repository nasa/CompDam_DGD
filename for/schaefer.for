Module schaefer_Mod
  ! Module for all Schaefer theory-related work

Contains

  subroutine schaefer(m, schaefer_a6, schaefer_b2, schaefer_n, schaefer_A, eps, eps_old, ndir, nshr, eps_vec_p_old, fp_old) 

    use matProp_Mod
    use matrixAlgUtil_Mod
    use stress_Mod
    Use forlog_Mod
    
    ! Inputs only
    Type(matProps), intent(IN) :: m
    Double PRECISION, intent(IN) :: schaefer_a6    ! Schaefer Method material property 
    Double PRECISION, intent(IN) :: schaefer_b2    ! Schaefer Method material property 
    Double PRECISION, intent(IN) :: schaefer_n ! fitting parameters for schaefre
    Double PRECISION, intent(IN) :: schaefer_A ! fitting parameters for schaefre
    Double PRECISION, intent(IN) :: eps(ndir, ndir)  ! Total strain (green-lagrange) 
    Double PRECISION, intent(IN) :: eps_old(ndir, ndir)  ! Total strain from previous iteration (green-lagrange) 
    Integer, intent(IN) :: ndir ! number direct components in stress and strain count(sig11,...)
    Integer, intent(IN) :: nshr ! number of shear components count(sig12, ....)
    ! Inputs that are updated in place
    Double PRECISION, intent(INOUT) :: eps_vec_p_old(ndir + nshr)   ! previously calculated total plastic strain. Updated at final step to new old
    Double PRECISION, intent(INOUT) :: fp_old ! yield function value

    Double PRECISION :: f ! yield criterion
    ! stress components
    Double PRECISION :: S22
    Double PRECISION :: S12
    ! additional terms
    Double PRECISION :: C(ndir + nshr, ndir + nshr) ! Stiffness Matrix
    Double PRECISION :: th ! theta is a variable defined for convenenience. The term f - schaefer_b2 * S22 shows up frequently. This term helps save typing
    Double PRECISION :: th2 ! th ** 2 (captured for convenience)
    Double PRECISION :: th3 ! th ** 3 (captured for convenience)
    Double PRECISION :: S(ndir + nshr) ! vector of stress matrix S 
    Double PRECISION :: dS(ndir + nshr) ! vector of incremental stress 
    Double PRECISION :: dfdS(ndir + nshr) ! derivative of scalar FUNCTION with respect to S
    Double PRECISION :: dfdSdfdS(ndir + nshr, ndir + nshr) ! outer product of dfdS 
    Double PRECISION:: dS22ddEp(ndir + nshr) ! derivative of stress component(S22) with respect to incremental plastic strain
    Double PRECISION:: dS12ddEp(ndir + nshr) ! derivative of stress compoentn (S12) with respect to incremental plastic strain
    Double PRECISION:: dfddEp(ndir + nshr) ! derivative of yield criterion with respect to incremental plastic strain
    Double PRECISION:: eye(ndir + nshr, ndir + nshr) ! identity matrix
    Double PRECISION:: dEp(ndir + nshr) ! incremental plastic strain
    Double PRECISION:: dEpUpdate(ndir + nshr) ! incremental plastic strain
    Double PRECISION:: eps_vec(ndir + nshr) ! strain matrix represented as vector (green-lagrange)
    Double PRECISION:: eps_old_vec(ndir + nshr) ! strain matrix represented as vector from previous iteration (green-lagrange)
    Double PRECISION:: dE(ndir + nshr) ! total incremental strain as vector
    Double PRECISION:: residual !! value keeps track of J0 == 0
    Double PRECISION:: dSddEp(ndir + nshr) ! derivative of stress with respect to incremntal plastic strain (dEp). Only holds ith derivative at a time (hence why its not a 6 by 6)
    Double PRECISION:: dfdSdfdSd22Ep ! Holds the values of hte ith derivative of [(dfdS)(dfdS)Transpose](2,2) with respect to dEpi
    Double PRECISION:: dfdSdfdSd24Ep! Holds the values of hte ith derivative of (dfdS)(dfdS)T24 with respect to dEpi
    Double PRECISION:: dfdSdfdSd44Ep! Holds the values of hte ith derivative of (dfdS)(dfdS)T44 with respect to dEpi
    Double PRECISION:: dfdSdfdSdEp(ndir + nshr, ndir + nshr)! Corresponds to the specific (ith derivative) and is the 6 by 6 matrix which is primarily zeros but incorporates the values of dfdSdfdSdijEp into itself
    Double PRECISION:: J1a(ndir + nshr) ! Component of J1 matrix used in newton raphson
    Double PRECISION:: J1b(ndir + nshr)! Component of J1 matrix used in newton raphson
    Double PRECISION:: J1c(ndir + nshr) ! Component of J1 matrix used in newton raphson
    Double PRECISION:: dS12 ! local variable used in ith column calc of J1 matrix. its a deriv val with respect to incremental plastic strain
    Double PRECISION:: dS22 ! local variable used in ith column calc of J1 matrix. its a deriv val with respect to incremental plastic strain
    Double PRECISION:: df! local variable used in ith column calc of J1 matrix. its a deriv val with respect to incremental plastic strain
    Double PRECISION:: schaefer_nr_tolerance ! tolerance value used to determine if convergence has occurred in newton raphson loop
    Integer :: counter_limit ! maximum number of Newton-Raphson loops before divergence is assumed
    Integer :: counter ! counter of iterations done thus far in newton raphson
    ! newton raphson equation updates dEp according to the following
    ! dEp = dEp - inv(J1) * J0
    Double PRECISION:: J0(ndir + nshr) ! 
    Double PRECISION:: J1(ndir + nshr, ndir + nshr)! Corresponds to the specific (ith derivative) and is the 6 by 6 matrix which is primarily zeros but incorporates the values of dfdSdfdSdijEp into itself
    Double PRECISION:: J1column(ndir + nshr)

    Double PRECISION, Parameter :: zero=0.d0, one=1.d0, two=2.d0, three=3.d0, four=4.d0, six=6.d0

    schaefer_nr_tolerance = 1e-6
    counter_limit = 10000000
    ! Build the stiffness matrix
    C = StiffFunc(ndir+nshr, m%E1, m%E2, m%E3, m%G12, m%G13, m%G23, m%v12, m%v13, m%v23, zero, zero, zero)
    eye = zero; Do I = 1,ndir+nshr; eye(I,I) = one; End Do ! Identity Matrix 

    ! Assume dEp = [0] (iniitally)
    dEp = zero

    ! convert matrices to vectors (3+nshrx1)
    eps_vec = Matrix2Vec(eps, nshr)
    eps_old_vec = Matrix2Vec(eps_old, nshr)
    dE = eps_vec - eps_old_vec

    counter = 0
    DO 
      counter = counter + 1
      S = Matrix2Vec(Hooke(C, Vec2Matrix(eps_vec - eps_vec_p_old - dEp), nshr), nshr);
      dS = Matrix2Vec(Hooke(C, Vec2Matrix(dE - dEp), nshr), nshr);

      ! get components of stress vector for convenience
      S22 = S(2)
      S12 = S(ndir + 1) 
      ! get yield FUNCTION from components of S and fitting parameters
      f = Getf(S, schaefer_a6, schaefer_b2, ndir, nshr)
      ! f is zero more than likely indicates a 0 stress state.
      ! In this case j predict that no new plastic strain (no update to eps_plastic)
      ! and break sub
      IF (f .eq. zero) THEN
        fp_old = f
        Return
      END IF



      ! get terms of f deriviates wtih respect to S 
      ! Partial of f with respect to Stress
      dfdS = zero
      dfdS(2) = S22 / (f - schaefer_b2 * S22) + schaefer_b2
      dfdS(ndir + 1) = schaefer_a6 * S12 / (f - schaefer_b2 * S22)

      ! dfdf * dfdsT (outer multiply the same vector with self)
      dfdSdfdS = OuterProduct(dfdS, dfdS)

      ! calc init J0
      J0 = GetJ0(schaefer_n, schaefer_A, f, dfdSdfdS, dS, dEp, ndir, nshr)

      ! derivative of stress with respect to incremental plastic strain
      ! because S = C: (eps_vec - Eoldp - dEp) and C, eps_vec, Eoldp are constants
      ! youre left with the negative of the stiffness matrix
      dS22ddEp = -C(2, :)
      dS12ddEp = -C(ndir + 1, :) 

      ! for convenience defined th (theta) to be f - schaefer_b2 * S22. This term shows up a bunch
      th = f- schaefer_b2 * S22
      th2 = th ** 2
      th3 = th ** 3

      ! Partial of f with respect to plastic strain. 
      dfddEp = (S22 / th + schaefer_b2) * dS22ddEp + ((schaefer_a6 * S12) / th) * dS12ddEp

      ! zero out J1
      J1 = zero

      ! J1 is a (ndir + nshr) x (ndir + nshr) matrix. It is the derivative of J0 with respect to dEpi.
      ! A column "i" represents the terms corresponding to the dEpi derivative
      DO i = 1, ndir + nshr, 1
        ! fill in dSdEp
        dSddEp = zero
        dSddEp(2) = dS22ddEp(i)
        dSddEp(ndir + 1) = dS12ddEp(i)

        ! local derivatives
        dS12 = dS12ddEp(i)
        dS22 = dS22ddEp(i)
        df = dfddEp(i)

        ! get the the term corresponding to the ith partial derivative of dfdSdfdS
        ! with respeoct to Epi
        ! this derivative corresponds to the (2, 2) 
        dfdSdfdSd22Ep = (schaefer_b2 + S22 / th) * ( (2 * dS22) / th + (2 * (schaefer_b2 * dS22 -df) * S22) / th2) 
        ! this derivative corresponds to the (2, 4) 
        dfdSdfdSd24Ep = (schaefer_a6 * dS12/th) * (schaefer_b2 + S22/th) +&
            (schaefer_a6 * S12 / th2) * (schaefer_b2 + S22 / th) * (schaefer_b2 * dS22 - df) +&
            (schaefer_a6 * S12 / th) * (dS22 / th + (schaefer_b2 * dS22 - df)* S22 / th2)
        ! this derivative corresponds to the (4, 4) 
        dfdSdfdSd44Ep =  (2 * schaefer_a6**2 * S12 * dS12) / th2 +&
            (schaefer_a6**2 * (2 * schaefer_b2 * dS22 - 2 * df ) * S12**2) / th3
        ! Fill with all zeros and then fill its componenets as theyre calculated
        dfdSdfdSdEp = zero
        dfdSdfdSdEp(2, 2) = dfdSdfdSd22Ep
        dfdSdfdSdEp(2, ndir + 1) = dfdSdfdSd24Ep
        dfdSdfdSdEp(ndir + 1, 2) = dfdSdfdSd24Ep
        dfdSdfdSdEp(ndir + 1, ndir + 1) = dfdSdfdSd44Ep

        J1a = (schaefer_n - 1) * f**(schaefer_n - 2) * df * MATMUL(dfdSdfdS, dS)
        J1b = f**(schaefer_n - 1) * MATMUL(dfdSdfdSdEp, dS) 
        J1c = f**(schaefer_n - 1) * MATMUL(dfdSdfdS, dSddEp)
        ! th is used upstream in calculation of J1 and because its in the denominator
        ! has the potenial to drive the result to inf (NaN)
        ! This occurs primarily because th is squared and cubed
        ! and thus reduces potenially small numbers to even smaller numbers.
        ! Whent this is the case, just predict no changed in plasticity
        ! if th^3 is the most likely to be zero compared to th^2 or th because its a higher exponent.
        ! and if its zero then there will be issues
        IF ((th3 .eq. zero)) THEN
          fp_old = f
          Return
        END IF
        J1(:, i) = eye(:, i) - schaefer_n * schaefer_A * (J1a + J1b + J1c)
        IF ((J1(1, i) .ne. J1(1, i)) .or. (J1(2, i) .ne. J1(2, i)) .or. (J1(3, i) .ne. J1(3, i)) .or. (J1(4, i) .ne. J1(4, i)) .or. (J1(5, i) .ne. J1(5, i)) .or. (J1(6, i) .ne. J1(6, i))) THEN
          PRINT *, '***************************************************'
          PRINT *, '***************************************************'
          PRINT *, 'th'
          PRINT *, th
          PRINT *, 'th**2'
          PRINT *, th**2
          PRINT *, 'th**3'
          PRINT *, th**3
          PRINT *, 'dfdSdfdSd22Ep'
          PRINT *, dfdSdfdSd22Ep
          PRINT *, 'dfdSdfdSd24Ep'
          PRINT *, dfdSdfdSd24Ep
          PRINT *, 'dfdSdfdSd44Ep'
          PRINT *, dfdSdfdSd44Ep
          PRINT *, 'dfdSdfdSdEp'
          PRINT *, dfdSdfdSdEp
          PRINT *, 'dSddEp(2)'
          PRINT *, dSddEp(2)
          PRINT *, 'dSddEp(ndir + 1)'
          PRINT *, dSddEp(ndir + 1)
          PRINT *, 'J1a'
          PRINT *, J1a
          PRINT *, 'J1b'
          PRINT *, J1b
          PRINT *, 'J1c'
          PRINT *, J1c
          PRINT *, 'eye(:, i)'
          PRINT *, eye(:, i)
          PRINT *, 'J1(:, i)'
          PRINT *, J1(:, i)
          PRINT *, '***************************************************'
          PRINT *, '***************************************************'
        END IF
      END DO
      ! update according dEp according to newton raphson
      IF (ndir + nshr .eq. 6) THEN
        dEpUpdate = MATMUL(MInverse6x6(J1), J0)
        ! check that numbers arent na
        IF (dEpUpdate(2) .ne. dEpUpdate(2) .or. dEpUpdate(4) .ne. dEpUpdate(4)) THEN
          Call log%debug_str('a NaN was calculated this is whack')
          PRINT *, 'counter'
          PRINT *, counter
          PRINT *, 'residual'
          PRINT *, residual
          PRINT *, 'eps'
          PRINT *, eps
          ! J0 = GetJ0(schaefer_n, schaefer_A, f, dfdSdfdS, dS, dEp, ndir, nshr)
          PRINT *, 'schaefer_n'
          PRINT *, schaefer_n
          PRINT *, 'schaefer_A'
          PRINT *, schaefer_A
          PRINT *, 'f'
          PRINT *, f
          PRINT *, 'dfdS'
          PRINT *, dfdS
          PRINT *, 'dfdSdfdS'
          PRINT *, dfdSdfdS
          PRINT *, 'dS'
          PRINT *, dS
          PRINT *, 'J0'
          PRINT *, J0
          PRINT *, 'J1'
          PRINT *, J1
          PRINT *, 'ndir'
          PRINT *, ndir
          PRINT *, 'nshr'
          PRINT *, nshr
          PRINT *, 'S '
          PRINT *, S 
          ! dfdS(2) = S22 / (f - schaefer_b2 * S22) + schaefer_b2
          ! dfdS(ndir + 1) = schaefer_a6 * S12 / (f - schaefer_b2 * S22)
          PRINT *, 'S22 '
          PRINT *, S22 
          PRINT *, 'S12 '
          PRINT *, S12 
          PRINT *, 'schaefer_b2'
          PRINT *, schaefer_b2
          PRINT *, 'schaefer_a6'
          PRINT *, schaefer_a6
          PRINT *, 'dEp'
          PRINT *, dEp
          PRINT *, 'dEpUpdate'
          PRINT *, dEpUpdate
          Call log%error('A Na was calculated. Break out')

        END IF
        dEp = dEp - MATMUL(MInverse6x6(J1), J0)
      ELSE IF (ndir + nshr .eq. 4) THEN
        dEp = dEp - MATMUL(MInverse4x4(J1), J0)
      ELSE 
        Call log%error('Schaefer material prepeak nonlinearity model only supports ndir + nshr=6 and 4 elements')
      END IF

      ! The newton raphson continually updates J0 (and as it converges J0) should
      ! approach 0. Determine its norm (and refer to this as the residual)
      residual = SQRT(DOT_PRODUCT(J0, J0))
      ! once we're within a certain tolerance break the loop
      if (residual .lt. schaefer_nr_tolerance) THEN
        ! If f is larger than previous step then additional plastic strain is predicted
        if (f > fp_old) THEN
          eps_vec_p_old = eps_vec_p_old + dEp
          fp_old = f
        ELSE 
          Call log%debug_str('No change predicted')
        END IF
        Return
      END IF
      IF (MODULO(counter, counter_limit) .eq. 0) THEN
        Call log%debug_str('Maybe too many iterations occurred during newton-raphson convergence in schaefer sub')
        Call log%debug_int('counter', counter)
      END IF
    END DO
  END subroutine schaefer

  FUNCTION Getf(S, schaefer_a6, schaefer_b2, ndir, nshr) result(f)
    Double PRECISION, intent(IN) :: schaefer_a6
    Double PRECISION, intent(IN) :: schaefer_b2
    Integer, intent(IN) :: ndir, nshr
    Double PRECISION, intent(IN) :: S(ndir + nshr, 1)
    
    Double PRECISION :: f
    Double PRECISION :: S12
    Double PRECISION :: S22

    S22 = S(2, 1)
    S12 = S(ndir + 1, 1)
    f = SQRT(S22**2 + schaefer_a6 * S12**2) + schaefer_b2 * S22
    Return
  End FUNCTION Getf

  ! Calculates J0 and returns it. J0 is a variable used in the newton Raphson
  ! determination of dEp
  FUNCTION GetJ0(schaefer_n, schaefer_A, f, dfdSdfdS, dS, dEp, ndir, nshr) result(J0)
    Double PRECISION, intent(IN) :: schaefer_n
    Double PRECISION, intent(IN) :: schaefer_A
    Double PRECISION, intent(IN) :: f
    Integer, intent(IN) :: ndir, nshr
    Double PRECISION, intent(IN) :: dfdSdfdS(ndir + nshr, ndir + nshr)
    Double PRECISION, intent(IN) :: dS(ndir + nshr)
    Double PRECISION, intent(IN) :: dEp(ndir + nshr)
    Double PRECISION :: J0(ndir + nshr)

    J0 = dEp - schaefer_n * schaefer_A * f**(schaefer_n - 1) * MATMUL(dfdSdfdS, dS)
    Return
  End FUNCTION GetJ0
End Module schaefer_mod
