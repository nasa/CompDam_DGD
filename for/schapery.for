Module schapery_mod
  ! Module for all Schapery theory-related work

Contains

  Function Schapery_damage(m, eps, Sr_old) result(Sr)
    ! Determine the current Schapery damage state variable based on the current strain state.

    Use matProp_Mod

    ! Arguments
    Type(matProps), intent(IN) :: m
    Double Precision, intent(IN) :: eps(3,3)            ! strain tensor
    Double Precision, intent(IN) :: Sr_old              ! previous converged value for Schapery damage
    Double Precision  :: Sr

    ! Locals
    Double Precision :: Sr_A, Sr_B, Sr_C
    Parameter (zero=0.d0, one=1.d0, two=2.d0, three=3.d0, four=4.d0, six=6.d0)
    ! -------------------------------------------------------------------- !

    ! Determine the reduce Schapery damage state variable, Sr
    Sr_A = ( eps(2,2)**2*m%E2*m%es(4) + (two*eps(1,2))**2*m%G12*m%gs(4) ) / two + one
    Sr_B = ( eps(2,2)**2*m%E2*m%es(3) + (two*eps(1,2))**2*m%G12*m%gs(3) ) / three
    Sr_C = ( eps(2,2)**2*m%E2*m%es(2) + (two*eps(1,2))**2*m%G12*m%gs(2) ) / six

    Sr = ( -Sr_B + SQRT( Sr_B*Sr_B - four*Sr_A*Sr_C ) ) / ( two*Sr_A )

    Sr = max(Sr, Sr_old)

    Return
  End Function Schapery_damage

  Function Schapery_reduction(Sr, c) result(reduction)
    ! Calculates the stiffness reduction according to Schapery theory.

    ! Arguments
    Double Precision, intent(IN) :: Sr           ! Reduced Schapery damage state variable
    Double Precision, intent(IN) :: c(4)         ! Schapery coefficients
    Double Precision :: reduction
    ! -------------------------------------------------------------------- !

    reduction = c(1) + c(2)*Sr + c(3)*Sr**2 + c(4)*Sr**3

    Return
  End Function Schapery_reduction

End Module schapery_mod
