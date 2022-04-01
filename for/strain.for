Module strain_mod
  ! Module for all strain-related calculations

Contains

  Subroutine Strains(F, m, DT, ndir, eps)
    ! The purpose of this subroutine is to calculate the strains (eps, plas12, inel12)
    ! for the given deformation gradient and temperature.
    !
    ! If a nonzero value for R_phi0 is provided, the strains are given in the fiber reference frame

    Use matProp_Mod
    Use forlog_Mod

    ! Arguments
    Type(matProps), intent(IN) :: m
    Double Precision, intent(IN) :: F(3,3)                    ! Deformation gradient
    Double Precision, intent(IN) :: DT                        ! Coefficients of thermal expansion and temperature change
    Integer, intent(IN) :: ndir
    Double Precision, intent(OUT) :: eps(ndir,ndir)           ! GL stain tensor, without plasticity

    ! Locals
    Parameter (zero=0.d0)
    ! -------------------------------------------------------------------- !

    ! Initialize strain tensor
    eps = zero

    ! Calculate the Green-Lagrange strain
    eps = GLStrain(F,ndir)

    ! Account for thermal strains
    If (DT /= zero) Then
      eps(1,1) = eps(1,1) - m%cte(1)*DT
      eps(2,2) = eps(2,2) - m%cte(2)*DT
      eps(3,3) = eps(3,3) - m%cte(3)*DT
    End If

    Return
  End Subroutine Strains


  Pure Function GLStrain(F, ndir)
    ! Computes the Green-Lagrange strain from the deformation gradient

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



End Module strain_mod
