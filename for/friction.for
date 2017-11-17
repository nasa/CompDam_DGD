Module friction_mod
  ! Determines whether the crack is slipping or sticking, and finds the resulting traction

Contains

  Pure Subroutine crack_traction_and_slip(delta, Pen, slide_old, slide, mu_L, mu_T, dmg, AdAe, stress, Sliding)
    ! Returns the traction on the crack and the amount of slip between the crack surfaces

    Include 'vaba_param.inc'

    ! Arguments
    Double Precision, intent(IN) :: delta(3)      ! cohesive displacement-jump [length]
    Double Precision, intent(IN) :: slide_old(2)  ! previous amount of displacement slip [length]
    Double Precision, intent(IN) :: Pen(3)        ! cohesive penalty stiffnesses [stress/length]
    Double Precision, intent(IN) :: mu_L, mu_T    ! coefficients of friction in longitudinal and transverse directions
    Double Precision, intent(IN) :: dmg, AdAe     ! cohesive damage variables
    Logical, intent(IN) :: Sliding

    ! Output
    Double Precision, intent(OUT) :: stress(3)    ! cohesive tractions [stress]
    Double Precision, intent(OUT) :: slide(2)     ! current amount of displacement slip [length]

    ! Local
    Double Precision :: stress_old(3), SD(3), shear, mu
    Double Precision, parameter :: zero=0.d0, one=1.d0
    ! -------------------------------------------------------------------- !

    If (delta(2) .LT. zero) Then  ! Closed crack

      ! Tractions using the old slip displacements
      stress_old(1) = Pen(1)*(delta(1) - slide_old(1))
      stress_old(2) = Pen(2)*delta(2)
      stress_old(3) = Pen(3)*(delta(3) - slide_old(2))

      If (.NOT. Sliding) Then  ! Sticking

        stress(1) = (one - dmg)*Pen(1)*delta(1) + AdAe*stress_old(1)
        stress(2) = stress_old(2)
        stress(3) = (one - dmg)*Pen(3)*delta(3) + AdAe*stress_old(3)

        slide(1) = slide_old(1)
        slide(2) = slide_old(2)

      Else If (Sliding) Then  ! Sliding

        ! Shear stress on the crack, using the old slip displacements
        shear = SQRT(stress_old(1)**2 + stress_old(3)**2)

        mu = (mu_L*stress_old(1)**2 + mu_T*stress_old(3)**2)/shear**2  ! Combined coefficient of friction

        SD(1) = -mu*stress_old(2)*stress_old(1)/shear
        SD(2) = stress_old(2)
        SD(3) = -mu*stress_old(2)*stress_old(3)/shear

        stress(1) = (one - dmg)*Pen(1)*delta(1) + AdAe*SD(1)
        stress(2) = SD(2)
        stress(3) = (one - dmg)*Pen(3)*delta(3) + AdAe*SD(3)

        slide(1) = (mu*stress_old(2) + shear)/Pen(1)*stress_old(1)/shear + slide_old(1)
        slide(2) = (mu*stress_old(2) + shear)/Pen(3)*stress_old(3)/shear + slide_old(2)

      End If

    Else  ! Open crack

      stress(:) = (one - dmg)*Pen(:)*delta(:)

      slide(1) = delta(1)
      slide(2) = delta(3)

    End If

    Return
  End Subroutine crack_traction_and_slip


  Pure Function crack_is_sliding(delta, Pen, slide, mu_L, mu_T) result(Sliding)
    ! Determines whether the crack surfaces are sliding or sticking

    Include 'vaba_param.inc'

    ! Input
    Double Precision, intent(IN) :: delta(3)      ! cohesive displacement-jump [length]
    Double Precision, intent(IN) :: slide(2)      ! previous amount of displacement slip [length]
    Double Precision, intent(IN) :: Pen(3)        ! cohesive penalty stiffnesses [stress/length]
    Double Precision, intent(IN) :: mu_L, mu_T    ! coefficients of friction in longitudinal and transverse directions

    ! Output
    Logical :: Sliding

    ! Locals
    Double Precision :: stress_old(3), shear_squared, mu
    Double Precision, parameter :: zero=0.d0
    ! -------------------------------------------------------------------- !

    If (delta(2) .LT. zero) Then  ! Closed crack

      ! Stresses using the old slip displacements
      stress_old(1) = Pen(1)*(delta(1) - slide(1))
      stress_old(2) = Pen(2)*delta(2)
      stress_old(3) = Pen(3)*(delta(3) - slide(2))

      ! Shear stress on the crack (squared), using the old slip displacements
      shear_squared = stress_old(1)**2 + stress_old(3)**2

      If (shear_squared .GT. zero) Then
        mu = (mu_L*stress_old(1)**2 + mu_T*stress_old(3)**2)/shear_squared  ! Combined coefficient of friction

        ! If the friction traction is greater than the shear traction,
        If ((mu*stress_old(2))**2 .GE. shear_squared) Then
          ! the crack surfaces are sticking.
          Sliding = .False.
        Else
          ! the crack surfaces are sliding.
          Sliding = .True.
        End IF

      Else
        Sliding = .False.

      End If

    Else  ! Open crack
      Sliding = .True.

    End If

    Return
  End Function crack_is_sliding

End Module friction_mod
