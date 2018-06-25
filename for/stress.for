Module stress_Mod
  ! Module for building stiffness tensors, tensor transformations, and converting different definitions of stress

Contains

  Pure Function StiffFunc(NTENS, E1, E2, E3, G12, G13, G23, v12, v13, v23, d1, d2, d3) result(Stiff)
    ! Constructs the damaged orthotropic stiffness tensor based on the calculated damage state variables

    ! Input
    Double Precision, intent(IN) :: E1,E2,E3,G12,G13,G23,v12,v23,v13,d1,d2,d3
    Integer, intent(IN) :: NTENS

    ! Output
    Double Precision :: Stiff(NTENS,NTENS)

    ! Locals
    Double Precision :: v21, v31, v32, delta
    Double Precision, parameter :: zero=0.d0, one=1.d0, two=2.d0
    ! -------------------------------------------------------------------- !

    v21 = v12*E2/E1
    If (NTENS > 4) Then
       v31 = v13*E3/E1
       v32 = v23*E3/E2
    Else
       v31 = zero
       v32 = zero
    End If

    delta = one/(one - v12*v21*(one-d1)*(one-d2) - v23*v32*(one-d2)*(one-d3) - v31*v13*(one-d1)*(one-d3) - two*v12*v23*v31*(one-d1)*(one-d2)*(one-d3))
    Stiff = zero

    Stiff(1,1) = delta*(one-d1)*(one - v23*v32*(one-d2)*(one-d3))*E1
    Stiff(2,2) = delta*(one-d2)*(one - v13*v31*(one-d1)*(one-d3))*E2

    Stiff(1,2) = delta*(one-d1)*(one-d2)*(v12 + v32*v13*(one-d3))*E2
    Stiff(2,1) = Stiff(1,2)

    If (NTENS == 3) Then  ! Shell element in Abaqus/Standard
      Stiff(3,3) = (one-d1)*(one-d2)*G12
    Else
      Stiff(3,3) = delta*(one-d3)*(one - v12*v21*(one-d1)*(one-d2))*E3
      Stiff(4,4) = (one-d1)*(one-d2)*G12
    End If

    If (NTENS > 4) Then ! Solid element
      Stiff(1,3) = delta*(one-d1)*(one-d3)*(v13+v12*v23*(one-d2))*E3
      Stiff(3,1) = Stiff(1,3)
      Stiff(2,3) = delta*(one-d2)*(one-d3)*(v23+v21*v13*(one-d1))*E3
      Stiff(3,2) = Stiff(2,3)
      Stiff(5,5) = (one-d2)*(one-d3)*G23
      Stiff(6,6) = (one-d1)*(one-d3)*G13
    End If

    Return
  End Function StiffFunc


  Subroutine StiffFuncNL(m, ndir, nshr, d1, d2, d3, eps, stiff, sr)
    ! Constructs the damaged orthotropic stiffness tensor accounting for damage state variables and recoverable 
    ! pre-peak nonlinearities including Schapery and elastic fiber nonlinearity

    Use matProp_Mod
    Use schapery_mod

    ! Arguments
    Type(matProps), intent(IN) :: m
    Integer, intent(IN) :: ndir, nshr
    Double Precision, intent(IN) :: d1, d2, d3                       ! Damage state variables
    Double Precision, intent(IN) :: eps(ndir,ndir)                   ! Green-Lagrange strain tensor
    Double Precision, intent(OUT) :: stiff(ndir+nshr,ndir+nshr)      ! Stiffness tensor
    Double Precision :: sr                                           ! Shapery state variable (no intent is specified so that zero can be used when schapery is N/A)

    ! Locals
    Double Precision :: E1
    Double Precision, parameter :: zero=0.d0, one=1.d0, two=2.d0
    ! -------------------------------------------------------------------- !

    If (m%Schapery) Then
      sr = Schapery_damage(m, eps, sr)
    Else
      sr = zero
    End If

    ! Elastic fiber nonlinearity
    E1 = m%E1*(1+m%cl*eps(1,1))

    ! Update stiffness matrix
    stiff = StiffFunc(ndir+nshr, E1, m%E2*Schapery_reduction(sr, m%es), m%E3, m%G12*Schapery_reduction(sr, m%gs), m%G13, m%G23, m%v12, m%v13, m%v23, d1, d2, d3)

    Return
  End Subroutine StiffFuncNL


  Function StiffRot(C, NTENS, theta) result(StiffOut)
    ! Rotates the stiffness matrix by the specified angle (radians)

    Use forlog_Mod

    ! Input
    Double Precision, intent(IN) :: C(NTENS,NTENS), theta
    Integer, intent(IN) :: NTENS

    ! Output
    Double Precision :: StiffOut(NTENS,NTENS)

    ! Locals
    Double Precision :: m, n
    Double Precision :: zero, one, two, four
    Parameter (zero=0.d0, one=1.d0, two=2.d0, four=4.d0)
    ! -------------------------------------------------------------------- !

    m = cos(theta)
    n = sin(theta)

    StiffOut = 0

    ! Depends on element type
    If (NTENS > 4) Then
      StiffOut(1,1) = m**four*C(1,1) + 2*m**two*n**two*(C(1,2) + 2*C(4,4)) + n**four*C(2,2)
      StiffOut(1,2) = n**two*m**two*(C(1,1) + C(2,2) - four*C(4,4)) + (n**four + m**four)*C(1,2)
      StiffOut(1,3) = m**two*C(1,3) + n**two*C(2,3)
      StiffOut(1,4) = n*m*(m**two*(C(1,1) - C(1,2) - two*C(4,4)) + n**two*(C(1,2) - C(2,2) + two*C(4,4)))
      StiffOut(2,1) = StiffOut(1,2)
      StiffOut(2,2) = n**four*C(1,1) + 2*m**two*n**two*(C(1,2) + 2*C(4,4)) + m**four*C(2,2)
      StiffOut(2,3) = n**two*C(1,3) + m**two*C(2,3)
      StiffOut(2,4) = n*m*(n**two*(C(1,1) - C(1,2) - two*C(4,4)) + m**two*(C(1,2) - C(2,2) + two*C(4,4)))
      StiffOut(3,1) = StiffOut(1,3)
      StiffOut(3,2) = StiffOut(2,3)
      StiffOut(3,3) = C(3,3)
      StiffOut(3,4) = m*n*(C(1,3) - C(2,3))
      StiffOut(4,1) = StiffOut(1,4)
      StiffOut(4,2) = StiffOut(2,4)
      StiffOut(4,3) = StiffOut(3,4)
      StiffOut(4,4) = n**two*m**two*(C(1,1) - two*C(1,2) + C(2,2)) + C(4,4)*(n**two - m**two)**two
      StiffOut(5,5) = m**two*C(5,5) + n**two*C(6,6)
      StiffOut(5,6) = m*n*(C(6,6) - C(5,5))
      StiffOut(6,5) = StiffOut(5,6)
      StiffOut(6,6) = n**two*C(5,5) + m**two*C(6,6)


    Else
      Call log%error("Plane stress stiffness transformation is not implemented")

    End If

    Return
  End Function StiffRot


  Pure Function convertToCauchy(stress, F) result(Cauchy)
    ! Converts the stress to Cauchy stress

    Use matrixAlgUtil_Mod

    ! Input
    Double Precision, intent(IN) :: stress(3,3)
    Double Precision, intent(IN) :: F(3,3)            ! Deformation gradient

    ! Output
    Double Precision :: Cauchy(3,3)
    ! -------------------------------------------------------------------- !

    Cauchy = MATMUL(F, MATMUL(stress, TRANSPOSE(F)))/MDet(F)

    Return
  End Function convertToCauchy


  Pure Function Hooke(C, strain, nshr) result(stress)
    ! Calculates the energy conjugate stress based on the input strain and stiffness tensor

    Use matrixAlgUtil_Mod

    ! Input
    Double Precision, intent(IN) :: C(3+nshr,3+nshr)         ! Stiffness
    Double Precision, intent(IN) :: strain(3,3)
    Integer, intent(IN) :: nshr

    ! Output
    Double Precision :: stress(3,3)

    ! Locals
    Double Precision :: strainVec(3+nshr), stressVec(3+nshr)
    Double Precision, parameter :: zero=0.d0, two=2.d0
    ! -------------------------------------------------------------------- !

    stress = zero

    ! Convert strain to vector format, engineering shear strain
    strainVec = Matrix2Vec(strain, nshr)
    Do I=4,3+nshr; strainVec(I) = two*strainVec(I); End Do

    ! Compute stress vector
    stressVec = MATMUL(C,strainVec)

    ! Convert back to matrix format
    stress = Vec2Matrix(stressVec)

    Return
  End Function Hooke

End Module stress_Mod
