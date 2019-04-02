
Subroutine gaussian_marsaglia(n, avg, sdv, dist, limits)
  ! Generates a gaussian distribution of n points

  ! Arguments
  Integer, intent(IN) :: n                                     ! Number of points
  Double Precision, intent(IN) :: avg                          ! mu
  Double Precision, intent(IN) :: sdv                          ! sigma
  Double Precision, intent(IN) :: limits(2)                    ! Assumes limit(1) is less than zero

  ! Output
  Double Precision, intent(OUT) :: dist(n)

  ! Locals
  Double Precision :: v1(n), v2(n), s(n), z(n)
  Double Precision :: a, b

  Double Precision, parameter :: zero=0.d0, half=0.5d0, one=1.d0, two=2.d0
  ! -------------------------------------------------------------------- !

  v1 = zero
  v2 = zero
  s = zero

  Do I = 1,n
    Do While (s(I) >= one .OR. s(I) == zero)
      Call RANDOM_NUMBER(a)
      v1(i) = two*a - one
      Call RANDOM_NUMBER(b)
      v2(i) = two*b - one
      s(i) = v1(i)**two + v2(i)**two
    End Do
    s(i) = SQRT((-two*LOG(s(i)))/s(i))
    z(i) = v1(i)*s(i)
    dist(i) = avg + z(i)*sdv
  End Do

  If (.NOT. (limits(1) == zero .AND. limits(2) == zero)) Then
    Do I = 1,n
      If (dist(I) .GT. limits(2)) Then
        dist(I) = dist(I) - two*(limits(2))
      Else If (dist(I) .LT. limits(1)) Then
        dist(I) = dist(I) - two*(limits(1))
      End If
    End Do
  End If

  Return
End Subroutine gaussian_marsaglia


Subroutine lognormal(n, avg, sdv, dist)
  ! Generates a gaussian distribution of n points

  ! Arguments
  Integer, intent(IN) :: n                                     ! Number of points
  Double Precision, intent(IN) :: avg                          ! mu
  Double Precision, intent(IN) :: sdv                          ! sigma

  ! Output
  Double Precision, intent(OUT) :: dist(n)

  ! Locals
  Double Precision :: normal(n)

  Double Precision, parameter :: zero=0.d0, half=0.5d0, one=1.d0, two=2.d0
  ! -------------------------------------------------------------------- !

  Call gaussian_marsaglia(n, LOG10(avg), sdv, normal, (/zero, zero/))
  Do I = 1,n
    dist(i) = EXP(avg+sdv*normal(i))
  End Do

  Return
End Subroutine lognormal


Subroutine vexternaldb(lOp, i_Array, niArray, r_Array, nrArray)

  Use parameters_Mod

  Implicit Double Precision (a-h, o-z)
  Integer, parameter :: j_sys_Dimension = 2, n_vec_Length = 136, maxblk = n_vec_Length

  ! Arguments
  Dimension i_Array(niArray), r_Array(nrArray)

  ! Contents of i_Array
  Integer, parameter :: i_int_nTotalNodes = 1, i_int_nTotalElements = 2, i_int_kStep = 3, i_int_kInc = 4, i_int_iStatus = 5, i_int_lWriteRestart = 6

  ! Possible values for the lOp argument
  Integer, parameter :: j_int_StartAnalysis = 0, j_int_StartStep = 1, j_int_SetupIncrement = 2, j_int_StartIncrement = 3, j_int_EndIncrement = 4, j_int_EndStep = 5, j_int_EndAnalysis = 6

  ! Possible values for i_Array(i_int_iStatus)
  Integer, parameter :: j_int_Continue = 0, j_int_TerminateStep = 1, j_int_TerminateAnalysis = 2

  ! Contents of r_Array
  Integer, parameter :: i_flt_TotalTime = 1, i_flt_StepTime = 2, i_flt_dTime = 3


  ! Locals
  Integer :: kStep                                     ! Current step number
  Integer :: kInc                                      ! Current increment number
  Integer :: analysis_status
  Double Precision :: randomNumbers1(10000), randomNumbers2(10000)
  Double Precision :: randomNumbers3(10000), randomNumbers4(10000)
  Double Precision :: limits(2)
  Integer :: n
  Integer :: values(8)
  Type(parameters) :: p

  Dimension INTV(1), REALV(1)    ! For abaqus warning messages
  Character(len=8) CHARV(1)      ! For Abaqus warning messages

  ! Common
  Common /analysis_termination/ analysis_status
  Common randomNumbers1, randomNumbers2, randomNumbers3, randomNumbers4
  ! -------------------------------------------------------------------- !

  kStep = i_Array(i_int_kStep)
  kInc  = i_Array(i_int_kInc)

  ! Start of the analysis
  If (lOp == j_int_StartAnalysis) Then

    ! Initialize common analysis_status variable
    analysis_status = 1

    ! Load the CompDam solution parameters
    p = loadParameters()

    ! Initialize random numbers for fiber misalignment here
    If (p%fkt_random_seed) Then
      Call DATE_AND_TIME(VALUES=values)
      Call RANDOM_SEED (PUT=(/values(8)/))
    End If
    Call RANDOM_NUMBER(randomNumbers1)
    Call RANDOM_NUMBER(randomNumbers2)

    ! Get a Lognormal distribution (polar)
    n = 10000
    Call lognormal(n, p%fkt_init_misalignment_polar_shape, p%fkt_init_misalignment_polar_scale, randomNumbers3)

    ! Get a Gaussian distribution (azi)
    limits = (/-180.d0, 180.d0/)
    Call gaussian_marsaglia(n, p%fkt_init_misalignment_azi_mu, p%fkt_init_misalignment_azi_sigma, randomNumbers4, limits)

  ! End of the increment
  Else If (lOp == j_int_EndIncrement) Then

    ! If analysis_status has been changed, cleanly terminate the analysis
    If (analysis_status == 0) i_Array(i_int_iStatus) = j_int_TerminateAnalysis

  End If

  Return
End Subroutine vexternaldb
