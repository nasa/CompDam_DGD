Subroutine vexternaldb(lOp, i_Array, niArray, r_Array, nrArray)

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
  Integer, parameter :: randomNumberCount = 10000
  Double Precision :: randomNumbers(randomNumberCount)

  ! Common
  Common /analysis_termination/ analysis_status
  Common randomNumbers
  ! -------------------------------------------------------------------- !

  kStep = i_Array(i_int_kStep)
  kInc  = i_Array(i_int_kInc)

  ! Start of the analysis
  If (lOp == j_int_StartAnalysis) Then

    ! Initialize common analysis_status variable
    analysis_status = 1

    ! Initialize random numbers for fiber misalignment here
    Call RANDOM_NUMBER(randomNumbers)

  ! End of the increment
  Else If (lOp == j_int_EndIncrement) Then

    ! If analysis_status has been changed, cleanly terminate the analysis
    If (analysis_status == 0) i_Array(i_int_iStatus) = j_int_TerminateAnalysis

  End If

  Return
End Subroutine vexternaldb
