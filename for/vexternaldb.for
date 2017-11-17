Subroutine vexternaldb(lOp, i_Array, niArray, r_Array, nrArray)

  include 'vaba_param.inc'

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
  Integer, parameter :: randomNumberCount = 1000
  Double Precision :: randomNumbers(randomNumberCount)

  ! Common
  Common randomNumbers
  ! -------------------------------------------------------------------- !


  kStep = i_Array(i_int_kStep)
  kInc  = i_Array(i_int_kInc)


  ! Note that you  can use the MPI communication between parallel Abaqus processes to gather
  ! and scatter the data.

  ! Start of the analysis
  If (lOp .EQ. j_int_StartAnalysis) Then

    ! User coding to set up the environment, open files, launch/connect to the external programs, etc.

    ! Initialize random numbers for fiber misalignment here
    Call RANDOM_NUMBER(randomNumbers)

  End If

  Return
End Subroutine vexternaldb