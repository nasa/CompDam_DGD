
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


Subroutine get_compdam_parameters(print_logs, p)
  ! Populates CompDam parameters with defaults, overwriting with any nondefault parameters
  ! provided in the input deck.

  Use parameters_Mod

  ! Arguments
  Logical, intent(IN) :: print_logs

  ! Output
  Type(parameters), intent(OUT) :: p

  ! Locals
  Character (len=80), allocatable :: cParams(:)
  Character (len=80) :: tablelabel
  Character (len=80) :: tcNames(10)
  Character (len=80) :: parameterName
  Double Precision, allocatable :: rParams(:)
  Double Precision :: parameterRealValue
  Integer, allocatable :: iParamsDataType(:), iParams(:)
  Integer :: parameterIntValue
  Integer :: JERROR, numParams, nRows, numRows
  Dimension jSize(4)

  Double Precision, parameter :: zero=0.d0, half=0.5d0, one=1.d0, two=2.d0
  ! -------------------------------------------------------------------- !

  Call queryTableCollectionSize(jSize, jError)
  If (JERROR /= 0) Then
    If (print_logs) Print *, "ERROR: Unable to query table collection size, using default parameters"
    Return
  End If
  If (jSize(1) == 0) Then
    If (print_logs) Print *, "INFO: Using CompDam default solution parameters"
    Return
  End If
  Call queryTableCollectionNames(tcNames, jSize(1))
  If (ANY(tcNames(1:jSize(1)) == "COMPDAM_PARAM")) Then
    If (print_logs) Print *, "INFO: COMPDAM_PARAM table collection found, reading parameters"
    Call setTableCollection("COMPDAM_PARAM", JERROR)
    If (JERROR /= 0) Then
      If (print_logs) Print *, "ERROR: setTableCollection() error, using default parameters"
      Return
    End If
    ! Process table CDP_FLOAT
    tablelabel = "CDP_FLOAT"
    Call queryParameterTable(tablelabel, numParams, numRows, JERROR)
    If (JERROR == 0 .AND. numParams == 2) Then
      Allocate(iParamsDataType(numParams))
      Allocate(iParams(numParams*numRows))
      Allocate(rParams(numParams*numRows))
      Allocate(cParams(numParams*numRows))
      Call getParameterTable(tablelabel, numParams, iParamsDataType, iParams, rParams, cParams, jError)
      Do I = 1, numRows
        parameterName = cParams(2*I-1)
        parameterRealValue = rParams(2*I)
        Select Case (trim(parameterName))

          Case ('TOL_DGD_F')
            Call verifyAndSaveProperty_dbl(trim(parameterName), parameterRealValue, p%tol_DGD_f_min, p%tol_DGD_f_max, p%tol_DGD_f, print_logs)

          Case ('DGDGC_MIN')
            Call verifyAndSaveProperty_dbl(trim(parameterName), parameterRealValue, p%dGdGc_min_min, p%dGdGc_min_max, p%dGdGc_min, print_logs)

          Case ('COMPLIMIT')
            Call verifyAndSaveProperty_dbl(trim(parameterName), parameterRealValue, p%compLimit_min, p%compLimit_max, p%compLimit, print_logs)

          Case ('PENSTIFFMULT')
            Call verifyAndSaveProperty_dbl(trim(parameterName), parameterRealValue, p%penStiffMult_min, p%penStiffMult_max, p%penStiffMult, print_logs)

          Case ('CUTBACK_AMOUNT')
            Call verifyAndSaveProperty_dbl(trim(parameterName), parameterRealValue, p%cutback_amount_min, p%cutback_amount_max, p%cutback_amount, print_logs)

          Case ('TOL_DIVERGENCE')
            Call verifyAndSaveProperty_dbl(trim(parameterName), parameterRealValue, p%tol_divergence_min, p%tol_divergence_max, p%tol_divergence, print_logs)

          Case ('GAMMA_MAX')
            Call verifyAndSaveProperty_dbl(trim(parameterName), parameterRealValue, p%gamma_max_min, p%gamma_max_max, p%gamma_max, print_logs)

          Case ('KB_DECOMPOSE_THRES')
            Call verifyAndSaveProperty_dbl(trim(parameterName), parameterRealValue, p%kb_decompose_thres_min, p%kb_decompose_thres_max, p%kb_decompose_thres, print_logs)

          Case ('FKT_FIBER_FAILURE_ANGLE')
            Call verifyAndSaveProperty_dbl(trim(parameterName), parameterRealValue, p%fkt_fiber_failure_angle_min, p%fkt_fiber_failure_angle_max, p%fkt_fiber_failure_angle, print_logs)

          Case ('SCHAEFER_NR_TOLERANCE')
            Call verifyAndSaveProperty_dbl(trim(parameterName), parameterRealValue, p%schaefer_nr_tolerance_min, p%schaefer_nr_tolerance_max, p%schaefer_nr_tolerance, print_logs)

          Case ('DEBUG_KILL_AT_TOTAL_TIME')
            Call verifyAndSaveProperty_dbl(trim(parameterName), parameterRealValue, p%debug_kill_at_total_time_min, p%debug_kill_at_total_time_max, p%debug_kill_at_total_time, print_logs)

          Case ('FKT_INIT_MISALIGNMENT_AZI_MU')
            Call verifyAndSaveProperty_dbl(trim(parameterName), parameterRealValue, p%fkt_init_misalignment_azi_mu_min, p%fkt_init_misalignment_azi_mu_max, p%fkt_init_misalignment_azi_mu, print_logs)

          Case ('FKT_INIT_MISALIGNMENT_AZI_SIGMA')
            Call verifyAndSaveProperty_dbl(trim(parameterName), parameterRealValue, p%fkt_init_misalignment_azi_sigma_min, p%fkt_init_misalignment_azi_sigma_max, p%fkt_init_misalignment_azi_sigma, print_logs)

          Case ('FKT_INIT_MISALIGNMENT_POLAR_SHAPE')
            Call verifyAndSaveProperty_dbl(trim(parameterName), parameterRealValue, p%fkt_init_misalignment_polar_shape_min, p%fkt_init_misalignment_polar_shape_max, p%fkt_init_misalignment_polar_shape, print_logs)

          Case ('FKT_INIT_MISALIGNMENT_POLAR_SCALE')
            Call verifyAndSaveProperty_dbl(trim(parameterName), parameterRealValue, p%fkt_init_misalignment_polar_scale_min, p%fkt_init_misalignment_polar_scale_max, p%fkt_init_misalignment_polar_scale, print_logs)

          Case ('FATIGUE_R_RATIO')
            Call verifyAndSaveProperty_dbl(trim(parameterName), parameterRealValue, p%fatigue_R_ratio_min, p%fatigue_R_ratio_max, p%fatigue_R_ratio, print_logs)

          Case ('CYCLES_PER_INCREMENT_INIT')
            Call verifyAndSaveProperty_dbl(trim(parameterName), parameterRealValue, p%cycles_per_increment_init_min, p%cycles_per_increment_init_max, p%cycles_per_increment_init, print_logs)

          Case ('CYCLES_PER_INCREMENT_MOD')
            Call verifyAndSaveProperty_dbl(trim(parameterName), parameterRealValue, p%cycles_per_increment_mod_min, p%cycles_per_increment_mod_max, p%cycles_per_increment_mod, print_logs)

          Case ('CYCLES_PER_INCREMENT_MAX')
            Call verifyAndSaveProperty_dbl(trim(parameterName), parameterRealValue, p%cycles_per_increment_max_min, p%cycles_per_increment_max_max, p%cycles_per_increment_max, print_logs)

          Case ('CYCLES_PER_INCREMENT_MIN')
            Call verifyAndSaveProperty_dbl(trim(parameterName), parameterRealValue, p%cycles_per_increment_min_min, p%cycles_per_increment_min_max, p%cycles_per_increment_min, print_logs)

          Case ('FATIGUE_DAMAGE_MIN_THRESHOLD')
            Call verifyAndSaveProperty_dbl(trim(parameterName), parameterRealValue, p%fatigue_damage_min_threshold_min, p%fatigue_damage_min_threshold_max, p%fatigue_damage_min_threshold, print_logs)

          Case ('FATIGUE_DAMAGE_MAX_THRESHOLD')
            Call verifyAndSaveProperty_dbl(trim(parameterName), parameterRealValue, p%fatigue_damage_max_threshold_min, p%fatigue_damage_max_threshold_max, p%fatigue_damage_max_threshold, print_logs)

          Case Default
            If (print_logs) Print *, "WARN: Parameter not recognized: " // trim(parameterName)
        End Select
      End Do
      Deallocate(iParamsDataType)
      Deallocate(iParams)
      Deallocate(rParams)
      Deallocate(cParams)
    Else
      If (print_logs) Print *, "INFO: Missing CDP_FLOAT, using default parameters for floats"
    End If

    ! Process table CDP_INT
    tablelabel = "CDP_INT"
    Call queryParameterTable(tablelabel, numParams, numRows, JERROR)
    If (JERROR == 0 .AND. numParams == 2) Then
      Allocate(iParamsDataType(numParams))
      Allocate(iParams(numParams*numRows))
      Allocate(rParams(numParams*numRows))
      Allocate(cParams(numParams*numRows))
      Call getParameterTable(tablelabel, numParams, iParamsDataType, iParams, rParams, cParams, jError)
      Do I = 1, numRows
        parameterName = cParams(2*I-1)
        parameterIntValue = iParams(2*I)
        Select Case (trim(parameterName))

          Case ('LOGLEVEL')
            Call verifyAndSaveProperty_int(trim(parameterName), parameterIntValue, p%logLevel_min, p%logLevel_max, p%logLevel, print_logs)

          Case ('LOGFORMAT')
            Call verifyAndSaveProperty_int(trim(parameterName), parameterIntValue, p%logFormat_min, p%logFormat_max, p%logFormat, print_logs)

          Case ('CUTBACKS_MAX')
            Call verifyAndSaveProperty_int(trim(parameterName), parameterIntValue, p%cutbacks_max_min, p%cutbacks_max_max, p%cutbacks_max, print_logs)

          Case ('MD_MAX')
            Call verifyAndSaveProperty_int(trim(parameterName), parameterIntValue, p%MD_max_min, p%MD_max_max, p%MD_max, print_logs)

          Case ('EQ_MAX')
            Call verifyAndSaveProperty_int(trim(parameterName), parameterIntValue, p%EQ_max_min, p%EQ_max_max, p%EQ_max, print_logs)

          Case ('ALPHA_SEARCH')
            Call verifyAndSaveProperty_int(trim(parameterName), parameterIntValue, p%alpha_search_min, p%alpha_search_max, p%alpha_search, print_logs)

          Case ('ALPHA_INC')
            Call verifyAndSaveProperty_int(trim(parameterName), parameterIntValue, p%alpha_inc_min, p%alpha_inc_max, p%alpha_inc, print_logs)

          Case ('ALPHA_MAX')
            Call verifyAndSaveProperty_int(trim(parameterName), parameterIntValue, p%alpha_max_min, p%alpha_max_max, p%alpha_max, print_logs)

          Case ('SCHAEFER_NR_COUNTER_LIMIT')
            Call verifyAndSaveProperty_int(trim(parameterName), parameterIntValue, p%schaefer_nr_counter_limit_min, p%schaefer_nr_counter_limit_max, p%schaefer_nr_counter_limit, print_logs)

          Case ('TERMINATE_ON_NO_CONVERGENCE')
            Call verifyAndSaveProperty_int(trim(parameterName), parameterIntValue, p%terminate_on_no_convergence_min, p%terminate_on_no_convergence_max, p%terminate_on_no_convergence, print_logs)

          Case ('FKT_RANDOM_SEED')
            Call verifyAndSaveProperty_int(trim(parameterName), parameterIntValue, p%fkt_random_seed_min, p%fkt_random_seed_max, p%fkt_random_seed, print_logs)

          Case ('FATIGUE_STEP')
            Call verifyAndSaveProperty_int(trim(parameterName), parameterIntValue, p%fatigue_step_min, p%fatigue_step_max, p%fatigue_step, print_logs)

          Case ('CHECK_FOR_SNAP_BACK')
            Call verifyAndSaveProperty_int(trim(parameterName), parameterIntValue, p%check_for_snap_back_min, p%check_for_snap_back_max, p%check_for_snap_back, print_logs)

          Case ('SET_STATUS_0_ON_D2')
            Call verifyAndSaveProperty_int(trim(parameterName), parameterIntValue, p%set_status_0_on_d2_min, p%set_status_0_on_d2_max, p%set_status_0_on_d2, print_logs)

          Case ('SET_STATUS_0_ON_D1T')
            Call verifyAndSaveProperty_int(trim(parameterName), parameterIntValue, p%set_status_0_on_d1T_min, p%set_status_0_on_d1T_max, p%set_status_0_on_d1T, print_logs)

          Case ('SET_STATUS_0_ON_D1C')
            Call verifyAndSaveProperty_int(trim(parameterName), parameterIntValue, p%set_status_0_on_d1C_min, p%set_status_0_on_d1C_max, p%set_status_0_on_d1C, print_logs)

          Case Default
            If (print_logs) Print *, "WARN: Parameter not recognized: " // trim(parameterName)
        End Select
      End Do
    Else
      If (print_logs) Print *, "INFO: Missing CDP_INT, using default parameters for integers"
    End If
  End If
  Return
End Subroutine get_compdam_parameters

Subroutine vexternaldb(lOp, i_Array, niArray, r_Array, nrArray)

  Use parameters_Mod
  Use version_Mod

  Implicit Double Precision (a-h, o-z)
  Integer, parameter :: j_sys_Dimension = 2, maxblk = 512

  Include 'mpif.h'

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
  Integer :: increment                                      ! Current increment number
  Integer :: analysis_status
  Double Precision :: randomNumbers1(10000), randomNumbers2(10000)
  Double Precision :: randomNumbers3(10000), randomNumbers4(10000)
  Double Precision :: limits(2)
  Integer :: n
  Integer :: values(8)
  Double Precision, Parameter :: zero=0.d0, one=1.d0

  Dimension INTV(1), REALV(1)    ! For Abaqus warning messages
  Character(len=8) CHARV(1)      ! For Abaqus warning messages

  ! Fatigue
  Integer :: fatigue_parameters(2), increment_old
  Double Precision :: cycles_per_increment(1), cycles_per_increment_old
  Double Precision :: cycles, cycles_old
  ! fatigue_parameters(1) is an integer flag for whether the current analysis step is a fatigue step (1 = .TRUE.)
  ! fatigue_parameters(2) describes the current rate of fatigue damage propagation (0 = slow, 1 = in-range, 2 = fast)
  Pointer(ptr_fatigue_int, fatigue_parameters)  ! pointer for fatigue parameters integer array
  Pointer(ptr_fatigue_dbl, cycles_per_increment)  ! pointer for fatigue parameters real array
  Logical :: update_inc2cycles_log  ! flag to update the inc2cycles log

  ! MPI
  Integer :: ABA_COMM_WORLD  ! communicator that Abaqus defines for its worker processes
  Integer :: process_rank  ! Process number or rank
  Integer :: ierr  ! Error status

  ! Job names, file names, and unit numbers
  Character(len=256) :: outputDir, jobName, inc2cycles_filename
  Integer :: lenOutputDir, lenJobName
  Integer, parameter :: inc2cycles_file_unit = 15

  ! Parameters storage
  Type(parameters) :: default_params
  Type(parameters) :: params               ! Storage for parameters loaded from inp deck
  Type(parameters) :: p                    ! Parameters for execution
  Pointer(ptr_params, params)

  ! Common
  Common /analysis_termination/ analysis_status
  Common randomNumbers1, randomNumbers2, randomNumbers3, randomNumbers4

  Interface UserUtilities

    Function GETCOMMUNICATOR() result(communicator)
      ! Utility function GETCOMMUNICATOR can be called from any Abaqus user subroutine. Returns a communicator (type INTEGER) that
      ! Abaqus defines for its worker processes, similar to MPI_COMM_WORLD. The communicator thus obtained can be used for
      ! subsequent MPI communication routines. In a nonparallel run, communicators do not exist and GET_COMMUNICATOR() returns 0.
      ! See 'Obtaining parallel processes information' in the Abaqus documentation for more details.
      Integer(kind=4) :: communicator
    End Function GETCOMMUNICATOR

  End Interface UserUtilities

  Interface
    ! See 'Allocatable arrays' in the Abaqus documentation for more details on the following functions.
    Function SMAIntArrayCreate(ID, SIZE, INITVAL)
      ! Create or resize a global integer array.
      Integer(kind=8) :: SMAIntArrayCreate  ! Returns an address that can be associated with a Fortran pointer
      Integer(kind=4) :: ID       ! Arbitrary integer chosen by the user, used later to access this array
      Integer(kind=4) :: SIZE     ! max value is INT_MAX ( 2,147,483,647 )
      Integer(kind=4) :: INITVAL  ! initialization value for each element in the array
    End Function SMAIntArrayCreate

    Function SMAIntArrayAccess(ID)
      ! Access an existing global integer array.
      Integer(kind=8) :: SMAIntArrayAccess  ! Returns an address that can be associated with a Fortran pointer
      Integer(kind=4) :: ID  ! Array ID
    End Function SMAIntArrayAccess

    Function SMARealArrayCreateDP(ID, SIZE, INITVAL)
      ! Create or resize a global real array.
      Integer(kind=8) :: SMARealArrayCreateDP  ! Returns a pointer to the newly allocated array
      Integer(kind=4), intent(IN) :: ID        ! Arbitrary integer chosen by the user, used later to locate this array
      Integer(kind=4), intent(IN) :: SIZE      ! max value is INT_MAX ( 2,147,483,647 )
      Double Precision, intent(IN) :: INITVAL  ! (optional) initial value for each element of the array
    End Function SMARealArrayCreateDP

    Function SMARealArrayAccess(ID)
      ! Access an existing global real array.
      Integer(kind=8) :: SMARealArrayAccess  ! Returns an address that can be associated with a Fortran pointer
      Integer(kind=4) :: ID   ! Array ID
    End Function SMARealArrayAccess

    FUNCTION SMAStructArrayAccess(ID)
          INTEGER(KIND=8) :: SMAStructArrayAccess  ! -- Returns an address that can be associated with a Fortran pointer
          INTEGER(KIND=4) :: ID                    ! Array ID
    END FUNCTION SMAStructArrayAccess

  End Interface

  INTERFACE SMAStructArrayCreate

    ! -- Creates an array with a given ID, length = NUM_ITEMS; no initialization
    integer(kind=8) FUNCTION SMAStructArrayCreateNoInit(ARRAY_ID, NUM_ITEMS, ITEM_SIZE)
        INTEGER(KIND=4),INTENT(IN) :: ARRAY_ID   ! arbitrary ID chosen by the user
        INTEGER(KIND=4),INTENT(IN) :: NUM_ITEMS  ! max value is INT_MAX ( 2,147,483,647 )
        INTEGER(KIND=8),INTENT(IN) :: ITEM_SIZE  ! size of one struct in bytes as returned by SIZEOF()
    END FUNCTION SMAStructArrayCreateNoInit

    ! -- Creates an array with a given ID and SIZE; each slot is initialized to INITVAL
    integer(kind=8) FUNCTION SMAStructArrayCreateInit(ARRAY_ID,NUM_ITEMS,ITEM_SIZE,INITVAL)
        INTEGER(KIND=4),INTENT(IN) :: ARRAY_ID   ! arbitrary ID chosen by the user
        INTEGER(KIND=4),INTENT(IN) :: NUM_ITEMS  ! max value is INT_MAX ( 2,147,483,647 )
        INTEGER(KIND=8),INTENT(IN) :: ITEM_SIZE  ! size of one struct in bytes as returned by SIZEOF()
        CLASS(*),INTENT(IN)        :: INITVAL    ! a struct used as initializer for each slot of the array
    END FUNCTION SMAStructArrayCreateInit

  END INTERFACE SMAStructArrayCreate
  ! -------------------------------------------------------------------- !

  kStep = i_Array(i_int_kStep)
  increment = i_Array(i_int_kInc)

  ! Start of the analysis
  If (lOp == j_int_StartAnalysis) Then

    ptr_params = SMAStructArrayCreate(1000, 1, sizeof(default_params), default_params)
    p = default_params

    ! Initialize common analysis_status variable
    analysis_status = 1

    Call VGETRANK(process_rank)
    If (process_rank == 0) Then
      print *, '================= CompDam ================='
      print *, 'Date: ' // trim(timestamp)
      print *, 'Version: ' // trim(hash)
      print *, '==========================================='

      ! Load the CompDam solution parameters, with logging
      Call get_compdam_parameters(.TRUE., p)
      params = p
    Else
      ! Load the CompDam solution parameters, without logging
      Call get_compdam_parameters(.FALSE., p)
      params = p
    End If

    ! Global integer array for fatigue parameters
    ptr_fatigue_int = SMAIntArrayCreate(1, 2, 0)  ! create the pointer for fatigue parameters array
    ptr_fatigue_dbl = SMARealArrayCreateDP(2, 1, p%cycles_per_increment_init)  ! create the pointer for fatigue parameters array

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

  ! Start of a step
  Else If (lOp == j_int_StartStep) Then

    ! Load the CompDam solution parameters
    ptr_params = SMAStructArrayAccess(1000)
    p = params
    ptr_fatigue_int = SMAIntArrayAccess(1)  ! Access the pointer for fatigue parameters

    ! Check to see if this is a fatigue step
    FatigueStepSetup: If (kStep == p%fatigue_step) Then  ! This is a fatigue step

      fatigue_parameters(1) = 1

      ! Create a log file for the fatigue cycles per solution increment
      Call VGETRANK(process_rank)
      If (process_rank == 0) Then
        ptr_fatigue_dbl = SMARealArrayAccess(2)
        Call VGETJOBNAME(jobName, lenJobName)
        Call VGETOUTDIR(outputDir, lenOutputDir)
        inc2cycles_filename = trim(outputDir) // '/' // trim(jobName) // '_inc2cycles.log'
        open(inc2cycles_file_unit, file=inc2cycles_filename, status='replace', action='write')
        write(inc2cycles_file_unit,"(A9, A22, A22)") 'increment', 'cycles_per_increment', 'cycles'
        write(inc2cycles_file_unit,"(I9, ES22.12E2, ES22.12E2)") 0, cycles_per_increment(1), zero
        close(inc2cycles_file_unit)
      End If

    Else FatigueStepSetup  ! This is not a fatigue step

      fatigue_parameters(1) = 0

    End If FatigueStepSetup

  ! End of an increment
  Else If (lOp == j_int_EndIncrement) Then

    ! If analysis_status has been changed, cleanly terminate the analysis
    If (analysis_status == 0) i_Array(i_int_iStatus) = j_int_TerminateAnalysis

    ptr_fatigue_int = SMAIntArrayAccess(1)  ! Access the pointer for fatigue_parameters

    ! Check the rate of fatigue damage progression and adjust solution parameters if necessary
    FatigueIncrement: If (fatigue_parameters(1) == 1) Then

      ptr_fatigue_dbl = SMARealArrayAccess(2)

      ! Synchronize fatigue damage progression parameter across processes
      ABA_COMM_WORLD = GETCOMMUNICATOR()
      If (ABA_COMM_WORLD /= 0) Then
        Call MPI_Allreduce(MPI_IN_PLACE, fatigue_parameters(2), 1, MPI_INTEGER, MPI_MAX, ABA_COMM_WORLD, ierr)
      End If

      FatigueParameterUpdate: If (fatigue_parameters(2) == 1) Then
        Continue
      Else FatigueParameterUpdate
        ptr_params = SMAStructArrayAccess(1000)
        p = params
        update_inc2cycles_log = .FALSE.
        ! If fatigue damage is progressing too slowly...
        If (fatigue_parameters(2) == 0 .AND. cycles_per_increment(1) < p%cycles_per_increment_max) Then
            ! ...increase the number of cycles represented by a solution increment.
          cycles_per_increment(1) = MIN(cycles_per_increment(1)*(one + p%cycles_per_increment_mod), p%cycles_per_increment_max)
          update_inc2cycles_log = .TRUE.
        ! If fatigue damage is progressing too quickly...
        Else If (fatigue_parameters(2) == 2 .AND. cycles_per_increment(1) > p%cycles_per_increment_min) Then
          ! ...decrease the number of cycles represented by a solution increment.
          cycles_per_increment(1) = MAX(cycles_per_increment(1)/(one + p%cycles_per_increment_mod), p%cycles_per_increment_min)
          update_inc2cycles_log = .TRUE.
        End If

        FatigueLog: If (update_inc2cycles_log) Then
          ! Update the inc2cycles log file
          Call VGETRANK(process_rank)
          If (process_rank == 0) Then
            Call VGETJOBNAME(jobName, lenJobName)
            Call VGETOUTDIR(outputDir, lenOutputDir)
            inc2cycles_filename = trim(outputDir) // '/' // trim(jobName) // '_inc2cycles.log'
            open(inc2cycles_file_unit, file=inc2cycles_filename, status='old', action='READWRITE', position='append')
            backspace(inc2cycles_file_unit)
            read(inc2cycles_file_unit,"(I9, ES22.12E2, ES22.12E2)") increment_old, cycles_per_increment_old, cycles_old
            cycles = cycles_old + (increment - increment_old)*cycles_per_increment_old
            write(inc2cycles_file_unit,"(I9, ES22.12E2, ES22.12E2)") increment, cycles_per_increment(1), cycles
            close(inc2cycles_file_unit)
          End If
        End If FatigueLog

      End If FatigueParameterUpdate

      ! Reset fatigue damage progression variable, i.e., fatigue_parameters(2)
      fatigue_parameters(2) = 0

    End If FatigueIncrement

  End If

  Return
End Subroutine vexternaldb
