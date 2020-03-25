#define type(x) TYPE(x), target

Module parameters_Mod
  ! Module for loading and validating DGD solution parameters

  Private

  ! Stores the set of parameters
  Type parameters

    ! DGDEvolve solution parameters
    Integer :: logLevel                                           ! Logger level
    Integer :: logFormat                                          ! 1=default (to log file), 2=csv
    Integer :: cutbacks_max                                       ! maximum number of cut-backs
    Integer :: MD_max                                             ! maximum number of damage increments per solution increment
    Integer :: EQ_max                                             ! maximum number of equilibrium iterations
    Integer :: alpha_inc                                          ! increment, in degrees, for which the matrix failure criterion is evaluated
    Integer :: schaefer_nr_counter_limit                          ! maximum number of Newton-Raphson loops before divergence is assumed
    Double Precision :: tol_DGD_f                                 ! tol_DGD_f = tol_DGD/YT
    Double Precision :: dGdGc_min                                 ! minimum amount of damage dissipated per MD increment
    Double Precision :: compLimit                                 ! minimum accepted det(F_bulk)
    Double Precision :: penStiffMult                              ! penalty stiffness multiplier
    Double Precision :: cutback_amount                            ! artificially slow rate of change of F_bulk
    Double Precision :: tol_divergence                            ! Tolerance for divergence of internal Newton Raphson loop
    Double Precision :: gamma_max                                 ! Maximum shear strain; when this value is exceeded, the element is deleted
    Double Precision :: kb_decompose_thres                        ! Ratio of kink band size to element length at which to decompose the element
    Double Precision :: fkt_fiber_failure_angle                   ! Angle at which fiber failure occurs; no further plastic shear deformation allowed
    Double Precision :: schaefer_nr_tolerance                     ! Tolerance value used to determine if convergence has occurred in newton raphson loop
    Logical :: terminate_on_no_convergence                        ! Set to false to trigger an error, delete the element, and allow the analysis to continue when no converged solution can be found
    Double Precision :: debug_kill_at_total_time                  ! Time at which to kill analysis (for debugging purposes)
    Logical :: fkt_random_seed                                    ! Set to true to use the time as the seed for initial fiber misalignment to get different results for each realization
    Double Precision :: fkt_init_misalignment_azi_mu              ! Initial fiber misalignment azimuthal average [degrees]
    Double Precision :: fkt_init_misalignment_azi_sigma           ! Initial fiber misalignment azimuthal variance [degrees]
    Double Precision :: fkt_init_misalignment_polar_shape         ! Initial fiber misalignment polar shape parameter [degrees]
    Double Precision :: fkt_init_misalignment_polar_scale         ! Initial fiber misalignment polar standard deviation [degrees]
    Double Precision :: fatigue_R_ratio                           ! R ratio for cohesive fatigue model, sigma_min / sigma_max
    Double Precision :: cycles_per_increment_init                 ! Fatigue cycles per solution increment
    Double Precision :: cycles_per_increment_max                  ! Maximum fatigue cycles per solution increment
    Double Precision :: cycles_per_increment_min                  ! Minimum fatigue cycles per solution increment
    Double Precision :: cycles_per_increment_mod                  ! Percent change in cycles_per_increment when out of range
    Double Precision :: fatigue_damage_min_threshold              ! Minimum required incremental fatigue damage to count as progression
    Double Precision :: fatigue_damage_max_threshold              ! Maximum allowable incremental fatigue damage per solution increment
    Integer :: fatigue_step                                       ! Step number for fatigue analysis

    ! min and max values for acceptable range
    Integer, private :: logLevel_min, logLevel_max
    Integer, private :: logFormat_min, logFormat_max
    Integer, private :: cutbacks_max_min, cutbacks_max_max
    Integer, private :: MD_max_min, MD_max_max
    Integer, private :: EQ_max_min, EQ_max_max
    Integer, private :: alpha_inc_min, alpha_inc_max
    Integer, private :: schaefer_nr_counter_limit_min, schaefer_nr_counter_limit_max
    Double Precision, private :: tol_DGD_f_min, tol_DGD_f_max
    Double Precision, private :: dGdGc_min_min, dGdGc_min_max
    Double Precision, private :: compLimit_min, compLimit_max
    Double Precision, private :: penStiffMult_min, penStiffMult_max
    Double Precision, private :: cutback_amount_min, cutback_amount_max
    Double Precision, private :: tol_divergence_min, tol_divergence_max
    Double Precision, private :: gamma_max_min, gamma_max_max
    Double Precision, private :: kb_decompose_thres_min, kb_decompose_thres_max
    Double Precision, private :: fkt_fiber_failure_angle_min, fkt_fiber_failure_angle_max
    Double Precision, private :: schaefer_nr_tolerance_min, schaefer_nr_tolerance_max
    Double Precision, private :: debug_kill_at_total_time_min, debug_kill_at_total_time_max
    Double Precision, private :: fkt_init_misalignment_azi_mu_min, fkt_init_misalignment_azi_mu_max
    Double Precision, private :: fkt_init_misalignment_azi_sigma_min, fkt_init_misalignment_azi_sigma_max
    Double Precision, private :: fkt_init_misalignment_polar_shape_min, fkt_init_misalignment_polar_shape_max
    Double Precision, private :: fkt_init_misalignment_polar_scale_min, fkt_init_misalignment_polar_scale_max
    Double Precision, private :: fatigue_R_ratio_min, fatigue_R_ratio_max
    Double Precision, private :: cycles_per_increment_init_min, cycles_per_increment_init_max
    Double Precision, private :: cycles_per_increment_max_min, cycles_per_increment_max_max
    Double Precision, private :: cycles_per_increment_min_min, cycles_per_increment_min_max
    Double Precision, private :: cycles_per_increment_mod_min, cycles_per_increment_mod_max
    Double Precision, private :: fatigue_damage_min_threshold_min, fatigue_damage_min_threshold_max
    Double Precision, private :: fatigue_damage_max_threshold_min, fatigue_damage_max_threshold_max
    Integer, private :: fatigue_step_min, fatigue_step_max

  End Type parameters

  ! Public interface
  Public :: parameters
  Public :: loadParameters
  Public :: writeParametersToFile

  ! Reference to object for singleton
  type(parameters), Save :: p

Contains

  Type(parameters) Function loadParameters()
    ! Loads parameters into module variable p

    Use forlog_Mod

    ! Locals
    Character(len=256) :: outputDir, jobName, fileName
    Integer :: lenOutputDir, lenJobName
    Logical :: fileExists
    ! -------------------------------------------------------------------- !

    ! Initializations
    Call initializeParameters()

#ifndef PYEXT
    ! Get the output directory (location to search for .parameters files)
    Call VGETOUTDIR(outputDir, lenOutputDir)

    ! Get the jobname
    Call VGETJOBNAME(jobName, lenJobName)

    ! Look to see if a parameters file exists
    ! First look for: jobName.parameters
    fileName = trim(outputDir) // '/' // trim(jobName) // '.parameters'
    Inquire(FILE=fileName, EXIST=fileExists)

    If (.NOT. fileExists) Then
      ! Try looking for: CompDam.parameters
      fileName = trim(outputDir) // '/CompDam.parameters'
      Inquire(FILE=fileName, EXIST=fileExists)
    End If

#else
    fileName = 'CompDam.parameters'
    Inquire(FILE=fileName, EXIST=fileExists)
#endif

    ! If the file is present, load parameters from file
    If (fileExists) Then
      Call log%info("loadParameters: Using parameters file " // trim(fileName))
      Call loadParametersFromFile(trim(fileName))
    Else
      Call log%info("loadParameters: Parameters file not found. Using default values.")
    End If

    ! Return a reference to the parameters object
    loadParameters = p

    Return
  End Function loadParameters


  Subroutine loadParametersFromFile(fileName)
    ! Populates p with the solution parameters in the specified file

    Use forlog_Mod

    ! Arguments
    Character(len=*), intent(IN) :: fileName

    ! Locals
    Integer, Parameter :: unit=109
    Integer :: iostat
    Character(len=256) :: line, key, value, tmp
    Character(len=30) :: featureFlags
    Integer :: commentTokenPos, equalTokenPos
    Double Precision, Parameter :: zero=0.d0, one=1.d0, two=2.d0
    ! -------------------------------------------------------------------- !

    ! Try to load CompDam parameters from a '.parameters' file
    Open(UNIT=unit, FILE=fileName, STATUS='old', ACTION='read', position='rewind', IOSTAT=iostat)
    If (iostat /= 0) Call log%error("loadParameters: Unable to access the .parameters file")
    ReadLines: Do

      ! Read the next line in the file
      Read(unit,'(A255)',IOSTAT=iostat) line

      If (iostat > 0) Then
        Call log%error("loadParameters: Unknown error reading file")
        Exit

      Else If (iostat < 0) Then
        Exit

      Else  ! Parse the line in the file
        ! Skip blank lines
        If (Len(Trim(line)) == 0) Cycle ReadLines

        ! Ignore comment lines (token: //)
        commentTokenPos = Index(line, '//')
        If (commentTokenPos == 1) Then
          Cycle ReadLines

        Else If (commentTokenPos > 1) Then
          ! This line has a trailing comment, retain everything before the comment
          line = line(1:commentTokenPos-1)
        End If

        ! Split the line into a key and value
        equalTokenPos = Index(line, '=')
        If (equalTokenPos == 0) Call log%error("loadParameters: Expected [name] = [value] format not found. Line image: " // trim(line))

        ! Parse key and value from the line
        key = line(1:equalTokenPos-1)
        value = line(equalTokenPos+1:)

        ! Check if the key is a known parameter. If it's a known parameter, check validity, and store it
        Select Case (trim(key))

          Case ('logLevel')
            Call verifyAndSaveProperty_int(trim(key), value, p%logLevel_min, p%logLevel_max, p%logLevel)

          Case ('logFormat')
            Call verifyAndSaveProperty_int(trim(key), value, p%logFormat_min, p%logFormat_max, p%logFormat)

          Case ('cutbacks_max')
            Call verifyAndSaveProperty_int(trim(key), value, p%cutbacks_max_min, p%cutbacks_max_max, p%cutbacks_max)

          Case ('MD_max')
            Call verifyAndSaveProperty_int(trim(key), value, p%MD_max_min, p%MD_max_max, p%MD_max)

          Case ('EQ_max')
            Call verifyAndSaveProperty_int(trim(key), value, p%EQ_max_min, p%EQ_max_max, p%EQ_max)

          Case ('alpha_inc')
            Call verifyAndSaveProperty_int(trim(key), value, p%alpha_inc_min, p%alpha_inc_max, p%alpha_inc)

          Case ('tol_DGD_f')
            Call verifyAndSaveProperty_dbl(trim(key), value, p%tol_DGD_f_min, p%tol_DGD_f_max, p%tol_DGD_f)

          Case ('dGdGc_min')
            Call verifyAndSaveProperty_dbl(trim(key), value, p%dGdGc_min_min, p%dGdGc_min_max, p%dGdGc_min)

          Case ('compLimit')
            Call verifyAndSaveProperty_dbl(trim(key), value, p%compLimit_min, p%compLimit_max, p%compLimit)

          Case ('penStiffMult')
            Call verifyAndSaveProperty_dbl(trim(key), value, p%penStiffMult_min, p%penStiffMult_max, p%penStiffMult)

          Case ('cutback_amount')
            Call verifyAndSaveProperty_dbl(trim(key), value, p%cutback_amount_min, p%cutback_amount_max, p%cutback_amount)

          Case ('tol_divergence')
            Call verifyAndSaveProperty_dbl(trim(key), value, p%tol_divergence_min, p%tol_divergence_max, p%tol_divergence)

          Case ('gamma_max')
            Call verifyAndSaveProperty_dbl(trim(key), value, p%gamma_max_min, p%gamma_max_max, p%gamma_max)

          Case ('kb_decompose_thres')
            Call verifyAndSaveProperty_dbl(trim(key), value, p%kb_decompose_thres_min, p%kb_decompose_thres_max, p%kb_decompose_thres)

          Case ('fkt_fiber_failure_angle')
            Call verifyAndSaveProperty_dbl(trim(key), value, p%fkt_fiber_failure_angle_min, p%fkt_fiber_failure_angle_max, p%fkt_fiber_failure_angle)

          Case ('schaefer_nr_tolerance')
            Call verifyAndSaveProperty_dbl(trim(key), value, p%schaefer_nr_tolerance_min, p%schaefer_nr_tolerance_max, p%schaefer_nr_tolerance)

          Case ('schaefer_nr_counter_limit')
            Call verifyAndSaveProperty_int(trim(key), value, p%schaefer_nr_counter_limit_min, p%schaefer_nr_counter_limit_max, p%schaefer_nr_counter_limit)

          Case ('terminate_on_no_convergence')
            Call verifyAndSaveProperty_logical(trim(key), adjustl(value), p%terminate_on_no_convergence)

          Case ('debug_kill_at_total_time')
            Call verifyAndSaveProperty_dbl(trim(key), value, p%debug_kill_at_total_time_min, p%debug_kill_at_total_time_max, p%debug_kill_at_total_time)

          Case ('fkt_random_seed')
            Call verifyAndSaveProperty_logical(trim(key), adjustl(value), p%fkt_random_seed)

          Case ('fkt_init_misalignment_azi_mu')
            Call verifyAndSaveProperty_dbl(trim(key), value, p%fkt_init_misalignment_azi_mu_min, p%fkt_init_misalignment_azi_mu_max, p%fkt_init_misalignment_azi_mu)

          Case ('fkt_init_misalignment_azi_sigma')
            Call verifyAndSaveProperty_dbl(trim(key), value, p%fkt_init_misalignment_azi_sigma_min, p%fkt_init_misalignment_azi_sigma_max, p%fkt_init_misalignment_azi_sigma)

          Case ('fkt_init_misalignment_polar_shape')
            Call verifyAndSaveProperty_dbl(trim(key), value, p%fkt_init_misalignment_polar_shape_min, p%fkt_init_misalignment_polar_shape_max, p%fkt_init_misalignment_polar_shape)

          Case ('fkt_init_misalignment_polar_scale')
            Call verifyAndSaveProperty_dbl(trim(key), value, p%fkt_init_misalignment_polar_scale_min, p%fkt_init_misalignment_polar_scale_max, p%fkt_init_misalignment_polar_scale)

          Case ('fatigue_R_ratio')
            Call verifyAndSaveProperty_dbl(trim(key), value, p%fatigue_R_ratio_min, p%fatigue_R_ratio_max, p%fatigue_R_ratio, .FALSE.)

          Case ('cycles_per_increment_init')
            Call verifyAndSaveProperty_dbl(trim(key), value, p%cycles_per_increment_init_min, p%cycles_per_increment_init_max, p%cycles_per_increment_init, .FALSE.)

          Case ('cycles_per_increment_mod')
            Call verifyAndSaveProperty_dbl(trim(key), value, p%cycles_per_increment_mod_min, p%cycles_per_increment_mod_max, p%cycles_per_increment_mod, .FALSE.)

          Case ('cycles_per_increment_max')
            Call verifyAndSaveProperty_dbl(trim(key), value, p%cycles_per_increment_max_min, p%cycles_per_increment_max_max, p%cycles_per_increment_max, .FALSE.)

          Case ('cycles_per_increment_min')
            Call verifyAndSaveProperty_dbl(trim(key), value, p%cycles_per_increment_min_min, p%cycles_per_increment_min_max, p%cycles_per_increment_min, .FALSE.)

          Case ('fatigue_damage_min_threshold')
            Call verifyAndSaveProperty_dbl(trim(key), value, p%fatigue_damage_min_threshold_min, p%fatigue_damage_min_threshold_max, p%fatigue_damage_min_threshold, .FALSE.)

          Case ('fatigue_damage_max_threshold')
            Call verifyAndSaveProperty_dbl(trim(key), value, p%fatigue_damage_max_threshold_min, p%fatigue_damage_max_threshold_max, p%fatigue_damage_max_threshold, .FALSE.)

          Case ('fatigue_step')
            Call verifyAndSaveProperty_int(trim(key), value, p%fatigue_step_min, p%fatigue_step_max, p%fatigue_step, .FALSE.)

          Case Default
            Call log%error("loadParameters: Parameter not recognized: " // trim(key))
        End Select
      End If
    End Do ReadLines

    ! Close the .parameters file
    Close(unit)

    Return
  End Subroutine loadParametersFromFile


  Subroutine initializeParameters()

    ! Locals
    Double Precision, Parameter :: zero=0.d0, one=1.d0
    ! -------------------------------------------------------------------- !

    ! Default values
    p%logLevel       = 2
    p%logFormat      = 1
    p%cutbacks_max   = 4
    p%MD_max         = 1000
    p%EQ_max         = 1000
    p%alpha_inc      = 10
    p%tol_DGD_f      = 1.d-4
    p%dGdGc_min      = 1.d-12
    p%compLimit      = 0.25d0
    p%penStiffMult   = 1.d4
    p%cutback_amount = 0.5d0
    p%tol_divergence = 1.d-6
    p%gamma_max      = 4.d0
    p%kb_decompose_thres = 0.99d0
    p%fkt_fiber_failure_angle = -1.d0
    p%schaefer_nr_tolerance = 1.d-6
    p%schaefer_nr_counter_limit = 10000000
    p%terminate_on_no_convergence = .TRUE.
    p%debug_kill_at_total_time = -1.0d0
    p%fkt_random_seed = .FALSE.
    p%fkt_init_misalignment_azi_mu = 0.d0
    p%fkt_init_misalignment_azi_sigma = 45.0d0
    p%fkt_init_misalignment_polar_shape = 0.676d0
    p%fkt_init_misalignment_polar_scale = 2.25d0
    p%fatigue_R_ratio = 0.1d0
    p%cycles_per_increment_init = 1.d-4  ! 10,000 solution increments per fatigue cycle
    p%cycles_per_increment_mod = 0.1d0  ! changes cycles_per_increment by 10% when out of range
    p%cycles_per_increment_max = 1.d5
    p%cycles_per_increment_min = 1.d-5
    p%fatigue_damage_min_threshold = 5.d-6  ! 200,000 solution increments to fail an element at this rate
    p%fatigue_damage_max_threshold = 1.d-4  ! 10,000 solution increments to fail an element at this rate
    p%fatigue_step = 1000000

    ! Maximum and minimum values for parameters to be read from CompDam.parameters file
    p%logLevel_min = 0
    p%logLevel_max = 4

    p%logFormat_min = 1
    p%logFormat_max = 2

    p%cutbacks_max_min = 0
    p%cutbacks_max_max = 10

    p%MD_max_min = 0
    p%MD_max_max = 100000

    p%EQ_max_min = 0
    p%EQ_max_max = 100000

    p%alpha_inc_min = 1
    p%alpha_inc_max = 90

    p%tol_DGD_f_min = Tiny(zero)
    p%tol_DGD_f_max = one

    p%dGdGc_min_min = Tiny(zero)
    p%dGdGc_min_max = one

    p%compLimit_min = Tiny(zero)
    p%compLimit_max = one

    p%penStiffMult_min = 1.d0
    p%penStiffMult_max = 1.d8

    p%cutback_amount_min = Tiny(zero)
    p%cutback_amount_max = one

    p%tol_divergence_min = zero
    p%tol_divergence_max = 100

    p%gamma_max_min = zero
    p%gamma_max_max = 100

    p%kb_decompose_thres_min = zero
    p%kb_decompose_thres_max = one

    p%fkt_fiber_failure_angle_min = -1*Huge(zero)
    p%fkt_fiber_failure_angle_max = 45.d0

    p%schaefer_nr_tolerance_min = Tiny(zero)
    p%schaefer_nr_tolerance_max = one

    p%schaefer_nr_counter_limit_min = 0
    p%schaefer_nr_counter_limit_max = Huge(0)

    p%debug_kill_at_total_time_min = -2.0d0
    p%debug_kill_at_total_time_max = 1000.d0

    p%fkt_init_misalignment_azi_mu_min = -180.d0
    p%fkt_init_misalignment_azi_mu_max = 180.d0

    p%fkt_init_misalignment_azi_sigma_min = Tiny(zero)
    p%fkt_init_misalignment_azi_sigma_max = Huge(zero)

    p%fkt_init_misalignment_polar_shape_min = Tiny(zero)
    p%fkt_init_misalignment_polar_shape_max = Huge(zero)

    p%fkt_init_misalignment_polar_scale_min = Tiny(zero)
    p%fkt_init_misalignment_polar_scale_max = Huge(zero)

    p%fatigue_R_ratio_min = -one
    p%fatigue_R_ratio_max = one

    p%cycles_per_increment_init_min = 1.d-6
    p%cycles_per_increment_init_max = 1.d+6

    p%cycles_per_increment_mod_min = Tiny(zero)
    p%cycles_per_increment_mod_max = 9.d0  ! corresponds to changing cycles_per_increment by a factor of 10

    p%cycles_per_increment_max_min = one
    p%cycles_per_increment_max_max = Huge(zero)

    p%cycles_per_increment_min_min = Tiny(zero)
    p%cycles_per_increment_min_max = one

    p%fatigue_damage_min_threshold_min = Tiny(zero)
    p%fatigue_damage_min_threshold_max = one

    p%fatigue_damage_max_threshold_min = Tiny(zero)
    p%fatigue_damage_max_threshold_max = one

    p%fatigue_step_min = 2
    p%fatigue_step_max = Huge(0)

    Return
  End Subroutine initializeParameters


  Subroutine verifyAndSaveProperty_dbl(key, value, min, max, saveTo, nondefaultWarn_arg)
    ! Checks if the value is within the specified bounds. Prints an error message
    ! which kills the analysis if a value is out of bounds.

    Use forlog_Mod

    !Arguments
    Character(len=*), intent(IN) :: key, value
    Double Precision, intent(IN) :: min, max
    Double Precision, intent(INOUT) :: saveTo
    Logical, intent(IN), optional :: nondefaultWarn_arg

    ! Locals
    Double Precision :: valueDbl
    Logical :: nondefaultWarn
    ! -------------------------------------------------------------------- !

    If (Present(nondefaultWarn_arg)) Then
      nondefaultWarn = nondefaultWarn_arg
    Else
      nondefaultWarn = .TRUE.
    End If

    ! Convert to double
    Read(value,*) valueDbl

    ! Verify that the value is within the specified bounds
    If (valueDbl < min) Then
      Call log%error(" PARAMETER ERROR " // trim(key) // " cannot be less than " // trim(str(min)) // ". Found value: " // trim(str(valueDbl)))

    Else If (valueDbl > max) Then
      Call log%error(" PARAMETER ERROR " // trim(key) // " cannot be greater than " // trim(str(max)) // ". Found value: " // trim(str(valueDbl)))

    End If

    ! Check for non-default
    If ((valueDbl .NE. saveTo) .AND. nondefaultWarn) Then
      Call log%warn("Non-default parameter: " // trim(key) // " = " // trim(str(valueDbl)) // ", default = " // trim(str(saveTo)))
    End If

    ! Save the value and set the flag
    saveTo = valueDbl

    Return
  End Subroutine verifyAndSaveProperty_dbl

  Subroutine verifyAndSaveProperty_int(key, value, min, max, saveTo, nondefaultWarn_arg)
    ! Checks if the value is within the specified bounds. Prints an error message
    ! which kills the analysis if a value is out of bounds.

    Use forlog_Mod

    !Arguments
    Character(len=*), intent(IN) :: key, value
    Integer, intent(IN) :: min, max
    Integer, intent(INOUT) :: saveTo
    Logical, intent(IN), optional :: nondefaultWarn_arg

    ! Locals
    Integer :: valueInt
    Logical :: nondefaultWarn
    ! -------------------------------------------------------------------- !

    If (Present(nondefaultWarn_arg)) Then
      nondefaultWarn = nondefaultWarn_arg
    Else
      nondefaultWarn = .TRUE.
    End If

    ! Convert to integer
    Read(value,*) valueInt

    ! Verify that the value is within the specified bounds
    If (valueInt < min) Then
      Call log%error(" PARAMETER ERROR " // trim(key) // " cannot be less than " // trim(str(min)) // ". Found value: " // trim(str(valueInt)))

    Else If (valueInt > max) Then
      Call log%error(" PARAMETER ERROR " // trim(key) // " cannot be greater than " // trim(str(max)) // ". Found value: " // trim(str(valueInt)))

    End If

    ! Check for non-default
    If ((valueInt .NE. saveTo) .AND. nondefaultWarn) Then
      Call log%warn("Non-default parameter: " // trim(key) // " = " // trim(str(valueInt)) // ", default = " // trim(str(saveTo)))
    End If

    ! Save the value and set the flag
    saveTo = valueInt

    Return
  End Subroutine verifyAndSaveProperty_int

  Subroutine verifyAndSaveProperty_logical(key, value, saveTo, nondefaultWarn_arg)
    ! Checks if the value is true or false. Prints an error message
    ! which kills the analysis if a value is not true or false.

    Use forlog_Mod

    !Arguments
    Character(len=*), intent(IN) :: key, value
    Logical, intent(INOUT) :: saveTo
    Logical, intent(IN), optional :: nondefaultWarn_arg

    ! Locals
    Logical :: valueLogical
    Logical :: nondefaultWarn
    ! -------------------------------------------------------------------- !

    If (Present(nondefaultWarn_arg)) Then
      nondefaultWarn = nondefaultWarn_arg
    Else
      nondefaultWarn = .TRUE.
    End If

    ! Check for T or F
    If ((value /= '.TRUE.') .AND. (value /= '.FALSE.')) Then
      Call log%error('Invalid entry found for parameter ' // trim(key) // ': ' // trim(value))
    End If

    ! Convert to logical
    Read(value,*) valueLogical

    ! Check for non-default
    If ((valueLogical .NE. saveTo) .AND. nondefaultWarn) Then
      Call log%warn("Non-default parameter: " // trim(key) // " = " // trim(str(valueLogical)))
    End If

    ! Save the value and set the flag
    saveTo = valueLogical

    Return
  End Subroutine verifyAndSaveProperty_logical


  Subroutine writeParametersToFile(fileUnit, p)
    ! Writes provided parameters to a file as a python dictionary
    ! Assumes that file opening and closing is handled elsewhere

    ! Arguments
    Integer, intent(IN) :: fileUnit
    Type(parameters), intent(IN) :: p

    ! Locals
    Character(len=32) :: nameValueFmt
    ! -------------------------------------------------------------------- !

    ! Defines the format for writing the floating point numbers
    nameValueFmt = "(A,E21.15E2,A)"

    ! Write the parameters
    write(fileUnit, "(A)") 'p = {'
    write(fileUnit, "(A,I1,A)")   '    "cutbacks_max": ', p%cutbacks_max, ','
    write(fileUnit, "(A,I5,A)")   '    "MD_max": ', p%MD_max, ','
    write(fileUnit, "(A,I5,A)")   '    "EQ_max": ', p%EQ_max, ','
    write(fileUnit, "(A,I2,A)")   '    "alpha_inc": ', p%alpha_inc, ','
    write(fileUnit, nameValueFmt) '    "tol_DGD_f": ', p%tol_DGD_f, ','
    write(fileUnit, nameValueFmt) '    "dGdGc_min": ', p%dGdGc_min, ','
    write(fileUnit, nameValueFmt) '    "compLimit": ', p%compLimit, ','
    write(fileUnit, nameValueFmt) '    "penStiffMult": ', p%penStiffMult, ','
    write(fileUnit, nameValueFmt) '    "cutback_amount": ', p%cutback_amount, ','
    write(fileUnit, nameValueFmt) '    "tol_divergence": ', p%tol_divergence, ','
    write(fileUnit, nameValueFmt) '    "gamma_max": ', p%gamma_max, ','
    write(fileUnit, nameValueFmt) '    "kb_decompose_thres": ', p%kb_decompose_thres, ', '
    write(fileUnit, nameValueFmt) '    "fkt_fiber_failure_angle": ', p%fkt_fiber_failure_angle, ', '
    write(fileUnit, nameValueFmt) '    "schaefer_nr_tolerance": ', p%schaefer_nr_tolerance, ', '
    write(fileUnit, "(A,I9,A)")   '    "schaefer_nr_counter_limit": ', p%schaefer_nr_counter_limit, ', '
    If (p%terminate_on_no_convergence) Then
      write(fileUnit,"(A)") '    "terminate_on_no_convergence": True,'
    Else
      write(fileUnit,"(A)") '    "terminate_on_no_convergence": False,'
    End If
    write(fileUnit, nameValueFmt) '    "debug_kill_at_total_time": ', p%debug_kill_at_total_time, ', '
    If (p%fkt_random_seed) Then
      write(fileUnit,"(A)") '    "fkt_random_seed": True,'
    Else
      write(fileUnit,"(A)") '    "fkt_random_seed": False,'
    End If
    write(fileUnit, nameValueFmt) '    "fkt_init_misalignment_azi_mu": ', p%fkt_init_misalignment_azi_mu, ', '
    write(fileUnit, nameValueFmt) '    "fkt_init_misalignment_azi_sigma": ', p%fkt_init_misalignment_azi_sigma, ', '
    write(fileUnit, nameValueFmt) '    "fkt_init_misalignment_polar_shape": ', p%fkt_init_misalignment_polar_shape, ', '
    write(fileUnit, nameValueFmt) '    "fkt_init_misalignment_polar_scale": ', p%fkt_init_misalignment_polar_scale
    write(fileUnit, "(A)") '}'

  End Subroutine writeParametersToFile


End Module parameters_Mod
