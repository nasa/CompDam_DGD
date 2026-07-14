#define type(x) TYPE(x), target

Module parameters_Mod
  Type parameters
    ! Solution parameters
    Integer :: logLevel = 2                                       ! Logger level
    Integer :: logFormat = 1                                      ! 1=default (to log file), 2=csv
    Integer :: cutbacks_max = 4                                   ! maximum number of cut-backs
    Integer :: MD_max = 1000                                      ! maximum number of damage increments per solution increment
    Integer :: EQ_max = 1000                                      ! maximum number of equilibrium iterations
    Integer :: alpha_search = 1                                   ! whether or not to search for the matrix crack angle with maximum failure index
    Integer :: alpha_inc = 10                                     ! increment, in degrees, for which the matrix failure criterion is evaluated
    Integer :: alpha_max = 60                                     ! maximum matrix crack angle alpha, in degrees
    Integer :: schaefer_nr_counter_limit = 10000000               ! maximum number of Newton-Raphson loops before divergence is assumed
    Double Precision :: tol_DGD_f = 1.d-4                         ! tol_DGD_f = tol_DGD/YT
    Double Precision :: dGdGc_min = 1.d-12                        ! minimum amount of damage dissipated per MD increment
    Double Precision :: compLimit = 0.25d0                        ! minimum accepted det(F_bulk)
    Double Precision :: penStiffMult = 1.d3                       ! penalty stiffness multiplier
    Double Precision :: cutback_amount = 0.5d0                    ! artificially slow rate of change of F_bulk
    Double Precision :: tol_divergence = 1.d-6                    ! Tolerance for divergence of internal Newton Raphson loop
    Double Precision :: gamma_max = 4.d0                          ! Maximum shear strain; when this value is exceeded, the element is deleted
    Double Precision :: kb_decompose_thres = 0.99d0               ! Ratio of kink band size to element length at which to decompose the element
    Double Precision :: fkt_fiber_failure_angle = -1.d0           ! Angle at which fiber failure occurs; no further plastic shear deformation allowed
    Double Precision :: schaefer_nr_tolerance = 1.d-6             ! Tolerance value used to determine if convergence has occurred in newton raphson loop
    Integer :: terminate_on_no_convergence = 1                    ! Set to false to trigger a warning, delete the element, and allow the analysis to continue when no converged solution can be found
    Double Precision :: debug_kill_at_total_time = -1.d0          ! Time at which to kill analysis (for debugging purposes)
    Integer :: fkt_random_seed = 0                                ! Set to true to use the time as the seed for initial fiber misalignment to get different results for each realization
    Double Precision :: fkt_init_misalignment_azi_mu = 0.d0       ! Initial fiber misalignment azimuthal average [degrees]
    Double Precision :: fkt_init_misalignment_azi_sigma = 45.d0   ! Initial fiber misalignment azimuthal variance [degrees]
    Double Precision :: fkt_init_misalignment_polar_shape = 0.676d0 ! Initial fiber misalignment polar shape parameter [degrees]
    Double Precision :: fkt_init_misalignment_polar_scale = 2.25d0  ! Initial fiber misalignment polar standard deviation [degrees]
    Double Precision :: fatigue_R_ratio = 0.1d0                   ! R ratio for cohesive fatigue model, sigma_min / sigma_max
    Double Precision :: cycles_per_increment_init = 1.d-4         ! Fatigue cycles per solution increment
    Double Precision :: cycles_per_increment_max = 1.d5           ! Maximum fatigue cycles per solution increment
    Double Precision :: cycles_per_increment_min = 1.d-5          ! Minimum fatigue cycles per solution increment
    Double Precision :: cycles_per_increment_mod = 0.1d0          ! Percent change in cycles_per_increment when out of range
    Double Precision :: fatigue_damage_min_threshold = 5.d-6      ! Minimum required incremental fatigue damage to count as progression
    Double Precision :: fatigue_damage_max_threshold = 1.d-4      ! Maximum allowable incremental fatigue damage per solution increment
    Integer :: fatigue_step = 100000                              ! Step number for fatigue analysis
    Integer :: check_for_snap_back = 0                            ! Whether or not to check for constitutive snap-back
    Integer :: set_status_0_on_d2 = 0                             ! Whether of not the status state variable is set to 0 when d2 = 1
    Integer :: set_status_0_on_d1T = 0                            ! Whether of not the status state variable is set to 0 when d1T = 1
    Integer :: set_status_0_on_d1C = 0                            ! Whether of not the status state variable is set to 0 when d1C = 1

    ! min and max values for acceptable range
    Integer :: logLevel_min = 0, logLevel_max = 4
    Integer :: logFormat_min = 1, logFormat_max = 2
    Integer :: cutbacks_max_min = 0, cutbacks_max_max = 10
    Integer :: MD_max_min = 0, MD_max_max = 100000
    Integer :: EQ_max_min = 0, EQ_max_max = 100000
    Integer :: alpha_search_min = 0, alpha_search_max = 1
    Integer :: alpha_inc_min = 1, alpha_inc_max = 90
    Integer :: alpha_max_min = 0, alpha_max_max = 90
    Integer :: schaefer_nr_counter_limit_min = 0, schaefer_nr_counter_limit_max = Huge(0)
    Integer :: terminate_on_no_convergence_min = 0, terminate_on_no_convergence_max = 1
    Integer :: fkt_random_seed_min = 0, fkt_random_seed_max = 1
    Double Precision :: tol_DGD_f_min = Tiny(0.d0), tol_DGD_f_max = 1.d0
    Double Precision :: dGdGc_min_min = Tiny(0.d0), dGdGc_min_max = 1.d0
    Double Precision :: compLimit_min = Tiny(0.d0), compLimit_max = 1.d0
    Double Precision :: penStiffMult_min = 1.d0, penStiffMult_max= 1.d8
    Double Precision :: cutback_amount_min = Tiny(0.d0), cutback_amount_max = 1.d0
    Double Precision :: tol_divergence_min = 0.d0, tol_divergence_max = 100.d0
    Double Precision :: gamma_max_min = 0.d0, gamma_max_max = 100.d0
    Double Precision :: kb_decompose_thres_min = 0.d0, kb_decompose_thres_max = 1.d0
    Double Precision :: fkt_fiber_failure_angle_min = -1*Huge(0.d0), fkt_fiber_failure_angle_max = 45.d0
    Double Precision :: schaefer_nr_tolerance_min = Tiny(0.d0), schaefer_nr_tolerance_max = 1.d0
    Double Precision :: debug_kill_at_total_time_min = -2.d0, debug_kill_at_total_time_max = 1000.d0
    Double Precision :: fkt_init_misalignment_azi_mu_min = -180.d0, fkt_init_misalignment_azi_mu_max = 180.d0
    Double Precision :: fkt_init_misalignment_azi_sigma_min = Tiny(0.d0), fkt_init_misalignment_azi_sigma_max = Huge(0.d0)
    Double Precision :: fkt_init_misalignment_polar_shape_min = Tiny(0.d0), fkt_init_misalignment_polar_shape_max = Huge(0.d0)
    Double Precision :: fkt_init_misalignment_polar_scale_min = Tiny(0.d0), fkt_init_misalignment_polar_scale_max = Huge(0.d0)
    Double Precision :: fatigue_R_ratio_min = -1.d0, fatigue_R_ratio_max = 1.d0
    Double Precision :: cycles_per_increment_init_min = Tiny(0.d0), cycles_per_increment_init_max = Huge(0.d0)
    Double Precision :: cycles_per_increment_max_min = Tiny(0.d0), cycles_per_increment_max_max = Huge(0.d0)
    Double Precision :: cycles_per_increment_min_min = Tiny(0.d0), cycles_per_increment_min_max = Huge(0.d0)
    Double Precision :: cycles_per_increment_mod_min = Tiny(0.d0), cycles_per_increment_mod_max = Huge(0.d0)
    Double Precision :: fatigue_damage_min_threshold_min = Tiny(0.d0), fatigue_damage_min_threshold_max = 1.d0
    Double Precision :: fatigue_damage_max_threshold_min = Tiny(0.d0), fatigue_damage_max_threshold_max = 1.d0
    Integer :: fatigue_step_min = 2, fatigue_step_max = Huge(0)
    Integer :: check_for_snap_back_min = 0, check_for_snap_back_max = 1
    Integer :: set_status_0_on_d2_min = 0, set_status_0_on_d2_max = 1
    Integer :: set_status_0_on_d1T_min = 0, set_status_0_on_d1T_max = 1
    Integer :: set_status_0_on_d1C_min = 0, set_status_0_on_d1C_max = 1

  End Type parameters

  ! Public interface
  Public :: verifyAndSaveProperty_dbl
  Public :: verifyAndSaveProperty_int
  Public :: writeParametersToFile

Contains
  Subroutine verifyAndSaveProperty_dbl(key, valueDbl, min, max, saveTo, print_logs)
    ! Checks if the value is within the specified bounds. Prints an error message
    ! which kills the analysis if a value is out of bounds.

    Use forlog_Mod

    ! Arguments
    Character(len=*), intent(IN) :: key
    Double Precision, intent(IN) :: valueDbl, min, max
    Double Precision, intent(INOUT) :: saveTo
    Logical, intent(IN), optional :: print_logs

    ! Locals
    Dimension INTV(1), REALV(1)    ! For Abaqus warning messages
    Character(len=8) CHARV(1)      ! For Abaqus warning messages
    ! -------------------------------------------------------------------- !

    ! Verify that the value is within the specified bounds
    If (valueDbl < min) Then
#ifndef PYEXT
      Call XPLB_ABQERR(-3, " PARAMETER ERROR " // trim(key) // " cannot be less than " // trim(str(min)) // ". Found value: " // trim(str(valueDbl)), INTV, REALV, CHARV)
#else
      Print *, " PARAMETER ERROR " // trim(key) // " cannot be less than " // trim(str(min)) // ". Found value: " // trim(str(valueDbl))
      stop 1
#endif

    Else If (valueDbl > max) Then
#ifndef PYEXT
      Call XPLB_ABQERR(-3, " PARAMETER ERROR " // trim(key) // " cannot be greater than " // trim(str(max)) // ". Found value: " // trim(str(valueDbl)), INTV, REALV, CHARV)
#else
      Print *, " PARAMETER ERROR " // trim(key) // " cannot be greater than " // trim(str(max)) // ". Found value: " // trim(str(valueDbl))
      stop 1
#endif

    End If

    ! Check for non-default
    If (valueDbl .NE. saveTo) Then
      If (print_logs) Print *, "  Non-default parameter: " // trim(key) // " = " // trim(str(valueDbl)) // ", default = " // trim(str(saveTo))
    Else
      If (print_logs)Print *, "  Default parameter: " // trim(key) // " = " // trim(str(valueDbl))
    End If

    ! Save the value and set the flag
    saveTo = valueDbl

    Return
  End Subroutine verifyAndSaveProperty_dbl

  Subroutine verifyAndSaveProperty_int(key, valueInt, min, max, saveTo, print_logs)
    ! Checks if the value is within the specified bounds. Prints an error message
    ! which kills the analysis if a value is out of bounds.

    Use forlog_Mod

    !Arguments
    Character(len=*), intent(IN) :: key
    Integer, intent(IN) :: valueInt, min, max
    Integer, intent(INOUT) :: saveTo
    Logical, intent(IN), optional :: print_logs

    ! Locals
    Dimension INTV(1), REALV(1)    ! For Abaqus warning messages
    Character(len=8) CHARV(1)      ! For Abaqus warning messages
    ! -------------------------------------------------------------------- !

    ! Verify that the value is within the specified bounds
    If (valueInt < min) Then
#ifndef PYEXT
      Call XPLB_ABQERR(-3," PARAMETER ERROR " // trim(key) // " cannot be less than " // trim(str(min)) // ". Found value: " // trim(str(valueInt)),INTV,REALV,CHARV)
#else
      Print *, " PARAMETER ERROR " // trim(key) // " cannot be less than " // trim(str(min)) // ". Found value: " // trim(str(valueInt))
      stop 1
#endif

    Else If (valueInt > max) Then
#ifndef PYEXT
      Call XPLB_ABQERR(-3," PARAMETER ERROR " // trim(key) // " cannot be greater than " // trim(str(max)) // ". Found value: " // trim(str(valueInt)),INTV,REALV,CHARV)
#else
      Print *, " PARAMETER ERROR " // trim(key) // " cannot be greater than " // trim(str(max)) // ". Found value: " // trim(str(valueInt))
      stop 1
#endif

    End If

    ! Check for non-default
    If (valueInt .NE. saveTo) Then
      If (print_logs) Print *, "  Non-default parameter: " // trim(key) // " = " // trim(str(valueInt)) // ", default = " // trim(str(saveTo))
    Else
      If (print_logs) Print *, "  Default parameter: " // trim(key) // " = " // trim(str(valueInt))
    End If

    ! Save the value and set the flag
    saveTo = valueInt

    Return
  End Subroutine verifyAndSaveProperty_int

  Subroutine writeParametersToFile(fileUnit, p)
    ! Writes provided parameters to a file as a python dictionary
    ! Assumes that file opening and closing is handled elsewhere

    ! Arguments
    Integer, intent(IN) :: fileUnit
    Type(parameters), intent(IN) :: p

    ! Locals
    Character(len=32) :: fmtFloat, fmtInt
    ! -------------------------------------------------------------------- !

    ! Defines the format for writing the floating point numbers
    fmtFloat = "(A,E21.15E2,A)"

    ! Write the parameters
    write(fileUnit, "(A)") 'p = {'
    write(fileUnit, "(A,I1,A)") '    "cutbacks_max": ', p%cutbacks_max, ','
    write(fileUnit, "(A,I5,A)") '    "MD_max": ', p%MD_max, ','
    write(fileUnit, "(A,I5,A)") '    "EQ_max": ', p%EQ_max, ','
    write(fileUnit, "(A,I2,A)") '    "alpha_search": ', p%alpha_search, ','
    write(fileUnit, "(A,I2,A)") '    "alpha_inc": ', p%alpha_inc, ','
    write(fileUnit, "(A,I2,A)") '    "alpha_max": ', p%alpha_max, ','
    write(fileUnit, fmtFloat)   '    "tol_DGD_f": ', p%tol_DGD_f, ','
    write(fileUnit, fmtFloat)   '    "dGdGc_min": ', p%dGdGc_min, ','
    write(fileUnit, fmtFloat)   '    "compLimit": ', p%compLimit, ','
    write(fileUnit, fmtFloat)   '    "penStiffMult": ', p%penStiffMult, ','
    write(fileUnit, fmtFloat)   '    "cutback_amount": ', p%cutback_amount, ','
    write(fileUnit, fmtFloat)   '    "tol_divergence": ', p%tol_divergence, ','
    write(fileUnit, fmtFloat)   '    "gamma_max": ', p%gamma_max, ','
    write(fileUnit, fmtFloat)   '    "kb_decompose_thres": ', p%kb_decompose_thres, ', '
    write(fileUnit, fmtFloat)   '    "fkt_fiber_failure_angle": ', p%fkt_fiber_failure_angle, ', '
    write(fileUnit, fmtFloat)   '    "schaefer_nr_tolerance": ', p%schaefer_nr_tolerance, ', '
    write(fileUnit, "(A,I9,A)") '    "schaefer_nr_counter_limit": ', p%schaefer_nr_counter_limit, ', '
    write(fileUnit, "(A,I2,A)") '    "terminate_on_no_convergence": ', p%terminate_on_no_convergence, ','
    write(fileUnit, "(A,I2,A)") '    "check_for_snap_back": ', p%check_for_snap_back, ','
    write(fileUnit, "(A,I2,A)") '    "fkt_random_seed": ', p%fkt_random_seed, ','
    write(fileUnit, fmtFloat)   '    "fkt_init_misalignment_azi_mu": ', p%fkt_init_misalignment_azi_mu, ', '
    write(fileUnit, fmtFloat)   '    "fkt_init_misalignment_azi_sigma": ', p%fkt_init_misalignment_azi_sigma, ', '
    write(fileUnit, fmtFloat)   '    "fkt_init_misalignment_polar_shape": ', p%fkt_init_misalignment_polar_shape, ', '
    write(fileUnit, fmtFloat)   '    "fkt_init_misalignment_polar_scale": ', p%fkt_init_misalignment_polar_scale, ', '
    write(fileUnit, "(A,I6,A)") '    "fatigue_step": ', p%fatigue_step, ','
    write(fileUnit, fmtFloat)   '    "fatigue_R_ratio": ', p%fatigue_R_ratio, ', '
    write(fileUnit, fmtFloat)   '    "cycles_per_increment_init": ', p%cycles_per_increment_init, ', '
    write(fileUnit, fmtFloat)   '    "cycles_per_increment_max": ', p%cycles_per_increment_max, ', '
    write(fileUnit, fmtFloat)   '    "cycles_per_increment_min": ', p%cycles_per_increment_min, ', '
    write(fileUnit, fmtFloat)   '    "cycles_per_increment_mod": ', p%cycles_per_increment_mod, ', '
    write(fileUnit, fmtFloat)   '    "fatigue_damage_min_threshold": ', p%fatigue_damage_min_threshold, ', '
    write(fileUnit, fmtFloat)   '    "fatigue_damage_max_threshold": ', p%fatigue_damage_max_threshold, ', '
    write(fileUnit, "(A,I2,A)") '    "set_status_0_on_d2": ', p%set_status_0_on_d2, ','
    write(fileUnit, "(A,I2,A)") '    "set_status_0_on_d1T": ', p%set_status_0_on_d1T, ','
    write(fileUnit, "(A,I2,A)") '    "set_status_0_on_d1C": ', p%set_status_0_on_d1C, ','
    write(fileUnit, "(A)") '}'

  End Subroutine writeParametersToFile
End Module parameters_Mod
