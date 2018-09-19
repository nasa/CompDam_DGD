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
    Double Precision :: tol_DGD_f                                 ! tol_DGD_f = tol_DGD/YT
    Double Precision :: dGdGc_min                                 ! minimum amount of damage dissipated per MD increment
    Double Precision :: compLimit                                 ! minimum accepted det(F_bulk)
    Double Precision :: penStiffMult                              ! penalty stiffness multiplier
    Double Precision :: cutback_amount                            ! artificially slow rate of change of F_bulk
    Double Precision :: tol_divergence                            ! Tolerance for divergence of internal Newton Raphson loop
    Double Precision :: gamma_max                                 ! Maximum shear strain; when this value is exceeded, the element is deleted
    Double Precision :: kb_decompose_thres                        ! Ratio of kink band size to element length at which to decompose the element

    ! min and max values for acceptable range
    Integer, private :: logLevel_min, logLevel_max
    Integer, private :: logFormat_min, logFormat_max
    Integer, private :: cutbacks_max_min, cutbacks_max_max
    Integer, private :: MD_max_min, MD_max_max
    Integer, private :: EQ_max_min, EQ_max_max
    Integer, private :: alpha_inc_min, alpha_inc_max
    Double Precision, private :: tol_DGD_f_min, tol_DGD_f_max
    Double Precision, private :: dGdGc_min_min, dGdGc_min_max
    Double Precision, private :: compLimit_min, compLimit_max
    Double Precision, private :: penStiffMult_min, penStiffMult_max
    Double Precision, private :: cutback_amount_min, cutback_amount_max
    Double Precision, private :: tol_divergence_min, tol_divergence_max
    Double Precision, private :: gamma_max_min, gamma_max_max
    Double Precision, private :: kb_decompose_thres_min, kb_decompose_thres_max

  End Type parameters

  ! Public interface
  Public :: parameters
  Public :: loadParameters

  ! Reference to object for singleton
  type(parameters), Save :: p

Contains

  Type(parameters) Function loadParameters()
    ! Loads parameters into module variable p

    Use forlog_Mod

    ! Locals
    Character(len=256) :: outputDir, fileName
    Integer :: lenOutputDir
    Logical :: fileExists
    ! -------------------------------------------------------------------- !

    ! Initializations
    Call initializeParameters()

#ifndef PYEXT
    ! Get the output directory (location to search for CompDam.parameters file)
    Call VGETOUTDIR(outputDir, lenOutputDir)

    ! Look to see if a parameters file exists
    ! First look for: CompDam.parameters
    fileName = trim(outputDir) // '/CompDam.parameters'
#else
    fileName = 'CompDam.parameters'
#endif
    Inquire(FILE=fileName, EXIST=fileExists)

    ! If the file is present, load parameters from file
    If (fileExists) Call loadParametersFromFile(trim(fileName))

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
    p%tol_divergence = 0.1d0
    p%gamma_max      = 4.d0
    p%kb_decompose_thres = 0.99d0


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

    Return
  End Subroutine initializeParameters


  Subroutine verifyAndSaveProperty_dbl(key, value, min, max, saveTo)
    ! Checks if the value is within the specified bounds. Prints an error message
    ! which kills the analysis if a value is out of bounds.

    Use forlog_Mod

    !Arguments
    Character(len=*), intent(IN) :: key, value
    Double Precision, intent(IN) :: min, max
    Double Precision, intent(OUT) :: saveTo

    ! Locals
    Double Precision :: valueDbl
    ! -------------------------------------------------------------------- !

    ! Convert to double
    Read(value,*) valueDbl

    ! Verify that the value is within the specified bounds
    If (valueDbl < min) Then
      Call log%error(" PARAMETER ERROR " // trim(key) // " cannot be less than " // trim(str(min)) // ". Found value: " // trim(str(valueDbl)))

    Else If (valueDbl > max) Then
      Call log%error(" PARAMETER ERROR " // trim(key) // " cannot be greater than " // trim(str(max)) // ". Found value: " // trim(str(valueDbl)))

    End If

    ! Save the value and set the flag
    saveTo = valueDbl

    Return
  End Subroutine verifyAndSaveProperty_dbl

  Subroutine verifyAndSaveProperty_int(key, value, min, max, saveTo)
    ! Checks if the value is within the specified bounds. Prints an error message
    ! which kills the analysis if a value is out of bounds.

    Use forlog_Mod

    !Arguments
    Character(len=*), intent(IN) :: key, value
    Integer, intent(IN) :: min, max
    Integer, intent(OUT) :: saveTo

    ! Locals
    Integer :: valueInt
    ! -------------------------------------------------------------------- !

    ! Convert to integer
    Read(value,*) valueInt

    ! Verify that the value is within the specified bounds
    If (valueInt < min) Then
      Call log%error(" PARAMETER ERROR " // trim(key) // " cannot be less than " // trim(str(min)) // ". Found value: " // trim(str(valueInt)))

    Else If (valueInt > max) Then
      Call log%error(" PARAMETER ERROR " // trim(key) // " cannot be greater than " // trim(str(max)) // ". Found value: " // trim(str(valueInt)))

    End If

    ! Save the value and set the flag
    saveTo = valueInt

    Return
  End Subroutine verifyAndSaveProperty_int


End Module parameters_Mod
