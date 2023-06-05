! #define type(x) TYPE(x), target

Module forlog_Mod
  ! General module for logging information in abaqus subroutines
  ! Using this module exports a variable log to access the module features
  !
  ! The log%level controls the verbosity
  !    level=0  All output is suppressed (including warnings and errors)
  !    level=1  Only issues that terminate the analysis are logged
  ! -> level=2  (Recommended) Issues that may impact the accuracy of the results are logged
  !    level=3  Verbose logging
  !    level=4  Extremely verbose logging for debugging. Generates large logs and increases analysis run time.
  !
  ! Info on constructions used here:
  ! https://software.intel.com/en-us/node/579826
  ! https://www.pgroup.com/lit/articles/insider/v3n1a3.htm
  ! http://fortranwiki.org/fortran/show/Object-oriented+programming


  Use vumatArg_Mod
#ifndef PYEXT
  Use version_Mod
#endif

  ! Definition of l
  Type forlog
    Integer :: level
    Integer :: fileUnit                              ! Fortran file unit number for the output generated by forlog
    Type(vumatArg) :: arg                            ! Enables access to all abaqus arguments
#ifndef PYEXT
    Character(len=255) :: fileName
    Integer :: format                                ! Specifies log file format (1=Single line, human readable, 10=JSON)
    Logical, Private :: hashPrinted = .FALSE.
#endif
  Contains
#ifndef PYEXT
    Procedure :: init                                ! Initializes the logger
    Procedure :: close
    Procedure :: location
#endif
    Procedure :: debug_str
    Procedure :: debug_int
    Procedure :: debug_dbl
    Procedure :: debug_vec
    Procedure :: debug_mat
    Generic :: debug => debug_str, debug_int, debug_dbl, debug_vec, debug_mat  ! level=4, Maximum amount of information
    Procedure :: info                                ! level=3, Verbose logging
    Procedure :: warn                                ! level=2, Log issues that may impact accuracy of results
    Procedure :: error                               ! level=1, Log issues that immediately terminate the analysis
    Procedure :: terminate                           ! level=1, log issues that cleanly terminate the anlaysis
  End Type forlog

  Interface str
    Module Procedure str_int, str_double, str_real, str_logical
  End Interface

  ! For public access to forlog
  Type(forlog), Save, Public :: log


Contains

#ifndef PYEXT
  ! To create a logger
  Subroutine init(this, level, vumatArgStruct, format)

    ! Arguments
    Class(forlog), intent(INOUT) :: this
    Integer, intent(IN) :: level                                 ! Any message <= to the level specified here is written
    Type(vumatArg), intent(IN) :: vumatArgStruct                 ! Enables access to all abaqus arguments
    Integer, intent(IN) :: format                                ! Specifies log file format (1=Single line, human readable, 2=JSON)

    ! Locals
    Integer :: stat, lenOutputDir, lenJobName
    Character(len=256) :: outputDir, jobName, cmd, line

    ! Parameters
    Double Precision :: zero
    Parameter (zero=0.d0)

    ! Initialize class variables
    this%fileUnit = 6
    this%level = level
    this%arg = vumatArgStruct
    this%format = format

    If (this%format == 2) Then
      ! Initialize output file
      ! Load the output directory
      CALL VGETOUTDIR(outputDir, lenOutputDir)
      ! Load the job name
      CALL VGETJOBNAME(jobName, lenJobName)

      ! Use a different file unit number
      this%fileUnit = 106

      ! Set the file name
      this%fileName = trim(outputDir) // '/' // trim(jobName) // '_debug.csv'

      ! Open the log file
      open(this%fileUnit, file=trim(this%fileName), position='append', recl=1000)
    End If

    ! Git SHA-1 hash
    If (this%arg%totalTime == zero .AND. .NOT. this%hashPrinted) Then
      print *, 'Git hash: ' // trim(hash)
      print *, 'Git commit timestamp: ' // trim(timestamp)
      this%hashPrinted = .TRUE.
    End If

    ! Store reference to instance
    log = this
  End Subroutine init


  ! Tidy up
  Subroutine close(this)
    ! Arguments
    Class(forlog) :: this

    close(this%fileUnit)
  End Subroutine close


  ! Sets log%level
  Subroutine setLogLevel(l)

    ! Arguments
    Integer, intent(IN) :: l

    log%level = l

    Return
  End Subroutine setLogLevel


  ! Print current location info
  Subroutine location(this, logLevel)

    ! Arguments
    Class(forlog), intent(IN) :: this
    Integer, intent(IN) :: logLevel

    ! Locals
    Dimension INTV(1), REALV(1)    ! For abaqus warning messages
    Character(len=8) CHARV(1)      ! For Abaqus warning messages

    If (this%level >= logLevel) Then
      write(this%fileUnit,*) trim(str_double(this%arg%totalTime)) // ", Element No.: " // trim(str_int(this%arg%nElement)) // ", Int pt: " // trim(str_int(this%arg%nMatPoint))
      ! Call XPLB_ABQERR(-3,msg,INTV,REALV,CHARV)
    End If
  End Subroutine location

  ! Private subroutine for writing records
  Subroutine writeToLog(this, msg)
    ! Arguments
    Class(forlog), intent(IN) :: this
    Character*(*), intent(IN) :: msg

    write(this%fileUnit,*) trim(str_double(this%arg%totalTime)) // msg
  End Subroutine writeToLog

#endif

  ! Most verbose logging
  Subroutine debug_str(this, msg)

    ! Arguments
    Class(forlog), intent(IN) :: this
    Character*(*), intent(IN) :: msg

    If (this%level >= 4) Then
#ifndef PYEXT
      Call writeToLog(this,  ", DEBUG, " // msg)
#else
      write(this%fileUnit,*) "DEBUG, " // msg
#endif
    End If
  End Subroutine debug_str


  Subroutine debug_int(this, msg, int)

    ! Arguments
    Class(forlog), intent(IN) :: this
    Character*(*), intent(IN) :: msg
    Integer, intent(IN) :: int

    If (this%level >= 4) Then
      write(this%fileUnit,*) msg // " = " // str_int(int)
    End If
  End Subroutine debug_int

  Subroutine debug_dbl(this, msg, dbl)

    ! Arguments
    Class(forlog), intent(IN) :: this
    Character*(*), intent(IN) :: msg
    Double Precision, intent(IN) :: dbl

    If (this%level >= 4) Then
      write(this%fileUnit,*) msg // " = " // str_double(dbl)
    End If
  End Subroutine debug_dbl

  Subroutine debug_vec(this, msg, vec)

    ! Arguments
    Class(forlog), intent(IN) :: this
    Character*(*), intent(IN) :: msg
    Double Precision, intent(IN) :: vec(:)

    If (this%level >= 4) Then
      write(this%fileUnit,*) msg // " = "
      Do i=1,size(vec)
        write(this%fileUnit,*) str_double(vec(i))
      End Do
    End If
  End Subroutine debug_vec

  Subroutine debug_mat(this, msg, mat)

    ! Arguments
    Class(forlog), intent(IN) :: this
    Character*(*), intent(IN) :: msg
    Double Precision, intent(IN) :: mat(:,:)

    Integer :: mat_shape(2)

    If (this%level >= 4) Then
      write(this%fileUnit,*) msg // " = "
      mat_shape = SHAPE(mat)
      Do i=1,mat_shape(1)
        If (mat_shape(2) == 2) Then
          write(this%fileUnit,*) str_double(mat(i,1)) // ", " // str_double(mat(i,2))
        Else If (mat_shape(2) == 3) Then
          write(this%fileUnit,*) str_double(mat(i,1)) // ", " // str_double(mat(i,2)) // ", " // str_double(mat(i,3))
        End If
      End Do
    End If
  End Subroutine debug_mat

  ! General information (verbose logging)
  Subroutine info(this, msg)

    ! Arguments
    Class(forlog), intent(IN) :: this
    Character*(*), intent(IN) :: msg

    If (this%level >= 3) Then
#ifndef PYEXT
      Call writeToLog(this, ", INFO, " // msg)
#else
      write(this%fileUnit,*) "INFO, " // msg
#endif
    End If
  End Subroutine info


  ! Issues that may impact accuracy of results
  Subroutine warn(this, msg)

    ! Arguments
    Class(forlog), intent(IN) :: this
    Character*(*), intent(IN) :: msg

    ! Locals
    Dimension INTV(1), REALV(1)    ! For abaqus warning messages
    Character(len=8) CHARV(1)      ! For Abaqus warning messages

    If (this%level >= 2) Then
#ifndef PYEXT
      Call writeToLog(this, ", WARNING, " // msg)
      Call XPLB_ABQERR(-1,msg,INTV,REALV,CHARV)
#else
      write(this%fileUnit,*) "WARN, " // msg
#endif
    End If
  End Subroutine warn


  ! Issues that require the analysis be immediately terminated, e.g., numerical error
  Subroutine error(this, msg)

    ! Arguments
    Class(forlog), intent(IN) :: this
    Character*(*), intent(IN) :: msg

    ! Locals
    Dimension INTV(1), REALV(1)    ! For Abaqus warning messages
    Character(len=8) CHARV(1)      ! For Abaqus warning messages

    If (this%level >= 1) Then
#ifndef PYEXT
      Call writeToLog(this, ", ERROR, " // msg)
      Call XPLB_ABQERR(-3,msg,INTV,REALV,CHARV)
#else
      write(this%fileUnit,*) "ERROR, " // msg
      stop 1
#endif
    End If
  End Subroutine error


  ! Issues that require the analysis be cleanly terminated, e.g., nonconvergence
  Subroutine terminate(this, msg)

    ! Arguments
    Class(forlog), intent(IN) :: this
    Character*(*), intent(IN) :: msg

    ! Locals
    Dimension INTV(1), REALV(1)    ! For abaqus warning messages
    Character(len=8) CHARV(1)      ! For Abaqus warning messages
    Common /analysis_termination/ analysis_status

    If (this%level >= 1) Then
#ifndef PYEXT
      Call writeToLog(this, ", ERROR, " // msg)
      analysis_status = 0
      Call XPLB_ABQERR(-2,msg,INTV,REALV,CHARV)
#else
      write(this%fileUnit,*) "ERROR, " // msg
      stop 1
#endif
    End If
  End Subroutine terminate


  ! Convert an integer to string
  Character(len=30) function str_int(k)
    integer, intent(in) :: k
    write (str_int, *) k
    str_int = adjustl(str_int)
  End Function str_int


  ! Convert a double to a string
  Character(len=30) Function str_double(k)
    Double Precision, intent(in) :: k
    write (str_double, *) k
    str_double = adjustl(str_double)
  End function str_double


  ! Convert a single to a string
  Character(len=30) Function str_real(k)
    Real, intent(in) :: k
    write (str_real, *) k
    str_real = adjustl(str_real)
  End Function str_real


  ! Convert a logical to a string
  Character(len=30) Function str_logical(k)
    Logical, intent(in) :: k
    write (str_logical, *) k
    str_logical = adjustl(str_logical)
  End Function str_logical


End Module forlog_Mod
