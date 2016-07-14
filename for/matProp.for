
Module matProp_Mod
  ! Module for loading and validating material properties

  Private

  Type matProps
    ! Stores the set of properties

    Character*80 :: name ! Material name, required

    ! Material properties, grouped by feature
    ! (Required inputs)
    Double Precision :: E1, E2, G12, v12, v23                            ! Req. Transversely isotropic.
    ! (Optional inputs)
    Double Precision :: YT, SL, GYT, GSL, eta_BK, alpha0, YC             ! Opt. Matrix failure
    Double Precision :: E3, G13, G23, v13                                ! Opt. Orthotropic. Defaults are calculated based on transverse isotropy
    Double Precision :: cte(3)                                           ! Opt. CTEs (Required if a temperature change is applied in the analysis). Defaults to zero
    Double Precision :: aPL, nPL                                         ! Opt. Shear nonlinearity
    Double Precision :: XT, fXT, GXT, fGXT                               ! Opt. Longitudinal tensile damage (max strain failure criterion, bilinear softening)
    Double Precision :: XC, fXC, GXC, fGXC                               ! Opt. Longitudinal compressive damage (max strain failure criterion, bilinear softening)
    Double Precision :: mu                                               ! Opt. Friction. Defaults to zero

    ! min and max values for acceptable range
    Double Precision, private :: modulus_min, modulus_max
    Double Precision, private :: poissons_min, poissons_max
    Double Precision, private :: strength_min, strength_max
    Double Precision, private :: toughness_min, toughness_max
    Double Precision, private :: eta_BK_min, eta_BK_max
    Double Precision, private :: cte_min, cte_max
    Double Precision, private :: aPL_min, aPL_max, nPL_min, nPL_max
    Double Precision, private :: alpha0_min, alpha0_max
    Double Precision, private :: mu_min, mu_max

    ! Flags that indicate if values have been set
    Logical, private :: E1_def, E2_def, G12_def, v12_def, v23_def
    Logical, private :: YT_def, SL_def, GYT_def, GSL_def, eta_BK_def
    Logical, private :: E3_def, G13_def, G23_def, v13_def
    Logical, private :: cte_def(3)
    Logical, private :: aPL_def, nPL_def
    Logical, private :: XT_def, fXT_def, GXT_def, fGXT_def
    Logical, private :: XC_def, fXC_def, GXC_def, fGXC_def
    Logical, private :: YC_def, alpha0_def
    Logical, private :: mu_def

    ! Calculated properties
    Double Precision :: v21, v31, v32
    Double Precision :: etaL, etaT, ST, alpha0_deg

    Integer :: strainDef
    Double Precision :: thickness

    ! Flags for features
    Logical :: matrixDam, shearNonlinearity, fiberTenDam, fiberCompDam, fiberCompDamNew, friction

  End Type matProps

  Type materialList
    ! An array of matProps

    Type(matProps), Allocatable :: materials(:)

  End Type materialList

  ! Public interface
  Public :: matProps
  Public :: materialList
  Public :: loadMatProps
  Public :: getCharElemLengths
  Public :: checkForSnapBack


  ! Reference to object for singleton
  Type(matProps) :: m
  Type(materialList), Save :: user


Contains

  Type(matProps) Function loadMatProps(materialName, issueWarnings, nprops, props)
    ! Loads material properties into module variable m

    Use forlog_Mod

    Include 'vaba_param.inc'

    ! Arguments
    Character(len=*), intent(IN) :: materialName
    Logical, intent(IN) :: issueWarnings
    Integer, intent(IN) :: nprops
    Double Precision :: props(nprops)

    ! Locals
    Character(len=256) :: outputDir, jobName, fileName
    Integer :: lenOutputDir, lenJobName
    Logical :: fileExists

    Double Precision :: valueDbl

    Double Precision, Parameter :: zero=0.d0, one=1.d0, two=2.d0

    ! -------------------------------------------------------------------- !

    ! Initializations
    Call initializeMinMaxValues()

    ! Record material name
    m%name = trim(materialName)

    ! Get the output directory (location to search for .props file)
    Call VGETOUTDIR(outputDir, lenOutputDir)

    ! Get the jobname
    Call VGETJOBNAME(jobName, lenJobName)

    ! Look to see if a material properties file exists
    ! First look for: jobName_materialName.props
    fileName = trim(outputDir) // '/' // trim(jobName) // '_' // trim(materialName) // '.props'
    Inquire(FILE=fileName, EXIST=fileExists)
    If (fileExists) Then
      Call loadPropertiesFromFile(trim(fileName), nprops, props)
    Else
      ! Try looking for: materialName.props
      fileName = trim(outputDir) // '/' // trim(materialName) // '.props'
      Inquire(FILE=fileName, EXIST=fileExists)
      If (fileExists) Then
        Call loadPropertiesFromFile(trim(fileName), nprops, props)
      Else
        Call log%info("loadMatProps: Material property file not found. Looking for: " // trim(fileName))
        Call log%info("loadMatProps: Reading material properties from the input deck instead.")
        Call loadPropertiesFromInp(nprops, props)
      End If
    End If

    ! Feature flags
    If (.NOT. m%friction) m%mu = zero

    ! Checks that a consistent set of properties has been defined
    Call consistencyChecks(issueWarnings)

    ! Calculated properties
    m%v21 = m%E2*m%v12/m%E1
    m%v31 = m%E3*m%v13/m%E1
    m%v32 = m%E3*m%v23/m%E2
    m%cte(3) = m%cte(2)

    ! Return a reference to the matProps object
    loadMatProps = m

    Return
  End Function loadMatProps


  Subroutine loadPropertiesFromFile(fileName, nprops, props)
    ! Populates m with the material properties in the specified .props file

    Use forlog_Mod

    Include 'vaba_param.inc'

    !Arguments
    Character(len=*), intent(IN) :: fileName
    Integer, intent(IN) :: nprops
    Double Precision :: props(nprops)

    ! Locals
    Integer, Parameter :: unit=106
    Integer :: iostat
    Character(len=256) :: line, key, value, tmp
    Character(len=30) :: featureFlags
    Integer :: commentTokenPos, equalTokenPos
    Double Precision, Parameter :: zero=0.d0, one=1.d0, two=2.d0
    ! -------------------------------------------------------------------- !

    Call log%debug("loadMatProps: props file name: " // fileName)

    ! Try to load material properties from a '.props' file
    Open(UNIT=unit, FILE=fileName, STATUS='old', ACTION='read', position='rewind', IOSTAT=iostat)
    If (iostat .NE. 0) Call log%error("loadMatProps: Unable to access the .props file")
    ReadLines: Do

      ! Read the next line in the file
      Read(unit,'(A255)',IOSTAT=iostat) line

      If (iostat .GT. 0) Then
        Call log%error("loadMatProps: Unknown error reading file")
        Exit

      Else If (iostat .LT. 0) Then
        Call log%debug("loadMatProps: Reached end of props file")
        Exit

      Else  ! Parse the line in the file
        ! Skip blank lines
        If (Len(Trim(line)) .EQ. 0) Cycle ReadLines

        ! Ignore comment lines (token: //)
        commentTokenPos = Index(line, '//')
        If (commentTokenPos .EQ. 1) Then
          Cycle ReadLines

        Else If (commentTokenPos .GT. 1) Then
          ! This line has a trailing comment, retain everything before the comment
          line = line(1:commentTokenPos-1)
        End If

        ! Split the line into a key and value
        equalTokenPos = Index(line, '=')
        If (equalTokenPos .EQ. 0) Call log%error("loadMatProps: Expected [name] = [value] format not found. Line image: " // trim(line))

        ! Parse key and value from the line
        key = line(1:equalTokenPos-1)
        value = line(equalTokenPos+1:)

        ! Check if the key is a known property. If it's a known property, check validity, and store it
        Select Case (trim(key))

          Case ('E1')
            Call verifyAndSaveProperty_str(trim(key), value, m%modulus_min, m%modulus_max, m%E1, m%E1_def)

          Case ('E2')
            Call verifyAndSaveProperty_str(trim(key), value, m%modulus_min, m%modulus_max, m%E2, m%E2_def)

          Case ('G12')
            Call verifyAndSaveProperty_str(trim(key), value, m%modulus_min, m%modulus_max, m%G12, m%G12_def)

          Case ('v12')
            Call verifyAndSaveProperty_str(trim(key), value, m%poissons_min, m%poissons_max, m%v12, m%v12_def)

          Case ('v23')
            Call verifyAndSaveProperty_str(trim(key), value, m%poissons_min, m%poissons_max, m%v23, m%v23_def)

          Case ('YT')
            Call verifyAndSaveProperty_str(trim(key), value, m%strength_min, m%strength_max, m%YT, m%YT_def)

          Case ('SL')
            Call verifyAndSaveProperty_str(trim(key), value, m%strength_min, m%strength_max, m%SL, m%SL_def)

          Case ('GYT')
            Call verifyAndSaveProperty_str(trim(key), value, m%toughness_min, m%toughness_max, m%GYT, m%GYT_def)

          Case ('GSL')
            Call verifyAndSaveProperty_str(trim(key), value, m%toughness_min, m%toughness_max, m%GSL, m%GSL_def)

          Case ('eta_BK')
            Call verifyAndSaveProperty_str(trim(key), value, m%eta_BK_min, m%eta_BK_max, m%eta_BK, m%eta_BK_def)

          Case ('E3')
            Call verifyAndSaveProperty_str(trim(key), value, m%modulus_min, m%modulus_max, m%E3, m%E3_def)

          Case ('G13')
            Call verifyAndSaveProperty_str(trim(key), value, m%modulus_min, m%modulus_max, m%G13, m%G13_def)

          Case ('G23')
            Call verifyAndSaveProperty_str(trim(key), value, m%modulus_min, m%modulus_max, m%G23, m%G23_def)

          Case ('v13')
            Call verifyAndSaveProperty_str(trim(key), value, m%poissons_min, m%poissons_max, m%v13, m%v13_def)

          Case ('alpha11')
            Call verifyAndSaveProperty_str(trim(key), value, m%cte_min, m%cte_max, m%cte(1), m%cte_def(1))

          Case ('alpha22')
            Call verifyAndSaveProperty_str(trim(key), value, m%cte_min, m%cte_max, m%cte(2), m%cte_def(2))

          Case ('alpha_PL')
            Call verifyAndSaveProperty_str(trim(key), value, m%aPL_min, m%aPL_max, m%aPL, m%aPL_def)

          Case ('n_PL')
            Call verifyAndSaveProperty_str(trim(key), value, m%nPL_min, m%nPL_max, m%nPL, m%nPL_def)

          Case ('XT')
            Call verifyAndSaveProperty_str(trim(key), value, m%strength_min, m%strength_max, m%XT, m%XT_def)

          Case ('fXT')
            Call verifyAndSaveProperty_str(trim(key), value, zero, one, m%fXT, m%fXT_def)

          Case ('GXT')
            Call verifyAndSaveProperty_str(trim(key), value, m%toughness_min, m%toughness_max, m%GXT, m%GXT_def)

          Case ('fGXT')
            Call verifyAndSaveProperty_str(trim(key), value, zero, one, m%fGXT, m%fGXT_def)

          Case ('XC')
            Call verifyAndSaveProperty_str(trim(key), value, m%strength_min, m%strength_max, m%XC, m%XC_def)

          Case ('fXC')
            Call verifyAndSaveProperty_str(trim(key), value, Tiny(zero), one, m%fXC, m%fXC_def)

          Case ('GXC')
            Call verifyAndSaveProperty_str(trim(key), value, m%toughness_min, m%toughness_max, m%GXC, m%GXC_def)

          Case ('fGXC')
            Call verifyAndSaveProperty_str(trim(key), value, zero, one, m%fGXC, m%fGXC_def)

          Case ('YC')
            Call verifyAndSaveProperty_str(trim(key), value, m%strength_min, m%strength_max, m%YC, m%YC_def)

          Case ('alpha0')
            Call verifyAndSaveProperty_str(trim(key), value, m%alpha0_min, m%alpha0_max, m%alpha0, m%alpha0_def)

          Case ('mu')
            Call verifyAndSaveProperty_str(trim(key), value, m%mu_min, m%mu_max, m%mu, m%mu_def)

          Case Default
            Call log%error("loadMatProps: Property not recognized: " // trim(key))
        End Select
      End If
    End Do ReadLines

    ! Close the .props file
    Close(unit)

    ! Require at least the first three material properties to be specified in the inp deck
    If (nprops .LT. 3) Call log%error("loadMatProps: Must define the first three properties in the inp deck")

    ! Load material properties from input deck (at most, first 8)
    Call loadPropertiesFromInp(nprops, props)

    Return
  End Subroutine loadPropertiesFromFile


  Subroutine loadPropertiesFromInp(nprops, props)
    Use forlog_Mod

    Include 'vaba_param.inc'

    !Arguments
    Integer, intent(IN) :: nprops
    Double Precision :: props(nprops)

    ! Locals
    Character(len=30) :: tmp, featureFlags
    Integer :: numMissingLeadingZeros
    Double Precision, Parameter :: zero=0.d0, one=1.d0, two=2.d0
    ! -------------------------------------------------------------------- !

    Do i=1, nprops
      Select Case (i)
        Case (1)
          ! Convert to string (remove decimal and trailing zeros)
          write (tmp, *) props(i)
          featureFlags = Adjustl(tmp(1:Index(tmp, '.')-1))
          numMissingLeadingZeros = 6 - Len_trim(featureFlags)
          If (numMissingLeadingZeros .GT. 0) Then
            Do j=1, numMissingLeadingZeros
              featureFlags = '0' // featureFlags
            End Do
          End If

          ! Parse each position into a feature flag
          ! Position 1: Matrix damage
          ! Position 2: Shear nonlinearity
          ! Position 3: Fiber tensile damage
          ! Position 4: Fiber compression damage
          ! Position 5: *reserved*
          ! Position 6: Friction

          ! Default friction to enabled
          If (Len_trim(featureFlags) .LT. 6) Call log%error("loadMatProps: All 6 feature flags must be specified, only found " // trim(str(Len_trim(featureFlags))))

          ! Set flags
          Do j=1,Len_trim(featureFlags)
            Select Case (j)
              Case (1)
                If (featureFlags(j:j) .EQ. '1') Then
                  m%matrixDam = .TRUE.
                  Call log%info("loadMatProps: Matrix damage is ENABLED")
                Else
                  m%matrixDam = .FALSE.
                  Call log%info("loadMatProps: Matrix damage is DISABLED")
                End If
              Case (2)
                If (featureFlags(j:j) .EQ. '1') Then
                  m%shearNonlinearity = .TRUE.
                  Call log%info("loadMatProps: shear-nonlinearity is ENABLED")
                Else
                  m%shearNonlinearity = .FALSE.
                  Call log%info("loadMatProps: shear-nonlinearity is DISABLED")
                End If
              Case (3)
                If (featureFlags(j:j) .EQ. '1') Then
                  m%fiberTenDam = .TRUE.
                  Call log%info("loadMatProps: fiber tensile damage is ENABLED")
                Else
                  m%fiberTenDam = .FALSE.
                  Call log%info("loadMatProps: fiber tensile damage is DISABLED")
                End If
              Case (4)
                If (featureFlags(j:j) .EQ. '1') Then
                  m%fiberCompDam = .TRUE.
                  Call log%info("loadMatProps: fiber compressive damage is ENABLED")
                Else
                  m%fiberCompDam = .FALSE.
                  Call log%info("loadMatProps: fiber compressive damage is DISABLED")
                End If
              Case (5)
              Case (6)
                If (featureFlags(j:j) .EQ. '1') Then
                  m%friction = .TRUE.
                  Call log%info("loadMatProps: friction is ENABLED")
                Else
                  m%friction = .FALSE.
                  Call log%info("loadMatProps: friction is DISABLED")
                End If
              Case Default
                Call log%error("loadMatProps: Unknown position found in feature flags")
            End Select
          End Do
        Case (2)
          ! m%strainDef = Int(props(i))
          m%strainDef = 2
        Case (3)
          m%thickness = props(i)
        ! props(4:8) are reserved for use in the future as needed
        Case (4)
        Case (5)
        Case (6)
        Case (7)
        Case (8)

        Case (9)
          Call verifyAndSaveProperty_double('E1', props(i), m%modulus_min, m%modulus_max, m%E1, m%E1_def)

        Case (10)
          Call verifyAndSaveProperty_double('E2', props(i), m%modulus_min, m%modulus_max, m%E2, m%E2_def)

        Case (11)
          Call verifyAndSaveProperty_double('G12', props(i), m%modulus_min, m%modulus_max, m%G12, m%G12_def)

        Case (12)
          Call verifyAndSaveProperty_double('v12', props(i), m%poissons_min, m%poissons_max, m%v12, m%v12_def)

        Case (13)
          Call verifyAndSaveProperty_double('v23', props(i), m%poissons_min, m%poissons_max, m%v23, m%v23_def)

        Case (14)
          Call verifyAndSaveProperty_double('YT', props(i), m%strength_min, m%strength_max, m%YT, m%YT_def)

        Case (15)
          Call verifyAndSaveProperty_double('SL', props(i), m%strength_min, m%strength_max, m%SL, m%SL_def)

        Case (16)
          Call verifyAndSaveProperty_double('GYT', props(i), m%toughness_min, m%toughness_max, m%GYT, m%GYT_def)

        Case (17)
          Call verifyAndSaveProperty_double('GSL', props(i), m%toughness_min, m%toughness_max, m%GSL, m%GSL_def)

        Case (18)
          Call verifyAndSaveProperty_double('eta_BK', props(i), m%eta_BK_min, m%eta_BK_max, m%eta_BK, m%eta_BK_def)

        Case (19)
          Call verifyAndSaveProperty_double('YC', props(i), m%strength_min, m%strength_max, m%YC, m%YC_def)

        Case (20)
          Call verifyAndSaveProperty_double('alpha0', props(i), m%alpha0_min, m%alpha0_max, m%alpha0, m%alpha0_def)

        Case (21)
          Call verifyAndSaveProperty_double('E3', props(i), m%modulus_min, m%modulus_max, m%E3, m%E3_def)

        Case (22)
          Call verifyAndSaveProperty_double('G13', props(i), m%modulus_min, m%modulus_max, m%G13, m%G13_def)

        Case (23)
          Call verifyAndSaveProperty_double('G23', props(i), m%modulus_min, m%modulus_max, m%G23, m%G23_def)

        Case (24)
          Call verifyAndSaveProperty_double('v13', props(i), m%poissons_min, m%poissons_max, m%v13, m%v13_def)

        Case (25)
          Call verifyAndSaveProperty_double('alpha11', props(i), m%cte_min, m%cte_max, m%cte(1), m%cte_def(1))

        Case (26)
          Call verifyAndSaveProperty_double('alpha22', props(i), m%cte_min, m%cte_max, m%cte(2), m%cte_def(2))

        Case (27)
          Call verifyAndSaveProperty_double('alpha_PL', props(i), m%aPL_min, m%aPL_max, m%aPL, m%aPL_def)

        Case (28)
          Call verifyAndSaveProperty_double('n_PL', props(i), zero, m%nPL_max, m%nPL, m%nPL_def)

        Case (29)
          Call verifyAndSaveProperty_double('XT', props(i), m%strength_min, m%strength_max, m%XT, m%XT_def)

        Case (30)
          Call verifyAndSaveProperty_double('fXT', props(i), zero, one, m%fXT, m%fXT_def)

        Case (31)
          Call verifyAndSaveProperty_double('GXT', props(i), m%toughness_min, m%toughness_max, m%GXT, m%GXT_def)

        Case (32)
          Call verifyAndSaveProperty_double('fGXT', props(i), zero, one, m%fGXT, m%fGXT_def)

        Case (33)
          Call verifyAndSaveProperty_double('XC', props(i), m%strength_min, m%strength_max, m%XC, m%XC_def)

        Case (34)
          Call verifyAndSaveProperty_double('fXC', props(i), Tiny(zero), one, m%fXC, m%fXC_def)

        Case (35)
          Call verifyAndSaveProperty_double('GXC', props(i), m%toughness_min, m%toughness_max, m%GXC, m%GXC_def)

        Case (36)
          Call verifyAndSaveProperty_double('fGXC', props(i), zero, one, m%fGXC, m%fGXC_def)

        Case (37)
        Case (38)
        Case (39)
        Case (40)
          Call verifyAndSaveProperty_double('mu', props(i), m%mu_min, m%mu_max, m%mu, m%mu_def)

        Case Default
          Call log%error("loadMatProps: Unknown property #"//trim(str(i))//" found. Expecting nprops <= 40")
      End Select

    End Do

    Return
  End Subroutine loadPropertiesFromInp


  Subroutine initializeMinMaxValues()

    ! Locals
    Double Precision, Parameter :: zero=0.d0, one=1.d0, two=2.d0
    Double Precision, parameter :: pi = 3.141592653589793239
    ! -------------------------------------------------------------------- !

    m%modulus_min = Tiny(zero)
    m%modulus_max = Huge(zero)

    m%poissons_min = zero
    m%poissons_max = one

    m%strength_min = Tiny(zero)
    m%strength_max = Huge(zero)

    m%toughness_min = Tiny(zero)
    m%toughness_max = Huge(zero)

    m%eta_BK_min = one
    m%eta_BK_max = Huge(zero)

    m%cte_min = -one
    m%cte_max = one

    m%aPL_min = zero
    m%aPL_max = Huge(zero)

    m%nPL_min = zero
    m%nPL_max = Huge(zero)

    m%alpha0_min = zero
    m%alpha0_max = pi/two

    m%mu_min = zero
    m%mu_max = one


    Return
  End Subroutine initializeMinMaxValues


  Subroutine verifyAndSaveProperty_str(key, value, min, max, saveTo, flag)
    ! Checks if the value is within the specifed bounds. Prints an error message
    ! which kills the analysis if a value is out of bounds. Otherwise, save the property
    ! and sets the flag to indicate that the property has been read in.

    Use forlog_Mod

    !Arguments
    Character(len=*), intent(IN) :: key, value
    Double Precision, intent(IN) :: min, max
    Double Precision, intent(OUT) :: saveTo
    Logical, intent(OUT) :: flag

    ! Locals
    Double Precision :: valueDbl
    ! -------------------------------------------------------------------- !

    ! Convert to double
    Read(value,*) valueDbl

    ! Verify that the value is within the specified bounds
    If (valueDbl .LT. min) Then
      Call log%error(" PROPERTY ERROR " // trim(key) // " cannot be less than " // trim(str(min)) // ". Found value: " // trim(str(valueDbl)))

    Else If (valueDbl .GT. max) Then
      Call log%error(" PROPERTY ERROR " // trim(key) // " cannot be greater than " // trim(str(max)) // ". Found value: " // trim(str(valueDbl)))

    Else
      Call log%debug(" loadMatProps: Loaded " // trim(key) // " = " // trim(str(valueDbl)))
    End If

    ! Save the value and set the flag
    saveTo = valueDbl
    flag = .TRUE.

    Return
  End Subroutine verifyAndSaveProperty_str


  Subroutine verifyAndSaveProperty_double(key, value, min, max, saveTo, flag)
    ! Checks if the value is within the specifed bounds. Prints an error message
    ! which kills the analysis if a value is out of bounds. Otherwise, save the property
    ! and sets the flag to indicate that the property has been read in.

    Use forlog_Mod

    !Arguments
    Character(len=*), intent(IN) :: key
    Double Precision, intent(IN) :: value, min, max
    Double Precision, intent(OUT) :: saveTo
    Logical, intent(OUT) :: flag

    ! Locals
    Double Precision, Parameter :: zero=0.d0
    ! -------------------------------------------------------------------- !

    ! Verify that the value is within the specified bounds
    If (value .LT. min) Then
      If (value .EQ. zero) Then
        flag = .FALSE.
        Call log%info(" loadMatProps: assuming " // trim(key) // " is not defined")
        Return
      Else
        Call log%error(" PROPERTY ERROR " // trim(key) // " cannot be less than " // trim(str(min)) // ". Found value: " // trim(str(value)))
      End If

    Else If (value .GT. max) Then
      Call log%error(" PROPERTY ERROR " // trim(key) // " cannot be greater than " // trim(str(max)) // ". Found value: " // trim(str(value)))

    Else
      Call log%debug(" loadMatProps: Loaded " // trim(key) // " = " // trim(str(value)))
    End If

    ! Save the value and set the flag
    saveTo = value
    flag = .TRUE.

    Return
  End Subroutine verifyAndSaveProperty_double


  Subroutine consistencyChecks(issueWarnings)
    ! Checks that a consistent set of properties has been defined

    Use forlog_Mod

    ! Arguments
    Logical, intent(IN) :: issueWarnings

    ! Locals
    Double Precision, Parameter :: zero=0.d0, one=1.d0, two=2.d0, four=4.d0
    ! -------------------------------------------------------------------- !

    ! Check that all elastic properties have been defined for transverse isotropy
    If (.NOT. m%E1_def) Call log%error('PROPERTY ERROR: Must define a value for E1')
    If (.NOT. m%E2_def) Call log%error('PROPERTY ERROR: Must define a value for E2')
    If (.NOT. m%G12_def) Call log%error('PROPERTY ERROR: Must define a value for G12')
    If (.NOT. m%v12_def) Call log%error('PROPERTY ERROR: Must define a value for v12')
    If (.NOT. m%v23_def) Call log%error('PROPERTY ERROR: Must define a value for v23')

    ! Check if orthotropic elastic properties have been defined
    If (m%E3_def .OR. m%G13_def .OR. m%G23_def) Then
      If (.NOT. m%E3_def) Call log%error('PROPERTY ERROR: Some orthotropic elastic properties are missing. Must define a value for E3.')
      If (.NOT. m%G13_def) Call log%error('PROPERTY ERROR: Some orthotropic elastic properties are missing. Must define a value for G13.')
      If (.NOT. m%G23_def) Call log%error('PROPERTY ERROR: Some orthotropic elastic properties are missing. Must define a value for G23.')
      If (.NOT. m%v13_def) Call log%error('PROPERTY ERROR: Some orthotropic elastic properties are missing. Must define a value for v13.')
      Call log%info('PROPERTY: Orthotropic constants are defined directly')
    Else
      ! Compute these properties assuming transverse isotropy
      Call log%info('PROPERTY: Assuming transverse isotropy')
      m%E3  = m%E2
      m%G13 = m%G12
      m%v13 = m%v12
      m%G23 = m%E2/two/(one + m%v23)
    End If

    ! Check if matrix damage properties have been specified
    If (m%matrixDam) Then
      If (.NOT. m%YT_def) Call log%error('PROPERTY ERROR: Some matrix damage properties are missing. Must define a value for YT.')
      If (.NOT. m%SL_def) Call log%error('PROPERTY ERROR: Some matrix damage properties are missing. Must define a value for SL_def.')
      If (.NOT. m%GYT_def) Call log%error('PROPERTY ERROR: Some matrix damage properties are missing. Must define a value for GYT.')
      If (.NOT. m%GSL_def) Call log%error('PROPERTY ERROR: Some matrix damage properties are missing. Must define a value for GSL.')
      If (.NOT. m%eta_BK_def) Call log%error('PROPERTY ERROR: Some matrix damage properties are missing. Must define a value for eta_BK.')
      If (.NOT. m%YC_def) Call log%error('PROPERTY ERROR: Some matrix damage properties are missing. Must define a value for YC.')
      If (.NOT. m%alpha0_def) Call log%error('PROPERTY ERROR: Some matrix damage properties are missing. Must define a value for alpha0.')

      m%alpha0_deg = NINT(m%alpha0*45.d0/ATAN(one))
      m%alpha0 = alpha0_DGD(m%alpha0, m%E1, m%E2, m%E3, m%G12, m%G13, m%G23, m%v12, m%v13, m%v23, m%YC)
      m%etaL = -m%SL*COS(two*m%alpha0)/(m%YC*COS(m%alpha0)*COS(m%alpha0))
      m%etaT = -one/TAN(two*m%alpha0)
      m%ST   = m%YC*COS(m%alpha0)*(SIN(m%alpha0) + COS(m%alpha0)/TAN(two*m%alpha0))

      Call log%info('PROPERTY: Matrix damage is enabled')
    Else
      Call log%info('PROPERTY: Matrix damage is disabled')
    End If

    ! Check if CTEs have been defined
    If (m%cte_def(1) .OR. m%cte_def(2)) Then
      If (.NOT. m%cte_def(1)) Call log%error('PROPERTY ERROR: Some CTE properties are missing. Must define a value for alpha11.')
      If (.NOT. m%cte_def(2)) Call log%error('PROPERTY ERROR: Some CTE properties are missing. Must define a value for alpha22.')
      Call log%info('PROPERTY: CTEs have been defined')
    Else
      If (issueWarnings) Call log%warn('PROPERTY: CTEs are being set to zero')
      m%cte(1) = zero
      m%cte(2) = zero
    End If

    ! Check if shear nonlinearity properties have been defined
    If (m%shearNonlinearity) Then
      If (.NOT. m%aPL_def) Call log%error('PROPERTY ERROR: Some shear-nonlinearity properties are missing. Must define a value for alpha_PL.')
      If (.NOT. m%nPL_def) Call log%error('PROPERTY ERROR: Some shear-nonlinearity properties are missing. Must define a value for n_PL.')
      Call log%info('PROPERTY: Shear-nonlinearity properties have been defined')
    Else
      Call log%info('PROPERTY: Shear-nonlinearity is disabled')
    End If

    ! Check if fiber tensile damage properties have been defined
    If (m%fiberTenDam) Then
      If (.NOT. m%XT_def) Call log%error('PROPERTY ERROR: Some fiber tensile damage properties are missing. Must define a value for XT.')
      If (.NOT. m%fXT_def) Call log%error('PROPERTY ERROR: Some fiber tensile damage properties are missing. Must define a value for fXT.')
      If (.NOT. m%GXT_def) Call log%error('PROPERTY ERROR: Some fiber tensile damage properties are missing. Must define a value for GXT.')
      If (.NOT. m%fGXT_def) Call log%error('PROPERTY ERROR: Some fiber tensile damage properties are missing. Must define a value for fGXT.')

      Call log%info('PROPERTY: fiber tensile damage properties have been defined')
    Else
      Call log%info('PROPERTY: fiber tensile damage is disabled')
    End If

    ! Check if the fiber compression damage properties have been defined
    If (m%fiberCompDam) Then
      If (.NOT. m%XC_def) Call log%error('PROPERTY ERROR: Some fiber compression damage properties are missing. Must define a value for XC.')
      If (.NOT. m%fXC_def) Call log%error('PROPERTY ERROR: Some fiber compression damage properties are missing. Must define a value for fXC.')
      If (.NOT. m%GXC_def) Call log%error('PROPERTY ERROR: Some fiber compression damage properties are missing. Must define a value for GXC.')
      If (.NOT. m%fGXC_def) Call log%error('PROPERTY ERROR: Some fiber compression damage properties are missing. Must define a value for fGXC.')

      Call log%info('PROPERTY: fiber compression damage properties have been defined')
    Else
      Call log%info('PROPERTY: fiber compression damage is disabled')
    End If

    Return
  End Subroutine consistencyChecks


  Pure Function getCharElemLengths(m, nshr, charLength, lcOld)
    ! Loads the characteristic lengths and checks for snap back

    ! Arguments
    Type(matProps), intent(IN) :: m
    Integer, intent(IN) :: nshr
    Double Precision, intent(IN) :: charLength       ! VUMAT input
    Double Precision, intent(IN) :: lcOld            ! State variable old that stores the characteristic length
    Double Precision :: getCharElemLengths(3)

    ! Locals
    Double Precision :: Lc(3)
    Double Precision, parameter :: zero=0.d0
    ! -------------------------------------------------------------------- !

    ! Load the thickness
    Lc(3) = m%thickness

    ! If the Lc has already been calculated
    If (lcOld .GT. zero) Then
      Lc(1) = lcOld
      Lc(2) = Lc(1)
      getCharElemLengths = Lc
      Return
    End If

    If (nshr .GT. 1) Then
      Lc(1) = SQRT(charLength**3/Lc(3))
    Else
      Lc(1) = charLength
    End If
    Lc(2) = Lc(1)

    getCharElemLengths = Lc

    Return
  End Function getCharElemLengths


  Subroutine checkForSnapBack(m, Lc)
    ! Issues a warning if an element greater than the snap back size

    Use forlog_Mod

    Include 'vaba_param.inc'

    ! Arguments
    Type(matProps), intent(IN) :: m
    Double Precision, intent(IN) :: Lc(3)                    ! Element length

    ! Locals
    Double Precision :: Lc_fT, Lc_fC, Lc_mT, Lc_SL, Lc_YC    ! Maximum element size for each damage mode
    Double Precision, Parameter :: zero=0.d0, one=1.d0, two=2.d0
    ! -------------------------------------------------------------------- !

    ! Fiber Tension
    If (m%fiberTenDam) Then
      If (two*m%GXT*m%E1 .LT. Lc(1)*m%XT**2) Then
        Lc_fT = two*m%GXT*m%E1/m%XT**2
        Call log%warn("Element size > fiber tension snap-back threshold. Set element size < " // trim(str(Lc_fT)))
      End If

      If (two*m%GXT*m%fGXT*m%E1 .LT. Lc(1)*(m%XT*m%fXT)**2) Then
        Lc_fT = two*m%GXT*m%fGXT*m%E1/(m%XT*m%fXT)**2
        Call log%warn("Snap-back in 1st part of fiber tensile softening law. Adjust bilinear ratios or decrease the element size.")
      End If

      If (two*m%GXT*(one - m%fGXT)*m%E1 .LT. Lc(1)*(m%XT*(one - m%fXT))**2) Then
        Lc_fT = two*m%GXT*(one - m%fGXT)*m%E1/(m%XT*(one - m%fXT))**2
        Call log%warn("Snap-back in 2nd part of fiber tensile softening law. Adjust bilinear ratios or decrease the element size.")
      End If
    End If

    ! Fiber Compression
    If (m%fiberCompDam) Then
      If (two*m%GXC*m%E1 .LT. Lc(1)*m%XC**2) Then
        Lc_fC = two*m%GXC*m%E1/m%XC**2
        Call log%warn("Element size > fiber compression snap-back threshold 1. Set element size < " // trim(str(Lc_fC)))
      End If

      If (two*m%GXC*m%fGXC*m%E1 .LT. Lc(1)*(m%XC*m%fXC)**2) Then
        Lc_fT = two*m%GXC*m%fGXC*m%E1/(m%XC*m%fXC)**2
        Call log%warn("Snap-back in 1st part of fiber compression softening law. Adjust bilinear ratios or decrease the element size.")
      End If

      If (two*m%GXC*(one - m%fGXC)*m%E1 .LT. Lc(1)*(m%XC*(one - m%fXC))**2) Then
        Lc_fT = two*m%GXC*(one - m%fGXC)*m%E1/(m%XC*(one - m%fXC))**2
        Call log%warn("Snap-back in 2nd part of fiber compression softening law. Adjust bilinear ratios or decrease the element size.")
      End If
    End If

    If (m%matrixDam) Then
      ! Matrix Tension
      If (two*m%GYT*m%E2 .LT. Lc(2)*m%YT**2) Then
        Lc_mT = two*m%GYT*m%E2/m%YT/m%YT
        Call log%warn("Element size > matrix mode I snap-back threshold. Set element size < " // trim(str(Lc_mT)))
      End If

      ! Matrix Shear
      If (two*m%GSL*m%G12 .LT. Lc(2)*m%SL**2) Then
        Lc_SL = two*m%GSL*m%G12/m%SL**2
        Call log%warn("Element size > matrix mode II snap-back threshold. Set element size < " // trim(str(Lc_SL)))
      End If

      ! Matrix Compression
      If (two*m%GSL*m%E2 .LT. Lc(2)*m%YC**2*COS(m%alpha0)) Then
        Lc_YC = two*m%GSL*m%E2/(m%YC**2*COS(m%alpha0))
        Call log%warn("Element size > matrix compression snap-back threshold. Set element size < " // trim(str(Lc_YC)))
      End If
    End If

    Return
  End Subroutine checkForSnapBack


  Function alpha0_DGD(alpha0, E1, E2, E3, G12, G13, G23, v12, v13, v23, Yc)
    ! Determines the orientation of the angle alpha0 when subject to sigma22 = -Yc

    Use forlog_Mod
    Use matrixAlgUtil_Mod
    Use CLT_Mod

    Include 'vaba_param.inc'

    ! Arguments
    Double Precision, Intent(IN) :: alpha0, E1, E2, E3, G12, G13, G23, v12, v13, v23, Yc

    ! Locals
    Double Precision :: C(6,6)         ! 3-D Stiffness
    Double Precision :: F(3)           ! Represents diagonal of deformation gradient tensor
    Double Precision :: Rs(3)          ! Residual vector
    Double Precision :: tolerance
    Double Precision :: err
    Double Precision :: Jac(3,3)       ! Jacobian
    Integer :: counter
    Double Precision, parameter :: zero=0.d0, one=1.d0, two=2.d0
    ! -------------------------------------------------------------------- !

    Call log%debug('Start of Function alpha0_DGD')

    tolerance = 1.d-4  ! alphaLoop tolerance

    ! Build the stiffness matrix
    C = StiffFunc(6, E1, E2, E3, G12, G13, G23, v12, v13, v23, zero, zero, zero)

    ! Make an initial guess
    F = (/ one, one, one /)

  ! -------------------------------------------------------------------- !
  !    alphaLoop Loop and solution controls definition                   !
  ! -------------------------------------------------------------------- !
    counter = 0  ! Counter for alphaLoop
    counter_max = 100

    alphaLoop: Do  ! Loop to determine the F which corresponds to sigma22 = -Yc
      counter = counter + 1
      ! -------------------------------------------------------------------- !
      !    Define the stress residual vector, R. R is equal to the           !
      !    difference in stress between the cohesive interface and the bulk  !
      !    stress projected onto the cohesive interface                      !
      ! -------------------------------------------------------------------- !
      Rs(1) = (C(1,1)*(F(1)*F(1) - one) + C(1,2)*(F(2)*F(2) - one) + C(1,3)*(F(3)*F(3) - one))/two
      Rs(2) = (C(2,1)*(F(1)*F(1) - one) + C(2,2)*(F(2)*F(2) - one) + C(2,3)*(F(3)*F(3) - one))/two + Yc*F(1)*F(3)/F(2)
      Rs(3) = (C(3,1)*(F(1)*F(1) - one) + C(3,2)*(F(2)*F(2) - one) + C(3,3)*(F(3)*F(3) - one))/two

      ! Check for convergence
      err = Length(Rs)

      ! If converged,
      If (err .LT. tolerance) Then
        alpha0_DGD = ATAN(F(2)/F(3)*TAN(alpha0))
        EXIT alphaLoop
      End If
      IF (counter .EQ. counter_max) Call log%error('Function alpha0_DGD failed to converge')

      ! Define the Jacobian matrix, J
      Jac = zero

      Jac(1,1) = C(1,1)*F(1)
      Jac(1,2) = C(1,2)*F(2)
      Jac(1,3) = C(1,3)*F(3)

      Jac(2,1) = C(2,1)*F(1) + two*Yc*F(3)/F(2)
      Jac(2,2) = C(2,2)*F(1) - two*Yc*F(3)*F(1)/(F(2)*F(2))
      Jac(2,3) = C(2,3)*F(1) + two*Yc*F(1)/F(2)

      Jac(3,1) = C(3,1)*F(1)
      Jac(3,2) = C(3,2)*F(2)
      Jac(3,3) = C(3,3)*F(3)

      ! Calculate the new diagonal deformation gradient
      F = F - MATMUL(MInverse(Jac), Rs)

    End Do alphaLoop

    Return
  End Function alpha0_DGD

End Module matProp_Mod
