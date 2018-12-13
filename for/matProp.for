#define type(x) TYPE(x), target

Module matProp_Mod
  ! Module for loading and validating material properties

  Private

  Type matProps
    ! Stores the set of properties

    Character(len=80) :: name ! Material name, required

    ! Material properties, grouped by feature
    ! (Required inputs)
    Double Precision :: E1, E2, G12, v12, v23                            ! Req. Transversely isotropic.
    ! (Optional inputs)
    Double Precision :: YT, SL, GYT, GSL, eta_BK, YC, alpha0             ! Opt. Matrix failure
    Double Precision :: E3, G13, G23, v13                                ! Opt. Orthotropic. Defaults are calculated based on transverse isotropy
    Double Precision :: cte(3)                                           ! Opt. CTEs (Required if a temperature change is applied in the analysis). Defaults to zero
    Double Precision :: aPL, nPL                                         ! Opt. Shear nonlinearity
    Double Precision :: XT, fXT, GXT, fGXT                               ! Opt. Longitudinal tensile damage (max strain failure criterion, bilinear softening)
    Double Precision :: XC, fXC, GXC, fGXC                               ! Opt. Longitudinal compressive damage (max strain failure criterion, bilinear softening)
    Double Precision :: w_kb                                             ! Opt. LaRC15 kink band model (Also requires XC, shear nonlinearity).
    Double Precision :: cl                                           ! Opt. fiber nonlinearity (defaults to zero, which is no fiber nonlinearity)
    Double Precision :: mu                                               ! Opt. Friction. Defaults to zero
    Double Precision :: es(4), gs(4)
    Double Precision :: schaefer_a6
    Double Precision :: schaefer_b2
    Double Precision :: schaefer_n
    Double Precision :: schaefer_A                                       ! Specified parameters for schaefer theory

    ! min and max values for acceptable range
    Double Precision, private :: modulus_min, modulus_max
    Double Precision, private :: poissons_min, poissons_max
    Double Precision, private :: strength_min, strength_max
    Double Precision, private :: toughness_min, toughness_max
    Double Precision, private :: eta_BK_min, eta_BK_max
    Double Precision, private :: cte_min, cte_max
    Double Precision, private :: aPL_min, aPL_max, nPL_min, nPL_max
    Double Precision, private :: w_kb_min, w_kb_max
    Double Precision, private :: cl_min, cl_max
    Double Precision, private :: alpha0_min, alpha0_max
    Double Precision, private :: mu_min, mu_max
    Double Precision, private :: schapery_min, schapery_max
    Double Precision, private :: schaefer_min, schaefer_max

    ! Flags that indicate if values have been set
    Logical, private :: E1_def, E2_def, G12_def, v12_def, v23_def
    Logical, private :: YT_def, SL_def, GYT_def, GSL_def, eta_BK_def
    Logical, private :: E3_def, G13_def, G23_def, v13_def
    Logical, private :: cte_def(3)
    Logical, private :: aPL_def, nPL_def
    Logical, private :: XT_def, fXT_def, GXT_def, fGXT_def
    Logical, private :: XC_def, fXC_def, GXC_def, fGXC_def
    Logical, private :: YC_def, w_kb_def, alpha0_def
    Logical, private :: cl_def
    Logical, private :: mu_def
    Logical, private :: es_def(4), gs_def(4)
    Logical, private :: schaefer_a6_def
    Logical, private :: schaefer_b2_def
    Logical, private :: schaefer_n_def
    Logical, private :: schaefer_A_def

    ! Calculated properties
    Double Precision :: v21, v31, v32
    Double Precision :: etaL, etaT, ST, phic, alpha0_deg

    Double Precision :: thickness

    ! Flags for features
    Logical :: matrixDam
    Logical :: shearNonlinearity12
    Logical :: shearNonlinearity13
    Logical :: schapery
    Logical :: schaefer
    Logical :: fiberTenDam
    Logical :: fiberCompDamBL
    Logical :: fiberCompDamFKT
    Logical :: friction

  End Type matProps

  Type materialList
    ! An array of matProps

    Type(matProps), Allocatable :: materials(:)

  End Type materialList

  ! Public interface
  Public :: matProps
  Public :: materialList
  Public :: loadMatProps
  Public :: checkForSnapBack
  Public :: initializePhi0
  Public :: consistencyChecks


  ! Reference to object for singleton
  type(matProps) :: m
  type(materialList), Save :: user


Contains

  Type(matProps) Function loadMatProps(materialName, nprops, props)
    ! Loads material properties into module variable m

    Use forlog_Mod

    ! Arguments
    Character(len=*), intent(IN) :: materialName
    ! Logical, intent(IN) :: issueWarnings
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

#ifndef PYEXT
    ! Get the output directory (location to search for .props file)
    Call VGETOUTDIR(outputDir, lenOutputDir)

    ! Get the jobname
    Call VGETJOBNAME(jobName, lenJobName)

    ! Look to see if a material properties file exists
    ! First look for: jobName_materialName.props
    fileName = trim(outputDir) // '/' // trim(jobName) // '_' // trim(materialName) // '.props'
    Inquire(FILE=fileName, EXIST=fileExists)
    ! Try looking for: materialName.props
    If (.NOT. fileExists) Then
      fileName = trim(outputDir) // '/' // trim(materialName) // '.props'
      Inquire(FILE=fileName, EXIST=fileExists)
    End If
#else
    fileName = trim(materialName) // '.props'
    Inquire(FILE=fileName, EXIST=fileExists)
#endif

    ! Handle cases
    If (fileExists .AND. nprops > 8) Then
      Call log%warn("loadMatProps: Found material property file " // trim(fileName) // " and properties in input deck; using properties in the input deck.")
      Call loadPropertiesFromInp(nprops, props)
    Else If (fileExists) Then
      Call log%info("loadMatProps: Using material property file " // trim(fileName))
      Call loadPropertiesFromFile(trim(fileName), nprops, props)
    Else
      Call loadPropertiesFromInp(nprops, props)
      Call log%info("loadMatProps: Material property file not found. Looking for: " // trim(fileName))
      Call log%info("loadMatProps: Reading material properties from the input deck instead.")
    End If

    ! Feature flags
    If (.NOT. m%friction) m%mu = zero

    ! ! Checks that a consistent set of properties has been defined
    ! Call consistencyChecks(issueWarnings)

    ! Calculated properties
    m%v21 = m%E2*m%v12/m%E1
    m%v31 = m%E3*m%v13/m%E1
    m%v32 = m%E3*m%v23/m%E2
    m%cte(3) = m%cte(2)

    ! TODO - additional admissibility checks

    ! Return a reference to the matProps object
    loadMatProps = m

    Return
  End Function loadMatProps


  Subroutine loadPropertiesFromFile(fileName, nprops, props)
    ! Populates m with the material properties in the specified .props file

    Use forlog_Mod

    !Arguments
    Character(len=*), intent(IN) :: fileName
    Integer, intent(IN) :: nprops
    Double Precision :: props(nprops)

    ! Locals
    Integer, Parameter :: unit=108
    Integer :: iostat
    Character(len=256) :: line, key, value, tmp
    Character(len=30) :: featureFlags
    Integer :: commentTokenPos, equalTokenPos
    Double Precision, Parameter :: zero=0.d0, one=1.d0, two=2.d0
    ! -------------------------------------------------------------------- !

    Call log%debug("loadMatProps: props file name: " // fileName)

    ! Try to load material properties from a '.props' file
    Open(UNIT=unit, FILE=fileName, STATUS='old', ACTION='read', position='rewind', IOSTAT=iostat)

    If (iostat /= 0) Call log%error("loadMatProps: Unable to access the .props file")
    ReadLines: Do

      ! Read the next line in the file
      Read(unit,'(A255)',IOSTAT=iostat) line

      If (iostat > 0) Then
        Call log%error("loadMatProps: Unknown error reading file")
        Return

      Else If (iostat < 0) Then
        Call log%debug("loadMatProps: Reached end of props file")
        Exit ReadLines

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

        If (equalTokenPos == 0) Call log%error("loadMatProps: Expected [name] = [value] format not found. Line image: " // trim(line))

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

          Case ('w_kb')
            Call verifyAndSaveProperty_str(trim(key), value, m%w_kb_min, m%w_kb_max, m%w_kb, m%w_kb_def)

          Case ('cl')
            Call verifyAndSaveProperty_str(trim(key), value, m%cl_min, m%cl_max, m%cl, m%cl_def)

          Case ('alpha0')
            Call verifyAndSaveProperty_str(trim(key), value, m%alpha0_min, m%alpha0_max, m%alpha0, m%alpha0_def)

          Case ('mu')
            Call verifyAndSaveProperty_str(trim(key), value, m%mu_min, m%mu_max, m%mu, m%mu_def)

          Case ('es0')
            Call verifyAndSaveProperty_str(trim(key), value, m%schapery_min, m%schapery_max, m%es(1), m%es_def(1))
          Case ('es1')
            Call verifyAndSaveProperty_str(trim(key), value, m%schapery_min, m%schapery_max, m%es(2), m%es_def(2))
          Case ('es2')
            Call verifyAndSaveProperty_str(trim(key), value, m%schapery_min, m%schapery_max, m%es(3), m%es_def(3))
          Case ('es3')
            Call verifyAndSaveProperty_str(trim(key), value, m%schapery_min, m%schapery_max, m%es(4), m%es_def(4))

          Case ('gs0')
            Call verifyAndSaveProperty_str(trim(key), value, m%schapery_min, m%schapery_max, m%gs(1), m%gs_def(1))
          Case ('gs1')
            Call verifyAndSaveProperty_str(trim(key), value, m%schapery_min, m%schapery_max, m%gs(2), m%gs_def(2))
          Case ('gs2')
            Call verifyAndSaveProperty_str(trim(key), value, m%schapery_min, m%schapery_max, m%gs(3), m%gs_def(3))
          Case ('gs3')
            Call verifyAndSaveProperty_str(trim(key), value, m%schapery_min, m%schapery_max, m%gs(4), m%gs_def(4))
          ! Schaefer Theory
          Case ('schaefer_a6')
            Call verifyAndSaveProperty_str(trim(key), value, m%schaefer_min, m%schaefer_max, m%schaefer_a6, m%schaefer_a6_def)
          Case ('schaefer_b2')
            Call verifyAndSaveProperty_str(trim(key), value, m%schaefer_min, m%schaefer_max, m%schaefer_b2, m%schaefer_b2_def)
          Case ('schaefer_n')
            Call verifyAndSaveProperty_str(trim(key), value, m%schaefer_min, m%schaefer_max, m%schaefer_n, m%schaefer_n_def)
          Case ('schaefer_A')
            Call verifyAndSaveProperty_str(trim(key), value, m%schaefer_min, m%schaefer_max, m%schaefer_A, m%schaefer_A_def)
          Case Default
            Call log%error("loadMatProps: Property not recognized: " // trim(key))
        End Select
      End If
    End Do ReadLines

    ! Close the .props file
    Close(unit)

    ! Require at least the first three material properties to be specified in the inp deck
    If (nprops < 3) Call log%error("loadMatProps: Must define the first three properties in the inp deck")

    ! Load material properties from input deck (at most, first 8)
    If (nprops > 8) Then
      Call loadPropertiesFromInp(8, props)
    Else
      Call loadPropertiesFromInp(nprops, props)
    End If

    Return
  End Subroutine loadPropertiesFromFile


  Subroutine loadPropertiesFromInp(nprops, props)
    Use forlog_Mod

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
          If (numMissingLeadingZeros > 0) Then
            Do j=1, numMissingLeadingZeros
              featureFlags = '0' // featureFlags
            End Do
          End If

          ! Parse each position into a feature flag
          ! Position 1: Matrix damage
          ! Position 2: Peak-peak nonlinearity (shear plasticity or schapery micro-damage)
          ! Position 3: Fiber tensile damage
          ! Position 4: Fiber compression damage, CDM [Mutually exclusive with Position 5]
          ! Position 5: Fiber compression damage, DGD [Mutually exclusive with Position 4]
          ! Position 6: Friction

          ! Default friction to enabled
          If (Len_trim(featureFlags) < 6) Call log%error("loadMatProps: All 6 feature flags must be specified, only found " // trim(str(Len_trim(featureFlags))))

          ! Set flags
          Do j=1,Len_trim(featureFlags)
            Select Case (j)
              Case (1)
                If (featureFlags(j:j) == '1') Then
                  m%matrixDam = .TRUE.
                  Call log%info("loadMatProps: Matrix damage ENABLED")
                Else
                  m%matrixDam = .FALSE.
                  Call log%info("loadMatProps: Matrix damage DISABLED")
                End If

              Case (2)
                If (featureFlags(j:j) == '1') Then
                  m%shearNonlinearity12 = .TRUE.
                  m%shearNonlinearity13 = .FALSE.
                  m%schapery = .FALSE.
                  m%schaefer = .FALSE.
                  Call log%info("loadMatProps: Shear nonlinearity ENABLED")
                Else If (featureFlags(j:j) == '2') Then
                  m%shearNonlinearity12 = .FALSE.
                  m%shearNonlinearity13 = .FALSE.
                  m%schapery = .TRUE.
                  m%schaefer = .FALSE.
                  Call log%info("loadMatProps: Schapery micro-damage ENABLED")
                Else If (featureFlags(j:j) == '3') Then
                  m%shearNonlinearity12 = .TRUE.
                  m%shearNonlinearity13 = .TRUE.
                  m%schapery = .FALSE.
                  Call log%info("loadMatProps: Shear nonlinearity (3-D) ENABLED")
                Else If (featureFlags(j:j) == '4') Then
                  m%shearNonlinearity12 = .FALSE.
                  m%shearNonlinearity13 = .TRUE.
                  m%schapery = .FALSE.
                  Call log%info("loadMatProps: Shear nonlinearity (1-3 plane) ENABLED")
                Else If (featureFlags(j:j) == '5') THEN 
                  m%shearNonlinearity12 = .FALSE.
                  m%shearNonlinearity13 = .FALSE.
                  m%schapery = .FALSE.
                  m%schaefer = .TRUE.		
                  Call log%info("loadMatProps: Schaefer non-linearity ENABLED")
                Else
                  m%shearNonlinearity12 = .FALSE.
                  m%shearNonlinearity13 = .FALSE.
                  m%schapery = .FALSE.
                  m%schaefer = .FALSE.		
                  Call log%info("loadMatProps: pre-peak nonlinearity DISABLED")
                End If

              Case (3)
                If (featureFlags(j:j) == '1') Then
                  m%fiberTenDam = .TRUE.
                  Call log%info("loadMatProps: fiber tensile damage ENABLED")
                Else
                  m%fiberTenDam = .FALSE.
                  Call log%info("loadMatProps: fiber tensile damage DISABLED")
                End If

              Case (4)
                If (featureFlags(j:j) == '1') Then
                  m%fiberCompDamBL = .TRUE.
                  m%fiberCompDamFKT = .FALSE.
                  Call log%info("loadMatProps: fiber comp. damage (CDM) ENABLED")
                Else If (featureFlags(j:j) == '2') Then
                  Call log%error("loadMatProps: fiber comp. model 2 not implemented. TODO")
                Else If (featureFlags(j:j) == '3') Then
                  m%fiberCompDamBL = .FALSE.
                  m%fiberCompDamFKT = .TRUE.
                  Call log%info("loadMatProps: fiber comp. damage (DGD) ENABLED")
                Else
                  m%fiberCompDamBL = .FALSE.
                  m%fiberCompDamFKT = .FALSE.
                  Call log%info("loadMatProps: fiber comp. damage DISABLED")
                End If

              Case (5)
                If (featureFlags(j:j) == '1') Then
                  Call log%info("loadMatProps: feature flag position 5 is not used.")
                End If

              Case (6)
                If (featureFlags(j:j) == '1') Then
                  m%friction = .TRUE.
                  Call log%info("loadMatProps: friction ENABLED")
                Else
                  m%friction = .FALSE.
                  Call log%info("loadMatProps: friction DISABLED")
                End If

              Case Default
                Call log%error("loadMatProps: Unknown position found in feature flags")
            End Select
          End Do
        Case (2)
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
          Call verifyAndSaveProperty_double('cl', props(i), m%cl_min, m%cl_max, m%cl, m%cl_def)

        Case (38)
          Call verifyAndSaveProperty_double('w_kb', props(i), m%w_kb_min, m%w_kb_max, m%w_kb, m%w_kb_def)

        Case (39)

        Case (40)
          Call verifyAndSaveProperty_double('mu', props(i), m%mu_min, m%mu_max, m%mu, m%mu_def)

        ! TODO: add schapery micro-damage input properties to the input file parser -- currently only works with
        ! external properties file.

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

    m%w_kb_min = Tiny(zero)
    m%w_kb_max = Huge(zero)

    m%cl_min = zero
    m%cl_max = 1000.d0

    m%alpha0_min = zero
    m%alpha0_max = pi/two

    m%mu_min = zero
    m%mu_max = one

    m%schapery_min = -Huge(zero)
    m%schapery_max = Huge(zero)

    m%schaefer_min = -Huge(zero)
    m%schaefer_max = Huge(zero)

    Return
  End Subroutine initializeMinMaxValues


  Subroutine verifyAndSaveProperty_str(key, value, min, max, saveTo, flag)
    ! Checks if the value is within the specified bounds. Prints an error message
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
    If (valueDbl < min) Then
      Call log%error(" PROPERTY ERROR " // trim(key) // " cannot be less than " // trim(str(min)) // ". Found value: " // trim(str(valueDbl)))

    Else If (valueDbl > max) Then
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
    ! Checks if the value is within the specified bounds. Prints an error message
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
    If (value < min) Then
      If (value == zero) Then
        flag = .FALSE.
        Call log%info(" loadMatProps: assuming " // trim(key) // " is not defined")
        Return
      Else
        Call log%error(" PROPERTY ERROR " // trim(key) // " cannot be less than " // trim(str(min)) // ". Found value: " // trim(str(value)))
      End If

    Else If (value > max) Then
      Call log%error(" PROPERTY ERROR " // trim(key) // " cannot be greater than " // trim(str(max)) // ". Found value: " // trim(str(value)))

    Else
      Call log%debug(" loadMatProps: Loaded " // trim(key) // " = " // trim(str(value)))
    End If

    ! Save the value and set the flag
    saveTo = value
    flag = .TRUE.

    Return
  End Subroutine verifyAndSaveProperty_double


  Subroutine consistencyChecks(m, issueWarnings)
    ! Checks that a consistent set of properties has been defined

    Use forlog_Mod

    ! Arguments
    Type(matProps), intent(INOUT) :: m
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
    If (m%shearNonlinearity12 .OR. m%shearNonlinearity13) Then
      If (.NOT. m%aPL_def) Call log%error('PROPERTY ERROR: Some shear-nonlinearity properties are missing. Must define a value for alpha_PL.')
      If (.NOT. m%nPL_def) Call log%error('PROPERTY ERROR: Some shear-nonlinearity properties are missing. Must define a value for n_PL.')
      Call log%info('PROPERTY: Shear-nonlinearity properties have been defined')
    Else
      Call log%info('PROPERTY: Shear-nonlinearity is disabled')
    End If

    ! Check if all Schapery micro-damage material properties have been defined
    If (m%schapery) Then
      If (.NOT. (all(m%es_def) .AND. all(m%gs_def))) Call log%error('PROPERTY ERROR: Some Schapery micro-damage are missing.')
      Call log%info('PROPERTY: Schapery micro-damage properties have been defined')
    Else
      Call log%info('PROPERTY: Schapery micro-damage is disabled')
      ! Default values that will cause no pre-peak nonlinearity
      m%es = zero
      m%es(1) = one
      m%gs = zero
      m%gs(1) = one
    End If

    ! Check is all the parameters for Schaefer theory have been defined
    If (m%schaefer) Then
      If (.NOT. m%schaefer_a6_def ) Call log%error('PROPERTY ERROR: Some schaefer theory properties are missing. Must define schaefer_a6')
      If (.NOT. m%schaefer_b2_def ) Call log%error('PROPERTY ERROR: Some schaefer theory properties are missing. Must define schaefer_b2')
      If (.NOT. m%schaefer_n_def ) Call log%error('PROPERTY ERROR: Some schaefer theory properties are missing. Must define schaefer_n')
      If (.NOT. m%schaefer_A_def ) Call log%error('PROPERTY ERROR: Some schaefer theory properties are missing. Must define schaefer_A')
      Call log%info('PROPERTY: Schaefer properties have been defined')
    Else
      Call log%info('PROPERTY: Schaefer is disabled')
      ! Default values that will cause no pre-peak nonlinearity
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

    ! Check if the CDM fiber compression damage properties have been defined
    If (m%fiberCompDamBL) Then
      If (.NOT. m%XC_def) Call log%error('PROPERTY ERROR: Some fiber compression damage properties are missing (BL model). Must define a value for XC.')
      If (.NOT. m%fXC_def) Call log%error('PROPERTY ERROR: Some fiber compression damage properties are missing (BL model). Must define a value for fXC.')
      If (.NOT. m%GXC_def) Call log%error('PROPERTY ERROR: Some fiber compression damage properties are missing (BL model). Must define a value for GXC.')
      If (.NOT. m%fGXC_def) Call log%error('PROPERTY ERROR: Some fiber compression damage properties are missing (BL model). Must define a value for fGXC.')

      Call log%info('PROPERTY: fiber compression damage BL (model 1) properties have been defined')
    Else
      Call log%info('PROPERTY: fiber compression damage BL (model 1) is disabled')
    End If

    ! Check if the DGD fiber compression damage properties have been defined
    If (m%fiberCompDamFKT) Then
      If (.NOT. m%XC_def) Call log%error('PROPERTY ERROR: Some fiber compression damage properties are missing (FKT model). Must define a value for XC.')
      If (.NOT. m%YC_def) Call log%error('PROPERTY ERROR: Some fiber compression damage properties are missing (FKT model). Must define a value for YC.')
      If (.NOT. m%w_kb_def) Call log%error('PROPERTY ERROR: Some fiber compression damage properties are missing (FKT model). Must define a value for w_kb.')
      If (.NOT. m%alpha0_def) Call log%error('PROPERTY ERROR: Some fiber compression damage properties are missing (FKT model). Must define a value for alpha0.')
      If (.NOT. m%SL_def) Call log%error('PROPERTY ERROR: Some fiber compression damage properties are missing (FKT model). Must define a value for SL.')

      ! Make sure shear nonlinearity is enabled
      If (.NOT. m%shearNonlinearity12) Call log%error('PROPERTY ERROR: In-plane shear-nonlinearity must be enabled with fiber compression damage FKT.')

      Call log%info('PROPERTY: fiber compression damage FKT (model 3) properties have been defined')
    Else
      Call log%info('PROPERTY: fiber compression damage FKT (model 3) is disabled')
    End If

    ! check if fiber nonlinearity has been defined
    If (m%cl_def) Then
      Call log%info('PROPERTY: fiber nonlinearity has been defined')
    Else
      m%cl = zero
    End If

    Return
  End Subroutine consistencyChecks


  Subroutine checkForSnapBack(m, Lc, elementNumber)
    ! Issues a warning if an element greater than the snap back size

    Use forlog_Mod

    ! Arguments
    Type(matProps), intent(IN) :: m
    Double Precision, intent(IN) :: Lc(3)                    ! Element length
    Integer, intent(IN) :: elementNumber                     ! Element number

    ! Locals
    Double Precision :: Lc_fT, Lc_fC, Lc_mT, Lc_SL, Lc_YC    ! Maximum element size for each damage mode
    Double Precision, Parameter :: zero=0.d0, one=1.d0, two=2.d0
    ! -------------------------------------------------------------------- !

    ! Fiber Tension
    If (m%fiberTenDam) Then
      If (two*m%GXT*m%E1 < Lc(1)*m%XT**2) Then
        Lc_fT = two*m%GXT*m%E1/m%XT**2
        Call log%warn("Element " // trim(str(elementNumber)) // " has size " // trim(str(Lc(1))) // " which is > fiber tension snap-back threshold. Set element size < " // trim(str(Lc_fT)))
      End If

      If (two*m%GXT*m%fGXT*m%E1 < Lc(1)*(m%XT*m%fXT)**2) Then
        Lc_fT = two*m%GXT*m%fGXT*m%E1/(m%XT*m%fXT)**2
        Call log%warn("Snap-back in 1st part of fiber tensile softening law. Adjust bilinear ratios or decrease the element size.")
      End If

      If (two*m%GXT*(one - m%fGXT)*m%E1 < Lc(1)*(m%XT*(one - m%fXT))**2) Then
        Lc_fT = two*m%GXT*(one - m%fGXT)*m%E1/(m%XT*(one - m%fXT))**2
        Call log%warn("Snap-back in 2nd part of fiber tensile softening law. Adjust bilinear ratios or decrease the element size.")
      End If
    End If

    ! Fiber Compression
    If (m%fiberCompDamBL) Then
      If (two*m%GXC*m%E1 < Lc(1)*m%XC**2) Then
        Lc_fC = two*m%GXC*m%E1/m%XC**2
        Call log%warn("Element " // trim(str(elementNumber)) // " has size " // trim(str(Lc(1))) // " which is > fiber compression snap-back threshold 1. Set element size < " // trim(str(Lc_fC)))
      End If

      If (two*m%GXC*m%fGXC*m%E1 < Lc(1)*(m%XC*m%fXC)**2) Then
        Lc_fT = two*m%GXC*m%fGXC*m%E1/(m%XC*m%fXC)**2
        Call log%warn("Snap-back in 1st part of fiber compression softening law. Adjust bilinear ratios or decrease the element size.")
      End If

      If (two*m%GXC*(one - m%fGXC)*m%E1 < Lc(1)*(m%XC*(one - m%fXC))**2) Then
        Lc_fT = two*m%GXC*(one - m%fGXC)*m%E1/(m%XC*(one - m%fXC))**2
        Call log%warn("Snap-back in 2nd part of fiber compression softening law. Adjust bilinear ratios or decrease the element size.")
      End If
    End If

    If (m%matrixDam) Then
      ! Matrix Tension
      If (two*m%GYT*m%E2 < Lc(2)*m%YT**2) Then
        Lc_mT = two*m%GYT*m%E2/m%YT/m%YT
        Call log%warn("Element " // trim(str(elementNumber)) // " has size " // trim(str(Lc(2))) // " which is > matrix mode I snap-back threshold. Set element size < " // trim(str(Lc_mT)))
      End If

      ! Matrix Shear
      If (two*m%GSL*m%G12 < Lc(2)*m%SL**2) Then
        Lc_SL = two*m%GSL*m%G12/m%SL**2
        Call log%warn("Element " // trim(str(elementNumber)) // " has size " // trim(str(Lc(2))) // " which is > matrix mode II snap-back threshold. Set element size < " // trim(str(Lc_SL)))
      End If

      ! Matrix Compression
      If (two*m%GSL*m%E2 < Lc(2)*m%YC**2*COS(m%alpha0)) Then
        Lc_YC = two*m%GSL*m%E2/(m%YC**2*COS(m%alpha0))
        Call log%warn("Element " // trim(str(elementNumber)) // " has size " // trim(str(Lc(2))) // " which is > matrix compression snap-back threshold. Set element size < " // trim(str(Lc_YC)))
      End If
    End If

    Return
  End Subroutine checkForSnapBack


  Function initializePhi0(phi0, G12, XC, aPL, nPL, Lc, position)
    ! Determines the phi0 for shear instability

    Use forlog_Mod

    ! Arguments
    Double Precision, intent(IN) :: phi0                         ! Guess for phi0
    Double Precision, intent(IN) :: G12, XC                      ! Material properties
    Double Precision, intent(IN) :: aPL, nPL                     ! Shear nonlinearity
    Double Precision, intent(IN) :: Lc(3)                        ! Element characteristic lengths
    Double Precision, intent(IN) :: position(3)                  ! Position
    Double Precision :: initializePhi0

    ! Locals
    Double Precision :: shearInstabilityPhi0
    Double Precision :: randomNumbers(10000)
    Double Precision :: rn
    Integer :: index, rnsize
    Double Precision, parameter :: zero=0.d0, half=0.5d0, one=1.d0, two=2.d0

    ! Common
    Common randomNumbers
    ! -------------------------------------------------------------------- !


    shearInstabilityPhi0 = (nPL-one)/G12*((G12-XC)/(XC*nPL*aPL**(one/nPL)))**(nPL/(nPL-one))

    If (phi0 == zero) Then  ! Special case to use phi0 corresponding to instability
      ! Calculate phi0
      initializePhi0 = shearInstabilityPhi0

    Else If (phi0 <= half) Then   ! Use initial condition value of phi0
      initializePhi0 = phi0

    Else If (phi0 >= one) Then   ! 1-D variation
      ! Set index to be integer count of row
      If (phi0 == one) Then ! 1-D variation
        index = NINT(position(1)/0.05d0) + offset
      Else If (phi0 == two) Then ! 2-D variation (different 1-D variation for each ply)
        index = NINT((position(1)/0.05d0) + (100.d0*position(3))) + offset
      Else
        Call log%error('Invalid initial condition received for phi0. Found ' // trim(str(phi0)))
      End If

      ! Wrap in a loop if there are more rows than randomNumbers is long
      rnsize = SIZE(randomNumbers, 1)
      If (index > rnsize) Then
        index = MOD(index, rnsize)
        Call log%warn('Random fiber misalignments are being wrapped in a loop. Increase randomNumberCount in vexternaldb.')
      End If

      ! In case the position is negative
      If (index < 0) Then
        index = rnsize+index
      End If

      ! Get the random number (0 to 1)
      rn = randomNumbers(index)

      ! Scale -phi to +phi
      initializePhi0 = rn*two*shearInstabilityPhi0 - shearInstabilityPhi0

    Else
      Call log%error('Invalid phi0. Found phi0 = '  // trim(str(phi0)) // ' See README.')

    End If

    Return
  End Function initializePhi0


End Module matProp_Mod
