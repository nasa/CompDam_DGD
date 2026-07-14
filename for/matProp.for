#define type(x) TYPE(x), target

Module matProp_Mod
  ! Module for loading and validating material properties

  Private

  Type matProps
    ! Stores the set of properties

    Character(len=80) :: name = "" ! Material name, required

    ! Material properties, grouped by feature
    ! (Required inputs)
    Double Precision :: E1, E2, G12, v12, v23                            ! Req. Transversely isotropic.
    ! (Optional inputs)
    Double Precision :: YT, SL, GYT, GSL, eta_BK, YC, alpha0             ! Opt. Matrix failure
    Double Precision :: E3, G13, G23, v13                                ! Opt. Orthotropic. Defaults are calculated based on transverse isotropy
    Double Precision :: cte(3)                                           ! Opt. Coefficients of thermal expansion (defaults to zero; if not provided, assumes cte(3) = cte(2))
    Double Precision :: T_sf                                             ! Opt. stress free temperature (defaults to zero)
    Double Precision :: aPL, nPL                                         ! Opt. Shear nonlinearity
    Double Precision :: XT, fXT, GXT, fGXT                               ! Opt. Longitudinal tensile damage (max strain failure criterion, bilinear softening)
    Double Precision :: XC, fXC, GXC, fGXC                               ! Opt. Longitudinal compressive damage (max strain failure criterion, bilinear softening)
    Double Precision :: w_kb                                             ! Opt. LaRC15 kink band model (Also requires XC, shear nonlinearity).
    Double Precision :: cl                                               ! Opt. fiber nonlinearity (defaults to zero, which is no fiber nonlinearity)
    Double Precision :: mu                                               ! Opt. Friction. Defaults to zero
    Double Precision :: es(4), gs(4)                                     ! Opt. Schapery theory inputs
    Double Precision :: schaefer_a6, schaefer_b2, schaefer_n, schaefer_A  ! Opt. Schaefer nonlinearity model inputs
    Double Precision :: fatigue_gamma, fatigue_epsilon, fatigue_eta, fatigue_p_mod  ! Opt. CF20 fatigue properties

    ! min and max values for acceptable range
    Double Precision, private :: modulus_min = Tiny(0.d0), modulus_max = Huge(0.d0)
    Double Precision, private :: poissons_min = 0.d0, poissons_max = 1.d0
    Double Precision, private :: strength_min = Tiny(0.d0), strength_max = Huge(0.d0)
    Double Precision, private :: toughness_min = Tiny(0.d0), toughness_max = Huge(0.d0)
    Double Precision, private :: eta_BK_min = Tiny(0.d0), eta_BK_max = Huge(0.d0)
    Double Precision, private :: cte_min = -1.d0, cte_max = 1.d0
    Double Precision, private :: T_sf_min = -Huge(0.d0), T_sf_max = Huge(0.d0)
    Double Precision, private :: aPL_min = 0.d0, aPL_max = Huge(0.d0), nPL_min = 0.d0, nPL_max = Huge(0.d0)
    Double Precision, private :: w_kb_min = Tiny(0.d0), w_kb_max = Huge(0.d0)
    Double Precision, private :: cl_min = 0.d0, cl_max = 33.d0  ! max ensures that the tangent stiffness is > 0 up to 1.5% strain in compression
    Double Precision, private :: alpha0_min = 0.d0, alpha0_max = 1.5707963267949d0  ! 0 to pi/2 radians
    Double Precision, private :: mu_min = 0.d0, mu_max = 1.d0
    Double Precision, private :: schapery_min = -Huge(0.d0), schapery_max = Huge(0.d0)
    Double Precision, private :: schaefer_min = -Huge(0.d0), schaefer_max = Huge(0.d0)
    Double Precision, private :: fatigue_gamma_min = Tiny(0.d0), fatigue_gamma_max = Huge(0.d0)
    Double Precision, private :: fatigue_epsilon_min = Tiny(0.d0), fatigue_epsilon_max = 1.d0
    Double Precision, private :: fatigue_eta_min = Tiny(0.d0), fatigue_eta_max = 1.d0
    Double Precision, private :: fatigue_p_mod_min = -Huge(0.d0), fatigue_p_mod_max = Huge(0.d0)

    ! Flags that indicate if values have been set
    Logical, private :: E1_def = .FALSE., E2_def = .FALSE., G12_def = .FALSE., v12_def = .FALSE., v23_def = .FALSE.
    Logical, private :: YT_def = .FALSE., SL_def = .FALSE., GYT_def = .FALSE., GSL_def = .FALSE., eta_BK_def = .FALSE.
    Logical, private :: E3_def = .FALSE., G13_def = .FALSE., G23_def = .FALSE., v13_def = .FALSE.
    Logical, private :: cte_def(3) = [.FALSE. ,.FALSE. ,.FALSE.]
    Logical, private :: T_sf_def = .FALSE.
    Logical, private :: aPL_def = .FALSE., nPL_def = .FALSE.
    Logical, private :: XT_def = .FALSE., fXT_def = .FALSE., GXT_def = .FALSE., fGXT_def = .FALSE.
    Logical, private :: XC_def = .FALSE., fXC_def = .FALSE., GXC_def = .FALSE., fGXC_def = .FALSE.
    Logical, private :: YC_def = .FALSE., w_kb_def = .FALSE., alpha0_def = .FALSE.
    Logical, private :: cl_def = .FALSE.
    Logical, private :: mu_def = .FALSE.
    Logical, private :: es_def(4) = [.FALSE. ,.FALSE. ,.FALSE. ,.FALSE.]
    Logical, private :: gs_def(4) = [.FALSE. ,.FALSE. ,.FALSE. ,.FALSE.]
    Logical, private :: schaefer_a6_def = .FALSE.
    Logical, private :: schaefer_b2_def = .FALSE.
    Logical, private :: schaefer_n_def = .FALSE.
    Logical, private :: schaefer_A_def = .FALSE.
    Logical, private :: fatigue_gamma_def = .FALSE., fatigue_epsilon_def = .FALSE.
    Logical, private :: fatigue_eta_def = .FALSE., fatigue_p_mod_def = .FALSE.

    ! Calculated properties
    Double Precision :: v21, v31, v32
    Double Precision :: etaL, etaT, ST, phic, alpha0_deg

    Double Precision :: thickness

    ! Flags for features
    Character(len=6) :: featureFlags
    Logical :: matrixDam
    Logical :: cohesive
    Logical :: embedded_cohesive
    Logical :: shearNonlinearity12
    Logical :: shearNonlinearity13
    Logical :: schapery
    Logical :: schaefer
    Logical :: fiberTenDam
    Logical :: fiberCompDamBL
    Logical :: fiberCompDamFKT12
    Logical :: fiberCompDamFKT13
    Logical :: accumulateDissipPlasEnergy
    Logical :: accumulateDissipFractEnergy
    Logical :: friction

  End Type matProps

  Type matNameList
    ! A list of material names for mapping purposes
    Character(len=80) :: name = ""
  End Type matNameList

  ! Public interface
  Public :: matProps
  Public :: matNameList
  Public :: loadMatProps
  Public :: checkForSnapBack
  Public :: getMaterialIndex
  Public :: initializePhi0
  Public :: writeMaterialPropertiesToFile


Contains

  Subroutine loadMatProps(m, cmname, nprops, props, first_call)
    Use forlog_Mod

    !Arguments
    Class(matProps), intent(INOUT) :: m
    Character(len=80), intent(IN) :: cmname
    Integer, intent(IN) :: nprops
    Double Precision, intent(IN) :: props(nprops)
    Logical, intent(IN) :: first_call

    ! Locals
    Character(len=30) :: tmp, featureFlags
    Integer :: numMissingLeadingZeros
    Double Precision, Parameter :: zero=0.d0, one=1.d0, two=2.d0
    ! -------------------------------------------------------------------- !

    m%name = trim(cmname)

    If (first_call) Call log%writeToLog("INFO:Start loadMatProps() for "// trim(m%name))

    Do i=1, nprops
      Select Case (i)
        Case (1)
          ! Convert to string (remove decimal and trailing zeros)
          If (first_call) Call log%writeToLog("INFO:  "// trim(m%name) //" featureFlag: " // str(props(i)))
          write (tmp, *) props(i)
          featureFlags = Adjustl(tmp(1:Index(tmp, '.')-1))
          numMissingLeadingZeros = 6 - Len_trim(featureFlags)
          If (numMissingLeadingZeros > 0) Then
            Do j=1, numMissingLeadingZeros
              featureFlags = '0' // featureFlags
            End Do
          End If

          m%featureFlags = trim(featureFlags)

          ! Parse each position into a feature flag
          ! Position 1: Matrix damage
          ! Position 2: Peak-peak nonlinearity (shear plasticity or schapery micro-damage)
          ! Position 3: Fiber tensile damage
          ! Position 4: Fiber compression damage
          ! Position 5: Energy accumulation
          ! Position 6: Friction

          ! Default friction to enabled
          If (Len_trim(featureFlags) < 6) Call log%error("loadMatProps:"// trim(m%name) //" All 6 feature flags must be specified, only found " // trim(str(Len_trim(featureFlags))))

          ! Set flags
          Do j=1,Len_trim(featureFlags)
            Select Case (j)
              Case (1)
                If (featureFlags(j:j) == '1') Then
                  m%matrixDam = .TRUE.
                  m%cohesive = .FALSE.
                  If (first_call) Call log%writeToLog("INFO:  Matrix damage ENABLED")
                Else
                  m%matrixDam = .FALSE.
                  If (featureFlags(j:j) == '2') Then
                    m%cohesive = .TRUE.
                    If (first_call) Call log%writeToLog("INFO:  Cohesive mode ENABLED")
                  Else
                    m%cohesive = .FALSE.
                    If (first_call) Call log%writeToLog("INFO:  Matrix damage DISABLED")
                  End If
                End If

              Case (2)
                If (featureFlags(j:j) == '1') Then
                  m%shearNonlinearity12 = .TRUE.
                  m%shearNonlinearity13 = .FALSE.
                  m%schapery = .FALSE.
                  m%schaefer = .FALSE.
                  If (first_call) Call log%writeToLog("INFO:  Shear nonlinearity (1-2 plane) ENABLED")
                Else If (featureFlags(j:j) == '2') Then
                  m%shearNonlinearity12 = .FALSE.
                  m%shearNonlinearity13 = .FALSE.
                  m%schapery = .TRUE.
                  m%schaefer = .FALSE.
                  If (first_call) Call log%writeToLog("INFO:  Schapery micro-damage ENABLED")
                Else If (featureFlags(j:j) == '3') Then
                  m%shearNonlinearity12 = .TRUE.
                  m%shearNonlinearity13 = .TRUE.
                  m%schapery = .FALSE.
                  If (first_call) Call log%writeToLog("INFO:  Shear nonlinearity (3-D) ENABLED")
                Else If (featureFlags(j:j) == '4') Then
                  m%shearNonlinearity12 = .FALSE.
                  m%shearNonlinearity13 = .TRUE.
                  m%schapery = .FALSE.
                  If (first_call) Call log%writeToLog("INFO:  Shear nonlinearity (1-3 plane) ENABLED")
                Else If (featureFlags(j:j) == '5') THEN
                  m%shearNonlinearity12 = .FALSE.
                  m%shearNonlinearity13 = .FALSE.
                  m%schapery = .FALSE.
                  m%schaefer = .TRUE.
                  If (first_call) Call log%writeToLog("INFO:  Schaefer nonlinearity ENABLED")
                Else
                  m%shearNonlinearity12 = .FALSE.
                  m%shearNonlinearity13 = .FALSE.
                  m%schapery = .FALSE.
                  m%schaefer = .FALSE.
                  If (first_call) Call log%writeToLog("INFO:  pre-peak nonlinearity DISABLED")
                End If

              Case (3)
                If (featureFlags(j:j) == '1') Then
                  m%fiberTenDam = .TRUE.
                  If (first_call) Call log%writeToLog("INFO:  fiber tensile damage ENABLED")
                Else
                  m%fiberTenDam = .FALSE.
                  If (first_call) Call log%writeToLog("INFO:  fiber tensile damage DISABLED")
                End If

              Case (4)
                If (featureFlags(j:j) == '1') Then
                  m%fiberCompDamBL = .TRUE.
                  m%fiberCompDamFKT12 = .FALSE.
                  m%fiberCompDamFKT13 = .FALSE.
                  If (first_call) Call log%writeToLog("INFO:  fiber comp. damage (CDM) ENABLED")
                Else If (featureFlags(j:j) == '2') Then
                  Call log%error("  fiber comp. model 2 not implemented. TODO")
                Else If (featureFlags(j:j) == '3') Then
                  m%fiberCompDamBL = .FALSE.
                  m%fiberCompDamFKT12 = .TRUE.
                  m%fiberCompDamFKT13 = .FALSE.
                  If (first_call) Call log%writeToLog("INFO:  fiber comp. damage (FKT) in-plane ENABLED")
                Else If (featureFlags(j:j) == '4') Then
                  m%fiberCompDamBL = .FALSE.
                  m%fiberCompDamFKT12 = .FALSE.
                  m%fiberCompDamFKT13 = .TRUE.
                  If (first_call) Call log%writeToLog("INFO:  fiber comp. damage (FKT) out-of-plane ENABLED")
                Else If (featureFlags(j:j) == '5') Then
                  m%fiberCompDamBL = .FALSE.
                  m%fiberCompDamFKT12 = .TRUE.
                  m%fiberCompDamFKT13 = .TRUE.
                  If (first_call) Call log%writeToLog("INFO:  fiber comp. damage (FKT) 3-D ENABLED")
                Else
                  m%fiberCompDamBL = .FALSE.
                  m%fiberCompDamFKT12 = .FALSE.
                  m%fiberCompDamFKT13 = .FALSE.
                  If (first_call) Call log%writeToLog("INFO:  fiber comp. damage DISABLED")
                End If

              Case (5)
                If (featureFlags(j:j) == '0') Then
                  m%accumulateDissipPlasEnergy = .TRUE.
                  m%accumulateDissipFractEnergy = .TRUE.
                  If (first_call) Call log%writeToLog("INFO:  ALLPD: plasticity and fracture will be summed")
                Else If (featureFlags(j:j) == '1') Then
                  m%accumulateDissipPlasEnergy = .FALSE.
                  m%accumulateDissipFractEnergy = .TRUE.
                  If (first_call) Call log%writeToLog("INFO:  ALLPD: plasticity is being ignored (only fracture energy)")
                Else If (featureFlags(j:j) == '2') Then
                  m%accumulateDissipPlasEnergy = .TRUE.
                  m%accumulateDissipFractEnergy = .FALSE.
                  If (first_call) Call log%writeToLog("INFO:  ALLPD: fracture is being ignored (only plastic energy)")
                Else
                  m%accumulateDissipPlasEnergy = .FALSE.
                  m%accumulateDissipFractEnergy = .FALSE.
                  If (first_call) Call log%writeToLog("INFO:  ALLPD: inelastic energy accumulation DISABLED")
                End If

              Case (6)
                If (featureFlags(j:j) == '1') Then
                  m%friction = .TRUE.
                  If (first_call) Call log%writeToLog("INFO:  friction ENABLED")
                Else
                  m%friction = .FALSE.
                  If (first_call) Call log%writeToLog("INFO:  friction DISABLED")
                End If

              Case Default
                Call log%error("  Unknown position found in feature flags")
            End Select
          End Do
        Case (2)  ! Reserved
        Case (3)
          m%thickness = props(i)

        Case (4)  ! Reserved
        Case (5)
          Call verifyAndSaveProperty_double(m, 'fatigue_gamma', props(i), m%fatigue_gamma_min, m%fatigue_gamma_max, m%fatigue_gamma, m%fatigue_gamma_def)

        Case (6)
          Call verifyAndSaveProperty_double(m, 'fatigue_epsilon', props(i), m%fatigue_epsilon_min, m%fatigue_epsilon_max, m%fatigue_epsilon, m%fatigue_epsilon_def)

        Case (7)
          Call verifyAndSaveProperty_double(m, 'fatigue_eta', props(i), m%fatigue_eta_min, m%fatigue_eta_max, m%fatigue_eta, m%fatigue_eta_def)

        Case (8)
          Call verifyAndSaveProperty_double(m, 'fatigue_p_mod', props(i), m%fatigue_p_mod_min, m%fatigue_p_mod_max, m%fatigue_p_mod, m%fatigue_p_mod_def)

        Case (9)
          Call verifyAndSaveProperty_double(m, 'E1', props(i), m%modulus_min, m%modulus_max, m%E1, m%E1_def)

        Case (10)
          Call verifyAndSaveProperty_double(m, 'E2', props(i), m%modulus_min, m%modulus_max, m%E2, m%E2_def)

        Case (11)
          Call verifyAndSaveProperty_double(m, 'G12', props(i), m%modulus_min, m%modulus_max, m%G12, m%G12_def)

        Case (12)
          Call verifyAndSaveProperty_double(m, 'v12', props(i), m%poissons_min, m%poissons_max, m%v12, m%v12_def)

        Case (13)
          Call verifyAndSaveProperty_double(m, 'v23', props(i), m%poissons_min, m%poissons_max, m%v23, m%v23_def)

        Case (14)
          Call verifyAndSaveProperty_double(m, 'YT', props(i), m%strength_min, m%strength_max, m%YT, m%YT_def)

        Case (15)
          Call verifyAndSaveProperty_double(m, 'SL', props(i), m%strength_min, m%strength_max, m%SL, m%SL_def)

        Case (16)
          Call verifyAndSaveProperty_double(m, 'GYT', props(i), m%toughness_min, m%toughness_max, m%GYT, m%GYT_def)

        Case (17)
          Call verifyAndSaveProperty_double(m, 'GSL', props(i), m%toughness_min, m%toughness_max, m%GSL, m%GSL_def)

        Case (18)
          Call verifyAndSaveProperty_double(m, 'eta_BK', props(i), m%eta_BK_min, m%eta_BK_max, m%eta_BK, m%eta_BK_def)

        Case (19)
          Call verifyAndSaveProperty_double(m, 'YC', props(i), m%strength_min, m%strength_max, m%YC, m%YC_def)

        Case (20)
          Call verifyAndSaveProperty_double(m, 'alpha0', props(i), m%alpha0_min, m%alpha0_max, m%alpha0, m%alpha0_def)

        Case (21)
          Call verifyAndSaveProperty_double(m, 'E3', props(i), m%modulus_min, m%modulus_max, m%E3, m%E3_def)

        Case (22)
          Call verifyAndSaveProperty_double(m, 'G13', props(i), m%modulus_min, m%modulus_max, m%G13, m%G13_def)

        Case (23)
          Call verifyAndSaveProperty_double(m, 'G23', props(i), m%modulus_min, m%modulus_max, m%G23, m%G23_def)

        Case (24)
          Call verifyAndSaveProperty_double(m, 'v13', props(i), m%poissons_min, m%poissons_max, m%v13, m%v13_def)

        Case (25)
          Call verifyAndSaveProperty_double(m, 'alpha11', props(i), m%cte_min, m%cte_max, m%cte(1), m%cte_def(1))

        Case (26)
          Call verifyAndSaveProperty_double(m, 'alpha22', props(i), m%cte_min, m%cte_max, m%cte(2), m%cte_def(2))

        Case (27)
          Call verifyAndSaveProperty_double(m, 'alpha_PL', props(i), m%aPL_min, m%aPL_max, m%aPL, m%aPL_def)

        Case (28)
          Call verifyAndSaveProperty_double(m, 'n_PL', props(i), zero, m%nPL_max, m%nPL, m%nPL_def)

        Case (29)
          Call verifyAndSaveProperty_double(m, 'XT', props(i), m%strength_min, m%strength_max, m%XT, m%XT_def)

        Case (30)
          Call verifyAndSaveProperty_double(m, 'fXT', props(i), zero, one, m%fXT, m%fXT_def)

        Case (31)
          Call verifyAndSaveProperty_double(m, 'GXT', props(i), m%toughness_min, m%toughness_max, m%GXT, m%GXT_def)

        Case (32)
          Call verifyAndSaveProperty_double(m, 'fGXT', props(i), zero, one, m%fGXT, m%fGXT_def)

        Case (33)
          Call verifyAndSaveProperty_double(m, 'XC', props(i), m%strength_min, m%strength_max, m%XC, m%XC_def)

        Case (34)
          Call verifyAndSaveProperty_double(m, 'fXC', props(i), Tiny(zero), one, m%fXC, m%fXC_def)

        Case (35)
          Call verifyAndSaveProperty_double(m, 'GXC', props(i), m%toughness_min, m%toughness_max, m%GXC, m%GXC_def)

        Case (36)
          Call verifyAndSaveProperty_double(m, 'fGXC', props(i), zero, one, m%fGXC, m%fGXC_def)

        Case (37)
          Call verifyAndSaveProperty_double(m, 'cl', props(i), m%cl_min, m%cl_max, m%cl, m%cl_def)

        Case (38)
          Call verifyAndSaveProperty_double(m, 'w_kb', props(i), m%w_kb_min, m%w_kb_max, m%w_kb, m%w_kb_def)

        Case (39)
          Call verifyAndSaveProperty_double(m, 'T_sf', props(i), m%T_sf_min, m%T_sf_max, m%T_sf, m%T_sf_def)

        Case (40)
          Call verifyAndSaveProperty_double(m, 'mu', props(i), m%mu_min, m%mu_max, m%mu, m%mu_def)

        Case (41)
          Call verifyAndSaveProperty_double(m, 'schaefer_a6', props(i), m%schaefer_min, m%schaefer_max, m%schaefer_a6, m%schaefer_a6_def)

        Case (42)
          Call verifyAndSaveProperty_double(m, 'schaefer_b2', props(i), m%schaefer_min, m%schaefer_max, m%schaefer_b2, m%schaefer_b2_def)

        Case (43)
          Call verifyAndSaveProperty_double(m, 'schaefer_n', props(i), m%schaefer_min, m%schaefer_max, m%schaefer_n, m%schaefer_n_def)

        Case (44)
          Call verifyAndSaveProperty_double(m, 'schaefer_A', props(i), m%schaefer_min, m%schaefer_max, m%schaefer_A, m%schaefer_A_def)

        Case (45)
          Call verifyAndSaveProperty_double(m, 'es0', props(i), m%schapery_min, m%schapery_max, m%es(1), m%es_def(1))

        Case (46)
          Call verifyAndSaveProperty_double(m, 'es1', props(i), m%schapery_min, m%schapery_max, m%es(2), m%es_def(2))

        Case (47)
          Call verifyAndSaveProperty_double(m, 'es2', props(i), m%schapery_min, m%schapery_max, m%es(3), m%es_def(3))

        Case (48)
          Call verifyAndSaveProperty_double(m, 'es3', props(i), m%schapery_min, m%schapery_max, m%es(4), m%es_def(4))

        Case (49)
          Call verifyAndSaveProperty_double(m, 'gs0', props(i), m%schapery_min, m%schapery_max, m%gs(1), m%gs_def(1))

        Case (50)
          Call verifyAndSaveProperty_double(m, 'gs1', props(i), m%schapery_min, m%schapery_max, m%gs(2), m%gs_def(2))

        Case (51)
          Call verifyAndSaveProperty_double(m, 'gs2', props(i), m%schapery_min, m%schapery_max, m%gs(3), m%gs_def(3))

        Case (52)
          Call verifyAndSaveProperty_double(m, 'gs3', props(i), m%schapery_min, m%schapery_max, m%gs(4), m%gs_def(4))

        Case (53)
          Call verifyAndSaveProperty_double(m, 'alpha33', props(i), m%cte_min, m%cte_max, m%cte(3), m%cte_def(3))

        Case Default
          Call log%error("loadMatProps:"// trim(m%name) //" Unknown property #"//trim(str(i))//" found. Expecting nprops <= 40")
      End Select

    End Do

    ! Feature flags
    If (.NOT. m%friction) m%mu = zero
    If (m%mu == zero) m%friction = .False.

    ! Calculate alpha0
    If (m%matrixDam) Then
      m%alpha0_deg = NINT(m%alpha0*45.d0/ATAN(one))
      m%alpha0 = alpha0_DGD(m)
    End If

    ! Checks that a consistent set of properties has been defined
    Call consistencyChecks(m, first_call)

    ! Calculated properties
    m%v21 = m%E2*m%v12/m%E1
    m%v31 = m%E3*m%v13/m%E1
    m%v32 = m%E3*m%v23/m%E2

    If (first_call) Call log%writeToLog("INFO:End loadMatProps() for "// trim(m%name))

    Return
  End Subroutine loadMatProps


  Subroutine verifyAndSaveProperty_str(m, key, value, min, max, saveTo, flag)
    ! Checks if the value is within the specified bounds. Prints an error message
    ! which kills the analysis if a value is out of bounds. Otherwise, save the property
    ! and sets the flag to indicate that the property has been read in.

    Use forlog_Mod

    !Arguments
    Class(matProps), intent(INOUT) :: m
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
      Call log%error(" PROPERTY ERROR " // trim(m%name) // ": " // trim(key) // " cannot be less than " // trim(str(min)) // ". Found value: " // trim(str(valueDbl)))

    Else If (valueDbl > max) Then
      Call log%error(" PROPERTY ERROR " // trim(m%name) // ": " // trim(key) // " cannot be greater than " // trim(str(max)) // ". Found value: " // trim(str(valueDbl)))

    Else
      Call log%debug(" loadMatProps: Loaded " // trim(key) // " = " // trim(str(valueDbl)))
    End If

    ! Save the value and set the flag
    saveTo = valueDbl
    flag = .TRUE.

    Return
  End Subroutine verifyAndSaveProperty_str


  Subroutine verifyAndSaveProperty_double(m, key, value, min, max, saveTo, flag)
    ! Checks if the value is within the specified bounds. Prints an error message
    ! which kills the analysis if a value is out of bounds. Otherwise, save the property
    ! and sets the flag to indicate that the property has been read in.

    Use forlog_Mod

    !Arguments
    Class(matProps), intent(INOUT) :: m
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
        Return
      Else
        Call log%error(" PROPERTY ERROR " // trim(m%name) // ": " // trim(key) // " cannot be less than " // trim(str(min)) // ". Found value: " // trim(str(value)))
      End If

    Else If (value > max) Then
      Call log%error(" PROPERTY ERROR " // trim(m%name) // ": " // trim(key) // " cannot be greater than " // trim(str(max)) // ". Found value: " // trim(str(value)))

    Else
      Call log%debug(" loadMatProps: Loaded " // trim(key) // " = " // trim(str(value)))
    End If

    ! Save the value and set the flag
    saveTo = value
    flag = .TRUE.

    Return
  End Subroutine verifyAndSaveProperty_double


  Subroutine consistencyChecks(m, write_info)
    ! Checks that a consistent set of properties has been defined

    Use forlog_Mod

    ! Arguments
    Type(matProps), intent(INOUT) :: m
    Logical, intent(IN) :: write_info

    ! Locals
    Double Precision :: eps_0tanstiff, eps_0_fn
    Double Precision, Parameter :: zero=0.d0, one=1.d0, two=2.d0, four=4.d0
    ! -------------------------------------------------------------------- !

    If (.NOT. m%cohesive) Then
      ! Check that all elastic properties have been defined for transverse isotropy
      If (.NOT. m%E1_def) Call log%error('PROPERTY ERROR:'// trim(m%name) // ', Must define a value for E1')
      If (.NOT. m%E2_def) Call log%error('PROPERTY ERROR:'// trim(m%name) // ', Must define a value for E2')
      If (.NOT. m%G12_def) Call log%error('PROPERTY ERROR:'// trim(m%name) // ', Must define a value for G12')
      If (.NOT. m%v12_def) Call log%error('PROPERTY ERROR:'// trim(m%name) // ', Must define a value for v12')
      If (.NOT. m%v23_def) Call log%error('PROPERTY ERROR:'// trim(m%name) // ', Must define a value for v23')

      ! Check if orthotropic elastic properties have been defined
      If (m%E3_def .OR. m%G13_def .OR. m%G23_def) Then
        If (.NOT. m%E3_def) Call log%error('PROPERTY ERROR:'// trim(m%name) // ', Some orthotropic elastic properties are missing. Must define a value for E3.')
        If (.NOT. m%G13_def) Call log%error('PROPERTY ERROR:'// trim(m%name) // ', Some orthotropic elastic properties are missing. Must define a value for G13.')
        If (.NOT. m%G23_def) Call log%error('PROPERTY ERROR:'// trim(m%name) // ', Some orthotropic elastic properties are missing. Must define a value for G23.')
        If (.NOT. m%v13_def) Call log%error('PROPERTY ERROR:'// trim(m%name) // ', Some orthotropic elastic properties are missing. Must define a value for v13.')
        If (write_info) Call log%writeToLog('INFO:  Orthotropic constants are defined directly')
      Else
        ! Compute these properties assuming transverse isotropy
        If (write_info) Call log%writeToLog('INFO:  Assuming transverse isotropy')
        m%E3  = m%E2
        m%G13 = m%G12
        m%v13 = m%v12
        m%G23 = m%E2/two/(one + m%v23)
      End If

    Else If (m%cohesive) Then  ! Check if cohesive element properties have been specified
      ! Check for cohesive stiffness material properties
      If (.NOT. m%E3_def) Call log%error('PROPERTY ERROR:'// trim(m%name) // ', Cohesive material laws require a definition for E3.')
      If (m%G13_def .AND. m%G23_def) Then
        m%embedded_cohesive = .TRUE.
        If (write_info) Call log%writeToLog('INFO:  G13 and G23 are defined for a cohesive element. Assuming finite thickness.')
      Else
        m%embedded_cohesive = .FALSE.
      End If
      ! Check for cohesive strength material properties
      If (.NOT. m%YT_def) Call log%error('PROPERTY ERROR:'// trim(m%name) // ', Some cohesive element properties are missing. Must define a value for YT.')
      If (.NOT. m%SL_def) Call log%error('PROPERTY ERROR:'// trim(m%name) // ', Some cohesive element properties are missing. Must define a value for SL.')
      If (m%YC_def .AND. m%alpha0_def) Then
        If (m%alpha0 .EQ. zero) Then
          Call log%warn('PROPERTY: alpha0 = 0, which leads to etaT = inf. Assuming etaL = etaT = zero and S_T = S_L.')
          m%etaL = zero
          m%etaT = zero
          m%ST   = m%SL
        Else
          ! Compute these properties for the transverse shear strength
          m%etaL = -m%SL*COS(two*m%alpha0)/(m%YC*COS(m%alpha0)*COS(m%alpha0))
          m%etaT = -one/TAN(two*m%alpha0)
          m%ST   = m%YC*COS(m%alpha0)*(SIN(m%alpha0) + COS(m%alpha0)/TAN(two*m%alpha0))
        End If
      Else
        If (write_info) Call log%writeToLog('INFO:  YC and/or alpha0 not defined. Assuming that S_T = S_L.')
        m%etaL = zero
        m%etaT = zero
        m%ST   = m%SL
      End If
      ! Check for cohesive fracture toughness material properties
      If (.NOT. m%GYT_def) Call log%error('PROPERTY ERROR:'// trim(m%name) // ', Some cohesive element properties are missing. Must define a value for GYT.')
      If (.NOT. m%GSL_def) Call log%error('PROPERTY ERROR:'// trim(m%name) // ', Some cohesive element properties are missing. Must define a value for GSL.')
      If (.NOT. m%eta_BK_def) Call log%error('PROPERTY ERROR:'// trim(m%name) // ', Some cohesive element properties are missing. Must define a value for eta_BK.')

      If (write_info) Call log%writeToLog('INFO:  All required cohesive element properties are defined.')
    End IF

    ! Check if matrix damage properties have been specified
    If (m%matrixDam) Then
      If (.NOT. m%YT_def) Call log%error('PROPERTY ERROR:'// trim(m%name) // ', Some matrix damage properties are missing. Must define a value for YT.')
      If (.NOT. m%SL_def) Call log%error('PROPERTY ERROR:'// trim(m%name) // ', Some matrix damage properties are missing. Must define a value for SL.')
      If (.NOT. m%GYT_def) Call log%error('PROPERTY ERROR:'// trim(m%name) // ', Some matrix damage properties are missing. Must define a value for GYT.')
      If (.NOT. m%GSL_def) Call log%error('PROPERTY ERROR:'// trim(m%name) // ', Some matrix damage properties are missing. Must define a value for GSL.')
      If (.NOT. m%eta_BK_def) Call log%error('PROPERTY ERROR:'// trim(m%name) // ', Some matrix damage properties are missing. Must define a value for eta_BK.')
      If (.NOT. m%YC_def) Call log%error('PROPERTY ERROR:'// trim(m%name) // ', Some matrix damage properties are missing. Must define a value for YC.')
      If (.NOT. m%alpha0_def) Call log%error('PROPERTY ERROR:'// trim(m%name) // ', Some matrix damage properties are missing. Must define a value for alpha0.')
      m%etaL = -m%SL*COS(two*m%alpha0)/(m%YC*COS(m%alpha0)*COS(m%alpha0))
      m%etaT = -one/TAN(two*m%alpha0)
      m%ST   = m%YC*COS(m%alpha0)*(SIN(m%alpha0) + COS(m%alpha0)/TAN(two*m%alpha0))

      If (write_info) Call log%writeToLog('INFO:  Matrix damage is enabled')
    Else
      If (write_info) Call log%writeToLog('INFO:  Matrix damage is disabled')
    End If

    ! Check if CTEs have been defined
    If (m%cte_def(1) .OR. m%cte_def(2) .OR. m%cte_def(3)) Then
      If (.NOT. m%cte_def(1)) Call log%error('PROPERTY ERROR:'// trim(m%name) // ', Some CTE properties are missing. Must define a value for alpha11.')
      If (.NOT. m%cte_def(2)) Call log%error('PROPERTY ERROR:'// trim(m%name) // ', Some CTE properties are missing. Must define a value for alpha22.')
      If (.NOT. m%cte_def(3)) Then
        If (write_info) Call log%writeToLog('INFO:  Assuming alpha33 = alpha22')
        m%cte(3) = m%cte(2)
      End If
      If (write_info) Call log%writeToLog('INFO:  CTEs have been defined')
    Else
      If (write_info) Call log%writeToLog('INFO:  CTEs are being set to zero')
      m%cte(1) = zero
      m%cte(2) = zero
      m%cte(3) = zero
    End If
    
    ! Check if stress free residual temperature has been defined
    If (m%T_sf_def) Then
      If (write_info) Call log%writeToLog('INFO:  residual stress free temperature has been defined')
    Else
      If (write_info) Call log%writeToLog('INFO:  the residual stress free temperature is being set to zero')
      m%T_sf = zero
    End If

    ! Check if shear nonlinearity properties have been defined
    If (m%shearNonlinearity12 .OR. m%shearNonlinearity13) Then
      If (.NOT. m%aPL_def) Call log%error('PROPERTY ERROR:'// trim(m%name) // ', Some shear-nonlinearity properties are missing. Must define a value for alpha_PL.')
      If (.NOT. m%nPL_def) Call log%error('PROPERTY ERROR:'// trim(m%name) // ', Some shear-nonlinearity properties are missing. Must define a value for n_PL.')
      If (m%aPL .EQ. zero) Call log%error('PROPERTY ERROR:'// trim(m%name) // ', Shear nonlinearity is enabled, but alpha_PL is zero')
      If (m%nPL .EQ. zero) Call log%error('PROPERTY ERROR:'// trim(m%name) // ', Shear nonlinearity is enabled, but n_PL is zero')
      If (write_info) Call log%writeToLog('INFO:  Shear-nonlinearity properties have been defined')
    Else
      If (write_info) Call log%writeToLog('INFO:  Shear-nonlinearity is disabled')
    End If

    ! Check if all Schapery micro-damage material properties have been defined
    If (m%schapery) Then
      If (.NOT. (all(m%es_def) .AND. all(m%gs_def))) Call log%error('PROPERTY ERROR:'// trim(m%name) // ', Some Schapery micro-damage are missing.')
      If (write_info) Call log%writeToLog('INFO:  Schapery micro-damage properties have been defined')
    Else
      If (write_info) Call log%writeToLog('INFO:  Schapery micro-damage is disabled')
      ! Default values that will cause no pre-peak nonlinearity
      m%es = zero
      m%es(1) = one
      m%gs = zero
      m%gs(1) = one
    End If

    ! Check is all the parameters for Schaefer nonlinearity model have been defined
    If (m%schaefer) Then
      If (.NOT. m%schaefer_a6_def ) Call log%error('PROPERTY ERROR:'// trim(m%name) // ', Some Schaefer nonlinearity properties are missing. Must define schaefer_a6')
      If (.NOT. m%schaefer_b2_def ) Call log%error('PROPERTY ERROR:'// trim(m%name) // ', Some Schaefer nonlinearity properties are missing. Must define schaefer_b2')
      If (.NOT. m%schaefer_n_def ) Call log%error('PROPERTY ERROR:'// trim(m%name) // ', Some Schaefer nonlinearity properties are missing. Must define schaefer_n')
      If (.NOT. m%schaefer_A_def ) Call log%error('PROPERTY ERROR:'// trim(m%name) // ', Some Schaefer nonlinearity properties are missing. Must define schaefer_A')
      If (write_info) Call log%writeToLog('INFO:  Schaefer properties have been defined')
    Else
      If (write_info) Call log%writeToLog('INFO:  Schaefer is disabled')
      ! Default values that will cause no pre-peak nonlinearity
    End If

    ! Check if fiber tensile damage properties have been defined
    If (m%fiberTenDam) Then
      If (.NOT. m%XT_def) Call log%error('PROPERTY ERROR:'// trim(m%name) // ', Some fiber tensile damage properties are missing. Must define a value for XT.')
      If (.NOT. m%fXT_def) Call log%error('PROPERTY ERROR:'// trim(m%name) // ', Some fiber tensile damage properties are missing. Must define a value for fXT.')
      If (.NOT. m%GXT_def) Call log%error('PROPERTY ERROR:'// trim(m%name) // ', Some fiber tensile damage properties are missing. Must define a value for GXT.')
      If (.NOT. m%fGXT_def) Call log%error('PROPERTY ERROR:'// trim(m%name) // ', Some fiber tensile damage properties are missing. Must define a value for fGXT.')

      If (write_info) Call log%writeToLog('INFO:  fiber tensile damage properties have been defined')
    Else
      If (write_info) Call log%writeToLog('INFO:  fiber tensile damage is disabled')
    End If

    ! Check if the CDM fiber compression damage properties have been defined
    If (m%fiberCompDamBL) Then
      If (.NOT. m%XC_def) Call log%error('PROPERTY ERROR:'// trim(m%name) // ', Some fiber compression damage properties are missing (BL model). Must define a value for XC.')
      If (.NOT. m%fXC_def) Call log%error('PROPERTY ERROR:'// trim(m%name) // ', Some fiber compression damage properties are missing (BL model). Must define a value for fXC.')
      If (.NOT. m%GXC_def) Call log%error('PROPERTY ERROR:'// trim(m%name) // ', Some fiber compression damage properties are missing (BL model). Must define a value for GXC.')
      If (.NOT. m%fGXC_def) Call log%error('PROPERTY ERROR:'// trim(m%name) // ', Some fiber compression damage properties are missing (BL model). Must define a value for fGXC.')

      If (write_info) Call log%writeToLog('INFO:  fiber compression damage BL (model 1) properties have been defined')
    Else
      If (write_info) Call log%writeToLog('INFO:  fiber compression damage BL (model 1) is disabled')
    End If

    ! Check if the DGD fiber compression damage properties have been defined
    If (m%fiberCompDamFKT12 .OR. m%fiberCompDamFKT13) Then
      If (.NOT. m%XC_def) Call log%error('PROPERTY ERROR:'// trim(m%name) // ', Some fiber compression damage properties are missing (FKT model). Must define a value for XC.')
      If (.NOT. m%YC_def) Call log%error('PROPERTY ERROR:'// trim(m%name) // ', Some fiber compression damage properties are missing (FKT model). Must define a value for YC.')
      If (.NOT. m%w_kb_def) Call log%error('PROPERTY ERROR:'// trim(m%name) // ', Some fiber compression damage properties are missing (FKT model). Must define a value for w_kb.')
      If (.NOT. m%alpha0_def) Call log%error('PROPERTY ERROR:'// trim(m%name) // ', Some fiber compression damage properties are missing (FKT model). Must define a value for alpha0.')
      If (.NOT. m%SL_def) Call log%error('PROPERTY ERROR:'// trim(m%name) // ', Some fiber compression damage properties are missing (FKT model). Must define a value for SL.')

      ! Make sure shear nonlinearity is enabled
      If (m%fiberCompDamFKT12 .AND. (.NOT. m%shearNonlinearity12)) Call log%error('PROPERTY ERROR:'// trim(m%name) // ', Shear-nonlinearity 1-2 must be enabled with fiber compression damage FKT 1-2.')
      If (m%fiberCompDamFKT13 .AND. (.NOT. m%shearNonlinearity13)) Call log%error('PROPERTY ERROR:'// trim(m%name) // ', Shear-nonlinearity 1-3 must be enabled with fiber compression damage FKT 1-3.')

      If (m%fiberCompDamFKT12 .AND. m%fiberCompDamFKT13 .AND. write_info) Call log%writeToLog('INFO:  fiber compression damage FKT 3-D (model 5) properties have been defined')
      If (m%fiberCompDamFKT12 .AND. write_info) Call log%writeToLog('INFO:  fiber compression damage FKT in-plane (model 3) properties have been defined')
      If (m%fiberCompDamFKT13 .AND. write_info) Call log%writeToLog('INFO:  fiber compression damage FKT out-of-plane (model 4) properties have been defined')
    Else
      If (write_info) Call log%writeToLog('INFO:  fiber compression damage FKT is disabled')
    End If

    ! check if fiber nonlinearity has been defined
    If (m%cl_def) Then
      ! Check that the tangent stiffness does not go to zero before fiber failure
      If (m%fiberCompDamBL) Then
        eps_0tanstiff = one/(two*m%cl)
        eps_0_fn = (-m%E1+SQRT(m%E1**two-four*m%E1*m%cl*m%XC))/(two*m%E1*m%cl)  ! Strain at fiber compression failure
        If (eps_0tanstiff < eps_0_fn) Then
          Call log%error('PROPERTY ERROR:'// trim(m%name) // ', Tangent stiffness goes (eps11 =  ' // trim(str(eps_0tanstiff)) // ')to zero before fiber compression failure strain is reached(eps11 = ' // &
            trim(str(eps_0_fn)) // '). Adjust cl, E1, and/or XC.')
        End If
      End If
      If (write_info) Call log%writeToLog('INFO:  fiber nonlinearity has been defined')
    Else
      If (write_info) Call log%writeToLog('INFO:  fiber nonlinearity has not been defined. Setting cl = 0')
      m%cl = zero
    End If

    ! check if fatigue properties have been defined
    If (m%fatigue_gamma_def .AND. m%fatigue_epsilon_def .AND. m%fatigue_eta_def .AND. m%fatigue_p_mod_def) Then
      If (write_info) Call log%writeToLog('INFO:  cohesive fatigue properties have been defined')
    Else
      If (write_info) Call log%writeToLog('INFO:  cohesive fatigue properties have not been defined. Using default values.')
      If (.NOT. m%fatigue_gamma_def) m%fatigue_gamma = 1.d7
      If (.NOT. m%fatigue_epsilon_def) m%fatigue_epsilon = 0.2d0
      If (.NOT. m%fatigue_eta_def) m%fatigue_eta = 0.95d0
      If (.NOT. m%fatigue_p_mod_def) m%fatigue_p_mod = zero
    End If

    Return
  End Subroutine consistencyChecks


  Subroutine checkForSnapBack(m, Lc, elementNumber, fiber_tension_snapback, fiber_compression_snapback)
    ! Issues an error if an element has a dimension larger than its snap-back size for active damage modes. For superposed softening
    ! laws, an attempt is made to switch to a linear softening response to avoid snap-back.

    Use forlog_Mod

    ! Arguments
    Type(matProps), intent(IN) :: m     ! Material properties
    Double Precision, intent(IN) :: Lc(3)  ! Element length
    Integer, intent(IN) :: elementNumber   ! Element number
    Logical, intent(OUT) :: fiber_tension_snapback, fiber_compression_snapback

    ! Locals
    Double Precision :: Lc_fT, Lc_fT_A, Lc_fT_B, Lc_fC, Lc_fC_A, Lc_fC_B  ! Maximum element sizes based on fiber properties
    Double Precision :: Lc_mT, Lc_SL, Lc_YC                               ! Maximum element sizes based on matrix properties
    Double Precision, Parameter :: zero=0.d0, one=1.d0, two=2.d0
    ! -------------------------------------------------------------------- !
    
    fiber_tension_snapback = .False.
    fiber_compression_snapback = .False.

    ! Fiber Tension
    If (m%fiberTenDam) Then
      If (m%fXT /= m%fGXT) Then  ! Superposed stress-strain responses
        Lc_fT_A = two*m%GXT*m%fGXT*m%E1/(m%XT**2*m%fXT)
        Lc_fT_B = two*m%GXT*(one - m%fGXT)*m%E1/(m%XT**2*(one - m%fXT))
        If (Lc(1) > Lc_fT_A .OR. Lc(1) > Lc_fT_B) Then
          Call log%warn("Snap-back in superposed fiber tensile softening laws of element "// trim(str(elementNumber)) //".")
          fiber_tension_snapback = .True.
        End If
      End If
      
      If (m%fGXT == m%fXT .OR. fiber_tension_snapback) Then  ! Linear softening response
        Lc_fT = two*m%GXT*m%E1/m%XT**2
        If (Lc(1) > Lc_fT) Then
          Call log%error("Element " // trim(str(elementNumber)) // " has size " // trim(str(Lc(1))) // " which is > fiber tension snap-back threshold. Set element size < " // trim(str(Lc_fT)))
        End If
      End If

    End If

    ! Fiber Compression
    If (m%fiberCompDamBL) Then
      If (m%fXC /= m%fGXC) Then  ! Superposed stress-strain responses
        Lc_fC_A = two*m%GXC*m%fGXC*m%E1/(m%XC**2*m%fXC)
        Lc_fC_B = two*m%GXC*(one - m%fGXC)*m%E1/(m%XC**2*(one - m%fXC))
        If (Lc(1) > Lc_fC_A .OR. Lc(1) > Lc_fC_B) Then
          Call log%warn("Snap-back in superposed fiber compression softening laws of element "// trim(str(elementNumber)) //".")
          fiber_compression_snapback = .True.
        End If
      End If
      
      If (m%fGXC == m%fXC .OR. fiber_compression_snapback) Then  ! Linear softening response
        Lc_fC = two*m%GXC*m%E1/m%XC**2
        If (Lc(1) > Lc_fC) Then
          Call log%error("Element " // trim(str(elementNumber)) // " has size " // trim(str(Lc(1))) // " which is > fiber compression snap-back threshold 1. Set element size < " // trim(str(Lc_fC)))
        End If
      End If

    End If

    If (m%matrixDam) Then
      ! Matrix Tension
      Lc_mT = two*m%GYT*m%E2/m%YT/m%YT
      If (Lc(2) > Lc_mT) Then
        Call log%error("Element " // trim(str(elementNumber)) // " has size " // trim(str(Lc(2))) // " which is > matrix mode I snap-back threshold. Set element size < " // trim(str(Lc_mT)))
      End If

      ! Matrix Shear
      Lc_SL = two*m%GSL*m%G12/m%SL**2
      If (Lc(2) > Lc_SL) Then
        Call log%error("Element " // trim(str(elementNumber)) // " has size " // trim(str(Lc(2))) // " which is > matrix mode II snap-back threshold. Set element size < " // trim(str(Lc_SL)))
      End If

      ! Matrix Compression
      Lc_YC = two*m%GSL*m%E2/(m%YC**2*COS(m%alpha0))
      If (Lc(2) > Lc_YC) Then
        Call log%error("Element " // trim(str(elementNumber)) // " has size " // trim(str(Lc(2))) // " which is > matrix compression snap-back threshold. Set element size < " // trim(str(Lc_YC)))
      End If
    End If

    Return
  End Subroutine checkForSnapBack

  Function alpha0_DGD(m)
    ! Determines the orientation of the angle alpha0 when subject to sigma22 = -Yc
    !
    ! The material input alpha0 is the orientation of the crack surface normal in the 2--3 plane
    ! with respect to the material 2 direction for a crack caused by pure matrix compression.
    ! alpha0 is the orientation of the resulting crack after the material has been unloaded, as it
    ! would be measured in a tested specimen. This function calculates the orientation of the
    ! crack normal defined by alpha0 when damage initiates due to pure compressive matrix loading.

    Use forlog_Mod
    Use matrixAlgUtil_Mod
    Use stress_Mod

    ! Arguments
    Type(matProps), intent(IN) :: m

    ! Locals
    Double Precision :: E3, G13, v13, G23    ! For transverse isotropy assumption
    Double Precision :: C(6,6)               ! 3-D Stiffness
    Double Precision :: F(3)                 ! Represents diagonal of deformation gradient tensor
    Double Precision :: Residual(3)          ! Residual vector
    Double Precision :: tolerance
    Double Precision :: err
    Double Precision :: Jac(3,3)             ! Jacobian
    Integer :: counter
    Double Precision, parameter :: zero=0.d0, one=1.d0, two=2.d0
    ! -------------------------------------------------------------------- !

    Call log%debug('Start of Function alpha0_DGD')

    tolerance = 1.d-4  ! alphaLoop tolerance

    ! Assume transverse isotropy if needed
    If (m%E3 < one) Then
      E3 = m%E2
    Else
      E3 = m%E3
    End If
    If (m%G13 < one) Then
      G13 = m%G12
    Else
      G13 = m%G13
    End If
    If (m%v13 == zero) Then
      v13 = m%v12
    Else
      v13 = m%v13
    End If
    If (m%G23 < one) Then
      G23 = m%E2/two/(one + m%v23)
    Else
      G23 = m%G23
    End If

    ! Build the stiffness matrix
    C = StiffFunc(6, m%E1, m%E2, E3, m%G12, G13, G23, m%v12, v13, m%v23, zero, zero, zero)

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
      Residual(1) = (C(1,1)*(F(1)*F(1) - one) + C(1,2)*(F(2)*F(2) - one) + C(1,3)*(F(3)*F(3) - one))/two
      Residual(2) = (C(2,1)*(F(1)*F(1) - one) + C(2,2)*(F(2)*F(2) - one) + C(2,3)*(F(3)*F(3) - one))/two + m%Yc*F(1)*F(3)/F(2)
      Residual(3) = (C(3,1)*(F(1)*F(1) - one) + C(3,2)*(F(2)*F(2) - one) + C(3,3)*(F(3)*F(3) - one))/two

      ! Check for convergence
      err = Length(Residual)

      ! If converged,
      If (err < tolerance) Then
        alpha0_DGD = ATAN(F(2)/F(3)*TAN(m%alpha0))
        EXIT alphaLoop
      End If
      IF (counter == counter_max) Call log%error('Function alpha0_DGD failed to converge, material: ' // trim(m%name))

      ! Define the Jacobian matrix, J
      Jac = zero

      Jac(1,1) = C(1,1)*F(1)
      Jac(1,2) = C(1,2)*F(2)
      Jac(1,3) = C(1,3)*F(3)

      Jac(2,1) = C(2,1)*F(1) + two*m%Yc*F(3)/F(2)
      Jac(2,2) = C(2,2)*F(1) - two*m%Yc*F(3)*F(1)/(F(2)*F(2))
      Jac(2,3) = C(2,3)*F(1) + two*m%Yc*F(1)/F(2)

      Jac(3,1) = C(3,1)*F(1)
      Jac(3,2) = C(3,2)*F(2)
      Jac(3,3) = C(3,3)*F(3)

      ! Calculate the new diagonal deformation gradient
      F = F - MATMUL(MInverse(Jac), Residual)

    End Do alphaLoop

    Call log%info('alpha0_DGD: ' // trim(str(alpha0_DGD)))

    Return
  End Function alpha0_DGD

  Function applyIdxRangeChecks(index, randomNumbers) result(output)
    ! Adjusts the index so that it is within bounds
    ! Handles large numbers and negative numbers

    Use forlog_Mod

    ! Arguments
    Integer, intent(IN) :: index
    Double Precision, intent(IN) :: randomNumbers(10000)

    ! Output
    Integer :: output

    ! Locals
    Integer :: rnsize
    Double Precision, parameter :: zero=0.d0
    ! -------------------------------------------------------------------- !

    ! Wrap in a loop if there are more rows than randomNumbers is long
    output = index
    rnsize = SIZE(randomNumbers, 1)
    If (output > rnsize) Then
      output = MOD(output, rnsize)
      Call log%warn('Random fiber misalignments are being wrapped in a loop. Increase randomNumberCount in vexternaldb.')
    End If

    ! In case the position is negative
    If (output < 0) Then
      output = rnsize+output
    End If

    Return
  End Function applyIdxRangeChecks


  Subroutine initializePhi0(m, Lc, position, phi012, phi013)
    ! Determines the phi0 for shear instability
    ! Calculates using polar and azimuth angles

    Use forlog_Mod

    ! Arguments
    Type(matProps), intent(IN) :: m
    ! Double Precision, intent(IN) :: G12, G13, XC                 ! Material properties
    ! Double Precision, intent(IN) :: aPL, nPL                     ! Shear nonlinearity
    Double Precision, intent(IN) :: Lc(3)                        ! Element characteristic lengths
    Double Precision, intent(IN) :: position(3)                  ! Position
    Double Precision, intent(INOUT) :: phi012, phi013            ! Return values


    ! Locals
    Double Precision :: shearInstabilityPhi0
    Double Precision :: randomNumbers1(10000), randomNumbers2(10000)
    Double Precision :: randomNumbers3(10000), randomNumbers4(10000)
    Double Precision :: rn1, rn2
    Double Precision :: phi0, azi
    Double Precision :: rad_to_deg, deg_to_rad
    Double Precision :: phi0_input
    Integer :: index, rnsize
    Double Precision, parameter :: zero=0.d0, half=0.5d0, one=1.d0, two=2.d0, three=3.d0

    ! Common
    Common randomNumbers1, randomNumbers2, randomNumbers3, randomNumbers4
    ! -------------------------------------------------------------------- !

    ! Initialization
    index = 0
    rad_to_deg = 45.d0/ATAN(one)  ! Converts radians to degrees when multiplied
    deg_to_rad = one/rad_to_deg   ! Converts degrees to radians when multiplied

    shearInstabilityPhi0 = (m%nPL-one)/m%G12*((m%G12-m%XC)/(m%XC*m%nPL*m%aPL**(one/m%nPL)))**(m%nPL/(m%nPL-one))

    ! Check for valid initial conditions for phi012 and phi013
    If (m%fiberCompDamFKT12 .AND. m%fiberCompDamFKT13) Then
      If (phi012 > half .AND. phi013 > half) Then
        If (ABS(phi012 - phi013) > 1d-8) Then
          Call log%error('Expect that the initial conditions for phi0_12 and phi0_13 are the same. Found phi012: ' // trim(str(phi012)) // ' and phi0_13: ' // trim(str(phi013)))
        End If
      Else If (phi012 < -1*half .OR. phi013 < -1*half) Then
        Call log%error('Unexpected initial conditions for phi0_12 and phi0_13. Phi012 and phi013 cannot be less than -0.5. Found phi012: ' // trim(str(phi012)) //   &
          ' and phi0_13: ' // trim(str(phi013)))
      Else If ((phi012 == zero .AND. phi013 /= zero) .OR. (phi012 /= zero .AND. phi013 == zero)) Then
        Call log%error('Unexpected initial conditions for phi0_12 and phi0_13. For 3-D kinking, phi012 and phi013 must be both zero or both nonzero. Found phi012: ' //  &
          trim(str(phi012)) // ' and phi0_13: ' // trim(str(phi013)))
      End If
      phi0_input = phi012
    Else If (m%fiberCompDamFKT12) Then
      phi0_input = phi012
    Else If (m%fiberCompDamFKT13) Then
      phi0_input = phi013
    Else
      phi0_input = zero
    End If

    If (phi0_input == zero) Then  ! Special case to use phi0 corresponding to instability
      If (m%fiberCompDamFKT12) Then
        phi012 = shearInstabilityPhi0
      End If
      If (m%fiberCompDamFKT13) Then
        phi013 = shearInstabilityPhi0
      End If

    Else If (phi0_input <= half) Then   ! Use initial condition value of phi0
      Return

    Else If (phi0_input == one) Then ! 1-D spatial variation
      index = NINT(position(1)/0.05d0)
      index = applyIdxRangeChecks(index, randomNumbers1)

      ! Get the random number (0 to 1)
      rn = randomNumbers1(index)
      rn2 = randomNumbers2(index)

      ! Polar angle, scale -phi to +phi
      phi0 = (rn*two-one)*shearInstabilityPhi0  ! polar angle

      ! Azimuth angle, scale -180 to 180
      azi = (rn2*two-one)*180.d0*deg_to_rad

      ! Phi012 and phi013
      phi012 = phi0*COS(azi)
      phi013 = phi0*SIN(azi)

    Else If (phi0_input == two) Then ! 2-D variation (different 1-D variation for each ply)
      index = NINT((position(1)/0.05d0) + (100.d0*position(3)/Lc(3)))
      index = applyIdxRangeChecks(index, randomNumbers1)

      ! Get the random number (0 to 1)
      rn = randomNumbers1(index)
      rn2 = randomNumbers2(index)

      ! Polar angle, scale -phi to +phi
      phi0 = (rn*two-one)*shearInstabilityPhi0  ! polar angle

      ! Azimuth angle, scale -180 to 180
      azi = (rn2*two-one)*180.d0*deg_to_rad

      ! Phi012 and phi013
      phi012 = phi0*COS(azi)
      phi013 = phi0*SIN(azi)

    ElseIf (phi0_input == three) Then
      index = NINT((position(1)/0.05d0) + (100.d0*position(3)/Lc(3)))
      index = applyIdxRangeChecks(index, randomNumbers1)

      phi0 = randomNumbers3(index)*deg_to_rad  ! polar angle
      azi = randomNumbers4(index)*deg_to_rad   ! Azimuth

      ! Phi012 and phi013
      phi012 = phi0*COS(azi)
      phi013 = phi0*SIN(azi)

    Else
      Call log%error('Invalid phi0. Found phi0 = '  // trim(str(phi0)) // ' See README.')

    End If

    Return
  End Subroutine initializePhi0


  Subroutine writeMaterialPropertiesToFile(fileUnit, m)
    ! Writes provided material properties to a file as a python dictionary
    ! Assumes that file opening and closing is handled elsewhere

    ! Arguments
    Integer, intent(IN) :: fileUnit
    Type(matProps), intent(IN) :: m

    ! Locals
    Character(len=32) :: nameValueFmt
    Double Precision, parameter :: zero=0.d0
    ! -------------------------------------------------------------------- !

    ! Defines the format for writing the floating point numbers
    nameValueFmt = "(A,E21.15E2,A)"

    ! Write the feature flags
    write(fileUnit,"(A)") 'featureFlags = {'
    write(fileUnit, "(A,A,A)") '    "integer": "', m%featureFlags, '",'
    If (m%matrixDam) Then
      write(fileUnit,"(A)") '    "matrixDam": True,'
    Else
      write(fileUnit,"(A)") '    "matrixDam": False,'
    End If
    If (m%shearNonlinearity12) Then
      write(fileUnit,"(A)") '    "shearNonlinearity12": True,'
    Else
      write(fileUnit,"(A)") '    "shearNonlinearity12": False,'
    End If
    If (m%shearNonlinearity13) Then
      write(fileUnit,"(A)") '    "shearNonlinearity13": True,'
    Else
      write(fileUnit,"(A)") '    "shearNonlinearity13": False,'
    End If
    If (m%schapery) Then
      write(fileUnit,"(A)") '    "schapery": True,'
    Else
      write(fileUnit,"(A)") '    "schapery": False,'
    End If
    If (m%schaefer) Then
      write(fileUnit,"(A)") '    "schaefer": True,'
    Else
      write(fileUnit,"(A)") '    "schaefer": False,'
    End If
    If (m%fiberTenDam) Then
      write(fileUnit,"(A)") '    "fiberTenDam": True,'
    Else
      write(fileUnit,"(A)") '    "fiberTenDam": False,'
    End If
    If (m%fiberCompDamBL) Then
      write(fileUnit,"(A)") '    "fiberCompDamBL": True,'
    Else
      write(fileUnit,"(A)") '    "fiberCompDamBL": False,'
    End If
    If (m%fiberCompDamFKT12) Then
      write(fileUnit,"(A)") '    "fiberCompDamFKT12": True,'
    Else
      write(fileUnit,"(A)") '    "fiberCompDamFKT12": False,'
    End If
    If (m%fiberCompDamFKT13) Then
      write(fileUnit,"(A)") '    "fiberCompDamFKT13": True,'
    Else
      write(fileUnit,"(A)") '    "fiberCompDamFKT13": False,'
    End If
    If (m%friction) Then
      write(fileUnit,"(A)") '    "friction": True'
    Else
      write(fileUnit,"(A)") '    "friction": False'
    End If
    write(fileUnit, "(A)") '}'

    ! Write the material properties
    write(fileUnit, "(A)") 'm = ['
    write(fileUnit, nameValueFmt) '    ', m%E1, ',  # 9: E1'
    write(fileUnit, nameValueFmt) '    ', m%E2, ',  # 10: E2'
    write(fileUnit, nameValueFmt) '    ', m%G12, ',  # 11: G12'
    write(fileUnit, nameValueFmt) '    ', m%v12, ',  # 12: v12'
    write(fileUnit, nameValueFmt) '    ', m%v23, ',  # 13: v23'
    write(fileUnit, nameValueFmt) '    ', m%YT, ',  # 14: YT'
    write(fileUnit, nameValueFmt) '    ', m%SL, ',  # 15: SL'
    write(fileUnit, nameValueFmt) '    ', m%GYT, ',  # 16: GYT'
    write(fileUnit, nameValueFmt) '    ', m%GSL, ',  # 17: GSL'
    write(fileUnit, nameValueFmt) '    ', m%eta_BK, ',  # 18: eta_BK'
    write(fileUnit, nameValueFmt) '    ', m%YC, ',  # 19: YC'
    write(fileUnit, nameValueFmt) '    ', m%alpha0, ',  # 20: alpha0'
    write(fileUnit, nameValueFmt) '    ', m%E3, ',  # 21: E3'
    write(fileUnit, nameValueFmt) '    ', m%G13, ',  # 22: G13'
    write(fileUnit, nameValueFmt) '    ', m%G23, ',  # 23: G23'
    write(fileUnit, nameValueFmt) '    ', m%v13, ',  # 24: v13'
    write(fileUnit, nameValueFmt) '    ', m%cte(1), ',  # 25: cte11'
    write(fileUnit, nameValueFmt) '    ', m%cte(2), ',  # 26: cte22'
    write(fileUnit, nameValueFmt) '    ', m%aPL, ',  # 27: aPL'
    write(fileUnit, nameValueFmt) '    ', m%nPL, ',  # 28: nPL'
    write(fileUnit, nameValueFmt) '    ', m%XT, ',  # 29: XT'
    write(fileUnit, nameValueFmt) '    ', m%fXT, ',  # 30: fXT'
    write(fileUnit, nameValueFmt) '    ', m%GXT, ',  # 31: GXT'
    write(fileUnit, nameValueFmt) '    ', m%fGXT, ',  # 32: fGXT'
    write(fileUnit, nameValueFmt) '    ', m%XC, ',  # 33: XC'
    write(fileUnit, nameValueFmt) '    ', m%fXC, ',  # 34: fXC'
    write(fileUnit, nameValueFmt) '    ', m%GXC, ',  # 35: GXC'
    write(fileUnit, nameValueFmt) '    ', m%fGXC, ',  # 36: fGXC'
    write(fileUnit, nameValueFmt) '    ', m%cl, ',  # 37: cl'
    write(fileUnit, nameValueFmt) '    ', m%w_kb, ',  # 38: w_kb'
    write(fileUnit, nameValueFmt) '    ', m%T_sf, ',  # 39: T_sf'
    write(fileUnit, nameValueFmt) '    ', m%mu, ',  # 40: mu'
    write(fileUnit, "(A)") ']'
    ! write(fileUnit, "(A,E22.15E2)") 'cte33 = ', m%cte(3)

  End Subroutine writeMaterialPropertiesToFile


End Module matProp_Mod
