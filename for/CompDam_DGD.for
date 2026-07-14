! This software may be used, reproduced, and provided to others only as permitted under the terms of the agreement
! under which it was acquired from the U.S. Government. Neither title to, nor ownership of, the software is hereby
! transferred. This notice shall remain on all copies of the software.
!
! Copyright 2016 United States Government as represented by the Administrator of the National Aeronautics and
! Space Administration. No copyright is claimed in the United States under Title 17, U.S. Code. All Other Rights
! Reserved.

#include "vumatArgs.for"
#include "version.for"
#include "forlog.for"
#include "matrixUtil.for"
#include "stress.for"
#include "matProp.for"
#include "stateVar.for"
#include "parameters.for"
#include "schapery.for"
#include "strain.for"
#include "schaefer.for"
#include "plasticity.for"
#include "fiberDamage.for"
#include "friction.for"
#include "cohesive.for"
#include "vucharlength.for"
#include "DGD.for"
#include "vexternaldb.for"


Subroutine VUMAT(  &
  ! Read only (unmodifiable) variables:
  jblock,ndir,nshr,nstatev,nfieldv,nprops,jInfoArray,stepTime,  &
  totalTime,dtArray,cmname,coordMp,charLength,props,density,    &
  strainInc,relSpinInc,tempOld,stretchOld,defgradOld,fieldOld,  &
  stressOld,stateOld,enerInternOld,enerInelasOld,tempNew,       &
  stretchNew,defgradNew,fieldNew,                               &
  ! Write only (modifiable) variables:
  stressNew,stateNew,enerInternNew,enerInelasNew)

  Include 'vaba_param.inc'

  Dimension jblock(*),props(*),density(*),coordMp(*),charLength(*),strainInc(*),relSpinInc(*),tempOld(*),tempNew(*), &
    stretchOld(*),defgradOld(*),fieldOld(*),stressOld(*),stateOld(*),enerInternOld(*),                               &
    enerInelasOld(*),stretchNew(*),defgradNew(*),fieldNew(*),stressNew(*),stateNew(*),                               &
    enerInternNew(*),enerInelasNew(*),dtArray(*),jInfoArray(*)

  Character(len=80) :: cmname

  Integer, parameter :: i_nblock=1, i_npt=2, i_layer=3, i_kspt=4, i_noel=5
  Integer, parameter :: i_info_AnnealFlag=1, i_info_effModDefn=5
  Dimension INTV(1), REALV(1)    ! For abaqus XPLB_ABQERR
  Character(len=8) CHARV(1)      ! For Abaqus XPLB_ABQERR
  iUpdateEffMod = jInfoArray(i_info_effModDefn)
  lanneal = jInfoArray(i_info_AnnealFlag)

  ! -------------------------------------------------------------------- !
  !    End VUMAT standard interface                                      !
  ! -------------------------------------------------------------------- !

  ! Throw error on single precision
  If (totalTime == 0 .AND. PRECISION(dt) < 15) Then
    Call XPLB_ABQERR(-3,'Found single precision data. CompDam must be run in double precision.',INTV,REALV,CHARV)

  Else
    Call CompDam(  &
      ! Read only (unmodifiable) variables:
      jblock(i_nblock),ndir,nshr,nstatev,nfieldv,nprops,lanneal,stepTime,     &
      totalTime,dtArray,cmname,coordMp,charLength,props,density,         &
      strainInc,relSpinInc,tempOld,stretchOld,defgradOld,fieldOld,  &
      stressOld,stateOld,enerInternOld,enerInelasOld,tempNew,       &
      stretchNew,defgradNew,fieldNew,                               &
      ! Write only (modifiable) variables:
      stressNew,stateNew,enerInternNew,enerInelasNew,               &
      ! Read only, extra arguments:
      jblock(i_noel), jblock(i_npt), jblock(i_layer), jblock(i_kspt), iUpdateEffMod)
  End If
End Subroutine VUMAT


Subroutine CompDam(  &
  ! Read only (unmodifiable) variables:
  nblock,ndir,nshr,nstatev,nfieldv,nprops,lanneal,stepTime,     &
  totalTime,dtArray,cmname,coordMp,charLength,props,density_abq,     &
  strainInc,relSpinInc,tempOld,stretchOld,defgradOld,fieldOld,  &
  stressOld,stateOld,enerInternOld,enerInelasOld,tempNew,       &
  stretchNew,defgradNew,fieldNew,                               &
  ! Write only (modifiable) variables:
  stressNew,stateNew,enerInternNew,enerInelasNew,               &
  ! Read only, extra arguments:
  nElement, nMatPoint, nLayer, nSecPoint, iUpdateEffMod)

  Use forlog_Mod
  Use matrixAlgUtil_Mod
  Use matProp_Mod
  Use stateVar_Mod
  Use parameters_Mod
  Use DGD_Mod
  Use plasticity_mod
  Use cohesive_mod
  Use friction_mod
  Use strain_mod  ! Only needed for packager
  Use stress_Mod  ! Only needed for packager

  Implicit Double Precision (a-h, o-z)
  Integer, parameter :: j_sys_Dimension = 2, maxblk = 512

  Double Precision :: props(nprops), density_abq(nblock), coordMp(nblock,*), charLength(nblock,*),       &
    dtArray(2*(nblock)+1),  &
    strainInc(nblock,ndir+nshr), relSpinInc(nblock,nshr), tempOld(nblock), tempNew(nblock),              &
    stretchOld(nblock,ndir+nshr), defgradOld(nblock,ndir+nshr+nshr), fieldOld(nblock,nfieldv),           &
    stressOld(nblock,ndir+nshr), stateOld(nblock,nstatev), enerInternOld(nblock), enerInelasOld(nblock), &
    stretchNew(nblock,ndir+nshr), defgradNew(nblock,ndir+nshr+nshr), fieldNew(nblock,nfieldv),           &
    stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev), enerInternNew(nblock), enerInelasNew(nblock)

  ! Extra arguments
  Integer :: nElement(nblock), nMatPoint, nLayer, nSecPoint, iUpdateEffMod

  Double Precision :: stepTime, totalTime, dt

  Character(len=80) :: cmname

  ! -------------------------------------------------------------------- !
  !    End VUMAT standard interface                                      !
  ! -------------------------------------------------------------------- !
  ! Damage variables, failure indices, damage thresholds
  Double Precision :: d1C_dummy, rfC_dummy

  ! Stress
  Double Precision :: Cauchy(3,3)
  Double Precision :: CauchyABQ(3,3)
  Double Precision :: T_coh(3), trial_traction_int(3), trial_traction_ele(3)  ! Cohesive tractions (shear, normal, shear)
  Double Precision :: Rot(3,3)
  Double Precision :: stiff(ndir+nshr,ndir+nshr), eps(3,3), stress(3,3)  ! Used for the packager

  ! Other
  Double Precision :: F(3,3)      ! Current deformation gradient tensor
  Double Precision :: F_old(3,3)  ! Previous deformation gradient tensor
  Double Precision :: U(3,3)      ! Current stretch tensor
  Double Precision :: density_current  ! Current density calculated from Fbulk (ignores the cohesive crack displacement)
  Logical :: Sliding
  
  ! Cohesive law
  Double Precision :: Pen(3), Pen_int(3)  ! Cohesive penalty stiffness
  Double Precision :: delta(3), delta_int(3)  ! Cohesive displacement-jump vector
  Double Precision :: dmg_penalty  ! Penalty-based damage variable
  Double Precision :: dGdGc, dGdGc_inc  ! Change in normalized energy dissipation
  Double Precision :: constit_thk   ! Cohesive element constitutive thickness
  Double Precision :: tar_mass, tar_dt  ! Target mass and time increment when effmod + constit_thk != 1

  ! Parameters
  Double Precision, parameter :: zero=0.d0, one=1.d0, two=2.d0, three=3.d0
  Double Precision, parameter :: COH_max_dmg_effmod=0.9999d0
  Integer, parameter :: max_num_materials=100
  Integer, parameter :: max_num_cpus=512
  Integer, parameter :: mat_name_struc_id=3000
  Integer, parameter :: mat_obj_struc_id=2000

  ! Structure of all VUMAT arguments
  Type(vumatArg) :: args

  ! For access to the logger
  Type(forlog) :: logger

  ! Material
  Logical :: fiber_tension_snapback, fiber_compression_snapback
  Logical :: is_new_material, found_mat
  Type(matProps) :: m, mat_prop, mat_props(max_num_materials*max_num_cpus)
  pointer(ptr_props, mat_props)

  Type(matNameList) :: mat_name, mat_names(max_num_materials)
  pointer(ptr_mat_names, mat_names)

  ! State variables
  Type(stateVars) :: sv

  ! Solution parameters
  Type(parameters) :: default_params    ! Default parameters
  Type(parameters) :: p                 ! Parameters for execution
  Type(parameters) :: params            ! Storage for parameters loaded from inp deck
  pointer(ptr_params, params)

  ! ! Compdam version printing
  ! Integer :: printed_compdam_version
  ! pointer(ptr_printed_compdam_version, printed_compdam_version)

  ! Fatigue
  Logical :: fatigue_step
  Integer :: fatigue_parameters(2)
  Pointer(ptr_fatigue_int, fatigue_parameters)

  ! Convenience
  Logical :: in_packager, in_1st_inc, in_packager_or_1st_inc, update_effmod
  Logical :: in_first_eval_cur_mat  ! First time the subroutine is called for the current material
  integer NUMPROCESSES, KPROCESSNUM
  ! Integer :: process_rank  ! Process number or rank

  Dimension INTV(1), REALV(1)    ! For abaqus warning messages
  Character(len=8) CHARV(1)      ! For Abaqus warning messages

  ! See: "SIMULIA/EstProducts/2024/SMAUsubs/PublicInterfaces/SMAAspUserArrays.hdr"
  !      "SIMULIA/EstProducts/2024/SMAUsubs/PublicInterfaces/SMAAspUserUtilities.hdr"
  Interface

      Function SMAIntArrayAccess(ID)
        ! Access an existing global integer array.
        Integer(kind=8) :: SMAIntArrayAccess  ! Returns an address that can be associated with a Fortran pointer
        Integer(kind=4) :: ID       ! Array ID
      End Function SMAIntArrayAccess

      FUNCTION SMAStructArrayAccess(ID)
        INTEGER(KIND=8) :: SMAStructArrayAccess  ! -- Returns an address that can be associated with a Fortran pointer
        INTEGER(KIND=4) :: ID                    ! Array ID
      END FUNCTION SMAStructArrayAccess

      SUBROUTINE MutexInit ( ID )
         INTEGER(KIND=4) :: ID
      END SUBROUTINE MutexInit

      ! Lock Mutex -- can be called anywhere after initialization
      SUBROUTINE MutexLock ( ID )
         INTEGER(KIND=4) :: ID
      END SUBROUTINE MutexLock

      ! Unlock Mutex -- can be called anywhere after initialization
      SUBROUTINE MutexUnlock ( ID )
         INTEGER(KIND=4) :: ID
      END SUBROUTINE MutexUnlock

      ! integer(kind=8) FUNCTION SMALocalIntArrayAccess(ID)
      !     INTEGER(KIND=4),INTENT(IN) :: ID
      ! END FUNCTION SMALocalIntArrayAccess

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

  ! INTERFACE SMALocalIntArrayCreate

  !   ! -- Creates an array with a given ID and SIZE; initialized to supplied INITVAL
  !   integer(kind=8) FUNCTION SMALocalIntArrayCreateInit(ID,SIZE,INITVAL)
  !       INTEGER(KIND=4),INTENT(IN) :: ID
  !       INTEGER(KIND=4),INTENT(IN) :: SIZE
  !       INTEGER(KIND=4),INTENT(IN) :: INITVAL
  !   END FUNCTION SMALocalIntArrayCreateInit

  !   ! -- Creates an array with a given ID and SIZE; initialized implicitly to INT_MAX
  !   integer(kind=8) FUNCTION SMALocalIntArrayCreateNoInit(ID,SIZE)
  !       INTEGER(KIND=4),INTENT(IN) :: ID
  !       INTEGER(KIND=4),INTENT(IN) :: SIZE
  !   END FUNCTION SMALocalIntArrayCreateNoInit

  ! END INTERFACE SMALocalIntArrayCreate
  ! -------------------------------------------------------------------- !

  ! Initializations for multiple processors
  Call MutexInit( 1 )
  Call MutexInit( 2 )
  Call MutexInit( 3 )
  CALL VGETNUMCPUS( NUMPROCESSES )
  CALL VGETRANK( KPROCESSNUM )

  ! EFFMOD
  ! The main use of EFFMOD is prevent damaged cohesive elements from making the stable time increment
  ! reduce thereby slowing the simulation results.
  ! EFFMOD is applicable to zero-thickness cohesive elements; the constitutive thickness can be modified
  ! from the default value of 1.
  ! Note that EFFMOD is currently only implemented for cohesive elements, implementation for
  ! continuum elements is future work. Such implementation can serve the same purpose, but must account
  ! for the fact that stiffness degradation is anisotropic.
  dt = dtArray(1)
  update_effmod = iUpdateEffMod .EQ. 1

  ! Initialize structure of VUMAT args
  Call args%init(nblock, ndir, nshr, nstatev, stepTime, totalTime, dt, cmname)

  in_packager = totalTime == 0
  in_1st_inc = stepTime <= dt .AND. .NOT. in_packager
  in_packager_or_1st_inc = stepTime <= dt .OR. in_packager
  in_first_eval_cur_mat = .FALSE.

  ! ! Print CompDam version
  ! If (in_packager) Then
  !   Call VGETRANK(process_rank)
  !   If (process_rank == 0) Then
  !     ptr_printed_compdam_version = SMALocalIntArrayAccess(1)
  !     If (ptr_printed_compdam_version .EQ. 0) ptr_printed_compdam_version = SMALocalIntArrayCreate(1, 1, 0)
  !     If (printed_compdam_version .EQ. 0) Then
  !       Call PrintVersion()
  !       printed_compdam_version = 1
  !     End If
  !   End If
  ! End If

  ! Load the solution parameters
  ! Use default parameters during the packager
  If (in_packager) Then
    p = default_params
  Else
    ptr_params = SMAStructArrayAccess(1000)
    If (ptr_params .EQ. 0) then
      Call XPLB_ABQERR(-2,'COMPDAM ERROR ptr_params does not exist',INTV,REALV,CHARV)
    End If
    p = params
  End If

  ! Initialize the logger
  Call log%init(p%logLevel, args, p%logFormat, NUMPROCESSES, KPROCESSNUM)

  ! Process material properties in the packager
  ! Keep track of which materials have been processed to avoid redundant logging
  If (in_packager) Then
    ! Initialize mat_names array to track the materials that have been processed
    Call MutexLock( 1 )  ! Lock to avoid multiple processors writing (is_new_material) at the same time
    ptr_mat_names = SMAStructArrayAccess(mat_name_struc_id)
    If (ptr_mat_names .EQ. 0) ptr_mat_names = SMAStructArrayCreate(mat_name_struc_id, max_num_materials, sizeof(mat_name), mat_name)

    ! Check if material has already been processed
    is_new_material = .TRUE.
    MatL1: Do II=1,max_num_materials
      If (mat_names(II)%name == cmname) Then
        is_new_material = .FALSE.
        Exit MatL1
      End If
    End Do MatL1

    ! Load the material as m
    Call loadMatProps(m, cmname, nprops, props, is_new_material)

    If (is_new_material) Then
      ! Add material name to mat_names array since it has now been processed
      MatL2: Do II=1,max_num_materials
        If (trim(mat_names(II)%name) == "") Then
          mat_names(II)%name = cmname
          Exit MatL2
        End If
      End Do MatL2
    End If
    Call MutexUnlock( 1 )
  End If

  ! Process material properties during the first increment of the solver
  ! Load the materials into a shared array to avoid redundant processing
  ! Each processor has its own copy of the materials in its domain
  If (in_1st_inc) Then
    Call MutexLock( 2 )
    ptr_props = SMAStructArrayAccess(mat_obj_struc_id)
    If (ptr_props .EQ. 0) Then
      ptr_props = SMAStructArrayCreate(mat_obj_struc_id, max_num_materials*max_num_cpus, sizeof(mat_prop), mat_prop)
    End If
    Call MutexUnlock( 2 )
    MatL3: Do II = max_num_materials*KPROCESSNUM+1, max_num_materials*(KPROCESSNUM+1)
      If (trim(mat_props(II)%name) == cmname) Then
        ! Material is already processed and stored, load
        m = mat_props(II)
        Exit MatL3
      End If
      If (trim(mat_props(II)%name) == "") Then
        ! Material is not stored, need to read, process, and store
        Call loadMatProps(m, cmname, nprops, props, .FALSE.)
        mat_props(II) = m

        Exit MatL3
      End If

    End Do MatL3
  End If

  ! Load the material from the shared array to avoid redundant processing
  If (.NOT. in_packager_or_1st_inc) Then
    found_mat = .FALSE.
    ptr_props = SMAStructArrayAccess(mat_obj_struc_id)
    MatL4: Do II = max_num_materials*KPROCESSNUM+1, max_num_materials*(KPROCESSNUM+1)
      If (trim(mat_props(II)%name) == cmname) Then
        m = mat_props(II)
        found_mat = .TRUE.
        Exit MatL4
      End If
    End Do MatL4
    If (.NOT. found_mat) Then
      Call log%error('Material properties for '//trim(cmname)//' not found in shared array.')
    End If
  End If

  ! This forces the material to behave elastically during the packager and initial increment. This
  ! was motivated by issues with Intel MKL in the packager.
  If (in_packager) Then
    m%matrixDam = .FALSE.
    m%shearNonlinearity12 = .FALSE.
    m%shearNonlinearity13 = .FALSE.
    m%schapery = .FALSE.
    m%schaefer = .FALSE.
    m%fiberTenDam = .FALSE.
    m%fiberCompDamBL = .FALSE.
    m%fiberCompDamFKT12 = .FALSE.
    m%fiberCompDamFKT13 = .FALSE.
    m%friction = .FALSE.
    fatigue_step = .FALSE.
  Else
    fatigue_step = .FALSE.
#ifndef PYEXT
    ! Is this a fatigue step?
    ptr_fatigue_int = SMAIntArrayAccess(1)
    If (fatigue_parameters(1) == 1) Then
      fatigue_step = .TRUE.
      Call log%debug("This is a fatigue step.")
    End If
#endif
  End If

  ! -------------------------------------------------------------------- !
  master: Do km = 1,nblock  ! Master Loop

  ! Update location information for debugging
  Call log%arg%update(nElement(km), nMatPoint, nLayer, nSecPoint)

  Call log%debug("Master loop. km = " // trim(str(km)))

  ! -------------------------------------------------------------------- !
  !    Recall previous elastic limits, plasticity, and damage states:    !
  ! -------------------------------------------------------------------- !
  sv = loadStateVars(nstatev, stateOld(km,:), m)

  ! -------------------------------------------------------------------- !
  !    When STATUS=0 integration point is slated for deletion            !
  ! -------------------------------------------------------------------- !
  If (sv%STATUS == 0) Then
    If (in_packager) Then
      Call log%warn('Found STATUS=0 for element '//trim(str(nElement(km)))//' in material '//trim(cmname)//'; CompDam should be initialized with STATUS=1 using initial conditions')
    Else
      Cycle master
    End If
  End If

  ! -------------------------------------------------------------------- !
  !    Cohesive elements:                                                !
  ! -------------------------------------------------------------------- !
  ElementType: If (m%cohesive) Then

    constit_thk = one  ! In most cases, constitutive thickness should be 1 (strain = displacement)
    If (m%embedded_cohesive) Then
      If (in_first_eval_cur_mat .AND. m%thickness > zero) Then
        Call log%info('Thick cohesive has a thickness specified in the material card; this values is ignored')
      End If
      If (in_first_eval_cur_mat .AND. charLength(km,1) .NE. one) Then
        Call log%error('The thickness specified in *cohesive section must be 1, found: ' // trim(str(charLength(km,1))))
      End If
      If (in_first_eval_cur_mat .AND. update_effmod) Then
        Call log%error('Finite thickness cohesive elements do not support EFFMOD')
      End If
    Else ! zero thickness
      If (update_effmod) Then
        ! Allows for proper scaling if constitutive thickness is not = 1,
        ! assuming user has provided a thickness in the *Material defintion ...
        If (m%thickness > zero) Then
          constit_thk = m%thickness
        End If
      Else
        If (in_first_eval_cur_mat .AND. charLength(km,1) .NE. one) Then
          Call log%error('The thickness specified in *cohesive section must be 1, found: ' // trim(str(charLength(km,1))))
        End If
        If (in_first_eval_cur_mat .AND. m%thickness > zero) Then
          Call log%warn('Zero thickness cohesive without effmod enabled has a thickness specified in the material card; this values is ignored')
        End If
      End If
    End If

    ! Normal direction
    delta(2) = sv%Fb2 + strainInc(km,1)*constit_thk 
    If (nshr == 2) Then
      delta(1) = sv%Fb1 + two*strainInc(km,3)*constit_thk  ! Longitudinal shear direction (1--3)
      delta(3) = sv%Fb3 + two*strainInc(km,2)*constit_thk  ! Transverse shear direction (2--3)
    Else
      delta(1) = sv%Fb1 + two*strainInc(km,2)*constit_thk  ! Shear direction (1--2)
      delta(3) = zero
    End If
    
    ! Store current cohesive displacement-jump
    sv%Fb1 = delta(1)  ! Longitudinal shear displacement-jump
    sv%Fb2 = delta(2)  ! Normal displacement-jump
    sv%Fb3 = delta(3)  ! Transverse shear displacement-jump
    
    ! Define the penalty stiffnesses
    If (m%embedded_cohesive) Then
      Pen(1) = m%G13
      Pen(2) = m%E3
      Pen(3) = m%G23
    Else
      ! Define penalty stiffnesses
      ! Normal (Mode I) penalty stiffness. This is a user-input numerical constant.
      Pen(2) = m%E3
      ! Longitudinal (Pen(1)) and transverse (Pen(3)) shear penalty stiffnesses
      !  Defined according to the material property relationships in Turon (2010). Also includes the normal
      !  compression load dependence of the shear strengths from LaRC04.
      If (in_packager) Then
        Pen(1) = Pen(2)*m%GYT*m%SL*m%SL/(m%GSL*m%YT*m%YT)  ! Mode II penalty stiffness, Turon (2010)
        Pen(3) = Pen(2)*m%GYT*m%ST*m%ST/(m%GSL*m%YT*m%YT)
      Else
        Pen(1) = Pen(2)*m%GYT*(m%SL - m%etaL*MAX(-m%YC, Pen(2)*MIN(zero, delta(2))))**2/(m%GSL*m%YT**2)
        Pen(3) = Pen(2)*m%GYT*(m%ST - m%etaT*MAX(-m%YC, Pen(2)*MIN(zero, delta(2))))**2/(m%GSL*m%YT**2)
      End If
    End If
    
    If (in_packager) Then
      ! Avoids calculating damage during the packager and the first solution increment
      dmg_penalty = zero
      dGdGc = zero
      T_coh = cohesive_traction(delta, Pen, dmg_penalty)
      ! User feedback when constit_thick != 1
      If (in_first_eval_cur_mat .AND. constit_thk .NE. one) Then
        If (ABS(m%thickness - charLength(km,1))/m%thickness > 1.d-2) Then
          Call log%error('m%thickness is more than 1% different from cohesive section thickness')
        End If
        tar_dt = m%thickness*SQRT(density_abq(km)/((Pen(1) + Pen(3))/2 + Pen(2)*0.75d0))
        Call log%info("User-provided constitutive thickness; target dt = " // trim(str(tar_dt)))
        Call log%warn("Use constitutive thickness != 1 with caution; numerical instability can occur as evident by reaction force oscillations starting at low loads")
      End If
    
    Else If (m%embedded_cohesive) Then  ! Finite-thickness cohesive law with physical stiffness inputs
      ! Separate the penalty stiffness into bulk and interface penalty stiffness, a la springs in series
      ! Together, Pen(2) and Pen_int(2) have the same stiffness as the user-defined normal penalty stiffness
      Pen_int(2) = Pen(2) * p%penStiffMult
      Pen(2)     = Pen(2) * p%penStiffMult / (p%penStiffMult - one)
      
      ! Penalty-based damage variables from previous increment
      dmg_penalty = damage_area_to_penalty(sv%d2, delta_initiation(m%YT, Pen_int(2)), delta_final(m%YT, m%GYT))
      
      dGdGc = zero  ! Initialize incremental dissipated plastic energy
      
      embedded_cohesive_tolerance = (1.d-8*m%YT)**2  ! Tolerance for convergence of the EmbeddedCohesive loop
      EmbeddedCohesive: Do
      
        delta_int(2) = interface_displacement_jumps(delta(2), Pen(2), dmg_penalty, Pen_int(2), .TRUE.)
          
        Pen_int(1) = Pen_int(2)*m%GYT*(m%SL - m%etaL*MAX(-m%YC, Pen_int(2)*MIN(zero, delta_int(2))))**2/(m%GSL*m%YT**2)
        Pen_int(3) = Pen_int(2)*m%GYT*(m%ST - m%etaT*MAX(-m%YC, Pen_int(2)*MIN(zero, delta_int(2))))**2/(m%GSL*m%YT**2)
        
        delta_int(1) = interface_displacement_jumps(delta(1), Pen(1), dmg_penalty, Pen_int(1), .FALSE.)
        delta_int(3) = interface_displacement_jumps(delta(3), Pen(3), dmg_penalty, Pen_int(3), .FALSE.)
      
        Call cohesive_damage(m, p, delta_int, Pen_int, delta_int(2), sv%B, sv%FIm, fatigue_step, sv%d2, dmg_penalty, dGdGc_inc)
        dGdGc = dGdGc + dGdGc_inc

        ! Calculate the error between the interface traction and the bulk traction, in units of traction squared
        embedded_cohesive_error = zero
        trial_traction_int = cohesive_traction(delta_int, Pen_int, dmg_penalty)
        trial_traction_ele = cohesive_traction(delta - delta_int, Pen, zero)
        Do E = 1,3
          embedded_cohesive_error = embedded_cohesive_error + (trial_traction_int(E) - trial_traction_ele(E))**2
        End Do

        If (embedded_cohesive_error < embedded_cohesive_tolerance) Exit EmbeddedCohesive
      End Do EmbeddedCohesive
      
      If (m%friction .AND. delta(2) <= zero) Then  ! Closed cracks with friction
        Sliding = crack_is_sliding(delta_int, Pen_int, sv%slide, m%mu, m%mu)
        Call crack_traction_and_slip(delta_int, Pen_int, sv%slide, sv%slide, m%mu, m%mu, dmg_penalty, sv%d2, T_coh, Sliding)
      Else  ! Closed cracks without friction and open cracks
        T_coh = cohesive_traction(delta_int, Pen_int, dmg_penalty)
        sv%slide(1) = delta_int(1)
        sv%slide(2) = delta_int(3)
      End If
      
    Else  ! Zero-thickness cohesive law with numerical penalty stiffnesses
      Call cohesive_damage(m, p, delta, Pen, delta(2), sv%B, sv%FIm, fatigue_step, sv%d2, dmg_penalty, dGdGc)
      
      If (m%friction .AND. delta(2) <= zero) Then  ! Closed cracks with friction
        Sliding = crack_is_sliding(delta, Pen, sv%slide, m%mu, m%mu)
        Call crack_traction_and_slip(delta, Pen, sv%slide, sv%slide, m%mu, m%mu, dmg_penalty, sv%d2, T_coh, Sliding)
      Else  ! Closed cracks without friction and open cracks
        T_coh = cohesive_traction(delta, Pen, dmg_penalty)
        sv%slide(1) = delta(1)
        sv%slide(2) = delta(3)
      End If
    End If

    ! Update stress values
    stressNew(km,1) = T_coh(2)  ! Normal stress
    If (nshr == 2) Then  ! 3-D
      stressNew(km,3) = T_coh(1)  ! Longitudinal shear stress (1--3)
      stressNew(km,2) = T_coh(3)  ! Transverse shear stress (2--3)
    Else  ! 2-D
      stressNew(km,2) = T_coh(1)  ! Shear stress
    End If

    ! Element deletion
    If (p%set_status_0_on_d2 == 1 .AND. sv%d2 >= 1) Then
      sv%STATUS = 0
    End If

    ! EFFMOD for time increment calculation
    If (update_effmod) Then
      dtArray(2:km+1) = (Pen(1) + Pen(3) + three/two*Pen(2))/5.d0*(1-min(dmg_penalty,COH_max_dmg_effmod))
      dtArray(km+2:2*(km)+1) = dtArray(2:km+1)
    End If

    ! Update internal and inelastic energy terms
    enerInelasNew(km) = dGdGc*(m%GYT + (m%GSL - m%GYT)*sv%B**m%eta_BK)
    enerInelasNew(km) = enerInelasOld(km) + enerInelasNew(km)/density_abq(km)/constit_thk
    enerInternNew(km) = DOT_PRODUCT(T_coh, delta/constit_thk)/two/density_abq(km)
    enerInternNew(km) = enerInternNew(km) + enerInelasNew(km)

  ! -------------------------------------------------------------------- !
  !    Solid elements:                                                   !
  ! -------------------------------------------------------------------- !
  Else ElementType
  
    ! -------------------------------------------------------------------- !
    !    Deformation Gradient Tensor and Right Stretch Tensor              !
    ! -------------------------------------------------------------------- !
    U = Vec2Matrix(stretchNew(km,:))
    F = Vec2Matrix(defgradNew(km,:))
    F_old = Vec2Matrix(defgradOld(km,:))

    ! -------------------------------------------------------------------- !
    ! As of Abaqus 6.16, the packager receives a defGradNew of (0.999, 0.999, 0.0, 0.001, 0.001)
    ! for S4R elements. The F(3,3) of 0.0 breaks the initial pass through the VUMAT and the model
    ! will not run. The following statement is a workaround to this problem.
    If (in_packager .AND. nshr == 1) F(3,3) = one

    ! -------------------------------------------------------------------- !
    !    Define the characteristic element lengths                         !
    ! -------------------------------------------------------------------- !
    If (sv%Lc(2) == zero) Then

      sv%Lc(1) = charLength(km, 1)
      sv%Lc(2) = charLength(km, 2)
      If (nshr == 1) Then
        sv%Lc(3) = m%thickness
      Else
        sv%Lc(3) = charLength(km, 3)
      End If
      
      If (in_1st_inc .AND. p%check_for_snap_back) Then
        Call checkForSnapBack(m, sv%Lc, nElement(km), fiber_tension_snapback, fiber_compression_snapback)
        
        If (fiber_tension_snapback .OR. fiber_compression_snapback) Then
          If (fiber_tension_snapback) m%fGXT = m%fXT
          If (fiber_compression_snapback) m%fGXC = m%fXC
        End If
      
      End If

      ! Check for strange values of Lc2 (may indicate that vucharlength did not run)
      If (in_1st_inc) Then ! Don't run in packager
        If (sv%Lc(2) .EQ. zero) Then
          Call log%error("Found Lc2 = 0 at element " // trim(str(nElement(km))) // ", material " // trim(m%name) // ". Perhaps *Characteristic Length is missing from the input deck?")
        Else If (sv%Lc(3) .EQ. zero) Then
          Call log%error("Found Lc3 = 0 at element " // trim(str(nElement(km))) // ", material " // trim(m%name) // ". Perhaps *Characteristic Length is missing from the input deck?")
        ! Else If (sv%Lc(2) > sv%Lc(1)*10.d0 .OR. sv%Lc(2) < sv%Lc(1)/10.d0) Then
        !   Call log%warn("Found Lc = [" // trim(str(sv%Lc(1))) //','// trim(str(sv%Lc(2))) //','// trim(str(sv%Lc(3))) // "]; perhaps *Characteristic Length is missing from the input deck?")
        End If
      End If

      Call log%debug("Characteristic element lengths:")
      Call log%debug(trim(str(sv%Lc(1)))//' '//trim(str(sv%Lc(2)))//' '//trim(str(sv%Lc(3))))

    End If

    If (in_packager_or_1st_inc) Then
      sv%debugpy_count = 0  ! Initialize
    End If

    ! -------------------------------------------------------------------- !
    !    Initialize phi0                                                   !
    ! -------------------------------------------------------------------- !
    If (in_packager_or_1st_inc .AND. (m%fiberCompDamFKT12 .OR. m%fiberCompDamFKT13)) Then
      Call initializePhi0(m, sv%Lc, charLength(km, 4:6), sv%phi0_12, sv%phi0_13)
      Call log%debug("Calculated initial misalignments, phi0_12: " // trim(str(sv%phi0_12)) // " and phi0_13: " // trim(str(sv%phi0_13)))
    End If

    ! -------------------------------------------------------------------- !
    !    Initialize inelastic energy to the previous value                 !
    ! -------------------------------------------------------------------- !
    enerInelasNew(km) = enerInelasOld(km)

    ! -------------------------------------------------------------------- !
    !    Initialize fiber failure                                          !
    ! -------------------------------------------------------------------- !
    If (m%fiberCompDamFKT12 .AND. p%fkt_fiber_failure_angle > zero) Then
      sv%Inel12c = intializeFiberFailure(sv%phi0_12, p%fkt_fiber_failure_angle, m%G12, m%aPL, m%nPL)
    Else
      sv%Inel12c = Huge(zero)   ! Turn off fiber failure by setting the associate inelastic strain to a very large number
    End If
    If (m%fiberCompDamFKT13 .AND. p%fkt_fiber_failure_angle > zero) Then
      sv%Inel13c = intializeFiberFailure(sv%phi0_13, p%fkt_fiber_failure_angle, m%G13, m%aPL, m%nPL)
    Else
      sv%Inel13c = Huge(zero)   ! Turn off fiber failure by setting the associate inelastic strain to a very large number
    End If

    ! -------------------------------------------------------------------- !
    !    Damage Calculations:                                              !
    ! -------------------------------------------------------------------- !
    ! Damage initiation prediction
    If (in_packager) Then
      Call MutexLock(3)
      stiff = StiffFunc(ndir+nshr, m%E1, m%E2, m%E3, m%G12, m%G13, m%G23, m%v12, m%v13, m%v23, zero, zero, zero)
      eps = GLStrain(F,ndir)
      stress = Hooke(stiff, eps, nshr)
      Cauchy = convertToCauchy(stress, F)

    Else If (.NOT. (m%matrixDam .AND. sv%d2 > zero) .AND. .NOT. ((m%fiberCompDamFKT12 .OR. m%fiberCompDamFKT13) .AND. sv%d1C > zero)) Then

      Call DGDInit(U,F,F_old,m,p,sv,ndir,nshr,tempNew(km)-m%T_sf,density_abq(km),Cauchy,enerInternNew(km),enerInelasNew(km),fatigue_step)

    End If

    ! Matrix crack damage evolution
    If (m%matrixDam .AND. sv%d2 > zero) Then

      Call DGDEvolve(U,F,F_old,m,p,sv,ndir,nshr,tempNew(km)-m%T_sf,density_abq(km),Cauchy,enerInternNew(km),enerInelasNew(km),fatigue_step)

    ! Fiber compression damage evolution (FKT decomposition)
    Else If ((m%fiberCompDamFKT12 .OR. m%fiberCompDamFKT13) .AND. sv%d1C > zero) Then

      Call DGDKinkband(U,F,F_old,m,p,sv,ndir,nshr,tempNew(km)-m%T_sf,density_abq(km),Cauchy,enerInternNew(km),enerInelasNew(km))

    End If

    ! -------------------------------------------------------------------- !
    !    Determine and store the stress tensor:                            !
    ! -------------------------------------------------------------------- !
    ! Rotation tensor
    Rot = MATMUL(F, MInverse(U))
    ! Cauchy stress in the current configuration
    CauchyABQ = MATMUL(TRANSPOSE(Rot), MATMUL(Cauchy, Rot))
    ! Convert to vector format
    stressNew(km,:) = Matrix2Vec(CauchyABQ, nshr)

    ! Check for zero or negative dilatational modulus
    If (in_packager) Then
      If (ABS(CauchyABQ(1,1) + CauchyABQ(2,2) + CauchyABQ(3,3)) < 1.d-8) Then
        Call log%warn('Dilatational modulus is close to zero during the packager.')
      End If
      If (CauchyABQ(1,1) > zero .AND. F(1,1) < one) Then
        Call log%warn('Dilatational modulus 11 is negative during the packager.')
      End If
      If (CauchyABQ(2,2) > zero .AND. F(2,2) < one) Then
        Call log%warn('Dilatational modulus 22 is negative during the packager.')
      End If
      If (CauchyABQ(3,3) > zero .AND. F(3,3) < one) Then
        Call log%warn('Dilatational modulus 33 is negative during the packager.')
      End If
      Call MutexUnlock(3)
    End If

    ! Element deletion
    If (p%set_status_0_on_d2 == 1 .AND. sv%d2 >= 1) Then
      sv%STATUS = 0
    End If
    If (p%set_status_0_on_d1C == 1 .AND. sv%d1C >= 1) Then
      sv%STATUS = 0
    End If
    If (p%set_status_0_on_d1T == 1 .AND. sv%d1T >= 1) Then
      sv%STATUS = 0
    End If

    ! -------------------------------------------------------------------- !
    !    Excessive shear strain errors                                      !
    ! -------------------------------------------------------------------- !
    If (m%shearNonlinearity12) Then
      If (sv%Inel12 < stateOld(km,13)) Then
        Call writeDGDArgsToFile(m,p,sv,U,F,F_old,ndir,nshr,tempNew(km)-m%T_sf,density_abq(km),log%arg,'CompDam_DGD')
        Call log%terminate('Decrease in inelastic strain 12.')
      End If
    End If
    If (m%shearNonlinearity13) Then
      If (sv%Inel13 < stateOld(km,21)) Then
        Call writeDGDArgsToFile(m,p,sv,U,F,F_old,ndir,nshr,tempNew(km)-m%T_sf,density_abq(km),log%arg,'CompDam_DGD')
        Call log%terminate('Decrease in inelastic strain 13.')
      End If
    End If
    If (m%shearNonlinearity12 .OR. m%shearNonlinearity13) Then
      If (sv%Inel12 > 1.d0 .OR. sv%Inel13 > 1.d0) Then
        Call writeDGDArgsToFile(m,p,sv,U,F,F_old,ndir,nshr,tempNew(km)-m%T_sf,density_abq(km),log%arg,'CompDam_DGD')
        Call log%terminate('Excessive inelastic shear strain.')
      End If
    End If

  End If ElementType

  ! -------------------------------------------------------------------- !
  !    Store the updated state variables:                                !
  ! -------------------------------------------------------------------- !
  stateNew(km,:) = storeStateVars(sv, nstatev, m)

  ! -------------------------------------------------------------------- !
  !    Kill job for debugging purposes                                   !
  ! -------------------------------------------------------------------- !
  If (p%logLevel > 2 .AND. p%debug_kill_at_total_time > zero) Then
    If (log%arg%totalTime >= p%debug_kill_at_total_time) Then
      Call writeDGDArgsToFile(m,p,sv,U,F,F_old,ndir,nshr,tempNew(km)-m%T_sf,density_abq(km),log%arg,'CompDam_DGD')
      Call log%error('debug_kill_at_total_time condition satisfied; terminating job.')
    End If
  End If

  End Do master ! End Master Loop

  Call log%debug('End of VUMAT')

  Return
End Subroutine CompDam

! Subroutine PrintVersion()
!   Use version_Mod
!   print *, '================= CompDam ================='
!   print *, 'Date: ' // trim(timestamp)
!   print *, 'Version: ' // trim(hash)
!   print *, '==========================================='
! End Subroutine PrintVersion
