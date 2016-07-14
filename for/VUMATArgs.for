
Module VUMATArg_Mod


  Type VUMATArg

    Integer :: nblock, ndir, nshr, nstatev
    Double Precision :: stepTime, totalTime, dt
    ! Double Precision :: strainInc(nblock,ndir+nshr)
    ! Double Precision :: tempOld(nblock), tempNew(nblock)
    ! Double Precision :: stretchOld(nblock,ndir+nshr), stretchNew(nblock,ndir+nshr)
    ! Double Precision :: defgradOld(nblock,ndir+nshr+nshr), defgradNew(nblock,ndir+nshr+nshr)
    ! Double Precision :: stressOld(nblock,ndir+nshr)
    ! Double Precision :: stateOld(nblock,nstatev)
    ! Double Precision :: enerInternOld(nblock)
    ! Double Precision :: enerInelasOld(nblock)
  Contains
    procedure :: init => VUMATArg_init

  End Type VUMATArg

Contains

  Pure Subroutine VUMATArg_init(a, nblock, ndir, nshr, nstatev, stepTime, totalTime, dt)

    ! Arguments
    Class(VUMATArg), intent(INOUT) :: a
    Integer, intent(IN) :: nblock, ndir, nshr, nstatev
    Double Precision, intent(IN) :: stepTime, totalTime, dt

    a%nblock = nblock
    a%ndir = ndir
    a%nshr = nshr
    a%nstatev = nstatev
    a%stepTime = stepTime
    a%totalTime = totalTime
    a%dt = dt

  End Subroutine VUMATArg_init


End Module VUMATArg_Mod
