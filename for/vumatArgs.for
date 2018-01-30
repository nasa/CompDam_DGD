Module vumatArg_Mod
  ! Module to load and store the general VUMAT parameters for passing of information to user modules

  Type vumatArg

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
    Integer :: nElement, nMatPoint, nLayer, nSecPoint
  Contains
    procedure :: init => vumatArg_init
    procedure :: update

  End Type vumatArg

Contains

  Pure Subroutine vumatArg_init(a, nblock, ndir, nshr, nstatev, stepTime, totalTime, dt)

    ! Arguments
    Class(vumatArg), intent(INOUT) :: a
    Integer, intent(IN) :: nblock, ndir, nshr, nstatev
    Double Precision, intent(IN) :: stepTime, totalTime, dt

    a%nblock = nblock
    a%ndir = ndir
    a%nshr = nshr
    a%nstatev = nstatev
    a%stepTime = stepTime
    a%totalTime = totalTime
    a%dt = dt

  End Subroutine vumatArg_init


  Pure Subroutine update(a, nElement, nMatPoint, nLayer, nSecPoint)

    ! Arguments
    Class(vumatArg), intent(INOUT) :: a
    Integer, intent(IN) :: nElement, nMatPoint, nLayer, nSecPoint

    a%nElement = nElement
    a%nMatPoint = nMatPoint
    a%nLayer = nLayer
    a%nSecPoint = nSecPoint

  End Subroutine update


End Module vumatArg_Mod
