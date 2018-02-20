Subroutine vucharlength(                                    &
  ! Read only (unmodifiable) variables:
  nblock, nfieldv, nprops, ncomp, ndim, nnode, nstatev,     &
  kSecPt, kLayer, kIntPt, jElType, jElem,                   &
  totalTime, stepTime, dt,                                  &
  cmname, coordMp, coordNode, direct, T, props,             &
  field, stateOld,                                          &
  ! Write only (modifiable) variables:
  charLength )

  Implicit Double Precision (a-h, o-z)
  Integer, parameter :: j_sys_Dimension = 2, n_vec_Length = 136, maxblk = n_vec_Length

  Dimension jElType(3), jElem(nblock), coordMp(nblock,ndim),  &
    coordNode(nblock,nnode,ndim),                             &
    direct(nblock,3,3), T(nblock,3,3), props(nprops),         &
    stateOld(nblock,nstatev), charLength(nblock,ncomp),       &
    field(nblock, nfieldv)

  Character(len=80) :: cmname

  Double Precision :: center(ndim)

  ! Parameters
  Double Precision, parameter :: zero=0.d0, one=1.d0, two=2.d0

  Master: Do k = 1, nblock

    ! Center is the geometric center of the nodes
    center = zero
    Do n = 1, nnode
      center(:) = center(:) + coordNode(k, n, :)
    End Do
    center(:) = center(:) / nnode

    ! The characteristic element lengths are calculated as the averaged length from the geometric center to each node,
    ! projected onto the material directions (direct)
    charLength(k, :) = zero
    Nodes: Do n = 1, nnode
      Dimensions: Do d = 1, ndim
        charLength(k, d) = charLength(k, d) + ABS(DOT_PRODUCT(coordNode(k, n, :) - center(:), direct(k, :, d)))
      End Do Dimensions
    End Do Nodes
    charLength(k, :) = charLength(k, :) * two / nnode

    ! TODO: incorporate the data in T(nblock,3,3) for use in shell and membrane elements. This version may only be
    ! working for shells where the material direction is set to theta = zero degrees.

    ! Material point coordinates in the material frame (for random initial fiber misalignments)
    If (ncomp == 6) Then
      charLength(k, 4) = DOT_PRODUCT(coordMp(k,:), direct(k,:,1))
      charLength(k, 5) = DOT_PRODUCT(coordMp(k,:), direct(k,:,2))
      charLength(k, 6) = DOT_PRODUCT(coordMp(k,:), direct(k,:,3))
    End If

  End Do Master

  Return
End Subroutine vucharlength
