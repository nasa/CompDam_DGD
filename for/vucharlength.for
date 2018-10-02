Subroutine vucharlength(  &
  ! Read only (unmodifiable) variables:
  nblock, nfieldv, nprops, ncomp, ndim, nnode, nstatev, kSecPt, kLayer, kIntPt, jElType, jElem,  &
  totalTime, stepTime, dt, cmname, coordMp, coordNode, direct, T, props, field, stateOld,        &
  ! Write only (modifiable) variables:
  charLength )

  Use matrixAlgUtil_Mod

  Implicit Double Precision (a-h, o-z)
  Integer, parameter :: j_sys_Dimension = 2, n_vec_Length = 136, maxblk = n_vec_Length

  Dimension jElType(3), jElem(nblock), coordMp(nblock,ndim), coordNode(nblock,nnode,ndim), direct(nblock,3,3),  &
    T(nblock,3,3), props(nprops), stateOld(nblock,nstatev), charLength(nblock,ncomp), field(nblock, nfieldv)

  Character(len=80) :: cmname

  Integer :: fiberDir, matrixDir, thickDir
  Double Precision :: edges(ndim, ndim), edges_n(ndim, ndim), edges_m(ndim, ndim)
  Double Precision :: thick_test, fiber_test, delta_a_max, delta_a_avg

  ! Parameters
  Double Precision, parameter :: zero=0.d0, Pi=ACOS(-1.d0)

  ! Evaluate vucharlength() only when element length state variables are undefined.
  If ( stateOld(1, 6) == zero ) Then

    Master: Do k = 1, nblock

      ! Material point coordinates in the material frame (for random initial fiber misalignments)
      If (ncomp == 6) Then
        charLength(k, 4) = DOT_PRODUCT(coordMp(k,:), direct(k,:,1))
        charLength(k, 5) = DOT_PRODUCT(coordMp(k,:), direct(k,:,2))
        charLength(k, 6) = DOT_PRODUCT(coordMp(k,:), direct(k,:,3))
      End If

      ! Determine the average length of the parallel element edges
      edges(1, :) = coordNode(k, 2, :) - coordNode(k, 1, :) + coordNode(k, 3, :) - coordNode(k, 4, :)
      edges(2, :) = coordNode(k, 3, :) - coordNode(k, 2, :) + coordNode(k, 4, :) - coordNode(k, 1, :)
      If (ndim == 3) Then
        edges(1, :) = edges(1, :) - coordNode(k, 5, :) + coordNode(k, 6, :) + coordNode(k, 7, :) - coordNode(k, 8, :)
        edges(2, :) = edges(2, :) - coordNode(k, 5, :) - coordNode(k, 6, :) + coordNode(k, 7, :) + coordNode(k, 8, :)
        edges(3, :) = coordNode(k, 5, :) - coordNode(k, 1, :) + coordNode(k, 6, :) - coordNode(k, 2, :) + &
                      coordNode(k, 7, :) - coordNode(k, 3, :) + coordNode(k, 8, :) - coordNode(k, 4, :)
      End If
      edges = 2 * edges / nnode

      ! Matrix of normalized edges, edges_n(ndim, ndim)
      Do n = 1, ndim; edges_n(n, :) = edges(n, :) / Length(edges(n, :)); End Do

      ! Matrix of edge components in the material coordinate system
      edges_m(1, 1) = ABS(DOT_PRODUCT(edges(1, :), direct(k, :, 1)))
      edges_m(1, 2) = ABS(DOT_PRODUCT(edges(1, :), direct(k, :, 2)))
      edges_m(2, 1) = ABS(DOT_PRODUCT(edges(2, :), direct(k, :, 1)))
      edges_m(2, 2) = ABS(DOT_PRODUCT(edges(2, :), direct(k, :, 2)))
      If (ndim == 3) Then
        edges_m(3, 1) = ABS(DOT_PRODUCT(edges(3, :), direct(k, :, 1)))
        edges_m(3, 2) = ABS(DOT_PRODUCT(edges(3, :), direct(k, :, 2)))

        edges_m(1, 3) = ABS(DOT_PRODUCT(edges(1, :), direct(k, :, 3)))
        edges_m(2, 3) = ABS(DOT_PRODUCT(edges(2, :), direct(k, :, 3)))
        edges_m(3, 3) = ABS(DOT_PRODUCT(edges(3, :), direct(k, :, 3)))
      End If

      ! Determine which sets of element edges the thickness, fiber, and matrix material directions are most closely aligned
      !  Evaluate the thickness direction
      thickDir = 3  ! Initial guess that the thickness material direction is most closely aligned with the global 3 direction
      If (ndim == 3) Then
        thick_test = MAX(edges_m(1,3)/Length(edges(1, :)), edges_m(2,3)/Length(edges(2, :)), edges_m(3,3)/Length(edges(3, :)))
        If (edges_m(1, 3)/Length(edges(1, :)) == thick_test) thickDir = 1
        If (edges_m(2, 3)/Length(edges(2, :)) == thick_test) thickDir = 2
      End If
      !  Evaluate the fiber direction
      fiberDir = MOD(thickDir, 3) + 1  ! Initial guess for the fiber material direction
      fiber_test = MAX(edges_m(1,1)/Length(edges(1, :)), edges_m(2,1)/Length(edges(2, :)), edges_m(3,1)/Length(edges(3, :)))
      If (edges_m(MOD(thickDir + 1, 3) + 1, 1)/Length(edges(MOD(thickDir + 1, 3) + 1, :)) == fiber_test) Then
        fiberDir = MOD(thickDir + 1, 3) + 1
      End If
      !  Evaluate the matrix direction
      matrixDir = 6 - fiberDir - thickDir  ! The matrix material direction is the only unspoken for set of edges

      ! delta_a_max is the maximum distance a crack may grow while passing through the element along the fiber direction
      delta_a_max = Length(edges(fiberDir, :)) * SIN(ACOS(DOT_PRODUCT(edges_n(fiberDir, :), edges_n(matrixDir, :)))) / &
                     SIN(ACOS(edges_m(matrixDir, 1) / Length(edges(matrixDir, :))))
      ! delta_a_avg is the average distance a crack will grow when passing through the element along the fiber direction when it's
      !  position along the height (i.e., matrix direction) of the element is unknown
      delta_a_avg = delta_a_max * edges_m(matrixDir, 2) / (edges_m(fiberDir, 2) + edges_m(matrixDir, 2))

      ! Define the characteristic element lengths in the fiber (1), matrix (2), and thickness (3) material directions
      charLength(k, 1) = edges_m(fiberDir, 1)
      charLength(k, 2) = Length(edges(matrixDir, :)) * SIN(ACOS(DOT_PRODUCT(-edges_n(fiberDir, :), edges_n(matrixDir, :)))) / &
                          SIN(ACOS(edges_m(fiberDir, 2) / Length(edges(fiberDir, :))))
      charLength(k, 2) = charLength(k, 2) * edges_m(fiberDir, 1) / delta_a_avg
      If (ndim == 3) charLength(k, 3) = edges_m(thickDir, 3)

    End Do Master

  End If

  Return
End Subroutine vucharlength
