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

  Integer :: fiberEdge, matrixEdge, thickEdge
  Double Precision :: center(ndim), edges(ndim, ndim), edges_n(ndim, ndim), edges_m(ndim, ndim), fiber_band(ndim), matrix_band(ndim)
  Double Precision :: thick_test, fiber_test, da_max_matrix, da_avg_matrix

  ! Parameters
  Double Precision, parameter :: zero=0.d0, Pi=ACOS(-1.d0), two=2.d0

  ! Evaluate vucharlength() only when element length state variables are undefined.
  runOnce: If ( stateOld(1, 6) == zero ) Then

    Master: Do k = 1, nblock

      ! Material point coordinates in the material frame (for random initial fiber misalignments)
      If (ncomp == 6) Then
        charLength(k, 4) = DOT_PRODUCT(coordMp(k,:), direct(k,:,1))
        charLength(k, 5) = DOT_PRODUCT(coordMp(k,:), direct(k,:,2))
        charLength(k, 6) = DOT_PRODUCT(coordMp(k,:), direct(k,:,3))
      End If

      elementShape: If (nnode == 6 .OR. nnode == 3) Then
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

      Else
        ! Determine the average length of the parallel element edges
        ! The below logic is for either 4-node or 8-node elements, and assumes that coordNode ordering is per element definition.
        ! TODO: Add capability for other element shapes
        edges(1, :) = coordNode(k, 2, :) - coordNode(k, 1, :) + coordNode(k, 3, :) - coordNode(k, 4, :)
        edges(2, :) = coordNode(k, 3, :) - coordNode(k, 2, :) + coordNode(k, 4, :) - coordNode(k, 1, :)
        If (ndim == 3) Then
          edges(1, :) = edges(1, :) - coordNode(k, 5, :) + coordNode(k, 6, :) + coordNode(k, 7, :) - coordNode(k, 8, :)
          edges(2, :) = edges(2, :) - coordNode(k, 5, :) - coordNode(k, 6, :) + coordNode(k, 7, :) + coordNode(k, 8, :)
          edges(3, :) = coordNode(k, 5, :) - coordNode(k, 1, :) + coordNode(k, 6, :) - coordNode(k, 2, :) + &
                        coordNode(k, 7, :) - coordNode(k, 3, :) + coordNode(k, 8, :) - coordNode(k, 4, :)
        End If
        edges = 2 * edges / nnode

        ! edges_n is a matrix of normalized edges
        ! edges_m is a matrix of edge components in the material coordinate system, all positive values
        Do n = 1, ndim
          edges_n(n, :) = edges(n, :) / Length(edges(n, :))
          Do m = 1, ndim
            edges_m(n, m) = ABS(DOT_PRODUCT(edges(n, :), direct(k, :, m)))
          End Do
        End Do

        ! Determine which sets of element edges the thickness, fiber, and matrix material directions are most closely aligned
        !  Determine the thickness-aligned edge
        thickEdge = 3  ! Initial guess that the thickness material direction is most closely aligned with the element 3 direction
        If (ndim == 3) Then
          thick_test = MAX(edges_m(1,3)/Length(edges(1, :)), edges_m(2,3)/Length(edges(2, :)), edges_m(3,3)/Length(edges(3, :)))
          If (edges_m(1, 3)/Length(edges(1, :)) == thick_test) thickEdge = 1
          If (edges_m(2, 3)/Length(edges(2, :)) == thick_test) thickEdge = 2
        End If
        !  Determine the fiber-aligned edge
        fiberEdge = MOD(thickEdge, 3) + 1  ! Initial guess for the fiber-aligned edge
        fiber_test = MAX(edges_m(1,1)/Length(edges(1, :)), edges_m(2,1)/Length(edges(2, :)), edges_m(3,1)/Length(edges(3, :)))
        If (edges_m(MOD(thickEdge + 1, 3) + 1, 1)/Length(edges(MOD(thickEdge + 1, 3) + 1, :)) == fiber_test) Then
          fiberEdge = MOD(thickEdge + 1, 3) + 1
        End If
        !  Determine the matrix-aligned edge
        matrixEdge = 6 - fiberEdge - thickEdge

        ! FIBER DIRECTION CHARACTERISTIC ELEMENT LENGTH
        If (edges_m(fiberEdge, 1) < edges_m(matrixEdge, 1)) Then
          longEdge = matrixEdge
          shortEdge = fiberEdge
        Else
          longEdge = fiberEdge
          shortEdge = matrixEdge
        End If
        ! da_max_fiber is the maximum distance fiber damage may propagate while passing through the element along the matrix direction
        da_max_fiber = ACOS(DOT_PRODUCT(edges_n(fiberEdge, :), -edges_n(matrixEdge, :)))
        da_max_fiber = SIN(da_max_fiber) * Length(edges(shortEdge, :)) / SIN(ACOS(edges_m(longEdge, 2) / Length(edges(longEdge, :))))
        ! da_avg_fiber is the average distance fiber damage will propagate when passing through the element along the matrix direction
        da_avg_fiber = da_max_fiber * MAX(edges_m(fiberEdge, 1), edges_m(matrixEdge, 1)) / (edges_m(fiberEdge, 1) + edges_m(matrixEdge, 1))

        ! Define the characteristic element length for the material fiber direction
        fiber_band = edges(fiberEdge, :) - DOT_PRODUCT(edges(fiberEdge, :), edges_n(matrixEdge, :)) * edges_n(matrixEdge, :)
        charLength(k, 1) = Length(fiber_band)**2 / ABS(DOT_PRODUCT(fiber_band, direct(k, :, 1)))
        charLength(k, 1) = charLength(k, 1) * edges_m(matrixEdge, 2) / da_avg_fiber

        ! MATRIX DIRECTION CHARACTERISTIC ELEMENT LENGTH
        If (edges_m(fiberEdge, 2) < edges_m(matrixEdge, 2)) Then
          longEdge = matrixEdge
          shortEdge = fiberEdge
        Else
          longEdge = fiberEdge
          shortEdge = matrixEdge
        End If
        ! da_max_matrix is the maximum distance a matrix crack may grow while passing through the element along the fiber direction
        da_max_matrix = ACOS(DOT_PRODUCT(edges_n(fiberEdge, :), edges_n(matrixEdge, :)))
        da_max_matrix = SIN(da_max_matrix) * Length(edges(shortEdge, :)) / SIN(ACOS(edges_m(longEdge, 1) / Length(edges(longEdge, :))))
        ! da_avg_matrix is the average distance a matrix crack will grow when passing through the element along the fiber direction
        da_avg_matrix = da_max_matrix * MAX(edges_m(fiberEdge, 2), edges_m(matrixEdge, 2)) / (edges_m(fiberEdge, 2) + edges_m(matrixEdge, 2))

        ! Define the characteristic element length for the material matrix direction
        matrix_band = edges(matrixEdge, :) - DOT_PRODUCT(edges(matrixEdge, :), edges_n(fiberEdge, :)) * edges_n(fiberEdge, :)
        charLength(k, 2) = Length(matrix_band)**2 / ABS(DOT_PRODUCT(matrix_band, direct(k, :, 2)))
        charLength(k, 2) = charLength(k, 2) * edges_m(fiberEdge, 1) / da_avg_matrix

        ! THICKNESS DIRECTION CHARACTERISTIC ELEMENT LENGTH
        If (ndim == 3) charLength(k, 3) = edges_m(thickEdge, 3)

      End If elementShape

    End Do Master

  End If runOnce

  Return
End Subroutine vucharlength
