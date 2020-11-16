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

  Integer :: fiberEdge, matrixEdge, thickEdge, matrixEdge_test, edge_pair_1, edge_pair_2
  Double Precision :: center(ndim), edges(ndim, ndim), tri_edges(3, ndim), wedge_edges(4, ndim)
  Double Precision :: fiber_test, matrix_test, thick_test
  Double Precision :: fiber_test_max, matrix_test_min, thick_test_max
  Double Precision :: da_max_matrix, da_avg_matrix, da_max_fiber, da_avg_fiber
  Double Precision :: psi, gamma, Lc_meshlines_fiber, Lc_meshlines_matrix

  Dimension INTV(1), REALV(1)    ! For Abaqus warning messages
  Character(len=8) CHARV(1)      ! For Abaqus warning messages

  ! Parameters
  Double Precision, parameter :: zero=0.d0, Pi=ACOS(-1.d0), two=2.d0, three=3.d0

  ! Evaluate vucharlength() only when element length state variables are undefined.
  runOnce: If (stateOld(1,6) == zero) Then

    ! Check for the appropriate number of components in *Characteristic Length, definition=USER, components=3
    If ((ncomp /= 2) .AND. (ncomp /= 3) .AND. (ncomp /= 6)) Then
      INTV(1) = ncomp
      Call XPLB_ABQERR(-3,"Invalid number of components in *Characteristic Length. Expecting 2, 3, or 6. Found %I",INTV,REALV,CHARV)
    End If

    Master: Do k = 1,nblock

      ! Material point coordinates in the material frame (for random initial fiber misalignments)
      If (ncomp == 6) Then
        charLength(k,4) = DOT_PRODUCT(coordMp(k,:), direct(k,:,1))
        charLength(k,5) = DOT_PRODUCT(coordMp(k,:), direct(k,:,2))
        charLength(k,6) = DOT_PRODUCT(coordMp(k,:), direct(k,:,3))
      End If

      ! The following assumes that nodes in coordNode() are numbered according to the Abqus documentation.

      elementShape: If (nnode == 3 .OR. nnode == 6) Then  ! For 3- and 6-node elements,

        ! Define edge vectors for each edge of the triangular face
        ! For two-dimensional triangular elements:
        tri_edges(1,:) = coordNode(k,2,:) - coordNode(k,1,:)
        tri_edges(2,:) = coordNode(k,3,:) - coordNode(k,2,:)
        tri_edges(3,:) = coordNode(k,1,:) - coordNode(k,3,:)
        ! For three-dimensional wedge elements:
        WedgeEdges: If (ndim == 3) Then
          wedge_edges(1,:) = (tri_edges(1,:) + coordNode(k,5,:) - coordNode(k,4,:)) / two
          wedge_edges(2,:) = (tri_edges(2,:) + coordNode(k,6,:) - coordNode(k,5,:)) / two
          wedge_edges(3,:) = (tri_edges(3,:) + coordNode(k,4,:) - coordNode(k,6,:)) / two
          wedge_edges(4,:) = (coordNode(k,4,:) - coordNode(k,1,:) + coordNode(k,5,:) - coordNode(k,2,:) + coordNode(k,6,:) - coordNode(k,3,:)) / three

        ! Check to ensure wedge element thickness is aligned with the material direction thickness
          thick_test_max = zero
          Do i=1,4
            thick_test = ABS(DOT_PRODUCT(wedge_edges(i,:), direct(k,:,3)))
            If (thick_test > thick_test_max) Then
              thick_test_max = thick_test
              thickEdge = i
            End If
          End Do
          If (thickEdge /= 4) Then
            Call XPLB_ABQERR(-3,"Wedge elements must have element thickness aligned with material thickness direction.",INTV,REALV,CHARV)
          Else
            tri_edges(1,:) = wedge_edges(1,:)
            tri_edges(2,:) = wedge_edges(2,:)
            tri_edges(3,:) = wedge_edges(3,:)
          End If
        End If WedgeEdges

        ! Determine the fiber-aligned edge vector
        ! The edge vector closest to parallel with the fiber material direction is selected.
        fiber_test_max = zero
        Do i=1,3
          fiber_test = ABS(DOT_PRODUCT(tri_edges(i,:)/Length(tri_edges(i,:)), direct(k,:,1)))  ! a value of 1.0 indicates perpendicular
          If (fiber_test > fiber_test_max) Then
            fiber_test_max = fiber_test
            fiberEdge = i
          End If
        End Do
        edges(1,:) = tri_edges(fiberEdge,:)

        ! Determine the matrix-aligned edge vector
        ! The edge vector most perpendicular to the fiber edge vector is selected.
        matrix_test_min = HUGE(zero)
        Do i=0,1
          matrixEdge_test = MOD(fiberEdge+i,3)+1
          matrix_test = ABS(DOT_PRODUCT(edges(1,:)/Length(edges(1,:)), -tri_edges(matrixEdge_test,:)/Length(tri_edges(matrixEdge_test,:))))
          If (matrix_test < matrix_test_min) Then
            matrix_test_min = matrix_test
            matrixEdge = matrixEdge_test
          End If
        End Do
        edges(2,:) = -tri_edges(matrixEdge,:)

        ! Flip the direction of the edge vectors if the fiber edge is opposite the fiber material direction
        If (DOT_PRODUCT(tri_edges(fiberEdge,:), direct(k,:,1)) < zero) Then
          edges(1,:) = -edges(1,:)
          edges(2,:) = -edges(2,:)
        End If

        If (ndim == 3) edges(3,:) = wedge_edges(thickEdge,:)
        
        ! Reset element edge vector integer variables to corresponding indices in edges()
        fiberEdge = 1
        matrixEdge = 2
        thickEdge = 3

      Else If (nnode == 4 .OR. nnode == 8) Then elementShape  ! For 4- and 8-node elements,
      
        ! Define edge vectors for the oppositely oriented pairs of element edge vectors
        edges(1,:) = coordNode(k,2,:) - coordNode(k,1,:) + coordNode(k,3,:) - coordNode(k,4,:)
        edges(2,:) = coordNode(k,3,:) - coordNode(k,2,:) + coordNode(k,4,:) - coordNode(k,1,:)
        If (ndim == 3) Then
          edges(1,:) = edges(1,:) - coordNode(k,5,:) + coordNode(k,6,:) + coordNode(k,7,:) - coordNode(k,8,:)
          edges(2,:) = edges(2,:) - coordNode(k,5,:) - coordNode(k,6,:) + coordNode(k,7,:) + coordNode(k,8,:)
          edges(3,:) = coordNode(k,5,:) - coordNode(k,1,:) + coordNode(k,6,:) - coordNode(k,2,:) + &
                       coordNode(k,7,:) - coordNode(k,3,:) + coordNode(k,8,:) - coordNode(k,4,:)
        End If
        edges = 2 * edges / nnode

        ! Determine the thickness-aligned edge vector
        ThicknessEdgeHex: If (ndim == 3) Then
          thick_test_max = zero
          Do i=1,3
            thick_test = ABS(DOT_PRODUCT(edges(i,:), direct(k,:,3)))
            If (thick_test > thick_test_max) Then
              thick_test_max = thick_test
              thickEdge = i
            End If
          End Do
        Else ThicknessEdgeHex
          thickEdge = 3
        End If ThicknessEdgeHex
        ! Determine the fiber-aligned and matrix-aligned edge vectors
        edge_pair_1 = MOD(thickEdge,3) + 1
        edge_pair_2 = MOD(edge_pair_1,3) + 1
        FiberEdgeTest: If (ABS(DOT_PRODUCT(edges(edge_pair_1,:)/Length(edges(edge_pair_1,:)), direct(k,:,1))) >= &
            ABS(DOT_PRODUCT(edges(edge_pair_2,:)/Length(edges(edge_pair_2,:)), direct(k,:,1)))) Then
          fiberEdge = edge_pair_1
          matrixEdge = edge_pair_2
        Else FiberEdgeTest
          fiberEdge = edge_pair_2
          matrixEdge = edge_pair_1
        End If FiberEdgeTest

      Else elementShape  ! Element with unsupported number of nodes

        INTV(1) = nnode
        Call XPLB_ABQERR(-3,"Unsupported element shape for VUCHARLENGTH. Expecting 3, 4, 6, or 8 nodes. Found %I",INTV,REALV,CHARV)

      End If elementShape


      ! Define attributes of the mesh related to its misalignment (psi) and skew (gamma)
      psi = ACOS(DOT_PRODUCT(edges(fiberEdge,:) / Length(edges(fiberEdge,:)), direct(k, :, 1)))
      gamma = ACOS(DOT_PRODUCT(edges(fiberEdge,:) / Length(edges(fiberEdge,:)), edges(matrixEdge,:) / Length(edges(matrixEdge,:))))

      ! Fiber-direction characteristic element length
      ! Lc_meshlines_fiber is the characteristic element length for fiber damage growing along the mesh lines
      Lc_meshlines_fiber = Length(edges(fiberEdge,:)) * SIN(gamma) / SIN(psi + gamma)

      If (ABS(DOT_PRODUCT(edges(fiberEdge,:), direct(k,:,1))) > ABS(DOT_PRODUCT(edges(matrixEdge,:), direct(k,:,1)))) Then
        da_max_fiber = Length(edges(matrixEdge,:)) * SIN(gamma) / COS(psi)
        da_avg_fiber = da_max_fiber * ABS(DOT_PRODUCT(edges(fiberEdge,:), direct(k,:,1)))
      Else
        da_max_fiber = Length(edges(fiberEdge,:)) * SIN(gamma) / -COS(psi + gamma)
        da_avg_fiber = da_max_fiber * ABS(DOT_PRODUCT(edges(matrixEdge,:), direct(k,:,1)))
      End If
      da_avg_fiber = da_avg_fiber / (ABS(DOT_PRODUCT(edges(fiberEdge,:), direct(k,:,1))) + &
                        ABS(DOT_PRODUCT(edges(matrixEdge,:), direct(k,:,1))))

      ! Lc_ML_fiber is the characteristic element length for fiber damage growing along the mesh lines
      charLength(k,1) = Lc_meshlines_fiber * ABS(DOT_PRODUCT(edges(matrixEdge,:), direct(k,:,2))) / da_avg_fiber

      ! Matrix-direction characteristic element length
      ! Lc_meshlines_matrix is the characteristic element length for matrix damage growing along the mesh lines
      Lc_meshlines_matrix = Length(edges(matrixEdge,:)) * SIN(gamma) / COS(psi)  ! Equation 16

      ! Equations 17 and 18
      If (ABS(DOT_PRODUCT(edges(fiberEdge,:), direct(k,:,2))) > ABS(DOT_PRODUCT(edges(matrixEdge,:), direct(k,:,2)))) Then
        da_max_matrix = Length(edges(matrixEdge,:)) * SIN(gamma) / SIN(psi)
        da_avg_matrix = da_max_matrix * ABS(DOT_PRODUCT(edges(fiberEdge,:), direct(k,:,2)))
      Else
        da_max_matrix = Length(edges(fiberEdge,:)) * SIN(gamma) / SIN(psi + gamma)
        da_avg_matrix = da_max_matrix * ABS(DOT_PRODUCT(edges(matrixEdge,:), direct(k,:,2)))
      End If
      da_avg_matrix = da_avg_matrix / (ABS(DOT_PRODUCT(edges(fiberEdge,:), direct(k,:,2))) + &
                        ABS(DOT_PRODUCT(edges(matrixEdge,:), direct(k,:,2))))

      charLength(k,2) = Lc_meshlines_matrix * ABS(DOT_PRODUCT(edges(fiberEdge,:), direct(k,:,1))) / da_avg_matrix  ! Equation 20

      ! Thickness-direction characteristic element length
      If (ndim == 3) charLength(k,3) = ABS(DOT_PRODUCT(edges(thickEdge,:), direct(k,:,3)))
      
    End Do Master

  End If runOnce

  Return
End Subroutine vucharlength
