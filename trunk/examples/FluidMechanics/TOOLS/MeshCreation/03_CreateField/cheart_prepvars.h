        TYPE ARRAY_MESH
          INTEGER Lx					! External Structure: Number of Tessellation Coefficients
          INTEGER Lt					! External Structure: Number of Tessellation elements
          INTEGER ID					! Internal Structure: Basis number to call letter map
          INTEGER Lb
          DOUBLE PRECISION, pointer :: X(:,:)		! External Structure: Tessellation coefficients
          INTEGER, pointer :: T(:,:)			! External Structure: Tessellation map
          INTEGER, pointer :: NBT(:,:)			! Internal Structure: Neighboring element map
          INTEGER, pointer :: B(:,:)
        END TYPE ARRAY_MESH

        TYPE ARRAY_INT
          DOUBLE PRECISION, pointer :: Y(:,:)		! Internal Structure: Basis @ given volume quadrature points
          DOUBLE PRECISION, pointer :: Y_f(:,:,:)	! Internal Structure: Basis @ given facet quadrature points
          DOUBLE PRECISION, pointer :: dY(:,:,:), TdY(:,:,:), dY_f(:,:,:,:)
        END TYPE ARRAY_INT

        TYPE ARRAY_BASE
          INTEGER n					! External Structure: Function Space dimension on T_e
          INTEGER nl					! Internal Structure: Breakdown for tensor input
          INTEGER DISCONT				! External Structure: Continuity constraint
          INTEGER DM					! External Structure: Field dimension
          TYPE(ARRAY_INT) I				! External Structure: Field integrations
          DOUBLE PRECISION, pointer :: Q(:,:)		! External Structure: Basis coefficient tensor
          DOUBLE PRECISION, pointer :: P(:,:)		! External Structure: Basis power tensor
          DOUBLE PRECISION, pointer :: XI(:,:)		! External Structure: Basis xi-point coordinates
          INTEGER, pointer :: B_ID(:,:)			! External Structure: Basis ordering tensor
          CHARACTER*2 CL				! External Structure: Basis call letter
        END TYPE ARRAY_BASE

        TYPE ARRAY_PROBLEM_BASE
          TYPE(ARRAY_BASE), pointer :: B(:)		! Internal Structure: Basis
          INTEGER n_pts					! External Structure: number of volume quadrature points
          INTEGER n_pts_f				! External Structure: number of facet quadrature points
          INTEGER n_ptsl				! Internal Structure: volume quadrature for tensor input
          INTEGER n_pts_fl				! Internal Structure: facet quadrature for tensor input
          INTEGER HEXA					! External Structure: tensor input parameter
          INTEGER FACES					! External Structure: number of master element facets
          INTEGER FNODES				! External Structure: number of nodes per facet (spatial map)
          INTEGER n_B					! External Structure: number of basis
          INTEGER DM					! External Structure: spatial dimension of master element
          INTEGER TRI_BASIS, QUAD_BASIS, TET_BASIS, HEX_BASIS
          DOUBLE PRECISION VL				! External Structure: volume of master element
          DOUBLE PRECISION VL_f(6)			! External Structure: facet volume of master element
          DOUBLE PRECISION nrm(3,6)			! Internal Structure: facet normals of master element
          DOUBLE PRECISION, pointer :: gpt(:,:)		! External Structure: volume quadrature points
          DOUBLE PRECISION, pointer :: gw(:)		! External Structure: volume quadrature weights
          DOUBLE PRECISION, pointer :: gpt_f(:,:,:)	! External Structure: facet quadrature points
          DOUBLE PRECISION, pointer :: gw_f(:)		! External Structure: facet quadrature weights
        END TYPE ARRAY_PROBLEM_BASE