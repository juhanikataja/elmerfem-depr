!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! *
! *  This program is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU General Public License
! *  as published by the Free Software Foundation; either version 2
! *  of the License, or (at your option) any later version.
! *
! *  This program is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! *  GNU General Public License for more details.
! *
! *  You should have received a copy of the GNU General Public License
! *  along with this program (in file fem/GPL-2); if not, write to the
! *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
! *  Boston, MA 02110-1301, USA.
! *
! *****************************************************************************/
!
!/******************************************************************************
! *
! *  Authors: Juhani Kataja
! *  Email:   juhani.kataja@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland
! *
! *  Original Date: 08 Jun 1997
! *
! *****************************************************************************/



MODULE AVSolverUtils
  USE DefUtils

  LOGICAL, PARAMETER :: DPIOLA=.TRUE.
  INTEGER, PARAMETER :: DEDGEBASISDEGREE=1

CONTAINS
  !-------------------------------------------------------------------------------
  FUNCTION ReallocVec(A, m, istat) result(reallocated)
    USE DefUtils
    IMPLICIT NONE
    REAL(KIND=dp), ALLOCATABLE, INTENT(INOUT) :: A(:)
    INTEGER, INTENT(IN) :: m
    INTEGER, INTENT(OUT) :: istat
    LOGICAL :: reallocated

    reallocated=.false.
    istat = 0

    IF(.NOT. ALLOCATED(A)) THEN
      ALLOCATE(A(m), stat=istat)
      reallocated = .true.
      RETURN
    ELSEIF(size(A,1)<m) THEN
      DEALLOCATE(A)
      ALLOCATE(A(m), stat=istat)
      reallocated = .true.
      RETURN
    ENDIF

  END FUNCTION ReallocVec

  !-------------------------------------------------------------------------------
  FUNCTION ReallocMat(A, m, n, istat) result(reallocated)
    USE DefUtils
    IMPLICIT NONE
    REAL(KIND=dp), ALLOCATABLE, INTENT(INOUT) :: A(:,:)
    INTEGER, INTENT(IN) :: m, n
    INTEGER, INTENT(OUT) :: istat
    LOGICAL :: reallocated

    reallocated = .false.

    IF(.NOT. ALLOCATED(A)) THEN
      ALLOCATE(A(m,n), stat=istat)
      reallocated = .true.
      RETURN
    ELSEIF(size(A, 1)<m .or. size(A,2)<n) THEN
      DEALLOCATE(A)
      ALLOCATE(A(m,n), stat=istat)
      reallocated = .true.
      RETURN
    ENDIF

  END FUNCTION ReallocMat

  SUBROUTINE LocalMatrix(STIFF, MASS, FORCE, LOAD, Acoef, Element, Solver, n_nodes, n_dofs, n_bubbles)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(INOUT) :: STIFF(:,:), FORCE(:), MASS(:,:), LOAD(:,:), &
         Acoef(:)
    TYPE(Solver_t) :: Solver
    INTEGER, INTENT(IN) :: n_nodes, n_dofs, n_bubbles

    REAL(KIND=dp), ALLOCATABLE, SAVE :: Basis(:), dBasisdx(:,:), &
         WBasis(:,:), RotWBasis(:,:)

    REAL(KIND=dp) :: L(3), DetJ, M(3), nu

    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Element_t), POINTER :: Element

    INTEGER :: ndtot
    LOGICAL :: reallocated, stat

    INTEGER :: t, i, j, np, p, q, istat, n

    n = n_nodes
    ndtot = n_dofs + n_bubbles

    CALL GetElementNodesVec( Nodes, UElement=Element )

    reallocated = ReallocVec(Basis, n, istat)
    IF(istat /= 0) CALL Fatal('AVSolver::LocalMatrix', &
         'Memory allocation error in basis storage reservation')
    reallocated = ReallocMat(dbasisdx, n, 3, istat)
    IF(istat /= 0) CALL Fatal('AVSolver::LocalMatrix', &
         'Memory allocation error in basis storage reservation')
    reallocated = ReallocMat(WBasis, ndtot, 3, istat)
    IF(istat /= 0) CALL Fatal('AVSolver::LocalMatrix', &
         'Memory allocation error in basis storage reservation')
    reallocated = ReallocMat(RotWBasis, ndtot, 3, istat)
    IF(istat /= 0) CALL Fatal('AVSolver::LocalMatrix', &
         'Memory allocation error in basis storage reservation')

    STIFF = 0.0_dp
    MASS = 0.0_dp
    FORCE = 0.0_dp

    IP = GaussPoints(Element, &
         EdgeBasis=.TRUE., &
         PReferenceElement=DPIOLA, &
         EdgeBasisDegree=DEDGEBASISDEGREE)

    np = n * Solver % Def_Dofs(GetElementFamily(Element), Element % BodyID, 1)

    GAUSS_LOOP: DO t=1,IP % n
      stat = EdgeElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
           IP % W(t), DetF = DetJ, Basis = Basis, EdgeBasis = WBasis, &
           RotBasis = RotWBasis, dBasisdx = dBasisdx, &
           BasisDegree = 1, ApplyPiolaTransform = .TRUE.)

      ! Material term
      nu = SUM( Basis(1:n) * Acoef(1:n) )

      ! Load term
      L = MATMUL(LOAD(1:3,1:n), Basis(1:n))
      M = MATMUL(LOAD(4:6,1:n), Basis(1:n))

      DO i = 1,ndtot-np
        p = i+np
        FORCE(p) = FORCE(p) + (SUM(L*WBasis(i,:)) + &
             SUM(M*RotWBasis(i,:)))*detJ*IP%s(t) 
        DO j = 1,ndtot-np
          q = j+np
          STIFF(p,q) = STIFF(p,q) + nu* &
               SUM(RotWBasis(i,:)*RotWBasis(j,:))*&
               detJ*IP%s(t)
        END DO
      END DO

    END DO GAUSS_LOOP

    ! DEBUG
#ifdef DEBUG
    print *, 'np = ', np
    DO i = 1, ndtot-np
      write (*,'(A, I1, A, *(E12.3))'),  'stiff(',i,') = ', stiff(i, 1:ndtot)
    END DO
    write (*,'(A, *(E12.3))'),  'force = ', force(1:ndtot)
#endif

  END SUBROUTINE LocalMatrix

  SUBROUTINE LocalMatrixThr(STIFF, MASS, FORCE, LOAD, Acoef, Element, Solver, n_nodes, n_dofs, n_bubbles)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT(INOUT) :: STIFF(:,:), FORCE(:), MASS(:,:), LOAD(:,:), &
         Acoef(:)
    TYPE(Solver_t) :: Solver
    INTEGER, INTENT(IN) :: n_nodes, n_dofs, n_bubbles

    REAL(KIND=dp), ALLOCATABLE, SAVE :: Basis(:), dBasisdx(:,:), &
         WBasis(:,:), RotWBasis(:,:)

    REAL(KIND=dp) :: L(3), DetJ, M(3), nu

    TYPE(Nodes_t), SAVE :: Nodes
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Element_t), POINTER :: Element

    !$OMP THREADPRIVATE(Nodes, Basis, dBasisdx, WBasis, RotWBasis)
    
    INTEGER :: ndtot
    LOGICAL :: reallocated, stat

    INTEGER :: t, i, j, np, p, q, istat, n

    n = n_nodes
    ndtot = n_dofs + n_bubbles

    CALL GetElementNodesVec( Nodes, UElement=Element )

    reallocated = ReallocVec(Basis, n, istat)
    IF(istat /= 0) CALL Fatal('AVSolver::LocalMatrix', &
         'Memory allocation error in basis storage reservation')
    reallocated = ReallocMat(dbasisdx, n, 3, istat)
    IF(istat /= 0) CALL Fatal('AVSolver::LocalMatrix', &
         'Memory allocation error in basis storage reservation')
    reallocated = ReallocMat(WBasis, ndtot, 3, istat)
    IF(istat /= 0) CALL Fatal('AVSolver::LocalMatrix', &
         'Memory allocation error in basis storage reservation')
    reallocated = ReallocMat(RotWBasis, ndtot, 3, istat)
    IF(istat /= 0) CALL Fatal('AVSolver::LocalMatrix', &
         'Memory allocation error in basis storage reservation')

    STIFF = 0.0_dp
    MASS = 0.0_dp
    FORCE = 0.0_dp

    IP = GaussPoints(Element, &
         EdgeBasis=.TRUE., &
         PReferenceElement=DPIOLA, &
         EdgeBasisDegree=DEDGEBASISDEGREE)

    np = n * Solver % Def_Dofs(GetElementFamily(Element), Element % BodyID, 1)

    GAUSS_LOOP: DO t=1,IP % n
      stat = EdgeElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
           IP % W(t), DetF = DetJ, Basis = Basis, EdgeBasis = WBasis, &
           RotBasis = RotWBasis, dBasisdx = dBasisdx, &
           BasisDegree = 1, ApplyPiolaTransform = .TRUE.)

      ! Material term
      nu = SUM( Basis(1:n) * Acoef(1:n) )

      ! Load term
      L = MATMUL(LOAD(1:3,1:n), Basis(1:n))
      M = MATMUL(LOAD(4:6,1:n), Basis(1:n))

      DO i = 1,ndtot-np
        p = i+np
        FORCE(p) = FORCE(p) + (SUM(L*WBasis(i,:)) + &
             SUM(M*RotWBasis(i,:)))*detJ*IP%s(t) 
        DO j = 1,ndtot-np
          q = j+np
          STIFF(p,q) = STIFF(p,q) + nu* &
               SUM(RotWBasis(i,:)*RotWBasis(j,:))*&
               detJ*IP%s(t)
        END DO
      END DO

    END DO GAUSS_LOOP
    
  END SUBROUTINE LocalMatrixThr

  SUBROUTINE GetReluctivity(Material,Acoef,n)
    !------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(ValueList_t), POINTER :: Material
    INTEGER :: n
    REAL(KIND=dp) :: Acoef(:)
    !------------------------------------------------------------------------------
    LOGICAL :: Found, FirstTime = .TRUE., Warned = .FALSE.
    REAL(KIND=dp) :: Avacuum

    SAVE Avacuum 

    IF ( FirstTime ) THEN
      Avacuum = GetConstReal( CurrentModel % Constants, &
           'Permeability of Vacuum', Found )
      IF(.NOT. Found ) Avacuum = PI * 4.0d-7
      FirstTime = .FALSE.
    END IF

    Acoef(1:n) = GetReal( Material, 'Relative Permeability', Found )
    IF ( Found ) THEN
      Acoef(1:n) = Avacuum * Acoef(1:n)
    ELSE
      Acoef(1:n) = GetReal( Material, 'Permeability', Found )
    END IF
    IF ( Found ) THEN
      Acoef(1:n) = 1._dp / Acoef(1:n)
    ELSE
      Acoef(1:n) = GetReal( Material, 'Reluctivity', Found )
    END IF
    IF( .NOT. Found .AND. .NOT. Warned .AND. &
         .NOT. ListCheckPresent(Material, 'H-B Curve') ) THEN
      CALL Warn('GetReluctivityR','Give > Relative Permeability < or > Reluctivity <  for material!')
      Warned = .TRUE.
    END IF

    !------------------------------------------------------------------------------
  END SUBROUTINE GetReluctivity

END MODULE AVSolverUtils

SUBROUTINE AVSolver( Model,Solver,dt,Transient )
  USE DefUtils
  USE AVSolverUtils

  IMPLICIT NONE
  !------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model

  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
  !------------------------------------------------------------------------------
  ! Local variables
  !------------------------------------------------------------------------------
  TYPE(Element_t), POINTER :: Element
  TYPE(ValueList_t), POINTER :: SolverParams, BodyForce, Material, BC, BodyParams
  TYPE(Mesh_t), POINTER :: Mesh

  LOGICAL :: Found

  REAL(KIND=dp), ALLOCATABLE :: LOAD(:,:), Acoef(:), &
       STIFF(:,:), FORCE(:), Tcoef(:,:,:), MASS(:,:)
  !$OMP THREADPRIVATE(LOAD, Acoef, Tcoef, STIFF, FORCE, MASS)
  SAVE LOAD, Acoef, Tcoef, STIFF, FORCE, MASS

  ! LOGICAL :: HasStabC, HasStab, HasStabM, SteadyGauge
  ! REAL(KIND=dp) :: gauge_penalize_c, gauge_penalize_m
  LOGICAL :: AllocationsDone=.FALSE.
  !$OMP THREADPRIVATE(AllocationsDone)
  
  INTEGER :: istat, nActive, n, nd, nb, t

  REAL(KIND=dp) :: Norm

  CALL Info('AVSolver','-------------------------------------------',Level=8 )
  CALL Info('AVSolver','Solving the AV equations with edge elements',Level=5 )

  SolverParams => GetSolverParams()

  Mesh => GetMesh()

  CALL ResetTimer('MGDynAssembly')
  CALL DefaultInitialize()

  !-Allocate storage for local matrices-------------------------------------------
  IF ( .NOT. AllocationsDone ) THEN
    N = Mesh % MaxElementDOFs
    ALLOCATE(FORCE(N), &
         LOAD(6,N),  &
         STIFF(N,N), &
         MASS(N,N), &
         Acoef(N), &
         STAT=istat)
    IF ( istat /= 0) THEN
      CALL Fatal('AVSolver', 'Memory allocation error.')
    END IF
    AllocationsDone = .TRUE.
  END IF
  !-------------------------------------------------------------------------------

  nActive = GetNOFActive()
  ELEMENT_LOOP : DO t=1,nActive
    Element => GetActiveElement(t)
    n  = GetElementNOFNodes() ! kulmat
    nd = GetElementNOFDOFs()  ! vapausasteet
    nb = GetElementNOFBDOFs()  ! sisäiset vapausasteet

    BodyForce => GetBodyForce(Element)
    LOAD = 0.0d0
    IF(ASSOCIATED(BodyForce)) THEN
      CALL GetRealVector( BodyForce, Load(1:3,1:n), 'Current Density', Found )
    END IF

    Material => GetMaterial( Element )
    Acoef = 0.0d0
    IF(ASSOCIATED(Material)) THEN
      Acoef(1:n) = GetReal( Material, 'Reluctivity', Found)
      CALL GetRealVector(Material, Load(4:6,1:n), 'Magnetization', Found )
    END IF

    CALL LocalMatrix(STIFF, MASS, FORCE, LOAD, Acoef, Element, Solver, n, nd, nb)

    CALL DefaultUpdateEquations( STIFF, FORCE, UElement=Element)

  END DO ELEMENT_LOOP

  CALL CheckTimer('MGDynAssembly', DELETE=.TRUE.)
  CALL DefaultFinishAssembly()
  CALL DefaultDirichletBCs()

  Norm = DefaultSolve()


  !-------------------------------------------------------------------------------
END SUBROUTINE AVSolver

SUBROUTINE AVSolver_init0(Model, Solver, dt, Transient)
  USE DefUtils

  IMPLICIT NONE

  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient

  LOGICAL :: Found
  TYPE(ValueList_t), POINTER :: SolverParams

  SolverParams => GetSolverParams()

  CALL Warn("AVSolver_init0", "Resetting element string to n:0 e:1 -brick b:3 -quad_face b:2")
  CALL ListAddString( SolverParams, "Element", "n:0 e:1 -brick b:3 -quad_face b:2")
END SUBROUTINE AVSolver_init0

SUBROUTINE AVSolver_threaded_init0(Model, Solver, dt, Transient)
  USE DefUtils

  IMPLICIT NONE

  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient

  LOGICAL :: Found
  TYPE(ValueList_t), POINTER :: SolverParams

  SolverParams => GetSolverParams()

  CALL Warn("AVSolver_threaded_init0", "Resetting element string to n:0 e:1 -brick b:3 -quad_face b:2")
  CALL ListAddString( SolverParams, "Element", "n:0 e:1 -brick b:3 -quad_face b:2")
END SUBROUTINE AVSolver_threaded_init0

SUBROUTINE AVSolver_threaded(Model, Solver, dt, Transient )
  USE DefUtils
  USE AVSolverUtils

  IMPLICIT NONE
  !------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model

  REAL(KIND=dp) :: dt
  LOGICAL :: Transient
  !------------------------------------------------------------------------------
  ! Local variables
  !------------------------------------------------------------------------------
  TYPE(Element_t), POINTER :: Element
  TYPE(ValueList_t), POINTER :: SolverParams, BodyForce, Material, BC, BodyParams
  TYPE(Mesh_t), POINTER :: Mesh

  LOGICAL :: Found

  REAL(KIND=dp), POINTER :: LOAD(:,:), Acoef(:), &
       STIFF(:,:), FORCE(:), Tcoef(:,:,:), MASS(:,:)
  !$OMP THREADPRIVATE(LOAD, Acoef, Tcoef, STIFF, FORCE, MASS)
  SAVE LOAD, Acoef, Tcoef, STIFF, FORCE, MASS

  ! LOGICAL :: HasStabC, HasStab, HasStabM, SteadyGauge
  ! REAL(KIND=dp) :: gauge_penalize_c, gauge_penalize_m
  LOGICAL :: AllocationsDone=.FALSE.
  !$OMP THREADPRIVATE(AllocationsDone)

  INTEGER :: istat, nActive, n, nd, nb, t

  REAL(KIND=dp) :: Norm

  CALL Info('AVSolver','-------------------------------------------',Level=8 )
  CALL Info('AVSolver','Solving the AV equations with edge elements',Level=5 )

  SolverParams => GetSolverParams()

  Mesh => GetMesh()

  CALL ResetTimer('MGDynAssembly')
  CALL DefaultInitialize()

  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP SHARED(Solver, Model, dt, Transient, SolverParams, Mesh, nActive) &
  !$OMP PRIVATE(Element, BodyForce, Material, BC, BodyParams, istat, &
  !$OMP n, nd, nb, t, Found)
  !-Allocate storage for local matrices-------------------------------------------
  IF ( .NOT. AllocationsDone ) THEN
    N = Mesh % MaxElementDOFs
    ALLOCATE(FORCE(N), &
         LOAD(6,N),  &
         STIFF(N,N), &
         MASS(N,N), &
         Acoef(N), &
         STAT=istat)
    IF ( istat /= 0) THEN
      CALL Fatal('AVSolver', 'Memory allocation error.')
    END IF
    AllocationsDone = .TRUE.
  END IF
  !-------------------------------------------------------------------------------

  !$OMP SINGLE
  nActive = GetNOFActive()
  !$OMP END SINGLE

  !$OMP DO SCHEDULE(STATIC)
  ELEMENT_LOOP : DO t=1,nActive
    Element => GetActiveElement(t)
    n  = GetElementNOFNodes() ! kulmat
    nd = GetElementNOFDOFs()  ! vapausasteet
    nb = GetElementNOFBDOFs()  ! sisäiset vapausasteet

    BodyForce => GetBodyForce(Element)
    LOAD = 0.0d0
    IF(ASSOCIATED(BodyForce)) THEN
      CALL GetRealVector( BodyForce, Load(1:3,1:n), 'Current Density', Found )
    END IF

    Material => GetMaterial( Element )
    Acoef = 0.0d0
    IF(ASSOCIATED(Material)) THEN
      Acoef(1:n) = GetReal( Material, 'Reluctivity', Found)
      CALL GetRealVector(Material, Load(4:6,1:n), 'Magnetization', Found )
    END IF

    CALL LocalMatrixThr(STIFF, MASS, FORCE, LOAD, Acoef, Element, Solver, n, nd, nb)

    CALL DefaultUpdateEquations( STIFF, FORCE, UElement=Element)

  END DO ELEMENT_LOOP
  !$OMP END DO
  !$OMP END PARALLEL

  CALL CheckTimer('MGDynAssembly', DELETE=.TRUE.)
  CALL DefaultFinishAssembly()
  CALL DefaultDirichletBCs()

  Norm = DefaultSolve()


  !-------------------------------------------------------------------------------

END SUBROUTINE AVSolver_threaded
