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

!-------------------------------------------------------------------------------
SUBROUTINE WhitneyAVSolver( Model,Solver,dt,Transient )
!-------------------------------------------------------------------------------
  USE DefUtils

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

  LOGICAL :: HasStabC, HasStab, HasStabM, Found, FoundMagnetization, &
      SteadyGauge

  REAL(KIND=dp), ALLOCATABLE :: LOAD(:,:), Acoef(:), &
      STIFF(:,:), FORCE(:), Tcoef(:,:,:), MASS(:,:)
  !$OMP THREADPRIVATE(LOAD, Acoef, Tcoef, STIFF, FORCE, MASS)
  SAVE LOAD, Acoef, Tcoef, STIFF, FORCE, MASS

  REAL(KIND=dp) :: gauge_penalize_c, gauge_penalize_m
  !$OMP THREADPRIVATE(AllocationsDone)
  LOGICAL :: AllocationsDone=.FALSE.

  INTEGER :: istat, nActive, n, nd, nb, t

  CALL Info('WhitneyAVSolver','-------------------------------------------',Level=8 )
  CALL Info('WhitneyAVSolver','Solving the AV equations with edge elements',Level=5 )

  SolverParams => GetSolverParams()

  !-Keyword kungfu related to weak gauge conditions-------------------------------
  SteadyGauge = GetLogical(GetSolverParams(), 'Use Lagrange Gauge', Found) .and. .not. Transient

  IF (SteadyGauge) THEN
    CALL Info("WhitneyAVSolver", "Utilizing Lagrange multipliers for gauge condition in steady state computation")
    IF(.not. ListCheckPresent( SolverParams, 'Linear System Refactorize') ) THEN
      CALL ListAddLogical( SolverParams, 'Linear System Refactorize', .TRUE. )
    END IF
  END IF

  gauge_penalize_c = GetCReal(GetSolverParams(), 'Lagrange Gauge Penalization coefficient', HasStabC)
  gauge_penalize_m = GetCReal(GetSolverParams(), 'Lagrange Gauge Penalization coefficient mass', HasStabM)
  HasStab = HasStabC .OR. HasStabM

  IF (HasStab .and. SteadyGauge) THEN
    WRITE (Message, *), 'Lagrange Gauge penalization coefficient', gauge_penalize_c
    CALL Info('WhitneyAVSolver', message)
    WRITE (Message, *), 'Lagrange Gauge penalization coefficient mass', gauge_penalize_m
    CALL Info('WhitneyAVSolver', message)
  END IF
  !-------------------------------------------------------------------------------

  Mesh => GetMesh()

  CALL DefaultInitialize()

  !-Allocate storage for local matrices-------------------------------------------

  IF ( .NOT. AllocationsDone ) THEN
    N = Mesh % MaxElementDOFs
    ALLOCATE(FORCE(N), &
        LOAD(6,N),  &
        STIFF(N,N), &
        MASS(N,N), &
       !  Tcoef(3,3,N), &
        Acoef(N), &
        STAT=istat)
    IF ( istat /= 0) THEN
      CALL Fatal('WhitneyAVSolver', 'Memory allocation error.')
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
    FoundMagnetization = .FALSE.
    IF(ASSOCIATED(BodyForce)) THEN
      CALL GetRealVector( BodyForce, Load(1:3,1:n), 'Current Density', Found )
      CALL GetRealVector( BodyForce, Load(4:6,1:n), 'Magnetization', FoundMagnetization )
    END IF

    Material => GetMaterial( Element )
    Acoef = 0.0d0
    IF(ASSOCIATED(Material)) THEN
      IF(.NOT. FoundMagnetization) THEN
        CALL GetRealVector( Material, Load(4:6,1:n), &
            'Magnetization', FoundMagnetization )
      END IF

      Acoef(1:n) = GetReal( Material, 'Reluctivity', Found)
    END IF

    CALL LocalMatrix(MASS, STIFF, FORCE, LOAD, Acoef)

  END DO ELEMENT_LOOP

  CONTAINS

    FUNCTION ReallocVec(A, m, istat) result(reallocated)
        IMPLICIT NONE
        REAL(KIND=dp), ALLOCATABLE, INTENT(INOUT) :: A(:)
        INTEGER, INTENT(IN) :: m
        INTEGER, INTENT(OUT) :: istat
        LOGICAL, INTENT(OUT) :: reallocated

        reallocated=.false.

        IF(NOT(ALLOCATED(A))) THEN
          ALLOCATE(A(m), istat=istat)
          reallocated = .true.
          RETURN
        ELSEIF(size(A,1)<m) THEN
          DEALLOCATE(A)
          ALLOCATE(A(m), istat=istat)
          reallocated = .true.
          RETURN
        ENDIF

    END FUNCTION ReallocVec

  FUNCTION ReallocMat(A, m, n, istat) result(reallocated)
      IMPLICIT NONE
      REAL(KIND=dp), ALLOCATABLE, INTENT(INOUT) :: A(:,:)
      INTEGER, INTENT(IN) :: m, n
      INTEGER, INTENT(OUT) :: istat
      LOGICAL, INTENT(OUT) :: reallocated

      reallocated = .false.

      IF(NOT(ALLOCATED(A))) THEN
        ALLOCATE(A(m), istat=istat)
        reallocated = .true.
        RETURN
      ELSEIF(size(A, 1)<m .or. size(A,2)<n) THEN
        DEALLOCATE(A)
        ALLOCATE(A(m,n), istat=istat)
        reallocated = .true.
        RETURN
      ENDIF

  END FUNCTION ReallocMat

    SUBROUTINE LocalMatrix(STIFF, MASS, FORCE, Acoef)
      IMPLICIT NONE
      REAL(KIND=dp) :: STIFF(:,:), FORCE(:), MASS(:,:)
      REAL(KIND=dp) :: Acoef(:,:)

      REAL(KIND=dp), ALLOCATABLE, SAVE :: Basis(:), dBasisdx(:,:), &
          WBasis(:,:), RotWBasis(:,:)
      TYPE(Nodes_t), SAVE :: Nodes
      TYPE(GaussIntegrationPoints_t) :: IP
      !$OMP THREADPRIVATE(Basis, dBasisdx, WBasis, RotWBasis, Nodes)

      INTEGER :: ndtot
      LOGICAL, PARAMETER :: PIOLA=.TRUE.
      INTEGER, PARAMETER :: EDGEBASISDEGREE=1

      INTEGER :: t

      ndtot = nd+nb

      CALL GetElementNodesVec( Nodes, UElement=Element )

      STIFF = 0.0_dp
      MASS = 0.0_dp
      FORCE = 0.0_dp

      CALL GaussPoints(Element, &
          EdgeBasis=.TRUE., &
          PReferenceElement=PIOLA, &
          EdgeBasisDegree=EDGEBASISDEGREE)

      np = n*Solver % DefDofs(GetElementFamily(Element), Element % BodyID, 1)

      DO t=1,IP % n
        stat = EdgeElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
            IP % W(t), DetF = DetJ, Basis = Basis, EdgeBasis = WBasis, &
            RotBasis = RotWBasis, dBasisdx = dBasisdx, &
            BasisDegree = EDGEBASISDEGREE, ApplyPiolaTransform = PIOLA)
        A = SUM( Basis(1:n) * Acoef(1:n) )
      END DO

    END FUNCTION LocalMatrix

!-------------------------------------------------------------------------------
END SUBROUTINE WhitneyAVSolver
!-------------------------------------------------------------------------------
