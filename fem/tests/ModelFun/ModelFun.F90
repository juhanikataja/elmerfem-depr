!-----------------------------------------------------------------------------
!> A prototype solver for advection-diffusion-reaction equation,
!> This equation is generic and intended for education purposes
!> but may also serve as a starting point for more complex solvers.
!------------------------------------------------------------------------------
SUBROUTINE AdvDiffSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  TYPE(Element_t),POINTER :: Element
  REAL(KIND=dp) :: Norm
  INTEGER :: n, nb, nd, t, active
  INTEGER :: iter, maxiter, dim
  LOGICAL :: Found
  TYPE(ValueList_t), POINTER :: Params
  INTEGER(KIND=AddrInt) :: Proc
  TYPE(ValueListEntry_t), POINTER :: ptr
  !------------------------------------------------------------------------------
  TYPE(ValueHandle_t) :: Load_h, DiffCoeff_h

  Params => GetSolverParams()
  
  CALL Info('ModelFun','Creating internal function >MySqrt<')
  Proc = GetProcAddr("ModelFun SquareRoot" )
  CALL ListAddConstReal(Params,"MySqrt",1.0_dp,Proc )
  ptr => ListFind(Params,"MySqrt") 
  ptr % Fdim = 1
  
  
  CALL ListInitElementKeyword( Load_h,'Body Force','Field Source')
  CALL ListInitElementKeyword( DiffCoeff_h,'Material','Diffusion Coefficient', EvaluateAtIp=.TRUE.)
  
  CALL DefaultStart()
  
  maxiter = ListGetInteger( GetSolverParams(),&
      'Nonlinear System Max Iterations',Found,minv=1)
  IF(.NOT. Found ) maxiter = 1

  dim = CoordinateSystemDimension()


  
  ! Nonlinear iteration loop:
  !--------------------------
  DO iter=1,maxiter
    
    ! System assembly:
    !----------------
    CALL DefaultInitialize()
    Active = GetNOFActive()
    DO t=1,Active
      Element => GetActiveElement(t)
      n  = GetElementNOFNodes()
      nd = GetElementNOFDOFs()
      nb = GetElementNOFBDOFs()

      CALL LocalMatrix(  Element, n, nd+nb )
    END DO

    CALL DefaultFinishBulkAssembly()
    CALL DefaultFinishBoundaryAssembly()
    CALL DefaultFinishAssembly()
    CALL DefaultDirichletBCs()

    ! And finally, solve:
    !--------------------
    Norm = DefaultSolve()

    IF( Solver % Variable % NonlinConverged > 0 ) EXIT

  END DO

  CALL DefaultFinish()
  
CONTAINS

! Assembly of the matrix entries arising from the bulk elements
!------------------------------------------------------------------------------
  SUBROUTINE LocalMatrix( Element, n, nd )
!------------------------------------------------------------------------------
    INTEGER :: n, nd
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: D,C,R, rho,Velo(3,n),a(3), Weight
    REAL(KIND=dp) :: Basis(nd),dBasisdx(nd,3),DetJ,LoadAtIP
    REAL(KIND=dp) :: MASS(nd,nd), STIFF(nd,nd), FORCE(nd)
    LOGICAL :: Stat,Found
    INTEGER :: i,t,p,q
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(ValueList_t), POINTER :: BodyForce, Material
    TYPE(Nodes_t) :: Nodes
    SAVE Nodes
!------------------------------------------------------------------------------


    CALL GetElementNodes( Nodes )
    MASS  = 0._dp
    STIFF = 0._dp
    FORCE = 0._dp
    a = 0.0_dp
    
    ! Numerical integration:
    !-----------------------
    IP = GaussPoints( Element )
    DO t=1,IP % n
      ! Basis function values & derivatives at the integration point:
      !--------------------------------------------------------------
      stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), &
              IP % W(t), detJ, Basis, dBasisdx )

      ! The source term at the integration point:
      !------------------------------------------
      D = ListGetElementReal( DiffCoeff_h, element=Element, found=Found , gausspoint = t)       
      ! print *, 'D=', D
      LoadAtIP = ListGetElementReal( Load_h, Basis, Element, Found ) 
      
      Weight = IP % s(t) * DetJ

      ! diffusion term (D*grad(u),grad(v)):
      ! -----------------------------------
      STIFF(1:nd,1:nd) = STIFF(1:nd,1:nd) + Weight * &
             D * MATMUL( dBasisdx, TRANSPOSE( dBasisdx ) )
      FORCE(1:nd) = FORCE(1:nd) + Weight * LoadAtIP * Basis(1:nd)
    END DO

    IF(TransientSimulation) CALL Default1stOrderTime(MASS,STIFF,FORCE)
    CALL LCondensate( nd-nb, nb, STIFF, FORCE )
    CALL DefaultUpdateEquations(STIFF,FORCE)
!------------------------------------------------------------------------------
  END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------



! Perform static condensation in case bubble dofs are present
!------------------------------------------------------------------------------
  SUBROUTINE LCondensate( N, Nb, K, F )
!------------------------------------------------------------------------------
    USE LinearAlgebra
    INTEGER :: N, Nb
    REAL(KIND=dp) :: K(:,:),F(:),Kbb(Nb,Nb), &
         Kbl(Nb,N), Klb(N,Nb), Fb(Nb)

    INTEGER :: m, i, j, l, p, Ldofs(N), Bdofs(Nb)

    IF ( Nb <= 0 ) RETURN

    Ldofs = (/ (i, i=1,n) /)
    Bdofs = (/ (i, i=n+1,n+nb) /)

    Kbb = K(Bdofs,Bdofs)
    Kbl = K(Bdofs,Ldofs)
    Klb = K(Ldofs,Bdofs)
    Fb  = F(Bdofs)

    CALL InvertMatrix( Kbb,nb )

    F(1:n) = F(1:n) - MATMUL( Klb, MATMUL( Kbb, Fb  ) )
    K(1:n,1:n) = K(1:n,1:n) - MATMUL( Klb, MATMUL( Kbb, Kbl ) )
!------------------------------------------------------------------------------
  END SUBROUTINE LCondensate
!------------------------------------------------------------------------------
  
!------------------------------------------------------------------------------
END SUBROUTINE AdvDiffSolver
!------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
FUNCTION SquareRoot(Variable, Element, IP, ipind) RESULT ( y )
!-------------------------------------------------------------------------------
  USE DefUtils
  IMPLICIT NONE
  TYPE(Variable_t), INTENT(IN) :: Variable
  TYPE(Element_t), INTENT(IN) :: Element
  TYPE(GaussIntegrationPoints_t), INTENT(IN) :: IP
  integer :: ipind
  REAL(KIND=dp) :: y
!-------------------------------------------------------------------------------
  TYPE(Nodes_t), SAVE :: Nodes
  REAL(KIND=dp), ALLOCATABLE, SAVE :: Basis(:,:), dBasisdx(:,:,:), DetJ(:)
  REAL(KIND=dp), ALLOCATABLE, SAVE :: local_phi(:), phi_cache(:)
  logical :: estat, recache
  integer :: elemind  = -1
  integer :: k, t, ngp, nd
  !$OMP THREADPRIVATE(Nodes, Basis, dBasisdx, DetJ, local_phi, elemind, phi_cache)

  print *, 'at SquareRoot'
  recache = .false.
  if (elemind /= Element % ElementIndex) recache = .true.
  elemind = Element % ElementIndex

  IF (recache) THEN
    call GetElementNodesVec(Nodes, Uelement=Element, &
        USolver = Variable % Solver, UMesh = Variable % Solver % Mesh)

    ngp = IP % n
    nd = GetElementNOFDOFs(Element, USolver=Variable % Solver)

    IF (.not. allocated(phi_cache) .or. &
        size(phi_cache, 1) < ngp) THEN
      ALLOCATE(phi_cache(ngp))
    END IF

    IF (.not. ALLOCATED(Basis) .or. size(basis,1) < nd .or. size(basis,2) < ngp) THEN
      ALLOCATE(basis(ngp, nd), detj(ngp))
    END IF

    IF(.not. ALLOCATED(local_phi) .or. size(local_phi,1) < nd) THEN
      ALLOCATE(local_phi(nd))
    END IF

    CALL GetScalarLocalSolution(local_phi, UVariable=Variable)

    estat = ElementInfoVec( Element, Nodes, ngp, IP % U, IP % V, &
        IP % W, detJ, size(basis,2), Basis )
    phi_cache = 0.0_dp
    
    DO t = 1, nd
      DO k = 1, ip % n
        phi_cache(k) = phi_cache(k) + basis(k,t) * local_phi(t)
      END DO
      DO k = 1, ip % n
        phi_cache(k) = phi_cache(k) + basis(k,t) * local_phi(t)
      END DO
    END DO
  end if

  y = sqrt(phi_cache(ipind))
  print *, 'y=', phi_cache, y



   
!-------------------------------------------------------------------------------
END FUNCTION SquareRoot
!-------------------------------------------------------------------------------
 
