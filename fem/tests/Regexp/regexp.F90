!------------------------------------------------------------------------------
  SUBROUTINE RegexpTester( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
!******************************************************************************
!
!  Test regular expressions
!
!  ARGUMENTS:
!
!  TYPE(Model_t) :: Model,  
!     INPUT: All model information (mesh, materials, BCs, etc...)
!
!  TYPE(Solver_t) :: Solver
!     INPUT: Linear equation solver options
!
!  DOUBLE PRECISION :: dt,
!     INPUT: Timestep size for time dependent simulations
!
!  LOGICAL :: TransientSimulation
!     INPUT: Steady state or transient simulation
!
!******************************************************************************
  USE DefUtils
  USE ISO_C_BINDING
  USE RegExpUtils

  IMPLICIT NONE
  !------------------------------------------------------------------------------
  TYPE(Model_t) :: Model
  TYPE(Solver_t), TARGET :: Solver
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
  !------------------------------------------------------------------------------ 
  type(regex_t) :: regex
  character(len=MAX_NAME_LEN), allocatable :: str_in
  character(len=:, kind=c_char), allocatable :: str, last_str
  character(len=:, kind=c_char), allocatable :: regstr
  logical :: found

  regstr = "(.+) \| (.+) -- (.+) \((\d+)\)"  // c_null_char
  str = ListGetString(GetSolverParams(), "str to parse", found)

  if (found) then
    call regex % match(str, regstr, 30)
  end if
  if (regex % getsub(regex % nummatch()) /= "3") stop

!------------------------------------------------------------------------------
   END SUBROUTINE RegexpTester
!------------------------------------------------------------------------------
