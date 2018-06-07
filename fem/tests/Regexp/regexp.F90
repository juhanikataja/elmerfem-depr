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
  integer(kind=c_int), parameter :: maxmatch=10
  integer(kind=c_int) :: ovector(3*maxmatch)
  type(reCompiled_t) :: comp
  character(len=MAX_NAME_LEN), allocatable :: str_in
  character(len=:, kind=c_char), allocatable :: str, last_str
  character(len=:, kind=c_char), allocatable :: regstr
  logical :: found

  regstr = "(.+) \| (.+) -- (.+) \((\d+)\)"  // c_null_char
  ovector = 0
  ! regstr = "\((\d+)\)"  // c_null_char
  str = ListGetString(GetSolverParams(), "str to parse", found)

  ! str = trim(str_in) // c_null_char
  if (found) then
    comp = RegCompile(regstr, 0)
    ovector = RegMatch(str, comp, maxmatch)
  end if
  print *, str
  print *, regstr

  block
    integer :: i
    do i = 1, maxmatch, 2
      associate (startpos => ovector(i)+1, &
            sublen => ovector(i+1)-ovector(i))
        if(sublen < 1) exit
        last_str = str(startpos:startpos+sublen-1)
        print *, 'startpos:', startpos, 'sublen: ', sublen
        print *, last_str
      end associate
    end do
    if(last_str /= "3") stop
  end block


!------------------------------------------------------------------------------
   END SUBROUTINE RegexpTester
!------------------------------------------------------------------------------
