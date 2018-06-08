module RegExpUtils ! {{{
  use iso_c_binding
  implicit none
  private

  INTEGER, PARAMETER :: DEFAULT_MAX_MATCH = 10

  type, bind(c) :: reCompiled_t
    type(c_ptr) :: reCompiled = C_NULL_PTR
    type(c_ptr) :: pcreExtra = C_NULL_PTR
  end type

  ! public :: RegCompile, RegMatch, RegGetSubstring, RegNumMatch

  type, public :: regex_t
    integer(kind=c_int), allocatable, public :: ovector(:)
    integer :: maxmatch
    type(reCompiled_t) :: reCompiled
    logical :: compiled = .false.
    logical :: matched = .false.
    character(kind=c_char, len=:), pointer :: match_str => NULL()
    contains 
    procedure, public :: Compile => RegCompile
    procedure, public :: Match => RegMatch
    procedure, public :: NumMatch => RegNumMatch
    procedure, public :: Init => RegInit
    procedure, public :: GetSub => RegGetSubstring
  end type
  
interface 
  subroutine pcre_f_compile(regstr, msglevel, recompiled) bind(C, name="pcre_f_compile") ! {{{
    import
    character(kind=c_char), intent(in) :: regstr(*)
    integer(kind=c_int), value :: msglevel
    type(reCompiled_t) :: recompiled
  end subroutine ! }}}

  subroutine pcre_f_match_compiled(Str, StrLen, ovector, maxmatch, reCompiled) bind(C) ! {{{
    import
    character(kind=c_char), intent(in) :: Str(*)
    integer(kind=c_int), value :: StrLen
    integer(kind=c_int) :: ovector(*)
    integer(kind=c_int), value :: maxmatch
    type(reCompiled_t) :: reCompiled
  end subroutine  ! }}}
end interface

contains

  subroutine RegInit(regex, maxmatch_in)  ! {{{
    class(regex_t) :: regex
    integer :: maxmatch
    integer, optional :: maxmatch_in

    maxmatch = DEFAULT_MAX_MATCH
    if(present(maxmatch_in)) maxmatch = maxmatch_in

    regex % maxmatch = maxmatch
    if(allocated(regex % ovector) .and. size(regex % ovector,1) /= 3*maxmatch) then
      deallocate(regex % ovector)
      allocate(regex % ovector(maxmatch*3))
    end if
    if(.not. allocated(regex % ovector)) allocate(regex % ovector(maxmatch*3))
  end subroutine ! }}}

  subroutine RegCompile(regex, RegStr, MsgLevel) ! {{{
    class(regex_t) :: regex
    character(kind=c_char), intent(in) :: regstr(*)
    integer(kind=c_int) :: MsgLevel
    type(reCompiled_t) :: reCompiled
    
    call pcre_f_compile(regstr, MsgLevel, regex % reCompiled)
  end subroutine !}}}

  subroutine RegMatch(regex, Str, RegStr, MaxMatch)  ! {{{
    class(regex_t) :: regex
    character(kind=c_char,len=:), allocatable, intent(in), target :: str
    integer, optional :: MaxMatch
    character(kind=c_char), intent(in), optional :: regstr(*)

    regex % match_str => str

    if(present(RegStr)) then
      call regex % Compile(RegStr, 0)
    end if

    if (.not. C_ASSOCIATED(regex % reCompiled % reCompiled)) then
        print *, 'RegMatcher not compiled'
        return
    end if

    if(present(maxmatch)) then
      if(.not.(allocated(regex % ovector)) .or. regex % maxmatch < maxmatch) call regex % init(maxmatch)
    else 
      if (.not.(allocated(regex % ovector))) call regex % init(maxmatch)
    end if
    regex % ovector = 0
    call pcre_f_match_compiled(str, len(str), regex % ovector, regex % maxmatch, regex % reCompiled )
  end subroutine !}}}

  function RegGetSubstring(regex, N) result(substring) ! {{{
    class(regex_t) :: regex
    integer  :: N
    character(len=:), pointer :: substring

    integer :: subs_len, subs_start, subs_end

    associate(ovector => regex % ovector)
      subs_start = ovector(2*(N-1)+1)+1
      subs_end = ovector(2*(N-1)+2)
    end associate
    subs_len = subs_end - subs_start
    substring => regex % match_str(subs_start:subs_end)

  end function ! }}}

   function RegNumMatch(regex)  result(NumMatch) ! {{{
    class(regex_t) :: regex
    integer :: k, m, first, last, NumMatch

    NumMatch = 0

    associate(ovector => regex % ovector)
      do k = 1,size(ovector,1),2
        first = ovector(k)
        last = ovector(k+1)
        if (last - first < 1) exit
        NumMatch = NumMatch + 1
      end do
    end associate
  end function ! }}}

end module ! }}}
