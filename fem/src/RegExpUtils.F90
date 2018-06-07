module RegExpUtils ! {{{
  use iso_c_binding
  implicit none
  private

  type, public, bind(c) :: reCompiled_t
    type(c_ptr) :: reCompiled = C_NULL_PTR
    type(c_ptr) :: pcreExtra = C_NULL_PTR
  end type

  public :: RegCompile, RegMatch
  
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

  function RegCompile(RegStr, MsgLevel) result(reCompiled) ! {{{
    character(kind=c_char), intent(in) :: regstr(*)
    integer(kind=c_int) :: MsgLevel
    type(reCompiled_t) :: reCompiled
   
    call pcre_f_compile(regstr, MsgLevel, reCompiled)
  end function !}}}

  function RegMatch(Str, reCompiled, MaxMatch) result(matchvec) ! {{{
    character(kind=c_char,len=:), allocatable, intent(in) :: str
    type(reCompiled_t) :: reCompiled
    integer(kind=c_int) :: matchvec(3*MaxMatch)
    integer :: MaxMatch

    if (.not. C_ASSOCIATED(reCompiled % reCompiled)) then
        print *, 'RegMatcher not compiled'
        return
    end if

    call pcre_f_match_compiled(str, len(str), matchvec, maxmatch, reCompiled )
  end function !}}}

end module ! }}}
