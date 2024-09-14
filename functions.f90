Module functions
  implicit none

  PUBLIC :: nextBlank, notBlank

Contains

Function nextBlank(line,offset)

  implicit none

! ========================================================================

  integer                         :: nextBlank
  character(*),intent(in)         :: line
  integer,intent(in)              :: offset
  integer                         :: i

  nextBlank=offset
  do i=offset,len_trim(line)
    if (line(i:i) == ' ') then
      nextBlank=i
      return
    endif

! End of line

    if (i == len_trim(line)) then
      nextBlank=i+1
    endif
  enddo

End Function nextBlank

Function notBlank(line,offset)

  implicit none

! ========================================================================

  integer                         :: notBlank
  character(*),intent(in)         :: line
  integer,intent(in)              :: offset
  integer                         :: i

  do i=offset,len_trim(line)
    if (line(i:i) /= ' ') then
      notBlank=i
      return
    endif
  enddo

End Function notBlank

Function check_end(line,start) Result(end_char)
  character(len=120),intent(in)                   :: line
  character(len=120)                              :: input
  integer,intent(in),optional                     :: start
  integer                                         :: end_char
 
  input=ADJUSTL(line)
  IF(present(start)) input=ADJUSTL(line(start:120))
  
  end_char=INDEX(input,' ')
  IF((end_char.eq.0).OR.(end_char.gt.INDEX(input,';')).AND.(INDEX(input,';').gt.0)) end_char=INDEX(input,';')

  RETURN
END Function check_end

Function uppercase(input) result(output)
  character(len=1)                                :: input
  integer                                         :: output

  output=ichar(input)-ichar('a')+ichar('A')

END Function uppercase
END Module functions
