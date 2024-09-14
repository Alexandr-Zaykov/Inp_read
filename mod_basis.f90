Module mod_basis

  use functions

  implicit none

! ========================================================================

  integer                         :: MaxBasis,MaxShells,MaxPrim,MaxContract
  integer                         :: nBasis
  integer                         :: nShells, nShells_A, nShells_B, nShells_C
  integer,allocatable             :: b_atom(:)
  integer,allocatable             :: nShell(:)
  integer,allocatable             :: angular(:,:)
  integer,allocatable             :: nprim(:,:)
  integer,allocatable             :: ncontract(:,:)
  real*8,allocatable              :: alpha(:,:,:)
  real*8,allocatable              :: contract(:,:,:,:)
  integer,allocatable             :: map(:)

  contains
  Function read_Atom()

  use mod_table

  implicit none

! ========================================================================

  logical                         :: read_Atom
  character(120)                  :: line
  integer                         :: offi, offj
  character(1)                    :: t1
  character(2)                    :: t2
  integer                         :: j

  read_Atom=.true.

1 continue
    read(4,'(a)',end=6,err=5)line
    if (len_trim(line) == 0) goto 1
    if (line(1:1) == '!') goto 1
    call tolower(line)

! Read atom symbol

  offi=notBlank(line,1)
  offj=nextBlank(line,offi)-1
  do j=1,103
    if (offi == offj) then
      if (Table(j)(1:1) == ' ') then
        t1=Table(j)(2:2)
        call tolower(t1)
        if (t1 == line(offi:offj)) then
          b_atom(nBasis)=j
          return
        endif
      endif
    else
      t2=Table(j)
      call tolower(t2)
      if (t2 == line(offi:offj)) then
        b_atom(nBasis)=j
        return
      endif
    endif
  enddo

5 continue
  write(6,*)'read_Atom: End of file while reading atomic symbol'
  call exit(8)
6 continue
  read_Atom=.false.

  return

End Function read_Atom

Function read_nAtom()

  use mod_table

  implicit none

! ========================================================================

  logical                         :: read_nAtom
! integer                         :: notBlank, nextBlank
  character(120)                  :: line
  integer                         :: offi, offj
  character(1)                    :: t1
  character(2)                    :: t2
  integer                         :: j

  read_nAtom=.true.

1 continue
    read(4,'(a)',end=6,err=5)line
    if (len_trim(line) == 0) goto 1
    if (line(1:1) == '!') goto 1
    call tolower(line)

! Read atom symbol

  offi=notBlank(line,1)
  offj=nextBlank(line,offi)-1
  do j=1,103
    if (offi == offj) then
      if (Table(j)(1:1) == ' ') then
        t1=Table(j)(2:2)
        call tolower(t1)
        if (t1 == line(offi:offj)) then
          return
        endif
      endif
    else
      t2=Table(j)
      call tolower(t2)
      if (t2 == line(offi:offj)) then
        return
      endif
    endif
  enddo

5 continue
  write(6,*)'read_nAtom: End of file while reading atomic symbol'
  call exit(8)
6 continue
  read_nAtom=.false.

  return

End Function read_nAtom

Function nextShell(endfil)

  implicit none

! ========================================================================

  logical                         :: nextShell
  logical                         :: endfil
  character(120)                  :: line

  nextShell=.true.
  endfil=.false.

1 continue
  read(4,'(a)',end=5)line
  if (len_trim(line) == 0) goto 1
  if (line(1:1) == '!') goto 1
  if (index(line,'****') /= 0) then
    nextShell=.false.
    return
  endif
  backspace(4)
  return

5 continue
  write(*,*) "End of file reached."
  endfil=.true.
  return

End Function nextShell


End Module mod_basis
