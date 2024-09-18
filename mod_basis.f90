Module mod_basis

  use functions

  implicit none

! ========================================================================

  integer                         :: MaxBasis,MaxShells,MaxPrim,MaxContract
  integer                         :: nBasis
  integer                         :: nShells, nAtoms
  integer,allocatable             :: nShells_X(:)
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
! FINDLOC !!
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

  END Function read_Atom

!=============================================================

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

  END Function read_nAtom

!=============================================================

  Function nextShell(endfil)

  implicit none

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

  END Function nextShell

!=============================================================

  Subroutine read_Momentum()

  implicit none

  character(120)                  :: line
  integer                         :: offi, offj
  character(1)                    :: azimut

1 continue
    read(4,'(a)',end=5)line
    if (len_trim(line) == 0) goto 1
    if (line(1:1) == '!') goto 1
    call tolower(line)

! Read angular momentum

2 continue
    if (len_trim(line) == 0) goto 2
    if (line(1:1) == '!') goto 2
    call tolower(line)
    offi=notBlank(line,1)
    offj=nextBlank(line,offi)-1
    read(line(offi:offj),'(a)',err=5)azimut
    select case (azimut)
      case ('s')
        angular(nBasis,nShell(nBasis))=0
      case ('p')
        angular(nBasis,nShell(nBasis))=1
      case ('d')
        angular(nBasis,nShell(nBasis))=2
      case ('f')
        angular(nBasis,nShell(nBasis))=3
      case ('g')
        angular(nBasis,nShell(nBasis))=4
      case ('h')
        angular(nBasis,nShell(nBasis))=5
      case ('i')
        angular(nBasis,nShell(nBasis))=6
      case default
        goto 5
    end select

! Read numebr of primitives

  offi=notBlank(line,offj+1)
  offj=nextBlank(line,offi)-1
  read(line(offi:offj),*,err=5)nprim(nBasis,nShell(nBasis))

! Read number of contractions

  offi=notBlank(line,offj+1)
  offj=nextBlank(line,offi)-1
  read(line(offi:offj),*,err=5)ncontract(nBasis,nShell(nBasis))

  return

5 continue
  write(6,*)'read_Momentum: Error reading first line of shell'
  call exit(8)

  END Subroutine read_Momentum

!=============================================================

  Subroutine read_nMomentum()

  implicit none

  character(120)                  :: line
  integer                         :: offi, offj
  integer                         :: p,c
  integer                         :: i

1 continue
    read(4,'(a)',end=5)line
    if (len_trim(line) == 0) goto 1
    if (line(1:1) == '!') goto 1

! Read angular momentum

2 continue
    if (len_trim(line) == 0) goto 2
    if (line(1:1) == '!') goto 2
    offi=notBlank(line,1)
    offj=nextBlank(line,offi)-1

! Read numebr of primitives

  offi=notBlank(line,offj+1)
  offj=nextBlank(line,offi)-1
  read(line(offi:offj),*,err=6)p
  if (p > MaxPrim) MaxPrim=p

! Read number of contractions

  offi=notBlank(line,offj+1)
  offj=nextBlank(line,offi)-1
  read(line(offi:offj),*,err=7)c
  if (c > MaxContract) MaxContract=c

! Skip p lines

  do i=1,p
    read(4,*)
  enddo

  return

5 continue
  write(6,*)'read_nMomentum: Error reading first line of shell'
  call exit(8)
6 continue
  write(6,*)'read_nMomentum: Error reading number of primitives'
  call exit(8)
7 continue
  write(6,*)'read_nMomentum: Error reading number of contractions'
  call exit(8)

  END Subroutine read_nMomentum

!=============================================================

  Subroutine read_Exponent()

  implicit none

  character(120)                  :: line
  integer                         :: i, j, offi, offj

  do i=1,nprim(nBasis,nShell(nBasis))
    read(4,'(a)',err=5)line
    offi=notBlank(line,1)
    offj=nextBlank(line,offi)-1
    read(line(offi:offj),*,end=5,err=5)alpha(nBasis,nShell(nBasis),i)
    do j=1,ncontract(nBasis,nShell(nBasis))
      offi=notBlank(line,offj+1)
      offj=nextBlank(line,offi)-1
      read(line(offi:offj),*,end=5,err=5)contract(nBasis,nShell(nBasis),i,j)
    enddo
  enddo

  return

5 continue
  write(6,*)'read_Exponent: Error reading primitives'
  call exit(8)

  END Subroutine read_Exponent

!=============================================================

  Subroutine parse_Basis()

  implicit none

  logical                         :: endfil

  rewind (4)
  nBasis=1
  nShell(nBasis)=1

! Reading atom symbol

1 continue
  if (.not. read_Atom()) then
    nBasis=nBasis-1
    return
  endif

! Reading shell

2 continue
  call read_Momentum()

! Read exponents, and contraction coefficients

  call read_Exponent()

! Next shell or basis ?
  if (nextShell(endfil)) then
    if (endfil) return
    nShell(nBasis)=nShell(nBasis)+1
    goto 2
  else
    if (endfil) return
    nBasis=nBasis+1
    nShell(nBasis)=1
    goto 1
  endif

  END Subroutine parse_Basis

!=============================================================

  Subroutine parse_nBasis()

  implicit none

  logical                         :: endfil
  integer                         :: s

  rewind (4)
  MaxBasis=1
  s=1

! Reading atom symbol

1 continue
  if (.not. read_nAtom()) then
    MaxBasis=MaxBasis-1
    return
  endif
! Reading shell

2 continue
  call read_nMomentum()

! Next shell or basis ?
  if (nextShell(endfil)) then
    if (endfil) return
    s=s+1
    goto 2
  else
    if (endfil) return
    MaxBasis=MaxBasis+1
    if (s > MaxShells) MaxShells=S
    s=1
    goto 1
  endif

  END Subroutine parse_nBasis


END Module mod_basis
