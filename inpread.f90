
Module functions
  implicit none
Contains
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

Module molecules
  implicit none
  
  Type xyz_data
    character(len=2),allocatable                  :: atom(:)
    real*8,allocatable                            :: geom(:,:)
    
  END Type xyz_data

  Type basis_data
  ! Here, I will add basis set data
  END Type basis_data

  Type molecule_attribs
    character(len=1)                              :: specifier
    integer                                       :: n_atoms=0
    real*8                                        :: Tr(3)=0.0d0,Ro(3)=0.0d0
    real*8                                        :: energy(4)=(/0.0d0,0.0d0,1.0d0,1.0d0/)
    type(xyz_data)                                :: xyz

  END Type molecule_attribs

  type(molecule_attribs),allocatable              :: molecule(:)


END Module molecules
  

Module input
  implicit none

  character(len=20),allocatable                   :: flag(:)
  character(len=120),allocatable                  :: comment(:)
  Type input_logic
    logical                                       :: is_critical=.FALSE., is_present=.FALSE.
    character(len=40)                             :: name,location="the specified location"

    Contains
      procedure, public  :: throw_error => logic_throw_error 
  END type input_logic
  Contains
    Subroutine logic_throw_error(this,assumption)
      class(input_logic), intent(inout)           :: this
      character(len=*), intent(in),optional       :: assumption
     
      IF (.NOT.this%is_present) then
        IF (this%is_critical) THEN
          PRINT*, "Error: ", this%name(1:len_trim(this%name)), " is not present in ",this%location(1:len_trim(this%location)),". ", assumption
          CALL EXIT(1)
        END IF
        PRINT*, "Warning: ", this%name(1:len_trim(this%name)), " is not present in ",this%location(1:len_trim(this%location)),". ", assumption
      END IF
      
    END subroutine logic_throw_error
END Module input

Program InpRead
  use functions
  use input
  use molecules
  use mod_table
  implicit none

  integer                                         :: start_position=0, end_position=0
  integer                                         :: n_molecules
  integer                                         :: i,j,k,l 
  character(len=2)                                :: element
  character(len=5)                                :: type
  character(len=6)                                :: signal
  character(len=8)                                :: inpfile
  character(len=120)                              :: line
  character(len=1),dimension(3)                   :: axis_letter 
  character(len=1),dimension(4)                   :: state_letter
  character(len=1),allocatable                    :: allowed_ends(:)
  character(len=4),dimension(5)                   :: allowed_types

  type(input_logic)                               :: basis_logic,flags_logic,inp_logic
  type(input_logic),allocatable                   :: geom_logic(:),trro_logic(:),energy_logic(:)

  character(len=40)                               :: debug(3)
  character(len=1)                                :: char_debug
 
  n_molecules=3
  allowed_types = (/'geom','ener','trro','basi','flag'/)
  state_letter = (/'s','t','c','a'/)
  axis_letter = (/'x','y','z'/)

  inpfile = 'test.inp' !this will be read
  inp_logic%name=inpfile
  inp_logic%is_critical=.TRUE.
  basis_logic%name='Basis set definition'
  basis_logic%location=inpfile
  flags_logic%name='Calculation definition'
  flags_logic%location=inpfile
  basis_logic%is_critical=.TRUE.
  
  OPEN(3,file=inpfile,status='old',err=8)

  DO WHILE(.NOT.(flags_logic%is_present))
    READ(3,'(a)',end=1) line
    line=ADJUSTL(line)
    call tolower(line)
    IF (line(1:6).ne.'flags:') CYCLE
    end_position=INDEX(line,'dimer')
    IF(end_position.eq.0) CYCLE
    PRINT*, "Assuming this is a dimer!"
    flags_logic%is_present=.TRUE.
    n_molecules=2
    EXIT
  ENDDO
 
1 CONTINUE    
  ALLOCATE(allowed_ends(n_molecules+1))
  DO i=97,96+n_molecules
    allowed_ends(i-96)=char(i)
  ENDDO
  allowed_ends(n_molecules+1)='s'

  REWIND(3)

  ALLOCATE(molecule(n_molecules),geom_logic(n_molecules),trro_logic(6*(n_molecules-1)),energy_logic(n_molecules*4))
  ! Important allocation. Must be 0 in length!
  ALLOCATE(comment(0),flag(0))

  DO i=1,n_molecules
    geom_logic(i)%name='Geometry of '//char(uppercase(allowed_ends(i)))
    geom_logic(i)%is_critical=.TRUE.
    geom_logic(i)%location=inpfile
    molecule(i)%specifier=allowed_ends(i)
    DO j=1,4
      energy_logic((i-1)*4+j)%name=char(uppercase(state_letter(j)))//' energy of '//char(uppercase(allowed_ends(i)))
      energy_logic((i-1)*4+j)%location=inpfile
    ENDDO    
  ENDDO
  DO i=1,n_molecules-1
   DO j=1,3
      trro_logic((i-1)*3+j)%name=  'T'//axis_letter(j)//' of '//char(uppercase(allowed_ends(i+1)))
      trro_logic((i-1)*3+j+6)%name='R'//axis_letter(j)//' of '//char(uppercase(allowed_ends(i+1)))
      trro_logic((i-1)*3+j)%location=inpfile
      trro_logic((i-1)*3+j+6)%location=inpfile
    ENDDO
  ENDDO
  

  ! READING BLOCK
  ! Prunes the comments from other stuff, sets types.
2 READ(3,'(a)',end=3) line
  line=ADJUSTL(line)
  signal=line(1:6)
  CALL tolower(signal)
  IF ((signal.eq.'print:').OR.(signal(1:1).eq.'!')) THEN
    i=7
    IF(signal(1:1).eq.'!') i=2
    ! Behold, Fortran 2003! 
    comment=[comment,line(i:120)]
    GOTO 2
  ENDIF
  IF (.NOT.( (signal(6:6).eq.':') .AND. (ANY(allowed_types.eq.signal(1:4))) .AND. (ANY(allowed_ends.eq.signal(5:5))) )) GOTO 2
  CALL tolower(line)
  type=signal(1:5)
  end_position=INDEX(line,';')
  line=ADJUSTL(line(7:end_position))

  ! PROCESSING BLOCK
  ! I can save some lines by making a function
  SELECT CASE (type(1:4))
  CASE ('geom')
    IF (.NOT.(ANY(allowed_ends(1:SIZE(allowed_ends)-1).eq.type(5:5)))) GOTO 2
    i=FINDLOC(allowed_ends(1:SIZE(allowed_ends)-1),type(5:5),dim=1)
    geom_logic(i)%is_present=.TRUE.
    end_position=check_end(line)
    line=ADJUSTL(line(1:end_position-1))
    ! Read the data here
    ! HANDLE WRONG INPUT HERE
    j=30+i
    OPEN(j,file=line,status='old',err=81)
    READ(j,'(i5)',err=82) molecule(i)%n_atoms
    IF(molecule(i)%n_atoms.eq.0) GOTO 82
    ALLOCATE(molecule(i)%xyz%atom(molecule(i)%n_atoms),molecule(i)%xyz%geom(3,molecule(i)%n_atoms))
    READ(j,'(i5)',err=82)
    DO k=1,molecule(i)%n_atoms
      READ(j,'(a2)',end=82,advance='no') molecule(i)%xyz%atom(k)
      element=ADJUSTR(molecule(i)%xyz%atom(k))
      IF(FINDLOC(table,element,dim=1).eq.0) GOTO 82
      READ(j,'(a)',end=82) line 
    !SAVE
      DO l=1,3
        line=ADJUSTL(line)
        end_position=check_end(line)
        READ(line(1:end_position),*,err=82,end=82) molecule(i)%xyz%geom(l,k)
        line=ADJUSTL(line(1:end_position))
      ENDDO
    ENDDO

    CLOSE(j)

  CASE ('ener')
    IF (.NOT.(ANY(allowed_ends(1:SIZE(allowed_ends)-1).eq.type(5:5)))) GOTO 2
    i=FINDLOC(allowed_ends(1:SIZE(allowed_ends)-1),type(5:5),dim=1)

    DO j=1,4
      energy_logic((i-1)*4+j)%is_present=.TRUE.
      start_position=INDEX(line,state_letter(j))+2
      IF(start_position.eq.2) THEN
        energy_logic((i-1)*4+j)%is_present=.FALSE.
        CYCLE
      ENDIF
      end_position=check_end(line,start_position)+start_position-2
      READ(line(start_position:end_position),*,err=83) molecule(i)%energy(j)
    ! READ the data here
    ENDDO
    
  CASE ('trro')
    IF (.NOT.(ANY(allowed_ends(2:SIZE(allowed_ends)-1).eq.type(5:5)))) GOTO 2
    i=FINDLOC(allowed_ends,type(5:5),dim=1)
    
    DO j=1,3
      trro_logic((i-2)*3+j)%is_present=.TRUE.
      trro_logic((i-2)*3+j+6)%is_present=.TRUE.

      start_position=INDEX(line,'t'//axis_letter(j))+3
      IF(start_position.eq.3) trro_logic((i-2)*3+j)%is_present=.FALSE.
      end_position=check_end(line,start_position)+start_position-2
    ! READ T data here
    
      start_position=INDEX(line,'r'//axis_letter(j))+3
      IF(start_position.eq.3) trro_logic((i-2)*3+j+6)%is_present=.FALSE.
      end_position=check_end(line,start_position)+start_position-2
    ! READ R data here

    ENDDO

  CASE ('basi')
    IF ('s'.ne.type(5:5)) GOTO 2
    end_position=check_end(line)
    IF(ANY((/0,1/).eq.end_position)) GOTO 2
    basis_logic%is_present=.TRUE.
    line=ADJUSTL(line(1:end_position-1))
    ! Read the basis_logic set here

  CASE ('flag')
    IF ('s'.ne.type(5:5)) GOTO 2
    end_position=check_end(line)
    IF(ANY((/0,1/).eq.end_position)) GOTO 2
    flags_logic%is_present=.TRUE.
    line=ADJUSTL(line(1:end_position-1))
    ! Read flags_logic here
    
  END SELECT

  GOTO 2
3 CONTINUE
  ! CHECKING BLOCK
  ! Check if the file contains the required data
  ! Throw errors (=> EXIT) after warnings
  CALL flags_logic%throw_error("Assuming this is a SP calculation with a trimer!")
  DO i=1,n_molecules
    DO j=1,4
      WRITE(line,'(f8.3)') molecule(i)%energy(j)
      line=ADJUSTL(line)
      CALL energy_logic((i-1)*4+j)%throw_error("Using the default value ("//line(1:4)//" eV)!")
    ENDDO
  ENDDO
  DO i=1,n_molecules-1
    DO j=1,3
      CALL trro_logic((i-1)*3+j)%throw_error("Assuming 0.0 Å!")
      CALL trro_logic((i-1)*3+j+3)%throw_error("Assuming 0.0°!")
    ENDDO
  ENDDO
  DO i=1,n_molecules
    CALL geom_logic(i)%throw_error()
  ENDDO
  CALL basis_logic%throw_error()
  DO i=1,SIZE(comment)
    PRINT'(x,2a)', "User comment:", comment(i)
  ENDDO


  CLOSE(3)
  GOTO 9

! ERROR BLOCK
8   call inp_logic%throw_error() 
81  geom_logic(i)%is_present=.FALSE.
    geom_logic(i)%location=line(1:len_trim(line))
    call geom_logic(i)%throw_error("File not found!")
82  geom_logic(i)%is_present=.FALSE.
    geom_logic(i)%location=line(1:len_trim(line))
    call geom_logic(i)%throw_error("Is the xyz format correct?")
83  energy_logic((i-1)*4+j)%is_present=.FALSE.
    energy_logic((i-1)*4+j)%is_critical=.TRUE.
    energy_logic((i-1)*4+j)%location=line(start_position:end_position)
    call energy_logic((i-1)*4+j)%throw_error("Is the energy format correct?")
 
9 CONTINUE  
End Program InpRead

Subroutine tolower(line)

  implicit none

! ========================================================================
! Convert upercase characters in line into lowercase

  character(*),intent(inout)      :: line
  integer                         :: jgap
  integer                         :: i

  jgap=ichar('a')-ichar('A')

  do i=1,len_trim(line)
    if (line(i:i) <= 'Z') then
      if (line(i:i) >= 'A') line(i:i)=char(ichar(line(i:i))+jgap)
    endif
  enddo

  return

End Subroutine tolower

