Module molecules
  use functions
  implicit none
 
  character(len=5),allocatable                    :: flag(:)
  integer                                         :: n_molecules=3   ! the goal is a trimer for now, and so it is a default value
  real*8                                          :: lamda=0.2 ! eV

  Type xyz_data
    character(len=2),allocatable                  :: atom(:)
    real*8,allocatable                            :: geom(:,:)
    integer,allocatable                           :: nuc_charge(:)
    
  END Type xyz_data

  Type basis_data
  ! Here, I will add basis set data
  END Type basis_data

  Type molecule_attribs
    character(len=1)                              :: specifier
    integer                                       :: n_atoms=0, orb_rot(2)=0
    real*8                                        :: Tr(3)=0.0d0,Ro(3)=0.0d0
    real*8                                        :: energy(4)=(/0.0d0,0.0d0,1.0d0,1.0d0/)
    type(xyz_data)                                :: xyz
    Contains
    procedure, public  :: show_geom => molecule_print_geometry 

  END Type molecule_attribs

  type(molecule_attribs),allocatable              :: molecule(:)

  Contains
  Subroutine molecule_print_geometry(this)
    class(molecule_attribs), intent(in)          :: this
    integer                                      :: i_priv

    PRINT*
    PRINT*, "Molecule ", uppercase(this%specifier) 
    PRINT '(x,a,i6)', "Number of atoms: ", this%n_atoms
    PRINT*, "Atomic coordinates [Å]:"
    DO i_priv=1,this%n_atoms
      IF(len_trim(adjustl(this%xyz%atom(i_priv))).eq.1) THEN
        PRINT '(x,a,3f12.6)', this%xyz%atom(i_priv), this%xyz%geom(i_priv,1),&
                             &this%xyz%geom(i_priv,2), this%xyz%geom(i_priv,3)
        CYCLE
      ENDIF
      PRINT '(2x,a,f11.6,2f12.6)', this%xyz%atom(i_priv), this%xyz%geom(i_priv,1),&
                                  &this%xyz%geom(i_priv,2), this%xyz%geom(i_priv,3)
    ENDDO

  END subroutine molecule_print_geometry


END Module molecules
  

Module input
  implicit none

  character(len=120),allocatable                  :: comment(:)
  character(len=1),dimension(3),parameter         :: axis_letter=(/'x','y','z'/)
  character(len=1),dimension(4),parameter         :: state_letter=(/'s','t','c','a'/)
  character(len=4),dimension(7),parameter         :: allowed_types=(/'geom','ener','trro','basi','flag','reor','lamd'/)
  character(len=1),allocatable                    :: allowed_ends(:)

  Type input_logic
    logical                                       :: is_critical=.FALSE., is_present=.FALSE.
    character(len=40)                             :: name="A thing"
    character(len=200)                            :: location="the specified location"

    Contains
      procedure, public  :: throw_error => logic_throw_error 
  END type input_logic

  type(input_logic)                               :: basis_logic,flags_logic,inp_logic,lamda_logic
  type(input_logic),allocatable                   :: geom_logic(:),trro_logic(:),energy_logic(:),orbrot_logic(:)

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

Module pipes

use,intrinsic :: iso_c_binding

implicit none

private
! Pipes module from DEGENERATE CONIC – https://degenerateconic.com/fortran-c-interoperability.html
! Cheers!
interface

    function popen(command, mode) bind(C,name='popen')
    import :: c_char, c_ptr
    character(kind=c_char),dimension(*) :: command
    character(kind=c_char),dimension(*) :: mode
    type(c_ptr) :: popen
    end function popen

    function fgets(s, siz, stream) bind(C,name='fgets')
    import :: c_char, c_ptr, c_int
    type (c_ptr) :: fgets
    character(kind=c_char),dimension(*) :: s
    integer(kind=c_int),value :: siz
    type(c_ptr),value :: stream
    end function fgets

    function pclose(stream) bind(C,name='pclose')
    import :: c_ptr, c_int
    integer(c_int) :: pclose
    type(c_ptr),value :: stream
    end function pclose

end interface

public :: c2f_string, get_command_as_string

contains

!**********************************************
! convert a C string to a Fortran string
!**********************************************
function c2f_string(c) result(f)

    implicit none

    character(len=*),intent(in) :: c
    character(len=:),allocatable :: f

    integer :: i

    i = index(c,c_null_char)

    if (i<=0) then
        f = c
    else if (i==1) then
        f = ''
    else if (i>1) then
        f = c(1:i-1)
    end if

end function c2f_string

!**********************************************
! return the result of the command as a string
!**********************************************
function get_command_as_string(command) result(str)

    implicit none

    character(len=*),intent(in) :: command
    character(len=:),allocatable :: str

    integer,parameter :: buffer_length = 1000

    type(c_ptr) :: h
    integer(c_int) :: istat
    character(kind=c_char,len=buffer_length) :: line

    str = ''
    h = c_null_ptr
    h = popen(command//c_null_char,'r'//c_null_char)

    if (c_associated(h)) then
        do while (c_associated(fgets(line,buffer_length,h)))
            str = str//c2f_string(line)
        end do
        istat = pclose(h)
    end if

end function get_command_as_string

END Module pipes

Program InpRead
  use functions
  use pipes
  use input
  use molecules
  use mod_table
  use mod_basis
  implicit none

  integer                                         :: start_position=0, end_position=0, line_num=0, num_com=0, num_flg=0
  integer                                         :: i=0,j=0,k=0,l=0 ! i:molecules j:other specifiers (state/axis) k,l:misc; iff i used=>move it to l
  real*8                                          :: time_start, time_end
  character(len=2)                                :: element
  character(len=5)                                :: type
  character(len=6)                                :: signal
  character(len=10)                               :: tmp
  character(len=120)                              :: line
  character(len=:),allocatable                    :: sys_message,inpfile,basis_name,library_location,full_path
!=========================DEBUG===============================

! NOTE:
! I am trying to interweave OOP FORTRAN with the "standard" procedural FORTRAN.
! OOP is great for its flexibility, but the procedural FORTRAN allows for nifty optimization.
 
  PRINT '(t20,a)',"╭━━━┳━━━┳━━━┳━━━┳━━━╮" 
  PRINT '(t20,a)',"┃╭━╮┃╭━╮┃╭━╮┃╭━╮┃╭━╮┃" 
  PRINT '(t20,a)',"┃╰━━┫╰━╯┃┃╱┃┃╰━╯┃┃╱╰╯"
  PRINT '(t20,a)',"╰━━╮┃╭━━┫╰━╯┃╭╮╭┫┃╱╭╮"
  PRINT '(t20,a)',"┃╰━╯┃┃╱╱┃╭━╮┃┃┃╰┫╰━╯┃"
  PRINT '(t20,a)',"╰━━━┻╯╱╱╰╯╱╰┻╯╰━┻━━━╯"
  !i = COMMAND_ARGUMENT_COUNT() => Fortran 2003, instead of IARGC, error handling for more than a single input file?, flags?
  ! GET_COMMAND_ARGUMENT does not allocate deferred type characters until the 2023 standard.  
  CALL GET_COMMAND_ARGUMENT(1,length=i)
  IF (i .lt. 1) GOTO 8
  ALLOCATE(character(i) :: inpfile)
  CALL GET_COMMAND_ARGUMENT(1,value=inpfile)

! INITIALIZE
  basis_name='NONE'
  full_path='NONE'
  inp_logic%name=inpfile
  basis_logic%name='Basis set definition'
  flags_logic%name='Calculation definition'
  lamda_logic%name='Reorganization energy value'
  basis_logic%location=inpfile
  flags_logic%location=inpfile
  lamda_logic%location=inpfile
  basis_logic%is_critical=.TRUE.
  inp_logic%is_critical=.TRUE.
  ! lambda is not critical, defaults (0.2 eV) initialized in "molecules" module
  
  OPEN(3,file=inpfile,status='old',err=80)

  DO WHILE(.NOT.(flags_logic%is_present))
    READ(3,'(a)',end=1) line
    line=ADJUSTL(line)
    CALL tolower(line)
    IF (line(1:6).ne.'flags:') CYCLE
    ! THIS WILL NEED TO BE EXPANDED n=2,3,4,...
    end_position=INDEX(line,'dimer')
    IF(end_position.eq.0) CYCLE
    n_molecules=2
    EXIT
  ENDDO
 
1 CONTINUE    
  ALLOCATE(allowed_ends(n_molecules+1))
  DO i=97,96+n_molecules
    ! a=97, z=122, allowed_ends=(/a,b,c,.../)
    allowed_ends(i-96)=char(i)
  ENDDO
  allowed_ends(n_molecules+1)='s'

  REWIND(3)

  ALLOCATE(molecule(n_molecules),geom_logic(n_molecules),orbrot_logic(n_molecules),&
          &trro_logic(6*(n_molecules-1)),energy_logic(n_molecules*4))

  DO i=1,n_molecules
    geom_logic(i)%name='Geometry of molecule '//uppercase(allowed_ends(i)) 
    geom_logic(i)%is_critical=.TRUE.
    geom_logic(i)%location=inpfile
    orbrot_logic(i)%name='Orbital rotation of molecule '//uppercase(allowed_ends(i)) 
    orbrot_logic(i)%is_critical=.FALSE.
    orbrot_logic(i)%location=inpfile
    molecule(i)%specifier=allowed_ends(i)
    DO j=1,4
      energy_logic((i-1)*4+j)%name=uppercase(state_letter(j))//' energy of '//uppercase(allowed_ends(i))
      energy_logic((i-1)*4+j)%location=inpfile
    ENDDO    
  ENDDO
  DO i=1,n_molecules-1
   DO j=1,3
      trro_logic((i-1)*n_molecules+j)%name=  'T'//axis_letter(j)//' of '//uppercase(allowed_ends(i+1))
      trro_logic((i-1)*n_molecules+j+3*(n_molecules-1))%name='R'//axis_letter(j)//' of '//uppercase(allowed_ends(i+1))
      trro_logic((i-1)*n_molecules+j)%location=inpfile
      trro_logic((i-1)*n_molecules+j+3*(n_molecules-1))%location=inpfile
    ENDDO
  ENDDO
  

  ! READING BLOCK
  ! Prunes the comments from other stuff, sets types.
2 READ(3,'(a)',end=3) line
  line_num=line_num+1
  line=ADJUSTL(line)
  signal=line(1:6)
  CALL tolower(signal)
  IF ((signal.eq.'print:').OR.(signal(1:1).eq.'!')) THEN
    IF(.NOT.ALLOCATED(comment)) ALLOCATE(comment(0))
    i=7
    IF(signal(1:1).eq.'!') i=2
    num_com=num_com+1
    line=ADJUSTL(line(i:120))
    comment=[comment,line]
    GOTO 2
  ENDIF
  IF (.NOT.( (signal(6:6).eq.':') .AND. (ANY(allowed_types.eq.signal(1:4))) .AND. (ANY(allowed_ends.eq.signal(5:5))) )) GOTO 2
  CALL tolower(line)
  type=signal(1:5)
  end_position=INDEX(line,';')
  IF(end_position.eq.0) GOTO 810
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
    ALLOCATE(molecule(i)%xyz%atom(molecule(i)%n_atoms),molecule(i)%xyz%geom(molecule(i)%n_atoms,3),&
            &molecule(i)%xyz%nuc_charge(molecule(i)%n_atoms))
    molecule(i)%xyz%geom=0.0d0
    READ(j,'(i5)',err=82)
    DO k=1,molecule(i)%n_atoms
      READ(j,'(a)',end=82) line 
      line=ADJUSTL(line)
      molecule(i)%xyz%atom(k)=line(1:2)
      line=ADJUSTL(line(3:120))
      element=ADJUSTR(molecule(i)%xyz%atom(k))
      molecule(i)%xyz%nuc_charge(k)=FINDLOC(table,element,dim=1)
      IF(molecule(i)%xyz%nuc_charge(k).eq.0) GOTO 82
    !SAVE
      DO l=1,3
        line=ADJUSTL(line)
        end_position=check_end(line)
        READ(line(1:end_position-1),*,err=82,end=82) molecule(i)%xyz%geom(k,l) ! This is not optimal, but not worth -> geom(l,k)
        line=ADJUSTL(line(end_position:120))
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
      ! TRANSLATIONS
      k=(i-2)*n_molecules+j
      start_position=INDEX(line,'t'//axis_letter(j))+3
      IF(start_position.ne.3) THEN
        trro_logic(k)%is_present=.TRUE.
        end_position=check_end(line,start_position)+start_position-2
        READ(line(start_position:end_position),*,err=85) molecule(i)%Tr(j)
      ENDIF

      !ROTATIONS
      k=k+3*(n_molecules-1)
      start_position=INDEX(line,'r'//axis_letter(j))+3
      IF(start_position.ne.3) THEN
        trro_logic(k)%is_present=.TRUE.
        end_position=check_end(line,start_position)+start_position-2
        READ(line(start_position:end_position),*,err=85) molecule(i)%Ro(j)
      ENDIF

    ENDDO

  CASE ('basi')
    IF ('s'.ne.type(5:5)) GOTO 2
    end_position=check_end(line)
    IF(ANY((/0,1/).eq.end_position)) GOTO 2
    basis_logic%is_present=.TRUE.
    basis_name=ADJUSTL(line(1:end_position-1))
    i=len_trim(get_command_as_string('printenv | grep "SPARC_LIB" | cut -c11-'))
    IF(i.eq.1) GOTO 86 ! curiously, this is not 0
    ALLOCATE(character(len=i) :: library_location)
    library_location=get_command_as_string('printenv | grep "SPARC_LIB" | cut -c11-')
    ! i-1 as it contains a newline character
    full_path=library_location(1:i-1)//'/basis_sets/'//basis_name
    OPEN(4,file=full_path,status='old',err=87)
    CALL parse_nbasis()
    ALLOCATE(b_atom(MaxBasis), nShell(MaxBasis+1),&
            &angular(MaxBasis, MaxShells), nprim(MaxBasis, MaxShells), ncontract(MaxBasis, MaxShells),&
            &alpha(MaxBasis, MaxShells, MaxPrim), contract(MaxBasis, MaxShells, MaxPrim, MaxContract))
    CALL parse_basis()
    CLOSE(4)
    

  CASE ('flag')
    IF ('s'.ne.type(5:5)) GOTO 2
    end_position=check_end(line)
    IF(11.lt.end_position) GOTO 84
    IF(ANY((/0,1/).eq.end_position)) GOTO 2
    flags_logic%is_present=.TRUE.
    DO WHILE(.TRUE.)
      num_flg=num_flg+1
      tmp=line(1:end_position-1)  
      IF(.NOT.(ANY(calc_flags.eq.tmp(1:5)))) GOTO 84
      IF(.NOT.ALLOCATED(flag)) ALLOCATE(flag(0))
      flag=[flag,tmp(1:5)]
      line=ADJUSTL(line(end_position:120))
      end_position=check_end(line)
      IF(11.lt.end_position) GOTO 84
      IF(ANY((/0,1/).eq.end_position)) GOTO 2
    ENDDO

  CASE ('reor')
    IF (.NOT.(ANY(allowed_ends(1:SIZE(allowed_ends)-1).eq.type(5:5)))) GOTO 2
    i=FINDLOC(allowed_ends(1:SIZE(allowed_ends)-1),type(5:5),dim=1)

    orbrot_logic(i)%is_present=.TRUE.
    end_position=check_end(line) 
    IF(ANY((/0,1/).eq.end_position)) GOTO 88

    DO j=1,2
      READ(line(1:end_position-1),*,end=88,err=88) molecule(i)%orb_rot(j)
      line=ADJUSTL(line(end_position:120))
      end_position=check_end(line) 
    ENDDO
  
  CASE ('lamd')
    IF ('a'.ne.type(5:5)) GOTO 2
    ! THIS IS AN INJECTION POINT FOR DIFFERENT LAMDAS, USE FINDLOC TO DEFINE LAMDA, LAMDB, LAMDC,...
    
    lamda_logic%is_present=.TRUE.
    end_position=check_end(line)
    IF(ANY((/0,1/).eq.end_position)) GOTO 89
    READ(line(1:end_position-1),*,err=89) lamda
    
  END SELECT
  GOTO 2

3 CONTINUE
  CLOSE(3)
  ! CHECKING BLOCK/PRINTING BLOCK
  ! Check if the file contains the required data
  ! Throw errors (=> EXIT) after warnings
  ! WARNINGS
  CALL flags_logic%throw_error("Assuming this is a SP calculation with a trimer!")
  DO i=1,n_molecules
    DO j=1,4
      WRITE(line,'(f8.3)') molecule(i)%energy(j)
      line=ADJUSTL(line)
      CALL energy_logic((i-1)*4+j)%throw_error("Default value: "//line(1:4)//" eV!")
    ENDDO
  ENDDO
  DO i=2,n_molecules
    DO j=1,3
      WRITE(line,'(f4.1)') molecule(i)%Tr(j)
      line=ADJUSTL(line)
      CALL trro_logic((i-2)*n_molecules+j)%throw_error("Default value: "//line(1:LEN_TRIM(line))//" Å!")
      WRITE(line,'(f4.1)') molecule(i)%Ro(j)
      line=ADJUSTL(line)
      CALL trro_logic((i-2)*n_molecules+j+3*(n_molecules-1))%throw_error("Default value: "//line(1:LEN_TRIM(line))//"°!")
    ENDDO
  ENDDO
  WRITE(line,'(f4.1)') lamda
  line=ADJUSTL(line)
  CALL lamda_logic%throw_error("Default value: "//line(1:LEN_TRIM(line))//" eV!")
  ! ERRORS
  CALL lib_struct()
  CALL aostruct()
  DO i=1,n_molecules
    CALL geom_logic(i)%throw_error()
  ENDDO
  CALL basis_logic%throw_error()
  DO i=1,SIZE(comment)
    PRINT'(x,2a)', "User comment:", comment(i)
  ENDDO
  
  DO i=1,n_molecules
    CALL molecule(i)%show_geom()
  ENDDO

  ! CALCULATION BLOCK
  CALL CPU_TIME(time_start)
  DO i=1,n_molecules
    CALL molecule_param(i)
    CALL nucrep
    CALL alloc_mem
    CALL scf(i)
    CALL dealloc_mem(i)
  ENDDO
  CALL CPU_TIME(time_end)
  PRINT '(36x,a,f14.2,a)', "*CPU time elapsed (total): ", time_end-time_start, " s"
  CALL exact_elements()

  GOTO 9

! ERROR BLOCK
8   PRINT*, 'Missing input file specification. Usage: SPARC ', '"','input name','"'
    sys_message=get_command_as_string('ls *.inp')
    PRINT *, 'Input (.inp) files in the directory: '//achar(13)//achar(10),sys_message
    CALL EXIT(1)
80  CALL inp_logic%throw_error() 
81  geom_logic(i)%is_present=.FALSE.
    geom_logic(i)%location=line(1:len_trim(line))
    CALL geom_logic(i)%throw_error("File not found!")
82  geom_logic(i)%is_present=.FALSE.
    geom_logic(i)%location=line(1:len_trim(line))
    CALL geom_logic(i)%throw_error("Is the xyz format correct?")
83  energy_logic((i-1)*4+j)%is_present=.FALSE.
    energy_logic((i-1)*4+j)%is_critical=.TRUE.
    energy_logic((i-1)*4+j)%location=line(start_position:end_position)
    CALL energy_logic((i-1)*4+j)%throw_error("Is the energy format correct?")
84  flags_logic%is_present=.FALSE.
    flags_logic%is_critical=.TRUE.
    flags_logic%location=line(1:end_position)
    CALL flags_logic%throw_error("Not a known calculation flag!")
85  trro_logic(k)%is_present=.FALSE.
    trro_logic(k)%is_critical=.TRUE.
    trro_logic(k)%location=line(start_position:end_position)
    CALL trro_logic(k)%throw_error("Is the format correct?")
86  basis_logic%is_present=.FALSE.
    basis_logic%name="Basis set library location"
    basis_logic%location="'$SPARC_LIB'"
    CALL basis_logic%throw_error("Please, EXPORT it to the location of the SPARC installation.")
87  basis_logic%is_present=.FALSE.
    basis_logic%name=basis_name
    basis_logic%location=full_path
    CALL basis_logic%throw_error("You may add custom basis sets in Psi3 (SPDF uncontracted) format there."//achar(13)//achar(10)//&
    &" List of currently available ones:"//achar(13)//achar(10)//&
    &get_command_as_string('ls "$SPARC_LIB/basis_sets/"'))
88  orbrot_logic(i)%is_present=.FALSE.
    orbrot_logic(i)%is_critical=.TRUE.
    orbrot_logic(i)%location=line(1:len_trim(line))
    CALL orbrot_logic(i)%throw_error("Is the format (integer) correct?")
89  lamda_logic%is_present=.FALSE.
    lamda_logic%is_critical=.TRUE.
    lamda_logic%location=line(1:end_position)
    CALL lamda_logic%throw_error("Is the format (real) correct?")


810 PRINT '(x,a,i4,3a)', "Error: Line number", line_num, " (",line(1:LEN_TRIM(line)),") is missing the ending ';'. Add it after the data input."
    CALL EXIT(1)


9 CONTINUE  
End Program InpRead

Subroutine tolower(string)
  implicit none

! ========================================================================
! Convert upercase characters in line into lowercase

  character(*),intent(inout)      :: string
  integer                         :: jgap
  integer                         :: i

  jgap=ichar('a')-ichar('A')

  do i=1,len_trim(string)
    if (string(i:i) <= 'Z') then
      if (string(i:i) >= 'A') string(i:i)=char(ichar(string(i:i))+jgap)
    endif
  enddo

  return

End Subroutine tolower

