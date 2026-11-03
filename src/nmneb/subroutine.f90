subroutine single_point(nvar,energy,gradient,nstate)
        use molegeom_module  , only: molegeom

    implicit none
        integer, parameter :: rk = kind(1.d0) 
        integer :: i,Natom,j,dim,nstate,nvar
        real(rk) :: energy
        real(rk) :: gradient(nvar)
        real(rk),dimension(:),allocatable :: Coords
        real(rk),dimension(:,:),allocatable :: Coordinates
        character(len=2),dimension(:),allocatable :: Symbol
        character(len=79) :: Keyword
        character(len=512) :: gauss_exe, formchk_exe, cmdline
        real(rk),parameter :: AtoB=5.2917720810086E-01

        if(nstate/=1)return
       
        Keyword=molegeom%Keyword

        Natom=molegeom%nat
        dim=3*Natom
        
        allocate(Symbol(Natom))
        allocate(Coordinates(3,Natom),Coords(dim))!,molegeom%ncoords(dim))
        

       !open(28,file="Coordnates")
       !do i=1,Natom
       !    read(28,*)Coordinates(:,i)
       !end do
       !close(28)

        Coords=molegeom%ncoords*AtoB
       
      !  write(*,"(3f16.9)")molegeom%ncoords
 
        do i=1,Natom
           do j=1,3
                Coordinates(j,i)=Coords(3*(i-1)+j)
           end do
        end do
        write(6,*)"subroutine single_point get coords"
        write(6,"(3f16.9)")Coordinates       
 
        Symbol=molegeom%sym
      
        Keyword=adjustl(Keyword)
 
        open(62,file="single_point.gjf")
        write(62,"('%nproc=16',/,'%mem=8GB',/,'%chk=single_point.chk')")
        write(62,*)Keyword
        write(62,"(/,'singlet of inout_1',/,/,'0 2')")
        do i=1,Natom
           write(62,"(a3,5x,3f15.9)")Symbol(i),Coordinates(:,i)
        end do
        write(62,*)
        write(62,*)
        close(62)
        
        ! Locate the Gaussian executables: environment variables GAUSS_EXE and
        ! FORMCHK_EXE take precedence; otherwise "g16"/"formchk" are used from PATH.
        call get_environment_variable("GAUSS_EXE", gauss_exe)
        if (len_trim(gauss_exe) == 0) gauss_exe = "g16"
        call get_environment_variable("FORMCHK_EXE", formchk_exe)
        if (len_trim(formchk_exe) == 0) formchk_exe = "formchk"
        write(6,"(a,a)") " Gaussian executable : ", trim(gauss_exe)
        write(6,"(a,a)") " formchk executable : ", trim(formchk_exe)

        ! Quote the executable if its path contains spaces (e.g. Windows "Program Files")
        cmdline = trim(gauss_exe)
        if (index(cmdline," ") > 0) cmdline = '"' // trim(cmdline) // '"'
        call system(trim(cmdline) // " single_point.gjf")
        cmdline = trim(formchk_exe)
        if (index(cmdline," ") > 0) cmdline = '"' // trim(cmdline) // '"'
        call system(trim(cmdline) // " single_point.chk")


        open(72,file="single_point.fchk")
        call located(72,"Cartesian Gradient")
        read(72,*)
        read(72,*) (gradient(i), i=1,dim)
        close(72)

        open(82,file="single_point.log")
        call located(82,"SCF Done")
        read(82,"(25x,f14.9)")energy
        close(82)


        open(77,file="storedEG")
        write(77,"('Energy=',f16.9)")energy
        write(77,*)'Gradient='
        write(77,"(3f16.9)")gradient
        close(77)


        deallocate(Symbol)
        deallocate(Coordinates,Coords)

        return
end subroutine

subroutine located(filed,label)
        character(len=*) :: label
        character(len=79) :: para
        integer :: find_result,error,filed
        rewind (filed)
          do while (.true.)
            read(filed, "(A79)",iostat=error)para
            if (error/=0) stop "routine located: label not found!"
            find_result=index( para,label)
            if ( find_result  > 0 ) exit
          end do
        backspace(filed)
        return
end subroutine

integer Function GetAtomicNumber(symbol)
  implicit none
  character(2),intent(in) :: symbol
  integer :: i
  character(2), parameter :: elements(111) = &
       (/ 'H ','He', &
          'Li','Be','B ','C ','N ','O ','F ','Ne', &
          'Na','Mg','Al','Si','P ','S ','Cl','Ar', &
          'K ','Ca','Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu', &
          'Zn','Ga','Ge','As','Se','Br','Kr', &
          'Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag', &
          'Cd','In','Sn','Sb','Te','I ','Xe', &
          'Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy', &
          'Ho','Er','Tm','Yb','Lu','Hf','Ta','W ','Re','Os','Ir','Pt', &
          'Au','Hg','Tl','Pb','Bi','Po','At','Rn', &
          'Fr','Ra','Ac','Th','Pa','U ','Np','Pu','Am','Cm','Bk','Cf', &
          'Es','Fm','Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds', &
          'Rg' /)
! **********************************************************************
  do i=1,111
       if (symbol==elements(i)) then
             GetAtomicNumber=i
       end if
  end do
end Function GetAtomicNumber



real*8 Function GetAtomicMass(symbol)
  implicit none
  integer, parameter :: rk = kind(1.d0)
  character(2),intent(in) :: symbol
  integer :: i
  character(2), parameter :: elements(30) = &
       (/ 'H ','He', &
          'Li','Be','B ','C ','N ','O ','F ','Ne', &
          'Na','Mg','Al','Si','P ','S ','Cl','Ar', &
          'K ','Ca','Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu', &
          'Zn'/)
  real(rk),parameter :: atomic_masses(30) = &
      (/1.00794, 4.002602, &
        6.941, 9.012182, 10.811, 12.0107, 14.0067, 15.9994, 18.9984032, 20.1797,&
        22.98976928, 24.305, 26.9815386, 28.0855, 30.973762, 32.065, 35.453, 39.948,&
        39.0983, 40.078,44.955912, 47.867, 50.9415, 51.9961, 54.938045, &
        55.845, 58.933195, 58.6934, 63.546, 65.38/)
! **********************************************************************
  do i=1,30
       if (symbol==elements(i)) then
            GetAtomicMass=atomic_masses(i)
       end if
  end do
end Function GetAtomicMass




                                                                                                                                                          
