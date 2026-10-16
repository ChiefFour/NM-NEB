subroutine get_input
      ! use dlf_parameter_module, only: rk
        use molegeom_module  , only: molegeom

        implicit none
        integer, parameter :: rk = kind(1.d0)
        
        character(len=79) :: Keyword
        character(len=2),dimension(:),allocatable :: sym
        character(len=2) :: symbol
        integer :: i,j,k,natom,iopt,icoord,nframe,nimage,maxcycle,&
                  &ncons,nconn,nspec,master,nvar1,nvar2,&
                  &nweight,nmass,n_po_scaling,printl,printf,n
        integer,external :: GetAtomicNumber
        integer,dimension(:),allocatable :: nz,spec
                                !(nuclear charge)
        real(rk)  :: tolerance,tolerance_e,maxstep,scalestep,nebk
        real(rk),dimension(:,:),allocatable :: Coord_1
        real(rk),dimension(:,:,:),allocatable :: Coord_2
        real(rk),dimension(:),allocatable :: coords_1,coords_2,mass
        real(rk),parameter :: AtoB=5.2917720810086E-01
        real(rk),external :: GetAtomicMass      
 
        open(100,file="input")

        call located(100,"Atom Number")
        read(100,*)
        read(100,*)natom
        molegeom%nat=natom

        call located(100,"nframe")
        read(100,*)
        read(100,*)nframe
        molegeom%nframe=nframe       
 
        allocate(nz(natom),sym(natom),mass(natom))
        allocate(Coord_1(3,natom),Coord_2(3,natom,nframe))
        allocate(coords_1(3*natom))
        
        !Coord_1,dimension:2  coords_1,dimension:1       
 
        call located(100,"printl")
        read(100,*)
        read(100,*)printl
        molegeom%printl=printl 

        call located(100,"printf")
        read(100,*)
        read(100,*)printf
        molegeom%printf=printf

        call located(100,"Keyword")
        read(100,*)
        read(100,"(A)")Keyword
        molegeom%Keyword=Keyword

        call located(100,"iop")
        read(100,*)
        read(100,*)iopt
        molegeom%iopt=iopt
 
        call located(100,"icoord")
        read(100,*)
        read(100,*)icoord
        molegeom%icoord=icoord

        call located(100,"Nimage")
        read(100,*)
        read(100,*)nimage
        molegeom%nimage=nimage       
 
        call located(100,"maxcycle")
        read(100,*)
        read(100,*)maxcycle
        molegeom%maxcycle=maxcycle
       

        call located(100,"maxstep")
        read(100,*)
        read(100,*)maxstep
        molegeom%maxstep=maxstep

        call located(100,"scalestep")
        read(100,*)
        read(100,*)scalestep
        molegeom%scalestep=scalestep
              
        call located(100,"tolerance")
        read(100,*)
        read(100,*)tolerance
        molegeom%tolerance=tolerance
       ! tolerance_e=tolerance/450.D0

        call located(100,"tolerance_e")
        read(100,*)
        read(100,*)tolerance_e
        molegeom%tolerance_e=tolerance_e

        call located(100,"nebk")
        read(100,*)
        read(100,*)nebk
        molegeom%nebk=nebk


        call located(100,"geometry 1")
        read(100,*)
        do i = 1,natom
                read(100,*)sym(i),Coord_1(:,i)
        end do
        
        molegeom%sym=sym

        do i =1,natom
           do j=1,3
                coords_1((i-1)*3+j)=Coord_1(j,i)
           end do
        end do

        molegeom%coords1=coords_1/AtoB


     if (nframe==1)then
        call located(100,"geometry 2")
        read(100,*)
        do k=1,nframe
                do j = 1,natom
                        read(100,*)sym(i),Coord_2(:,j,k)
                end do
        end do
     else if(nframe>1) then
        call located(100,"geometry 2")
        read(100,*)
        read(100,*)
        do k=1,nframe
                do j = 1,natom
                        read(100,*)sym(i),Coord_2(:,j,k)
                end do
                read(100,*)
        end do
     else
        write(*,*)"The value of nframe is error"
        stop
     end if

        molegeom%Coords_2=Coord_2/AtoB



        do i =1,natom
           symbol=sym(i)
           nz(i)=GetAtomicNumber(symbol)
        end do
        molegeom%nz=nz



        do i =1,natom
           symbol=sym(i)
           mass(i)= GetAtomicMass(symbol)
                write(*,*)GetAtomicMass(symbol)
        end do
        !molegeom%mass=mass
        molegeom%mass=mass! *(1.660538782E-27)/9.109E-31
        write(*,"(f18.14)") molegeom%mass

        close(100)

        deallocate(nz,sym,mass,Coord_1,Coord_2,coords_1)

        print*,"Getting information from input is finshed!"

        return
end subroutine get_input

