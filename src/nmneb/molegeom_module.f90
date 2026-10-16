module molegeom_module
   !     use dlf_parameter_module, only: rk
        implicit none
        
        integer, parameter :: rk = kind(1.d0)   
      type molegeom_type 
        integer :: nat
        integer :: nvar1
        integer :: nvar2
        integer :: iopt,icoord,iline,nframe,nimage,maxcycle
        integer :: ncons,nconn,nspec,master
        integer :: nweight,nmass,massweight,n_po_scaling
        integer :: inithessian
        integer :: printl,printf        


        real(rk) :: tolerance_e
        real(rk) :: tolerance    
        real(rk) :: maxstep
        real(rk) :: scalestep
        real(rk) :: update
        real(rk) :: delta
        real(rk) :: nebk
        real(rk) :: energy

        character(len=79) :: Keyword !Keword of Gaussian

!array

        integer,allocatable :: nz(:)

        real(rk),allocatable :: Coords_1(:,:)
        real(rk),allocatable :: Coords1(:) !(3*natom)  the coord of gerometry1(initio)
       ! real(rk),allocatable :: coord2(:,:,:) !(3,natom,nfram)     
        real(rk),allocatable :: Coords_2(:,:,:)
        real(rk),allocatable :: Coords2(:)
        real(rk),allocatable :: ncoords(:) !(3*natom) the coord of iimage
        real(rk),allocatable :: spec(:)
        real(rk),allocatable :: gradient(:)
        real(rk),allocatable :: mass(:)

        character(len=2),allocatable :: sym(:)
        

      end type molegeom_type
      type(molegeom_type) :: molegeom
end module
