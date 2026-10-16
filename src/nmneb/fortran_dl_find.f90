subroutine fortranb_dl_find
        use molegeom_module  , only: molegeom

        implicit none
        integer, parameter :: rk = kind(1.d0)
        
        integer :: i,j,k,natom,iopt,icoord,nframe,nimage,maxcycle,&
                  &ncons,nconn,nspec,master,nvar1,nvar2,&
                  &nweight,nmass,n_po_scaling,printl,printf
        integer :: nat,n
        integer,dimension(:),allocatable :: nz,spec
                                !(nuclear charge)
        real(rk),dimension(:,:),allocatable :: Coord_1
        real(rk),dimension(:,:,:),allocatable :: Coord_2
        real(rk),dimension(:),allocatable :: coords_1,coords_2,mass
        real(rk),parameter :: AtoB=5.2917720810086E-01
        
        natom=molegeom%nat
        ncons=molegeom%ncons
        nconn=molegeom%nconn
        nspec=natom + natom + 5*ncons + 2*nconn + nat
        ! spec(nspec) : integer specification array
  !    nat     entries of freezing/residue number
  !    nz      entries of nuclear charges (same order as coords)
  !    5*ncons entries of constraints (typ, atom1,atom2, atom3, atom4)
  !    2*nconn entries of connections (atom1 atom2)
  !    nat     entries of microiterative region specification
  !    i.e. nspec= nat + nz + 5*ncons + 2*nconn + nat
        
        master=1
        nframe=molegeom%nframe
        nweight=0
        nmass=molegeom%nmass
        n_po_scaling=0
 
        nvar1=3*natom
        nvar2 = nframe*nvar1 + nweight + nmass + n_po_scaling
! tmpcoords2(nvarin2) : a multi-purpose real array
  !    nframe*nat*3    entries of coordinates of nframe structures
  !    nweight         entries of weights (nat or 0)
  !    nmass           entries of atomic masses (nat or 0)
  !    n_po_scaling    entries of radii scaling factors in the parallel 
  !                    optimization (0 [meaning all radii set to the base 
  !                    value], or a pre-known nivar)
  !    i.e. nvarin2= nframe*nat*3 + nweight + nmass + n_po_scaling
       
        
        allocate(coords_2(nvar2),Coord_2(3,nat,nframe),mass(natom))

        Coord_2=molegeom%Coords_2
      
        mass=molegeom%mass

        coords_2=0.0
 
        n=0
        do k =1,nframe
            do j = 1,natom
                do i=1,3
                     n=n+1
                     coords_2(n)=Coord_2(i,j,k)
                  end do
             end do
        end do


        if (nmass==0)then
            molegeom%coords2=coords_2
        else
           coords_2(natom*nframe*3+nweight+1:natom*nframe*3+nweight+nmass)=mass(1:natom)
           molegeom%coords2=coords_2
        end if



        molegeom%nvar1=nvar1  
        molegeom%nvar2=nvar2
        molegeom%nat=natom 
        molegeom%ncons=ncons
        molegeom%nconn=nconn 
        
        nspec =molegeom%nat + natom + 5*molegeom%ncons +2*molegeom%nconn+molegeom%nat
        
        allocate(spec(nspec))

        do i=1,nspec
                spec(i)=0
        end do

        nz=molegeom%nz
        spec(natom+1:2*natom)=nz(:)
        molegeom%spec=spec
        nvar1=3*natom

        print*,"Routine Main is finshied!"
        
        open(6,file="output")
        write(6,*)        " Check parameter",natom
        write(6,"(a,i8)") " nvar=" ,nvar1
        write(6,"(a,i8)") " nvar2=",nvar2
        write(6,"(a,i8)") " nspec=", nspec
        write(6,"(a,i8)") " master=",master
        write(6,"(a,i8)") " nmass=",nmass
        write(6,*)"mass="
        write(6,"(7x,f18.14)")mass

        call dl_find(nvar1, nvar2, nspec, master)

        write(6,*)'DL-FIND COMPLETED'
        close(6)

        deallocate(coords_2,Coord_2,mass) 
        deallocate(spec)

        return
end subroutine fortranb_dl_find

