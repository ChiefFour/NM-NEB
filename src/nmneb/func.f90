!> Oct. 12, 2014 -bsuo- Interface between BDF and DL-FIND
!                       Modified from dl-find crystal interface
subroutine dlf_get_params(nvar,nvar2,nspec,coords,coords2,spec,ierr,&
   &tolerance,printl,maxcycle,maxene,tatoms,icoord,&
   &iopt,iline,maxstep,scalestep,lbfgs_mem,nimage,nebk,dump,restart,&
   &nz,ncons,nconn,update,maxupd,delta,soft,inithessian,carthessian,tsrel,&
   &maxrot,tolrot,nframe,nmass,nweight,timestep,fric0,fricfac,fricp,&
   &imultistate, state_i,state_j,pf_c1,pf_c2,gp_c3,gp_c4,ln_t1,ln_t2,&
   &printf,tolerance_e,distort,massweight,minstep,maxdump,task,temperature,&
   &po_pop_size,po_radius,po_contraction,po_tolerance_r,po_tolerance_g,&
   &po_distribution,po_maxcycle,po_init_pop_size,po_reset,po_mutation_rate,&
   &po_death_rate,po_scalefac,po_nsave,ntasks,tdlf_farm,n_po_scaling,&
   &neb_climb_test,neb_freeze_test,nzero, coupled_states, qtsflag,&
!   they have been deleted in the updated dlf_global_module
!   tolerance_s,tolerance_rmsg,tolerance_rmss )
    imicroiter, maxmicrocycle, micro_esp_fit)

!  use dlf_parameter_module, only: rk
  use molegeom_module  , only: molegeom
  implicit none
  integer, parameter :: rk = kind(1.d0)
  integer   ,intent(in)      :: nvar 
  integer   ,intent(in)      :: nvar2
  integer   ,intent(in)      :: nspec
  real(rk), intent(inout)    :: coords(nvar)  ! start coordinates
  real(rk), intent(inout)    :: coords2(max(nvar2,1)) ! a real array that can be 
                                              ! used depending on the
                                              ! calculation
                                              ! e.g. a second set of
                                              ! coordinates
  integer   ,intent(inout)   :: spec(nspec)   ! specifications like fragment or frozen
  integer   ,intent(out)     :: ierr
  real(rk)  ,intent(inout)   :: tolerance
  real(rk)  ,intent(inout)   :: tolerance_e
  integer   ,intent(inout)   :: printl
  integer   ,intent(inout)   :: maxcycle
  integer   ,intent(inout)   :: maxene
  integer   ,intent(inout)   :: tatoms
  integer   ,intent(inout)   :: icoord
  integer   ,intent(inout)   :: iopt
  integer   ,intent(inout)   :: iline
  real(rk)  ,intent(inout)   :: maxstep
  real(rk)  ,intent(inout)   :: scalestep
  integer   ,intent(inout)   :: lbfgs_mem
  integer   ,intent(inout)   :: nimage
  real(rk)  ,intent(inout)   :: nebk
  integer   ,intent(inout)   :: dump
  integer   ,intent(inout)   :: restart
  integer   ,intent(inout)   :: nz
  integer   ,intent(inout)   :: ncons
  integer   ,intent(inout)   :: nconn
  integer   ,intent(inout)   :: update
  integer   ,intent(inout)   :: maxupd
  real(rk)  ,intent(inout)   :: delta
  real(rk)  ,intent(inout)   :: soft
  integer   ,intent(inout)   :: inithessian
  integer   ,intent(inout)   :: carthessian
  integer   ,intent(inout)   :: tsrel
  integer   ,intent(inout)   :: maxrot
  real(rk)  ,intent(inout)   :: tolrot
  integer   ,intent(inout)   :: nframe
  integer   ,intent(inout)   :: nmass
  integer   ,intent(inout)   :: nweight
  real(rk)  ,intent(inout)   :: timestep
  real(rk)  ,intent(inout)   :: fric0
  real(rk)  ,intent(inout)   :: fricfac
  real(rk)  ,intent(inout)   :: fricp
  integer   ,intent(inout)   :: imultistate
  integer   ,intent(inout)   :: state_i
  integer   ,intent(inout)   :: state_j
  real(rk)  ,intent(inout)   :: pf_c1  
  real(rk)  ,intent(inout)   :: pf_c2  
  real(rk)  ,intent(inout)   :: gp_c3  
  real(rk)  ,intent(inout)   :: gp_c4
  real(rk)  ,intent(inout)   :: ln_t1  
  real(rk)  ,intent(inout)   :: ln_t2  
  integer   ,intent(inout)   :: printf
  real(rk)  ,intent(inout)   :: distort
  integer   ,intent(inout)   :: massweight
  real(rk)  ,intent(inout)   :: minstep
  integer   ,intent(inout)   :: maxdump
  integer   ,intent(inout)   :: task
  real(rk)  ,intent(inout)   :: temperature
  integer   ,intent(inout)   :: po_pop_size
  real(rk)  ,intent(inout)   :: po_radius
  real(rk)  ,intent(inout)   :: po_contraction
  real(rk)  ,intent(inout)   :: po_tolerance_r
  real(rk)  ,intent(inout)   :: po_tolerance_g
  integer   ,intent(inout)   :: po_distribution
  integer   ,intent(inout)   :: po_maxcycle
  integer   ,intent(inout)   :: po_init_pop_size
  integer   ,intent(inout)   :: po_reset
  real(rk)  ,intent(inout)   :: po_mutation_rate
  real(rk)  ,intent(inout)   :: po_death_rate
  real(rk)  ,intent(inout)   :: po_scalefac
  integer   ,intent(inout)   :: po_nsave
  integer   ,intent(inout)   :: ntasks
  integer   ,intent(inout)   :: tdlf_farm
  integer   ,intent(inout)   :: n_po_scaling
  real(rk)  ,intent(inout)   :: neb_climb_test
  real(rk)  ,intent(inout)   :: neb_freeze_test
  integer   ,intent(inout)   :: nzero
  integer   ,intent(inout)   :: coupled_states
  integer   ,intent(inout)   :: qtsflag
  integer   ,intent(inout)   :: imicroiter
  integer   ,intent(inout)   :: maxmicrocycle
  integer   ,intent(inout)   :: micro_esp_fit
  character*1  :: chdum
  real(rk),parameter :: bohr=5.2917720810086E-01
  print*,"subroutine dlf_get_params is active on"


! local variables
! **********************************************************************

! Strategy: BDF's INPOPT reads in the DLFIND section of the input deck
!           and stores the corresponding variables in dlfind_module,
!           with a
!           bdf_glob%.
!           When dl_find is called, dlf_default_init sets initial values
!           to 
!           junk in most cases, to recognise user input.
!           Then dlf_get_params is called to fill arguments with input
!           values.
!           Then dlf_default_set assigns default values to stuff not
!           assigned
!           via input.
!
!           So, here in dlf_get_params, dummy arguments corresponding to 
!           variables that will never be under user control can simply
!           not be 
!           assigned if the defaults are what we want.
!           Also need to give initial values to the dlf_ variables that
!           match 
!           the default_init ones, to safely cover an absent input
!           directive.
!           Inconsistent options will be dealt with in or just after the
!           input 
!           reading routine.

! A note on units: the array section containing the masses has them in
! atomic 
! mass units.  Coordinates are in bohr (atomic units) and energies are
! returned 
! from BDF or other package in hartree (atomic units), with
! corresponding gradients.
! The temperature should be input in kelvin.
! The coords in dlf_put_coords are in bohr, so that crystal can do its
! own usual 
! conversions.
! The routine write_xyz in dlf_util.f90 does a conversion to Angstrom
! before writing.
  write(6,*) " Transfer DL-FIND parameter from Input to DL-FIND kernel routine ..."
  ! Things hardwired
  ierr = 0   ! not used here for error-checking
  tatoms = 1 ! could possibly use tatoms = 0 to fudge input of coords in terms 
             ! of symmetry allowed directions.  Need to investigate
             ! further...
  tdlf_farm = 0
        
                
        
  write(6,*) "nvar,nvar2",nvar,nvar2
  nz = molegeom%nat  ! cywu20210514
! the second set of coordinates
! if(dlf_nframe > 0 .and. .not. allocated(glob%xcoords2)) then
!   if(nvar2 < 3*bdf_dlf_glob%nat*dlf_nframe) then
!     write(6,*) 'Error! nvar2 in sub. dlf_get_params is not large
!     enough!'
!     call qexit()
!   end if
!   allocate(bdf_dlf_glob%xcoords2(1,nvar2,1))
!   call
!   rw_molegeom_in_punch(2,1,bdf_dlf_glob%nat*dlf_nframe,chdum,bdf_dlf_glob%xcoords2)
! end if

  ! Things defined by reading the dl-find section of the CRYSTAL input
  ! deck
! if (nvar2 > 0 .and. allocated(bdf_dlf_glob%xcoords2)) then
!   coords2(:) = 0.0D0
!   call dcopy(nvar2,bdf_dlf_glob%xcoords2,1,coords2,1)
! endif
  
  nweight = molegeom%nweight
  nmass   = molegeom%nmass
  massweight=molegeom%massweight
  nframe  = molegeom%nframe
  maxcycle= molegeom%maxcycle
  !maxene  = bdf_dlf_glob%maxene
  spec(:) = molegeom%spec(:) ! cywu20210514

  printl  = molegeom%printl
  printf  = molegeom%printf

  tolerance_e      = molegeom%tolerance_e
  tolerance        = molegeom%tolerance

  write(6,"(a,f12.8)") " dl-find tolerance_e    : ",tolerance_e
  write(6,"(a,f12.8)") " dl-find tolerance      : ",tolerance

  iopt =molegeom%iopt
  iline=molegeom%iline
  maxstep=molegeom%maxstep
  scalestep=molegeom%scalestep
!  lbfgs_mem= dlf_lbfgs_mem
  update= molegeom%update
!  maxupd= dlf_maxupd
!  delta= dlf_delta
  delta= molegeom%delta
  inithessian = molegeom%inithessian 
  nimage= molegeom%nimage
  icoord= molegeom%icoord
  
  coords=molegeom%coords1
  !write(*,*)"----------------------------"
  !write(6,*)"Inoput Coords1:"
  !write(6,"(3f16.6)")coords*bohr
  !write(*,*)"----------------------------"

  coords2=molegeom%coords2
  !write(*,*)"----------------------------"
  !write(6,*)"Inoput Coords2:"
  !write(6,"(3f16.6)")coords2*bohr
  !write(*,*)"----------------------------"
! iline must be 0 for NEB and Dimer
  if( icoord > 99) iline = 0
  nebk= molegeom%nebk

! MULTI-STATE
!  imultistate = bdf_dlf_glob%imultistate
!  state_i = bdf_dlf_glob%state_i
!  state_j = bdf_dlf_glob%state_j
!  pf_c1   = bdf_dlf_glob%pf_c1
!  pf_c2   = bdf_dlf_glob%pf_c2
!  gp_c3   = bdf_dlf_glob%gp_c3
!  gp_c4   = bdf_dlf_glob%gp_c4
!  ln_t1   = bdf_dlf_glob%ln_t1
!  ln_t2   = bdf_dlf_glob%ln_t2
!  if(bdf_dlf_glob%needcoupling.eq.2) coupled_states=0

!  ncons= bdf_dlf_glob%ncons ! cywu20210515
!  nconn= bdf_dlf_glob%nconn ! cywu20210515
!  dump= dlf_dump
!  restart= dlf_restart
!
!  timestep= dlf_timestep
!  fric0= dlf_fric0
!  fricfac= dlf_fricfac
!  fricp= dlf_fricp
!
!  task = dlf_task
!
!  temperature = dlf_temperature
!
!  po_pop_size = dlf_po_pop_size
!  po_radius = dlf_po_radius
!  po_contraction = dlf_po_contraction
!  po_tolerance_r = dlf_po_tolerance_r
!  po_tolerance_g = dlf_po_tolerance_g
!  po_distribution = dlf_po_distribution
!  po_maxcycle = dlf_po_maxcycle
!  po_init_pop_size = dlf_po_init_pop_size
!  po_reset = dlf_po_reset
!  po_mutation_rate = dlf_po_mutation_rate
!  po_death_rate = dlf_po_death_rate
!  po_scalefac = dlf_po_scalefac
!  po_nsave = dlf_po_nsave
!  n_po_scaling = dlf_n_po_scaling

  ! serial run
  ntasks = 1 !dlf_ntasks

  tdlf_farm = 0
  nzero=0 ! new option, not yet implemented in crystal
print*,"Subroutine dlf_get_params is finished!"
return
end subroutine dlf_get_params


subroutine dlf_get_gradient(nvar,coords,energy,gradient,iimage,status)
! get the gradient at a given geometry
!  use dlf_parameter_module, only: rk
!  use dlf_stat,only:stat   
! use dlfind_module, only : dlf_inf,dlf_forze,dlf_counter,bdf_dlf_glob
! use molegeom,      only : moleobj
! use units, only : bohr
! use control, only : toptask ! wzk20190704
! use bdfopt_mod, only : optinfo
  use molegeom_module, only: molegeom
  implicit none
  integer, parameter :: rk = kind(1.d0)
  integer   ,intent(in)    :: nvar
  real(rk)  ,intent(inout) :: coords(nvar)   ! also needs to be feedback
  real(rk)  ,intent(out)   :: energy
  real(rk)  ,intent(out)   :: gradient(nvar)
  integer   ,intent(in)    :: iimage
  integer   ,intent(out)   :: status
  integer   ,save          :: call_counter
  integer                  :: initio_counter
  !--- local ---
  integer,parameter :: nstate=1
  logical                  :: err1, err2
  integer :: fn,n,i ! scratch file
  real(rk),parameter :: bohr=5.2917720810086E-01
  character*128 :: bdftask
        
  
  if (initio_counter==0)then
     call_counter=0
     initio_counter=1
  end if

  call_counter=call_counter+1
        
  write(6,*) 
  write(6,*) "[dl_find_get_grad]"
  write(6,*) "The times of Subroutine dlf_get_gradient is called",call_counter
  write(6,*) 

  status = 0 ! status /= 0 means error
  err1 = .false.
  err2 = .false.

  ! write original coord into scratch file
  
  molegeom%ncoords=coords
 !open(18,file="Coordnates")
 !      write(18,"(3f16.6)")coords
 !close(18)
  !write(6, *)"-------------------------------------------"
  !write(6, *)"subroutine get_gradient get coords"
  !write(6,"(3f16.6)")coords*bohr
  !write(6, *)"-------------------------------------------"
  call single_point(nvar,energy,gradient,nstate)
  write(6,"(a)") " Gradient from gauss="
  n=1
  do i=1,molegeom%nat
    write(6,"(a,3(2x,f14.8))") molegeom%sym(i),gradient(n:n+2)
    n=n+3
  enddo

  ! read rotated coord from scratch file
 !call
 !rw_molegeom_in_punch(2,0,moleobj%natom,moleobj%label,moleobj%coord)
 !call dcopy(nvar,moleobj%coord,1,coords,1)

 !! read energy/grad and pass them to dl-find
 !if(optinfo%multistate == 0) then
 !  call rw_egrad_in_punch(2,moleobj%natom,energy,gradient,nstate)
 !else
 !  call multi_state_mix(moleobj%natom,energy,gradient)
 !end if
 !call
 !print_molegeom_opt_grad_info(moleobj%natom,moleobj%label,coords,gradient,energy,nstate)

   ! wzk20180120 adapted from bdfopt_bdf.F90
   write(*,*) 'Projecting out translational and rotational componentsin gradient...'
   ! wzk20190704: only do projection if there are no external charges
!  if(toptask%task(52).eq.0)then
!   call bdfProjTransRot(1, moleobj%natom, moleobj%coord,&
!       reshape(gradient, (/3, moleobj%natom/))) !
!       libmath_util/bdf3DVecUtil.F90
!  endif

  !==========================================================
  ! append geom,gradient and energy into pes.txt
  ! BDFTASK
! bdftask=" "
! call envvar("BDFTASK",bdftask)
! bdftask=trim(adjustl(bdftask))//".pes1"
  

  fn=1031
  write(fn,*)
!  write(fn,"(a,i8)") " Optimization cycle ",stat%ccycle 
  write(fn,"(a,f18.8)") " Energy=",energy
  write(fn,"(a)") " Geometry="
  do i=1,molegeom%nat
   !  write(fn,"(a,3(2x,f14.8))")molegeom%sym(i),molegeom%ncoords((i-1)*3+1:(i-1)*3+3)*bohr
        write(fn,"(a,3(2x,f14.8))")molegeom%sym(i),coords((i-1)*3+1:(i-1)*3+3)*bohr
  enddo
  write(fn,"(a)") " Gradient="
  n=1
  do i=1,molegeom%nat
    write(fn,"(a,3(2x,f14.8))") molegeom%sym(i),gradient(n:n+2)
    n=n+3
  enddo
  close(fn)
  !==========================================================
return
end subroutine dlf_get_gradient


subroutine dlf_get_multistate_gradients(nvar,coords,energy,gradient,&
                                        coupling,needcoupling,iimage,status)

! ! only a dummy routine up to now
! ! for conical intersection search
! use dlf_parameter_module, only: rk
! use dlf_stat, only : stat
! use dlfind_module, only : dlf_inf,dlf_forze,dlf_counter,bdf_dlf_glob
! use molegeom,      only : moleobj
! use units, only : bohr
! implicit none
! integer   ,intent(in)    :: nvar
! integer   ,intent(inout) :: coords(nvar)
! real(rk)  ,intent(in)    :: energy(2)
! real(rk)  ,intent(out)   :: gradient(nvar,2)
! real(rk)  ,intent(out)   :: coupling(nvar)
! logical   ,intent(in)    :: needcoupling
! integer   ,intent(in)    :: iimage
! integer   ,intent(out)   :: status
! !--- local ---
! integer,parameter :: nstate=2
! integer :: fn,i,n ! scratch file
! character*128 :: bdftask
!
! write(6,*) 
! write(6,*) "[dlf_get_multistate_gradients] needcoupling=",needcoupling
! write(6,*) 

! status = 0 ! status /= 0 means error

! ! write original coord into scratch file
! call dcopy(nvar,coords,1,moleobj%coord,1)
! call
! rw_molegeom_in_punch(1,0,moleobj%natom,moleobj%label,moleobj%coord)

! call bdf_single_point(nstate,.true.,1)

! ! read rotated coord from scratch file
! call
! rw_molegeom_in_punch(2,0,moleobj%natom,moleobj%label,moleobj%coord)
! call dcopy(nvar,moleobj%coord,1,coords,1)
! 
! !>>> multi-state version <<<
! ! read energy/grad and pass them to dl-find 
! call read_multi_egrad_in_punch(moleobj%natom,energy,gradient,nstate)
! call
! print_molegeom_opt_grad_info(moleobj%natom,moleobj%label,coords,gradient,energy,nstate)

! !>>> LOAD NAC - Z.Li@2014.11.14 <<<
! if(needcoupling)then
!    write(*,*) " LOAD NAC ... "
!    call rw_nac_in_punch(2,moleobj%natom,coupling)
!    call
!    print_molegeom_opt_nac_info(moleobj%natom,moleobj%label,coupling)
! endif
! !==========================================================
! ! append geom,gradient and energy into pes.txt
! ! BDFTASK
! bdftask=" "
! call envvar("BDFTASK",bdftask)
! bdftask=trim(adjustl(bdftask))//".pes2"
! fn=1031
! dlf_counter=dlf_counter+1
! if(dlf_counter.eq.1)then
!    open(fn,file=bdftask,form="formatted")
! else
!    open(fn,file=bdftask,form="formatted",access="append") 
! endif 
! write(fn,*)
! write(fn,"(a,i8)") " Optimization cycle ",stat%ccycle 
! write(fn,"(a,2f18.8)") " Energy=",energy
! write(fn,"(a)") " Geometry="
! do i=1,moleobj%natom
!    write(fn,"(a,3(2x,f14.8))")
!    moleobj%label(i),moleobj%coord(1:3,i)*bohr
! enddo
! write(fn,"(a)") " Gradient1="
! n=1
! do i=1,moleobj%natom
!   write(fn,"(a,3(2x,f14.8))") moleobj%label(i),gradient(n:n+2,1)
!   n=n+3
! enddo
! write(fn,"(a)") " Gradient2="
! n=1
! do i=1,moleobj%natom
!   write(fn,"(a,3(2x,f14.8))") moleobj%label(i),gradient(n:n+2,2)
!   n=n+3
! enddo
! ! NAC
! if(needcoupling)then
!    write(fn,"(a)") " NAC12="
!    n=1
!    do i=1,moleobj%natom
!      write(fn,"(a,3(2x,f30.8))") moleobj%label(i),coupling(n:n+2)
!      n=n+3
!    enddo
! endif
! close(fn)
! !==========================================================

return
end subroutine dlf_get_multistate_gradients


subroutine dlf_get_hessian(nvar,coords,hessian,status)
! get the hessian at a given geometry
! use dlf_parameter_module
! implicit none
! integer   ,intent(in)    :: nvar
! real(rk)  ,intent(in)    :: coords(nvar)
! real(rk)  ,intent(out)   :: hessian(nvar,nvar)
! integer   ,intent(out)   :: status

! Dummy routine currently; numerical Hessian could be provided here if
! desired.
! call dlf_fail("dlf_get_hessian called but not coded currently")
return
end subroutine dlf_get_hessian


subroutine dlf_put_coords(nvar,mode,energy,coords,iam)
! copy coords to molebobj%coord and write into file *.optgeom 
!  use dlf_parameter_module
!  use dlf_global, only: stdout
  use molegeom_module,   only: molegeom
  implicit none
  integer, parameter :: rk = kind(1.d0)
  integer   ,intent(in)    :: nvar
  integer   ,intent(in)    :: mode
  integer   ,intent(in)    :: iam
  real(rk)  ,intent(in)    :: energy
  real(rk)  ,intent(in)    :: coords(nvar)

  integer                  :: i, j, ii, lpr
  logical                  :: od
  character(len=80)        :: fname
  character(len=60)        :: label
  real(rk),parameter :: bohr=5.2917720810086E-01

  write(6, *)"-------------------------------------------"
  write(6,*) "dl-find put coord",iam
  ! Only do this writing of files if I am the rank-zero processor
 !if (iam /= 0) return
 !  molegeom%ncoords=coords
  write(6, *)"-------------------------------------------"
  write(6,"(3f16.6)")coords*bohr
  write(6, *)"-------------------------------------------"
 !open(31,file="Coordnates.put")
     !  write(31,"(3f16.6)")coords
 !close(31)    
        
  ! Bsuo: write new coord to scratch file
 ! call
 ! rw_molegeom_in_punch(1,0,moleobj%natom,moleobj%label,moleobj%coord)
 
end subroutine dlf_put_coords


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! NOT QUITE USEFUL
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine dlf_error()
  implicit none
 ! call error_exit("dlf_error, Failure in DL-FIND",100)
end subroutine dlf_error


subroutine dlf_update()
  implicit none
  ! only a dummy routine here.
end subroutine dlf_update


!subroutine dlf_put_procinfo(dlf_nprocs, dlf_iam, dlf_global_comm)
! implicit none
! integer, intent(in) :: dlf_nprocs ! total number of processors
! integer, intent(in) :: dlf_iam ! my rank, from 0, in mpi_comm_world
! integer, intent(in) :: dlf_global_comm ! world-wide communicator
! ! not needed as CRYSTAL will sort out the MPI stuff
! !call dlf_fail("dlf_put_procinfo called in error")
!end subroutine dlf_put_procinfo


!subroutine dlf_get_procinfo(dlf_nprocs, dlf_iam, dlf_global_comm)
! implicit none
! integer :: dlf_nprocs ! total number of processors
! integer :: dlf_iam ! my rank, from 0, in mpi_comm_world
! integer :: dlf_global_comm ! world-wide communicator
!!! dummy argument = corresponding variable in the calling program
!  dlf_nprocs = global_nproc
!  dlf_iam = global_iam
!  dlf_global_comm = global_comm
! end subroutine dlf_get_procinfo


!subroutine dlf_put_taskfarm(dlf_ntasks, dlf_nprocs_per_task,
!dlf_iam_in_task, &
!                       dlf_mytask, dlf_task_comm, dlf_ax_tasks_comm)

! implicit none

! integer, intent(in) :: dlf_ntasks          ! number of taskfarms
! integer, intent(in) :: dlf_nprocs_per_task ! no of procs per farm
! integer, intent(in) :: dlf_iam_in_task     ! my rank, from 0, in my
! farm
! integer, intent(in) :: dlf_mytask          ! rank of my farm, from 0
! integer, intent(in) :: dlf_task_comm       ! communicator within each
! farm
! integer, intent(in) :: dlf_ax_tasks_comm   ! communicator involving
! the 
                                             ! i-th proc from each farm
! **********************************************************************

!!! variable in the calling program = corresponding dummy argument 

! not needed as CRYSTAL will sort out the MPI stuff
  !call dlf_fail("dlf_put_taskfarm called in error")

!end subroutine dlf_put_taskfarm


!subroutine dlf_get_taskfarm(dlf_ntasks, dlf_nprocs_per_task,
!dlf_iam_in_task, &
!                       dlf_mytask, dlf_task_comm, dlf_ax_tasks_comm)
! implicit none

! integer :: dlf_ntasks          ! number of taskfarms
! integer :: dlf_nprocs_per_task ! no of procs per farm
! integer :: dlf_iam_in_task     ! my rank, from 0, in my farm
! integer :: dlf_mytask          ! rank of my farm, from 0
! integer :: dlf_task_comm       ! communicator within each farm
! integer :: dlf_ax_tasks_comm   ! communicator involving the
!                                ! i-th proc from each farm
! **********************************************************************
!!! dummy argument = corresponding variable in the calling program
! we should only be in this routine AFTER make_task_farm in libxMPI.f 
! has been called (with its replacement of the iam, nproc for
! MPI_COMM_WORLD 
! with values for the taskfarms)
!  dlf_ntasks = ntasks
!  dlf_nprocs_per_task = nproc
!  dlf_iam_in_task = iam
!  dlf_mytask = mytask 
!  dlf_task_comm = task_comms
!  dlf_ax_tasks_comm = ax_tasks_comms
! end subroutine dlf_get_taskfarm


!subroutine dlf_output(dum_stdout, dum_stderr)
! use dlf_parameter_module, only: rk
!  use dlf_global, only: glob,stderr,stdout,keep_alloutput
! implicit none
! integer :: dum_stdout
! integer :: dum_stderr
! integer :: ierr
! sort out output units; particularly important on multiple processors
! set unit numbers for main output and error messages
! if (dum_stdout >= 0) stdout = dum_stdout 
! if (dum_stderr >= 0) stderr = dum_stderr
! if (glob%nprocs > 1) then
!    ! write some info on the parallelization
!    write(stdout,'(1x,a,i10,a)')"I have rank ",glob%iam," in
!    mpi_comm_world"
!    write(stdout,'(1x,a,i10)')"Total number of processors =
!    ",glob%nprocs
!    if (keep_alloutput) then
!       write(stdout,'(1x,a)')"Keeping output from all processors"
!    else
!       write(stdout,'(1x,a)')"Not keeping output from processors /= 0"
!    end if
! end if
!end subroutine dlf_output

