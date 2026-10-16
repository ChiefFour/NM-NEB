program main
        implicit none
        integer :: Natom,DOF
        integer :: i,j,k,find_result
        real*8,allocatable,dimension(:) :: coord_r,coord_p
        real*8,allocatable,dimension(:,:) :: mode_r,mode_p
        real*8,allocatable,dimension(:) :: omega_r,omega_p
        real*8,allocatable,dimension(:,:) :: D,F
        integer :: row
        integer,allocatable,dimension(:) :: selected,best_selection
        real*8,allocatable,dimension(:) :: row_sums
        real*8 :: current_sum,max_sum

        open(12,file="reactant.log")
        call located(12,"NAtoms",find_result)
        read(12,"(14x,i2)")Natom
        print*,Natom
        close(12)

        DOF=3*Natom-5

        allocate(coord_r(3*Natom),coord_p(3*Natom))
        allocate(mode_r(DOF,3*Natom),mode_p(DOF,3*Natom))
        allocate(omega_r(DOF),omega_p(DOF)) 
        allocate(D(DOF,DOF),F(DOF,DOF))
        allocate(selected(DOF),best_selection(DOF))
        allocate(row_sums(DOF))
        
        coord_r=0.0d0
        coord_p=0.0d0
        mode_r=0.0d0
        mode_p=0.0d0
        omega_r=0.0d0
        omega_p=0.0d0
 
        call ReadInput(Natom,DOF,coord_r,coord_p,mode_r,mode_p,omega_r,omega_p)
        print*,"Input read finsh !"

 
        call Displacement(Natom,DOF,coord_r,coord_p,mode_r,mode_p,D) 
        print*,"Computation of Displacement finsh !"
        
        call FC_factor(DOF,Omega_R,Omega_P,D,F)
        print*,"Frank-Condon factor finsh !"


        selected = 0
        best_selection = 0
        max_sum = 0.0

        call Back_tracking(DOF,F,row,current_sum,selected,&
                           best_selection,max_sum,row_sums)

        print*,"The max value of the overlap:", max_sum
        print*,"The combination of the largest overlap:"
        
        do i=1,DOF
           print*,"row:",i,"column:",best_selection(i)
        end do
        print*,"The value of every rows:",row_sums


        print*,"================================="


        do i=1,DOF
           print*,F(i,:)
        end do



        deallocate(coord_r,coord_p)
        deallocate(mode_r,mode_p)
        deallocate(omega_r,omega_p)
        deallocate(D,F)
        deallocate(selected,best_selection)
        deallocate(row_sums)
end 

