
subroutine located(filed,label,find_result)
        implicit none
        character(len=*) :: label
        character(len=79) :: para
        integer :: find_result,error,filed
        find_result=0
        rewind (filed)
          do while (.true.)
            read(filed, "(A79)",iostat=error)para
            if (error/=0) exit
            find_result=index( para,label)
            if ( find_result  > 0 ) exit
          end do
        backspace(filed)
        return
end subroutine



subroutine ReadInput(Natom,DOF,R,P,mode_r_two,mode_P_two,omega_r,omega_p)
        implicit none
        integer :: Natom,DOF,n,i,j,k,l,m,find_result,n_mod
        character(len=79) :: line,order_1,order_2,symbol(Natom),a,b
        real*8 :: n1,n2,n3,n4,n5,n6
        real*8 :: mode_read(Natom,3*DOF)
        real*8 :: mode_R(DOF,Natom,3),mode_P(DOF,Natom,3)
        real*8 :: mode_R_two(DOF,3*Natom),mode_P_two(DOF,3*Natom)
        real*8 :: Omega_R(DOF),Omega_P(DOF)
        real*8 :: Coord_R(Natom,3),Coord_P(Natom,3),Delta_Coord(Natom,3)
        real*8 :: R(3*Natom),P(3*Natom)
        real*8,dimension(:,:,:),allocatable :: mode_n
!!!========================================================!!!
!!!========================================================!!!
!!!          Only suit for Nonlinear molecules             !!!
!!!========================================================!!!
!!!========================================================!!!

        write(line,*)Natom
        write(*,*)"Natom=",line
        write(*,*)"DOF=",DOF

        order_1="grep -A "//trim(adjustl(line))//" 'Atom  AN' reactant.log > mode_reactant"
        print*,order_1
        call system(order_1)

        order_1="grep 'Frequencies' reactant.log > omega_r"
        print*,order_1
        call system(order_1)

        order_2="grep -A "//trim(adjustl(line))//" 'Atom  AN' product.log > mode_product"
        print*,order_2
        call system(order_2)

        order_2="grep 'Frequencies' product.log > omega_p"
        print*,order_2
        call system(order_2)

!!=====================================================================
!!         Read the Coordinates of Reactant and Product              !!
!!
!!=====================================================================
 !! read coordinates
        open(90,file="reactant.com")
        do i=1,7
           read(90,*)
        end do
        do i=1,natom
           read(90,*)Symbol(i),Coord_R(i,:)
        end do
        close(90)
        open(95,file="product.com")
        do i=1,7
           read(95,*)
        end do
        do i=1,natom
           read(95,*)Symbol(i),Coord_P(i,:)
        end do
        close(95)

        Delta_Coord=Coord_P-Coord_R
        !The coordinate vector difference

 ! translate two dimension to one dimension
        k=0
        do i=1,Natom
           do j=1,3
                k=k+1
                R(k)=Coord_R(i,j)
                P(k)=Coord_P(i,j)
           end do
        end do

!!================================================================
!!       Read the mode in the reactants and products            !!
!!================================================================
        print*,"Coordination finish!"


        n_mod=mod(DOF,3)
        n=(DOF-n_mod)/3

        print*,"n_mod:",n_mod
        print*,"n:",n

        allocate(mode_n(n,Natom,9))

        open(100,file="mode_reactant",position="append")
        write(100,*)"---"
        rewind(100)
        do i=1,n
        read(100,*)
             do j=1,natom
                 read(100,*)a,b,mode_n(i,j,:)
             end do
        read(100,*)
        end do

        do i=1,n
            do j=1,Natom
                mode_read(j,(i-1)*9+1:i*9)=mode_n(i,j,:)
            end do
        end do

     if (n_mod/=0)then
        read(100,*)a
       ! print*,a
        if (n_mod==1)then
            do j=1,Natom
                read(100,*)a,b,n1,n2,n3
                mode_read(j,3*DOF)=n3
                mode_read(j,3*DOF-1)=n2
                mode_read(j,3*DOF-2)=n1
             end do
        else
            do j=1,Natom
                read(100,*)a,b,n1,n2,n3,n4,n5,n6
                mode_read(j,3*DOF)=n6
                mode_read(j,3*DOF-1)=n5
                mode_read(j,3*DOF-2)=n4
                mode_read(j,3*DOF-3)=n3
                mode_read(j,3*DOF-4)=n2
                mode_read(j,3*DOF-5)=n1
             end do
        end if
      end if

        close(100)
       print*,"mode_r:1"

        do i=1,DOF
           do j=1,Natom
                mode_R(i,j,:)=mode_read(j,3*(i-1)+1:3*i)
           end do
        end do

        do i=1,DOF
             l=1
             do j=1,Natom
                  do k=1,3
                        mode_R_two(i,l)=mode_R(i,j,k)
                        l=l+1
                  end do
            end do
        end do

        print*,"mode_p"

        mode_n=0.0d0
        mode_read=0.0d0
        open(100,file="mode_product",position="append")
        write(100,*)"---"
        rewind(100)
        do i=1,n
        read(100,*)
             do j=1,natom
                 read(100,*)a,b,mode_n(i,j,:)
             end do
        read(100,*)
        end do

        do i=1,n
            do j=1,Natom
                mode_read(j,(i-1)*9+1:i*9)=mode_n(i,j,:)
            end do
        end do

     if (n_mod/=0)then
        read(100,*)
        if (n_mod==1)then
            do j=1,Natom
                read(100,*)a,b,n1,n2,n3
                mode_read(j,3*DOF)=n3
                mode_read(j,3*DOF-1)=n2
                mode_read(j,3*DOF-2)=n1
             end do
        else
            do j=1,Natom
                read(100,*)a,b,n1,n2,n3,n4,n5,n6
                mode_read(j,3*DOF)=n6
                mode_read(j,3*DOF-1)=n5
                mode_read(j,3*DOF-2)=n4
                mode_read(j,3*DOF-3)=n3
                mode_read(j,3*DOF-4)=n2
                mode_read(j,3*DOF-5)=n1
             end do
        end if
      end if

        close(100)

        do i=1,DOF
           do j=1,Natom
                mode_P(i,j,:)=mode_read(j,3*(i-1)+1:3*i)
           end do
        end do

        do i=1,DOF
             l=1
             do j=1,Natom
                  do k=1,3
                        mode_P_two(i,l)=mode_P(i,j,k)
                        l=l+1
                  end do
            end do
        end do

        deallocate(mode_n)

        print*,"Mode finish!"
 !!================================================!!!
 !!               Read the Frequences              !!!
 !!================================================!!!
        Omega_R=0.0d0
        Omega_P=0.0d0
        open(110,file="omega_r")
        do i=1,n
            read(110,*)a,b,n1,n2,n3
            Omega_R((i-1)*3+1)=n1
            Omega_R((i-1)*3+2)=n2
            Omega_R((i-1)*3+3)=n3
        end do
        if (n_mod==1)then
           read(110,*)a,b,Omega_R(DOF)
        else if(n_mod==2)then
           read(110,*)a,b,Omega_R(DOF-1),Omega_R(DOF)
        end if
        close(110)

        open(120,file="omega_p")
        do i=1,n
            read(120,*)a,b,n1,n2,n3
            Omega_P((i-1)*3+1)=n1
            Omega_P((i-1)*3+2)=n2
            Omega_P((i-1)*3+3)=n3
        end do
        if (n_mod==1)then
           read(120,*)a,b,Omega_P(DOF)
        else if(n_mod==2)then
           read(120,*)a,b,Omega_P(DOF-1),Omega_P(DOF)
        end if
        close(120)

        print*,"Omega Reactant:",Omega_R
        print*,"Omega Product:",Omega_P
end subroutine ReadInput




subroutine Displacement(Natom,DOF,coord_r,coord_p,mode_r,mode_p,D)
        implicit none
        integer :: Natom,DOF,i,j
        real*8 :: coord_r(3*Natom),coord_p(3*Natom),delta_coord(3*Natom)
        real*8 :: proj_r(3*Natom)
        real*8 :: mode_r(DOF,3*Natom),mode_p(DOF,3*Natom)
        real*8 :: D(DOF,DOF)
        real*8 :: atmwgt(Natom),atmwgt_2(3*Natom)
        character(len=79) :: order,a
        real*8,parameter :: amuTokg=1.6606d-27
        real*8,parameter :: aTom=1.0d-10

        order="grep 'AtmWgt' reactant.log > atmwgt"
        call system(order)

        open(200,file="atmwgt")
        if (Natom<=10)then
           read(200,*)a,atmwgt(:)
        else
           read(200,*)a,atmwgt(1:10)
           read(200,*)a,atmwgt(11:Natom)
        end if
        close(200)

        print*,"The weight of atom:",atmwgt
        
        do i=1,Natom
           atmwgt_2(3*(i-1)+1:3*(i-1)+3)=atmwgt(i)
        end do

        atmwgt_2=atmwgt_2*amuTokg!! dimension:kg
        atmwgt_2=dsqrt(atmwgt_2)

        coord_r=coord_r*aTom
        coord_p=coord_p*aTom !!! dimension:m

        coord_r=coord_r*atmwgt_2
        coord_p=coord_p*atmwgt_2
        
        delta_coord=coord_p-coord_r
       
        proj_r=dot_product(delta_coord,coord_r)*delta_coord/&
               dot_product(delta_coord,delta_coord)
 
        do i=1,DOF
           do j=1,DOF
                D(i,j)=dot_product(mode_p(j,:),delta_coord)+&
                       dot_product(proj_r,mode_p(j,:)-mode_r(i,:))
           end do
        end do

end subroutine Displacement



subroutine FC_factor(DOF,Omega_R,Omega_P,D,F)
        implicit none
        integer :: DOF,i,j
        real*8 :: Omega_R(DOF),Omega_P(DOF) !J^-1 * s^-2
        real*8 :: D(DOF,DOF)  !kg^0.5 * m
        real*8 :: Alpha_R(DOF),Alpha_P(DOF) !kg-1 * m^-2
        real*8 :: F(DOF,DOF)  !Dimensionless
        real*8,parameter :: hbar=1.05457266d-34 !(J*s)
        real*8,parameter :: cmToHZ=2.997925d10  !(HZ,s^-1)


        Alpha_R=Omega_R*cmToHZ/hbar
        Alpha_P=Omega_P*cmToHZ/hbar

        do i=1,DOF
            do j=1,DOF
                F(i,j)=dsqrt((4*Alpha_R(i)*Alpha_P(j))/(Alpha_R(i)+Alpha_P(j))**2)*&
                       exp(0.0d0-Alpha_R(i)*Alpha_P(j)*D(i,j)**2/(Alpha_R(i)+Alpha_P(j)))
            end do
        end do
end subroutine FC_factor



recursive subroutine Back_tracking(DOF, F, row, current_sum, selected,&
                                   best_selection, max_sum, row_sums)
    implicit none
    integer :: DOF,i
    real*8 :: F(DOF, DOF)
    integer, intent(in) :: row
    real*8, intent(in) :: current_sum
    integer, intent(inout) :: selected(DOF), best_selection(DOF)
    real*8, intent(inout) :: max_sum
    real*8, intent(inout) :: row_sums(DOF)  ! 新增：存储每一行的和
    integer :: col

    if (row > DOF) then
        if (current_sum > max_sum) then
            max_sum = current_sum
            best_selection = selected  ! 保存当前最优选择
            ! 计算每一行的和
            do i = 1, DOF
                row_sums(i) = F(i, selected(i))
            end do
        end if
        return
    end if

    do col = 1, DOF
        if (all(selected(1:row-1) /= col)) then
            selected(row) = col
            call Back_tracking(DOF, F, row+1, current_sum + F(row, col),&
                 selected, best_selection, max_sum, row_sums)
            selected(row) = 0  ! 回溯
        end if
    end do
end subroutine Back_tracking


