program main
        implicit none

        integer,parameter :: Natom=8
        integer :: DOF,n,i,j,k,l,m,find_result,Z(1),n_mod
        character(len=79) :: line,order_1,order_2,symbol(Natom),a,b
        real*8 :: E,DeffCon,n1,n2,n3,n4,n5,n6
        real*8,dimension(:,:),allocatable :: mode_read
        real*8,dimension(:,:,:),allocatable :: mode_R,mode_P,mode_P_New
        real*8,dimension(:),allocatable :: mode_R_one,mode_P_one
        real*8,dimension(:),allocatable :: mode_R_one_square,mode_R_Normal
        real*8,dimension(:),allocatable :: mode_P_one_square,mode_P_Normal
        real*8,dimension(:,:),allocatable :: U
        real*8,dimension(:),allocatable :: mode_test
        real*8,dimension(:,:),allocatable :: mode_multi
        real*8,dimension(:),allocatable :: overlap
        real*8 :: Rsum,Psum
        real*8,dimension(:,:),allocatable :: Coord_R,Coord_P,Delta_Coord
        real*8,dimension(:,:),allocatable :: ts_guess
        real*8,dimension(:),allocatable :: R,P
        real*8,dimension(:,:),allocatable :: Coor_change_positive
        real*8,dimension(:,:),allocatable :: Coor_change_negative
        real*8,dimension(:),allocatable :: projection_R,projection_P
        real*8,dimension(:,:,:),allocatable :: mode_n

!!!========================================================!!!
!!!========================================================!!!
!!!          Only suit for Nonlinear molecules             !!!
!!!========================================================!!!
!!!========================================================!!!


        DOF=3*Natom-6     !Dgree of Freedom 
 
        write(line,*)Natom
        write(*,*)"Natom=",line
        write(*,*)"DOF=",DOF

        order_1="grep -A "//trim(adjustl(line))//" 'Atom  AN' reactant.log > mode_reactant"
        
       print*,order_1
        call system(order_1)

        order_2="grep -A "//trim(adjustl(line))//" 'Atom  AN' product.log > mode_product"

        print*,order_2
!!!!    grep -A 3 'Atom  AN' reactant.log > mode_reactant                              
!!!!     grep -A 3 'Atom  AN' product.log > mode_product
     call system(order_2)

!!=====================================================================
!!         Read the Coordinates of Reactant and Product              !!                                                            !
!!=====================================================================

        allocate(Coord_R(natom,3),Coord_P(natom,3),Delta_Coord(natom,3))
        allocate(ts_guess(natom,3))
        allocate(R(3*Natom),P(3*Natom))
        allocate(mode_read(Natom,3*DOF),mode_R(DOF,Natom,3),mode_P(DOF,Natom,3)) 
        allocate(mode_R_one(DOF*Natom*3),mode_P_one(DOF*Natom*3))
        allocate(mode_R_one_square(DOF*Natom*3),mode_R_Normal(DOF*Natom*3))
        allocate(mode_P_one_square(DOF*Natom*3),mode_P_Normal(DOF*Natom*3))
        allocate(U(DOF*Natom*3,DOF*Natom*3))
        allocate(mode_test(DOF*Natom*3))
        allocate(mode_multi(natom,3),mode_P_New(DOF,Natom,3))
        allocate(Coor_change_positive(Natom,3),Coor_change_negative(Natom,3))
        allocate(projection_R(DOF),projection_P(DOF))
        allocate(overlap(DOF))
        

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

        k=0
        do i=1,Natom
           do j=1,3
                k=k+1
                R(k)=Coord_R(i,j)
                P(k)=Coord_P(i,j)
           end do
        end do
        DeffCon=dsqrt(dot_product(P,P))-dsqrt(dot_product(R,R))
        write(*,*)"The difference of construction:",DeffCon


!!================================================================
!!       Read the mode in the reactants and products            !!
!!================================================================ 

        n_mod=mod(DOF,3)
        n=(DOF-n_mod)/3

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
                mode_R(i,j,:)=mode_read(j,3*(i-1)+1:3*i)
           end do
        end do
 

    !   do i=1,DOF 
    !      write(*,*)"The mode of reactant:"
    !      do j=1,Natom
    !         write(*,"(3f10.6)")mode_R(i,j,:)
    !      end do
    !         write(*,*)
    !   end do


        l=1
        do i=1,DOF
             do j=1,Natom
                  do k=1,3
                        mode_R_one(l)=mode_R(i,j,k)       
                        l=l+1
                  end do 
            end do
        end do  



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


    !   do i=1,DOF
    !      write(*,*)"The mode of product:"
    !      do j=1,Natom
    !         write(*,"(3f10.6)")mode_P(i,j,:)
    !      end do
    !         write(*,*)
    !   end do

        l=1
        do i=1,DOF
             do j=1,Natom
                  do k=1,3
                        mode_P_one(l)=mode_P(i,j,k)
                        l=l+1
                  end do
            end do
        end do

        deallocate(mode_n)
!!==============================================================
!!             Normalized the mode vector                     !!
!!          Compute the translational matrix                  !!
!!==============================================================
    !   do i=1,3*DOF*Natom
    !        mode_R_one_square(i)=mode_R_one(i)**2
    !   end do
    !   mode_R_Normal=mode_R_one/(sqrt(sum(mode_R_one_square)))
    !   
    !   do i=1,3*DOF*Natom
    !        mode_P_one_square(i)=mode_P_one(i)**2
    !   end do
    !   mode_P_Normal=mode_P_one/(sqrt(sum(mode_P_one_square)))

    !   Rsum=0.0d0
    !   Psum=0.0d0
    !   do i=1,3*DOF*Natom
    !         Rsum=Rsum+mode_R_Normal(i)**2
    !         Psum=Psum+mode_P_Normal(i)**2
    !   end do
    !   print*,"mode_R_Normal:"
    !   print*,mode_R_Normal
    !   print*,"mode_P_Normal:"
    !   print*,mode_P_Normal
    !   print*,"Rsum-Normalize:",Rsum,"Psum-Normalize:",Psum
!!==================================!!
!!   the first dimension is column  !! 
!!   the second dimension is line   !!
!!==================================!!
!!Uij=R_j*P_i

!       do i=1,DOF*Natom*3
!               do j=1,DOF*Natom*3
!                       U(j,i)=mode_R_Normal(j)*mode_P_Normal(i)
!               end do
!       end do

!       mode_test=matmul(mode_R_Normal,U)
!! The latter multiplies the former !!
!      
!       print*,"Taranslational Matrix:"
!  !     print*,U 
!       print*,"mode_test:"
!       print*,mode_test
!       print*,"mode_P_Normal:"
!       print*,mode_P_Normal


!===============================================================
!   Compare the vibration mode and calculate the overlap
!                (dot product and sum)
!===============================================================
        
        mode_multi=0.0d0
        do i=1,DOF
            do j=1,DOF
               do k=1,Natom
                    mode_multi(k,:)=mode_R(i,k,:)*mode_P(j,k,:)
                end do
            overlap(j)=abs(sum(mode_multi(:,:)))
            end do
            Z= maxloc(overlap)
          ! write(*,*)"Reactant mode",i,"Product mode",Z,"is the &
          !              same vibration mode"
          ! write(*,*)"The value of Overlap is",overlap(:)
          ! write(*,*)
        end do

!===============================================================
!   The coordinate vector difference is projected to mode
!===============================================================
        do i=1,DOF
           do j=1,Natom
            projection_R(i)=projection_R(i)+sum(mode_R(i,j,:)*Delta_Coord(j,:))
            projection_P(i)=projection_P(i)+sum(mode_P(i,j,:)*Delta_Coord(j,:))
           end do
        end do
        
        print*,"projection_R:",projection_R
        print*,"projection_P:",projection_P
        
        ts_guess(:,:)=Coord_R(:,:)&
        +0.4585*mode_R(4,:,:)-0.4082*mode_R(2,:,:)+0.2208*mode_R(8,:,:)
        
        write(*,*)"TS_guess:"
        do i=1,Natom
           write(*,"(a4,3f15.8)")Symbol(i),ts_guess(i,:)
        end do
        
 !!!    call modePES(Coord_R,Coord_P,mode_R,mode_P,Symbol,DOF,Natom)


!!===============================================================!!
!!                  Deallocate Arrays                            !!
!!===============================================================!!

        deallocate(Coord_R,Coord_P,Delta_Coord,ts_guess)
        deallocate(mode_read,mode_R,mode_P,mode_multi,mode_P_New)
        deallocate(mode_R_one,mode_P_one)
        deallocate(mode_R_one_square,mode_R_Normal)
        deallocate(mode_P_one_square,mode_P_Normal)
        deallocate(U)
        deallocate(mode_test)
        deallocate(Coor_change_positive,Coor_change_negative)
        deallocate(projection_R,projection_P)
        deallocate(overlap)


end


subroutine modePES(Coord_R,Coord_P,mode_R,mode_P,Symbol,DOF,Natom)

!!========================================================================
!!   Calculate the harmonic oscillator model relaxation near mode       !!
!!========================================================================

!!========================================================================
!!   Equilibrium position vibration relaxation to obtain a mode PES     !!
!!========================================================================

        implicit none

        integer :: Natom
        integer :: DOF,n,i,j,k,find_result,point
        character(len=79) :: line,order_1,order_2,symbol(Natom),a,b
        real*8 :: E
        real*8 :: mode_R(DOF,Natom,3),mode_P(DOF,Natom,3)
        real*8 :: overlap
        real*8 :: Coord_R(natom,3),Coord_P(natom,3),Delta_Coord(natom,3)
        real*8 :: Coor_change_positive(natom,3)
        real*8 :: Coor_change_negative(natom,3)
        character(len=79) :: data_(DOF,30)
!!!!!!!!!!!!!!Reactant!!!!!!!!!!!!!!

point=15

       do i=1,DOF
          do j=1,point
               Coor_change_positive(:,:)=Coord_R(:,:)+j*0.1*mode_R(i,:,:)         
               Coor_change_negative(:,:)=Coord_R(:,:)-(16-j)*0.1*mode_R(i,:,:)
               
               open(110,file="coor_positive")
               do k=1,natom
                  write(110,"(a4,3f12.7)")Symbol(k),Coor_change_positive(k,:)
               end do
                  write(110,*)
                  write(110,*)
               close(110) 
               call system("cat head.dat coor_positive > positive.com")

               open(120,file="coor_negative")
               do k=1,natom
                  write(120,"(a4,3f12.7)")Symbol(k),Coor_change_negative(k,:)
               end do
                  write(120,*)
                  write(120,*)
               close(120)
               call system("cat head.dat coor_negative > negative.com")

               call system("g16 positive.com")
               call system("g16 negative.com")

               open(130,file="positive.log")
               call located(130,"SCF Done:",find_result)
               open(140,file="E-positive-reactant",position="append")
               if (find_result/=0) then
                  read(130,"(25x,f16.10)")E
                  write(140,*)E
               else
                  write(140,*)"NO label!"
               end if
               close(130)
               close(140)
               
               open(150,file="negative.log")
               call located(150,"SCF Done:",find_result)
               open(160,file="E-negative-reactant",position="append")
               if (find_result/=0) then
                  read(150,"(25x,f16.10)")E
                  write(160,*)E
               else
                  write(160,*)"NO label!"
               end if
               close(150)
               close(160)


          end do
          open(140,file="E-positive-reactant",position="append")
          write(140,*)
          close(140)
          open(160,file="E-negative-reactant",position="append")
          write(160,*)
          close(160)
       end do

        open(160,file="E-negative-reactant")
        do i=1,DOF
           do j=1,point
                read(160,"(a)")data_(i,j)
           end do
                read(160,*)
        end do
        close(160,status="delete")


        open(140,file="E-positive-reactant")
        do i=1,DOF
           do j=1,point
                read(140,"(a)")data_(i,j+point)
           end do
                read(140,*)
        end do
        close(140,status="delete")
        
        open(141,file="E-Reactant")
        do i=1,DOF
            write(141,*)i
            write(141,*)"=========================="
            do j=1,point*2
                write(141,*)data_(i,j)
            end do
                write(141,*)
        end do
        close(141)


!!! ===============================================================
!!!!!!!!!!!!!!!!Product!!!!!!!!!!!!!!

        do i=1,DOF
          do j=1,point
               Coor_change_positive(:,:)=Coord_P(:,:)+(j)*0.1*mode_P(i,:,:)
               Coor_change_negative(:,:)=Coord_P(:,:)-(16-j)*0.1*mode_P(i,:,:)

               open(115,file="coor_positive")
               do k=1,natom
                  write(115,"(a4,3f12.7)")Symbol(k),Coor_change_positive(k,:)
               end do
                  write(115,*)
                  write(115,*)
               close(115)
               call system("cat head.dat coor_positive > positive.com")

               open(125,file="coor_negative")
               do k=1,natom
                  write(125,"(a4,3f12.7)")Symbol(k),Coor_change_negative(k,:)
               end do
                  write(125,*)
                  write(125,*)
               close(125)
               call system("cat head.dat coor_negative > negative.com")

               call system("g16 positive.com")
               call system("g16 negative.com")

               open(135,file="positive.log")
               call located(135,"SCF Done:",find_result)
               open(145,file="E-positive-product",position="append")
               if (find_result/=0) then
                  read(135,"(25x,f16.10)")E
                  write(145,*)E
               else
                  write(145,*)"NO label!"
               end if
               close(135)
               close(145)

               open(155,file="negative.log")
               call located(155,"SCF Done:",find_result)
               open(165,file="E-negative-product",position="append")
               if (find_result/=0) then
                  read(155,"(25x,f16.10)")E
                  write(165,*)E
               else
                  write(165,*)"NO label!"
               end if
               close(155)
               close(165)

          end do
          open(145,file="E-positive-product",position="append")
          write(145,*)
          close(145)
          open(165,file="E-negative-product",position="append")
          write(165,*)
          close(165)
       end do

        open(165,file="E-negative-product")
        do i=1,DOF
           do j=1,point
                read(165,"(a)")data_(i,j)
           end do
                read(165,*)
        end do
        close(165,status="delete")


        open(145,file="E-positive-product")
        do i=1,DOF
           do j=1,point
                read(145,"(a)")data_(i,j+point)
           end do
                read(145,*)
        end do
        close(145,status="delete")

        open(146,file="E-Product")
        do i=1,DOF
            write(146,*)i
            write(146,*)"=========================="
            do j=1,point*2
                write(146,*)data_(i,j)
            end do
                write(146,*)
        end do
        close(146)



end subroutine modePES

!subroutine Deallocate_
!!!===============================================================!!
!!!===============================================================!!
!       deallocate(Coord_R,Coord_P,Delta_Coord)
!       deallocate(mode_read,mode_R,mode_P,mode_multi,mode_P_New)
!       deallocate(Coor_change_positive,Coor_change_negative)
!       deallocate(projection)


!end subroutine Deallocate_


subroutine located(filed,label,find_result)
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
