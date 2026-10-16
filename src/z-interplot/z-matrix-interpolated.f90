program main
        implicit none
        integer :: natom,nimage,line,i,j,k
        integer,allocatable,dimension(:,:) :: info_connection
        integer,allocatable,dimension(:) :: info_connection_1
        character(len=4) :: a,b,c,d
        character(len=4),allocatable,dimension(:,:) :: info
        character(len=6),allocatable,dimension(:) :: natom_symbol,connect_symbol   
        character(len=79) :: test
        character(len=79),allocatable,dimension(:) :: Cartesian
        real*8,allocatable,dimension(:) :: para_R,para_P,delta_para,para_I
        real*8,allocatable,dimension(:) :: para_P_3n,para_R_3n

        natom=8
        nimage=7
        line=3*natom-6
       
        allocate(info(natom,8))
        allocate(natom_symbol(natom))
        allocate(info_connection(natom,3))
        allocate(info_connection_1(natom*3))
        allocate(connect_symbol(line),para_R(line),para_P(line))
        allocate(para_R_3n(3*natom),para_P_3n(3*natom))
        allocate(delta_para(3*natom),para_I(3*natom))
        
        allocate(Cartesian(natom))

!get information of head.dat(connection information)
        call system("rm coor-z")

        info_connection=0.0d0
        info="O"

        if (natom>=3)then
        open(50,file="info_connection")
           read(50,"(a)")natom_symbol(1)
           info_connection(1,:)=0   
        
           read(50,*)info(2,1:3)
           info_connection(2,:)=0
           natom_symbol(2)=info(2,1)
           a=info(2,2)
           read(a,*)info_connection(2,1)
           

           read(50,*)info(3,1:5)
           info_connection(3,:)=0
           natom_symbol(3)=info(3,1)
           a=info(3,2)
           read(a,*)info_connection(3,1)
           b=info(3,4)
           read(b,*)info_connection(3,2)

        do i=4,natom,1
           info_connection(i,:)=0
           read(50,*)info(i,:)
           natom_symbol(i)=info(i,1)
           a=info(i,2)
           read(a,*)info_connection(i,1)
           b=info(i,4)
           read(b,*)info_connection(i,2)
           c=info(i,6)
           read(c,*)info_connection(i,3)
        end do
        endif
        close(50)
        
        open(100,file="initio.dat")
        do i = 1,line   
            read(100,*)connect_symbol(i),para_R(i)
        end do
        close(100)

        open(200,file="final.dat")
        do i = 1,line
           read(200,*)connect_symbol(i),para_P(i)
        end do
        close(200)

        k=1
        do i=1,natom
                do j=1,3
                  info_connection_1(k)=info_connection(i,j)
                  k=k+1
                end do
        end do
        
        ! 0  0  0               B1
        ! B1 0  0               ...
        ! B2 A1 0               B5
        ! B3 A2 D1  ========    A1
        ! B4 A3 D2   TRANS      ...      
        ! B5 A4 D3              A5
        ! B6 A5 D4              D1
        !                       ...
        !                       D4
                
        para_R_3n=0.0d0
        do i=1,3*natom-6 
            if (i<=natom-1)then
                para_R_3n(3*i+1)=para_R(i)
            else if ((natom-1) < i .and. i <= (2*natom-3)) then
                para_R_3n(3*(i-Natom+2)+2) = para_R(i)
                 !!  3*(i-(N-1)+1) + 2
                 !!        Bond
            else
                para_R_3n(3*(i-2*Natom+5)+3) = para_R(i)
                 !!  3*(i-(N-2)-(N-1)+2) +3
                 !!        Ager   Bond
            end if
        end do

        para_P_3n=0.0d0
        do i=1,3*natom-6
            if (i<=natom-1)then
                para_P_3n(3*i+1)=para_P(i)
            else if (natom-1 < i .and. i <= 2*natom-3) then
                para_P_3n(3*(i-Natom+2)+2) = para_P(i)
            else
                para_P_3n(3*(i-2*Natom+5)+3) = para_P(i)
            end if
        end do
        

         do i=1,nimage  !!!! two endpoints + interpolated point
                        !!!! Nimage from input_DL-FIND is 
                        !!!! two endpoints + interpolated point +  climbing image
            print*,i            
            delta_para=(para_P_3n-para_R_3n)/(nimage-1)
            para_I=delta_para*(i-1)+para_R_3n
      
            open(131,file="coor-z",position="append")
            write(131,*)"==========================="
            write(131,*)"the Nimage :",i

        do j=1,3*natom-6
            if (j<=natom-1)then
                write(131,*)connect_symbol(j),para_I(3*j+1)
            else if (natom-1 < j .and. j <= 2*natom-3) then
                write(131,*)connect_symbol(j),para_I(3*(j-Natom+2)+2)
            else
                write(131,*)connect_symbol(j),para_I(3*(j-2*Natom+5)+3)
            end if
        end do
            close(131)
       
        !get information of head.dat(connection information)
            open(300,file="input")
            write(300,*)natom
            write(300,*)info_connection_1
            write(300,*)para_I 
            write(300, '(100(a2,","))',advance='no') (natom_symbol(j), j =1,natom-1)
            write(300, '(a2)')natom_symbol(natom)
            close(300)
            call system("./zmat_to_cart")
             
        !read the Cartesian from Z-mat-to-Cartesian
            open(400,file="output.txt")
            do j=1,natom
               read(400,"(a)")Cartesian(j)        
            end do
            close(400)!,status="delete")

        !write the initio guess Cartesian of 
            open(500,file="xyz.dat",position="append")
            if(i==1)then
                write(500,*)"geometry 1"
                do j=1,natom
                   write(500,"(a)")Cartesian(j)
                end do
                write(500,*)"end geometry 1"
                write(500,*)
            else if (i==2)then
                write(500,*)"geometry 2"
                write(500,*)"#geom of nimage",i
                do j=1,natom
                   write(500,"(a)")Cartesian(j)
                end do
            else
                write(500,*)"#geom of nimage",i
                do j=1,natom
                   write(500,"(a)")Cartesian(j)
                end do
            end if
            close(500)
            
        end do
            call system("cat 'head.dat' 'xyz.dat' 'tail.dat' > input_DL-FIND") 

        deallocate(info_connection)
        deallocate(info_connection_1)
        deallocate(natom_symbol)
        deallocate(connect_symbol)
        deallocate(para_R,para_P,delta_para,para_I)
        deallocate(Cartesian)





end


subroutine located(filed,label,line)
        character(len=*) :: label
        character(len=79) :: para
        integer :: find_result,error,filed,line
        line = 0
        rewind (filed)
          do while (.true.)
            line=line+1
            read(filed, "(A79)",iostat=error)para
            if (error/=0) stop "routine located: label not found!"
            find_result=index( para,label)
            if ( find_result  > 0 ) exit
          end do
        backspace(filed)
        return
end subroutine

