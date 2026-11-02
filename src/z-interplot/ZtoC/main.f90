        program main
        implicit none
        integer :: natom
        integer, allocatable :: mzmat(:,:)
        real*8, allocatable :: molzmat(:,:)
        real*8, allocatable :: cartcoord(:,:)
        character(len=1), allocatable :: atom(:) 
        integer :: i, j , pos
        integer, allocatable :: mz(:)
        real*8, allocatable :: molz(:)
        character(10000) :: line
        character(100) :: atom_names
        open(11, file='input', status='old') ! 打开输入文件

        ! 读取原子数
        read(11, *) natom

        !分配数组大小
        allocate(mzmat(3,natom))
        allocate(molzmat(3,natom))
        allocate(cartcoord(3,natom))
        allocate(atom(natom))
        allocate(mz(3*natom))
        allocate(molz(3*natom))
        ! 读取 Z 矩阵坐标
         read(11,*)mz(:) 
         mzmat = reshape(mz, [3, natom])
        ! 读取分子坐标
        read(11,*)molz(:)
         molzmat = reshape(molz, [3, natom])
        ! 读取原子符号
        read(11,*)line
        do i = 1, natom
        ! 寻找下一个逗号的位置
        pos = index(line, ',')
        if (pos == 0) then
            atom_names = trim(line)  !如果找不到逗号，说明是最后一个原子名称
        else
            atom_names = trim(adjustl(line(:pos-1)))  !从开头到逗号前的部分是一个原子名称
            line = line(pos+1:)  ! 更新字符串，去掉已经读取的部分
        end if
        atom(i) = atom_names  ! 将原子名称存储到数组中
        end do

         ! 关闭输入文件
       close(11) 

      ! 调用子程序进行坐标转换
        call zmat_to_cart(natom, mzmat, molzmat, cartcoord, atom)

       open(12, file='output.txt', status='replace') ! 打开输出文件

       ! 写入笛卡尔坐标到输出文件
       do j = 1, natom
        write(12, '(3F10.5)') cartcoord(1:3, j)
       end do
       close(12) ! 关闭输出文件
       call system('bash format_out.sh')
       end program main

