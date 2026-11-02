subroutine zmat_to_cart(natom, mzmat, molzmat, cartcoord, atom)
    implicit none
    integer, intent(inout) :: natom
    integer, dimension(3,natom), intent(in) :: mzmat
    real*8, dimension(3,natom), intent(in) :: molzmat
    real*8, dimension(3,natom), intent(inout) :: cartcoord
    character(len=1), dimension(natom), intent(inout) :: atom

    ! 定义常量 pi
    real*8, parameter :: pi = 3.14159265358979323846d0

    ! 声明局部变量
    integer :: iatom, jatom
    real*8 :: x(3), r, theta, phi, a(3), b(3), c(3)
    real*8 :: ab(3), bc(3), n(3), ncbc(3), norm
    real*8 :: M(3,3)


    do iatom = 1, natom
        if(any(mzmat(1:3,iatom).ge.iatom))goto 999
        select case(iatom)
        case(1)
            cartcoord(1:3,iatom)=0.d0
        case(2)
            cartcoord(1,iatom)=molzmat(1,iatom)
            cartcoord(2:3,iatom)=0.d0
        case(3)
            jatom = mzmat(1,iatom)
            if(jatom.eq.1)then
                cartcoord(1,iatom)=molzmat(1,iatom)*&
                    cosdeg(molzmat(2,iatom))
            else
                cartcoord(1,iatom)=cartcoord(1,jatom)-&
                    molzmat(1,iatom)*cosdeg(molzmat(2,iatom))
            endif
            cartcoord(2,iatom)=molzmat(1,iatom)*&
                sindeg(molzmat(2,iatom))
            cartcoord(3,iatom)=0.d0
        case default
            ! the NeRF algorithm, see Reference
            a = cartcoord(1:3,mzmat(3,iatom))
            b = cartcoord(1:3,mzmat(2,iatom))
            c = cartcoord(1:3,mzmat(1,iatom))
            r = molzmat(1,iatom)
            theta = molzmat(2,iatom)
            phi = molzmat(3,iatom)

            x(1) = r*cosdeg(theta)
            x(2) = r*sindeg(theta)*cosdeg(phi)
            x(3) = r*sindeg(theta)*sindeg(phi)

            ab = b-a
            bc = c-b
            call norm3d(bc,norm)
            bc = bc/norm
            call cross3d(ab,bc,n)
            call norm3d(n,norm)
            n = n/norm
            call cross3d(n,bc,ncbc)

            M(1:3,1) = bc
            M(1:3,2) = ncbc
            M(1:3,3) = n

            x(1) = -x(1)
            cartcoord(1,iatom) = c(1) + sum(M(1,1:3)*x)
            cartcoord(2,iatom) = c(2) + sum(M(2,1:3)*x)
            cartcoord(3,iatom) = c(3) + sum(M(3,1:3)*x)
        end select
    enddo

    ! 剥离掉任何虚拟原子（X）
    iatom = natom
    do while(iatom.ge.1)
        if(atom(iatom).eq.'X')then
            do jatom = iatom,natom-1
                cartcoord(1:3,jatom) = cartcoord(1:3,jatom+1)
                atom(jatom) = atom(jatom+1)
            enddo
            natom = natom-1
        endif
        iatom = iatom-1
    enddo

    ! 您需要在调用子程序时自行进行Angstrom到Bohr的转换
    return

999 continue
    write(*,*) 'zmat_to_cart: Error reading Z-matrix at atom ', iatom
!    call qexit()

    contains
        real*8 function cosdeg(ang)
            real*8, intent(in) :: ang
            cosdeg = dcos(ang*pi/180.d0)
        end function cosdeg
        real*8 function sindeg(ang)
            real*8, intent(in) :: ang
            sindeg = dsin(ang*pi/180.d0)
        end function sindeg

end subroutine zmat_to_cart

subroutine norm3d(vec, vecnorm)
implicit none
! dummy
real*8, intent(in)  :: vec(3)
real*8, intent(out) :: vecnorm

vecnorm = dsqrt(vec(1)*vec(1) + vec(2)*vec(2) + vec(3)*vec(3))

end subroutine norm3d


! Dot product of two 3-D vectors: not implemented, because internal function
! dot_product already suffices


! Cross product of two 3-D vectors
! Adapted from icoord_mod/rdic_common_mod.F90, subroutine cross_product
subroutine cross3d(vec1, vec2, vecout)
implicit none
! dummy
real*8,intent(in)   :: vec1(3),vec2(3)
real*8,intent(out)  :: vecout(3)

vecout(1) = vec1(2)*vec2(3) - vec1(3)*vec2(2)
vecout(2) = vec1(3)*vec2(1) - vec1(1)*vec2(3)
vecout(3) = vec1(1)*vec2(2) - vec1(2)*vec2(1)

end subroutine cross3d

