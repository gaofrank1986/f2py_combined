
module gen_matrix
    
    use mesh!mesh data
    use wave 
    use mod_func!supplement funcs
    !use matrix_funcs
    implicit none
    
    real(8),allocatable :: amata(:,:,:),bmata(:,:)
    real(8),allocatable :: param_local(:,:),param(:,:)
    real(8),parameter :: ex(4) = (/1.,1.,-1.,-1./)
    real(8),parameter :: ey(4) = (/1.,-1.,-1.,1./)
    real(8),parameter :: rsn(4,4) = reshape((&
                        &/1.,  1.,  1.,  1.,& 
                        & 1., -1.,  1., -1.,&
                        & 1.,  1., -1., -1.,& 
                        & 1., -1., -1.,  1./),(/4,4/)) 
contains
    include './mfunc/solid.f90'    
    subroutine test()
            print *, rsn
    end subroutine
    
    subroutine calc_matrix()
        implicit none
       
        integer :: inode,ii,ielem,i
        real(8) :: xp,yp,zp,s_angle,BMATRIX(4,8),AMATRIX(4,8) 
        real(8) :: angle(nnode),phi,dpox,dpoy,dpoz
    
        allocate(amata(nnode,nnode,nsys))
        allocate(bmata(nnode,nsys))
        allocate(param_local(nsys,4),param(inode,4)) 
        amata = 0
        bmata = 0        
        
        do inode = 1,nnode
            XP=XYZ(1,INODE)
            YP=XYZ(2,INODE)
            ZP=XYZ(3,INODE)
            CALL SOLIDANGLE(INODE,NNODE,NELEM,&
            &NCN,NCON,NODQUA,H,XYZ,DXYZE,S_ANGLE)    
            S_ANGLE=1.0d0-S_ANGLE
            angle(inode) = s_angle
            amata(inode,inode,1:nsys) = s_angle
        end do

        do inode = 1,nnf
        write(110,*) "bmata",bmata(1,1)
            do ielem = 1,nelemf

                call comp_link(ielem,inode,ii)
   !             call wrapper_func(0,ielem,inode,xp,yp,zp,amatrix,bmatrix,ii) 
                call intg_branch(0,ielem,inode,xp,yp,zp,amatrix,bmatrix,ii) 
                call common_block(0,0,ielem,inode,amatrix,bmatrix)
            end do ! ielem = 1,nelemf

            do ielem = nelemf+1,nelem

                call comp_link(ielem,inode,ii) 
  !              call wrapper_func(0,ielem,inode,xp,yp,zp,amatrix,bmatrix,ii) 
                call  intg_branch(1,ielem,inode,xp,yp,zp,amatrix,bmatrix,ii) 
                call common_block(1,0,ielem,inode,amatrix,bmatrix)
            end do!ielem = nelemf+1,nelem

            param(inode,:) = param_local(1,:)
            !1-> A3 2->C31 3>C32 4->C33
            
            PHI=POXY(XP,YP,ZP)

            CALL DINP(XP,YP,ZP,DPOX,DPOY,DPOZ)       
            BMATA(INODE,1)=BMATA(INODE,1)-param(inode,1)*PHI-&
                            &param(INODE,2)*DPOX-param(INODE,2)*DPOY



 
        end do! inode   
        do inode = nnf+1,nnode
        write(110,*) "bmata",bmata(1,1)
            do ielem = 1,nelemf
                ii = 0
!                call wrapper_func(1,ielem,inode,xp,yp,zp,amatrix,bmatrix,ii) 
                call  intg_branch(1,ielem,inode,xp,yp,zp,amatrix,bmatrix,ii) 
                call common_block(1,1,ielem,inode,amatrix,bmatrix)
            end do!ielem = 1,nelemf

            do ielem = nelemf+1,nelem
                call comp_link(ielem,inode,ii)
 !               call wrapper_func(1,ielem,inode,xp,yp,zp,amatrix,bmatrix,ii)
                call  intg_branch(1,ielem,inode,xp,yp,zp,amatrix,bmatrix,ii) 
                call common_block(0,1,ielem,inode,amatrix,bmatrix)
            end do
        enddo
    end subroutine calc_matrix        
    !==========================================================
    !=========================================================
    !==
    !==
    !==
    !=========================================================
    
    subroutine intg_branch(flag,ielem,inode,xp,yp,zp,amatrix,bmatrix,ii) 
                    
        implicit none

        integer,intent(in) :: flag,ielem,inode,ii
        real(8),intent(in) :: xp,yp,zp
        real(8),intent(out) :: amatrix(4,8),bmatrix(4,8)

        if (flag.eq.0) then
            IF (II .EQ. 0)   THEN 
            !CALL NORM_ELE1(IELEM,XP,YP,ZP,AMATRIX,BMATRIX)
                print *,"call norm_elem1"
            ELSE 
            !CALL SING_ELE1(INODE,IELEM,NODQUA(INODE),XP,YP,ZP,&
            !                  & AMATRIX,BMATRIX)
                print *,"call sing_elem1"
            END IF                
        else
            IF (II .EQ. 0)   THEN 
            !CALL NORM_ELE0(IELEM,XP,YP,ZP,AMATRIX,BMATRIX)
                print *,"call norm_elem0"
            ELSE 
            !CALL SING_ELE0(INODE,IELEM,NODQUA(INODE),XP,YP,ZP,&
            !                  & AMATRIX,BMATRIX)
                print *,"call sing_elem0"
            END IF                
        end if
    end subroutine 

    subroutine comp_link(ielem,inode,ii)
        implicit none
        integer,intent(in) :: ielem,inode
        integer,intent(out):: ii
        integer :: i
        II=0
        DO  I=1, NODNOE(INODE)
            IF(IELEM .EQ. NODELE(INODE,I)) THEN
                II=II+1
            ENDIF
        ENDDO
    end subroutine 

    !==========================================================
    !=========================================================
    !==
    !==
    !==
    !=========================================================

    subroutine common_block(r_flag,s_flag,ielem,inode,amatrix,bmatrix)

        implicit none
        
        integer,intent(in) :: r_flag,s_flag,ielem,inode        
        integer :: j,cur_node,cur_nrml,ip,i,i_param,is
        real(8) :: xsb,ysb,zsb,nx,ny,nz,dpdn,dpox,dpoy,dpoz,phi
        real(8),intent(in) :: amatrix(4,8),bmatrix(4,8)
        do J=1,  NCN(IELEM) 

            cur_node = NCON(IELEM,J)!current node id
            cur_nrml = NCOND(IELEM,J)!current nrml id

            do  ip=1, nsys 
                xsb=ex(ip)*xyz(1,cur_node)
                ysb=ey(ip)*xyz(2,cur_node)
                zsb=       xyz(3,cur_node)

                nx=ex(ip)*dxyz(1,cur_nrml)
                ny=ey(ip)*dxyz(2,cur_nrml)
                nz=       dxyz(3,cur_nrml)
                
                call dinp(xsb,ysb,zsb,dpox,dpoy,dpoz)       
                write(113,*) "xyzb",xsb,ysb,zsb
                write(113,*) "nxyz",nx,ny,nz
                dpdn=dpox*nx+dpoy
                write(113,*) "dpoxyzb",dpox,dpoy,dpoz
                dpdn=dpox*nx+dpoy**ny+dpoz*nz 
                write(111,*) "dpdn",dpdn
! (0,0) (1,0)
                do  is=1, nsys    
                    if (r_flag .eq. 0) then
                        amata(inode,cur_node,ip)=amata(inode,cur_node,ip)+&
                                         &rsn(is,ip)*bmatrix(is,j)
                        bmata(inode,ip)=bmata(inode,ip)+rsn(is,ip)*&
                                    &amatrix(is,j)*poxy(xsb,ysb,zsb)
                    else 
                        if(cur_node.gt. nnf)  then        
                            amata(inode,cur_node,ip)=amata(inode,cur_node,ip)-&
                                     &rsn(is,ip)*amatrix(is,j)
                        else
                                phi=poxy(xsb,ysb,zsb)
                                write(111,*) "phi",phi
                                bmata(inode,ip)=bmata(inode,ip)+rsn(is,ip)*amatrix(is,j)*phi
                        endif
                        bmata(inode,ip)=bmata(inode,ip)-rsn(is,ip)*bmatrix(is,j)*dpdn 
                    end if
                enddo!is
            
                if (s_flag .eq.0) then
                    do i_param = 0,3                
                        call dinp0(i,xsb,ysb,zsb,phi,dpox,dpoy,dpoz)       
                        dpdn=dpox*nx+dpoy*ny+dpoz*nz
                        write(111,*) "dpdn,phi",dpdn,phi
                        param_local = 0
                        do  is=1, nsys             
                            param_local(ip,i_param+1)=param_local(ip,i_param+1)-rsn(is,ip)&
                                                        &*bmatrix(is,j)*dpdn 
                            param_local(ip,i_param+1)=param_local(ip,i_param+1)+rsn(is,ip)&
                                                        &*amatrix(is,j)*phi
                        enddo!is
                    end do!i_param
                end if

        end do;end do !ip/j
    end subroutine
end module
