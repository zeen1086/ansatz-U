program rearrange_HR

        implicit none
        character(len=1024) :: the_input

        complex(kind(0d0)),allocatable  ::H_wan(:,:,:),H_k(:,:),H_buff(:,:),eig_vec(:,:,:)
        complex(kind(0d0)),allocatable  ::AAR(:,:,:,:)
        complex(kind(0d0)),allocatable  :: super_H_wan(:,:,:),super_AAR(:,:,:,:),Dh(:,:,:)
        complex(kind(0d0)),allocatable  :: H_check(:,:),u_check(:,:),delHH(:,:,:),M_tmp(:,:)
        real(kind(0d0)),allocatable::eig_check(:)
        integer :: ik,jk,nk,nkx,nky,nkz,irept,nrept,io,read3_int(3)
        integer :: super_nrept,super_nx,super_ny,super_nz,H_dim
        integer :: i,j,wan_x,wan_y,wan_z,nstep,ix,iy,iz,iwan,jwan
        integer :: i_dim,j_dim
        integer :: jx,jy,jz,period_x,period_y,period_z,i_shift,n1,n2,n3
        integer :: max_x,max_y,max_z,min_x,min_y,min_z,super_i,super_j
        integer :: super_x,super_y,super_z,super_rept,super_size
        real(kind(0d0)),allocatable     ::k_list(:,:),eigen(:,:)
        integer :: n_cell,nbnd,nwan,i_cell
        integer,allocatable     :: ir_vec(:,:),super_ir_vec(:,:)
        real(kind(0d0)),allocatable     :: nb(:,:)
        character(len=20)       :: irept_char
        complex(kind(0d0))      :: cI,phase1,phase2
        real(kind(0d0))         :: pi,dE
        logical, allocatable    :: rept_TF(:)
        integer,allocatable     :: deg(:),inv_k(:)
        character(len=200)      :: useless_string,nscf_dir,THE_dir,path



CALL execute_command_line('mkdir -p super_H' ) 
 



call get_command_argument(1,the_input)
if (len_trim(the_input)==0) error stop 'Supply explicit THE.in path as argument 1'
open(unit=101,file=trim(the_input),action='read',access='sequential',status='old')
read(101,*)!useless_string,func_dir
read(101,'(A9,A200)')useless_string,nscf_dir
read(101,'(A8,A200)')useless_string,THE_dir
read(101,*)!'(A7,I3)')  useless_string,n_func
read(101,*)!'(A13,I2)') useless_string,convert_size_x
read(101,*)!'(A13,I2)') useless_string,convert_size_y
read(101,*)!'(A13,I2)') useless_string,convert_size_z
read(101,'(A5,I3)') useless_string, nwan
read(101,*)! useless_string, nspin
read(101,*)!'(A3,I3)') useless_string, nkx
read(101,*)!'(A3,I3)') useless_string, nky
read(101,*)!'(A3,I3)') useless_string, nkz
read(101,*) ! A6 axial
read(101,'(A9,I3)')  useless_string, period_x
read(101,'(A9,I3)')  useless_string, period_y
read(101,'(A9,I3)')  useless_string, period_z
!read(101, useless_string, n_electron
!read(101,'(A4,F8.4)') useless_string, eta
!read(101,'(A4,F16.8)') useless_string, fac
!read(101,'(A6,I3)') useless_string,b_min
!read(101,'(A6,I3)') useless_string,b_max
close(101)



cI=dcmplx(0.0d0,1.0d0)
pi=3.1415926535d0

nbnd=nwan


irept=0
open(unit=101,file=trim(nscf_dir)//'./ir_vec_HYC',action='read',access='sequential')
do
        read(101,'(3I8)',iostat=io) read3_int(:)
                if(io/=0) exit
                irept=irept+1
        enddo
close(101)
        nrept=irept
        allocate(ir_vec(nrept,3))
        allocate(H_wan(nwan,nwan,nrept))
        allocate(deg(nrept))
allocate(H_check(nwan,nwan))
allocate(eig_check(nwan))
allocate(u_check(nwan,nwan))
allocate(delHH(nwan,nwan,3))
allocate(M_tmp(nwan,nwan))
allocate(Dh(nwan,nwan,3))
        irept=1
        open(unit=101,file=trim(nscf_dir)//'./ir_vec_HYC',action='read',access='sequential')
        do
                read(101,'(3I8)',iostat=io) read3_int(:)
                if(io/=0) exit
                ir_vec(irept,:)=read3_int(:)
                irept=irept+1
        enddo
        close(101)


open(unit=101,file=trim(nscf_dir)//'./H_deg_HYC',action='read',access='sequential')
do
        read(101,'(2I8)',iostat=io) i,j
        deg(i)=j
        if (io/=0)  exit
enddo
close(101)



do irept=1,nrept
        call int2char(irept, irept_char)
        open(unit=101,file=trim(nscf_dir)//'./hr_HYC_'//trim(irept_char),action='read'&
                        &,access='sequential',form='unformatted')
                read(101) H_wan(:,:,irept)
           ! print*,abs(H_wan(1,1:5,irept))    
        close(101)
 !       stop

enddo







allocate(AAR(nwan,nwan,nrept,3))



!! TO compute Hall effect!!
!! need postw90  !!

!open(unit=101,file=trim(nscf_dir)//'./AA_R_HYC',action='read',form='unformatted',access='sequential')
!read(101) read3_int(:)
!print*,read_int3
!read(101) AAR
!close(101)

!!!!!!!!!!!!!!!!!!!!!!


max_x=-99
max_y=-99
max_z=-99
min_x=99
min_y=99
min_z=99

do irept=1,nrept
        !if (ir_vec(irept,1) .gt. max_x) max_x=ir_vec(irept,1)
        !if (ir_vec(irept,2) .gt. max_y) max_y=ir_vec(irept,2)
        !if (ir_vec(irept,3) .gt. max_z) max_z=ir_vec(irept,3)
        if (ir_vec(irept,1)+period_x-1 .gt. max_x) max_x=ir_vec(irept,1)+period_x-1
        if (ir_vec(irept,2)+period_y-1 .gt. max_y) max_y=ir_vec(irept,2)+period_y-1
        if (ir_vec(irept,3)+period_z-1 .gt. max_z) max_z=ir_vec(irept,3)+period_z-1
        if (ir_vec(irept,1) .lt. min_x) min_x=ir_vec(irept,1)
        if (ir_vec(irept,2) .lt. min_y) min_y=ir_vec(irept,2)
        if (ir_vec(irept,3) .lt. min_z) min_z=ir_vec(irept,3)
enddo
!print*,max_x,max_y,max_z
!print*,min_x,min_y,min_z

max_x=max_x/period_x
max_y=max_y/period_y
max_z=max_z/period_z
min_x=(min_x+1)/period_x-1
min_y=(min_y+1)/period_y-1
min_z=(min_z+1)/period_z-1


super_nx=(max_x-min_x+1)
super_ny=(max_y-min_y+1)
super_nz=(max_z-min_z+1)

!print*,max_x,max_y,max_z
!print*,min_x,min_y,min_z
!stop

super_nrept=(max_x-min_x+1)*(max_y-min_y+1)*(max_z-min_z+1)
allocate(super_ir_vec(super_nrept,3))
allocate(super_H_wan(nwan*period_x*period_y*period_z,nwan*period_x*period_y*period_z&
                        &,super_nrept))

super_size=period_x*period_y*period_z

H_dim=nwan*period_x*period_y*period_z
allocate(super_AAR(H_dim,H_dim,super_nrept,3))


allocate(H_k(H_dim,H_dim))
allocate(H_buff(H_dim,H_dim))
allocate(eig_vec(1,H_dim,H_dim))
allocate(eigen(1,H_dim))
allocate(nb(period_x*period_y*period_z,3))

do ix=min_x,max_x
do iy=min_y,max_y
do iz=min_z,max_z
        irept=1+(iz-min_z)+super_nz*(iy-min_y)+&
                &super_nz*super_ny*(ix-min_x)
        super_ir_vec(irept,:)=(/ix,iy,iz/)
        !print*,super_ir_vec(irept,:)
enddo
enddo
enddo


do ix=0,period_x-1
do iy=0,period_y-1
do iz=0,period_z-1
        i=1+iz+period_z*iy+period_z*period_y*ix
        nb(i,:)=(/ix*1.0d0/period_x,iy*1.0d0/period_y,iz*1.0d0/period_z/)
enddo
enddo
enddo



super_H_wan(:,:,:)=dcmplx(0.0d0,0.0d0)
super_AAR(:,:,:,:)=dcmplx(0.0d0,0.0d0)

do ix=0,period_x-1
do iy=0,period_y-1
do iz=0,period_z-1
        print*,ix,iy,iz
        do iwan=1,nwan
                i=1+iz+period_z*iy+period_z*period_y*ix
                super_i=1+(iwan-1)+nwan*(i-1)

                do irept=1,nrept
                        jx=ix+ir_vec(irept,1)
                        jy=iy+ir_vec(irept,2)
                        jz=iz+ir_vec(irept,3)
                        if (jx .lt. 0) then
                                super_x=(jx+1)/period_x-1
                        else
                           super_x=jx/period_x
                        endif
                        if (jy .lt. 0) then
                                super_y=(jy+1)/period_y-1
                        else
                           super_y=jy/period_y
                        endif
                        if (jz .lt. 0) then
                                super_z=(jz+1)/period_z-1
                        else
                           super_z=jz/period_z
                        endif

                        super_rept=1+(super_z-min_z)+&
                                &super_nz*(super_y-min_y)+&
                                &super_nz*super_ny*(super_x-min_x)
                        
                        jx=jx-super_x*period_x
                        jy=jy-super_y*period_y
                        jz=jz-super_z*period_z
                        
                        j=1+jz+period_z*jy+period_z*period_y*jx

                        do jwan=1,nwan
                                super_j=1+(jwan-1)+nwan*(j-1)
                                super_H_wan(super_i,super_j,super_rept)=&
                                &H_wan(iwan,jwan,irept)/deg(irept)
                                super_AAR(super_i,super_j,super_rept,:)=&
                                &AAR(iwan,jwan,irept,:)/deg(irept)


                        enddo


                enddo
        enddo
enddo
enddo
enddo




do irept=1,super_nrept
        call int2char(irept,irept_char)
        open(unit=102,file='./super_H/hr_HYC_'//trim(irept_char),&
                &action='write',access='sequential',position='append',status='replace',form='unformatted')
                write(102) super_H_wan(:,:,irept)
        close(102)
enddo


open(unit=102,file='./super_H/AAR',&
                &action='write',access='sequential',position='append',status='replace',form='unformatted')
                write(102) super_AAR(:,:,:,:)
close(102)


open(unit=102,file='./super_H/ir_vec_HYC',action='write',access='sequential',position='append',status='replace')
do irept=1,super_nrept
        write(102,'(3I8)') super_ir_vec(irept,:)
enddo
close(102)



H_k(:,:)=0.0d0
do iwan=1,nwan
print*,iwan,nwan
do ix=0,period_x-1
do iy=0,period_y-1
do iz=0,period_z-1
        i=1+iz+period_z*iy+period_z*period_y*ix
        super_i=1+(iwan-1)+nwan*(i-1)
do jwan=1,nwan
do jx=0,period_x-1
do jy=0,period_y-1
do jz=0,period_z-1
        j=1+jz+period_z*jy+period_z*period_y*jx
        super_j=1+(jwan-1)+nwan*(j-1)

do n1=1,super_size
        i_dim=1+(iwan-1)+nwan*(n1-1)
        phase1=exp(-ci*2*pi*dot_product(1.0d0*(/ix,iy,iz/),nb(n1,:)))
do n2=1,super_size
        j_dim=1+(jwan-1)+nwan*(n2-1)
        phase2=phase1*exp(ci*2*pi*dot_product(1.0d0*(/jx,jy,jz/),nb(n2,:)))
do irept=1,super_nrept
        H_k(i_dim,j_dim)=H_k(i_dim,j_dim)+phase2*super_H_wan(super_i,super_j,irept)
enddo
enddo
enddo
enddo
enddo
enddo
enddo
enddo
enddo
enddo
enddo
H_buff=H_k/super_size
        call  diasym(H_buff,eigen(1,:),H_dim,eig_vec(1,:,:))

open(unit=102,file='./cmap_Gamma',action='write',access='sequential',position='append',status='replace')
do i_dim=1,H_dim
        write(102,*) eigen(1,i_dim)
enddo
close(102)

!print*,eigen(1,:)

contains

subroutine int2char(in_int,out_char)
        integer,intent(in)              :: in_int
        character(len=20),intent(out)    :: out_char

IF (in_int .LT. 10 )    WRITE(out_char,'(I1)')in_int
IF ((in_int .LT. 100)  .and. (in_int .GE. 10)   )WRITE(out_char,'(I2)')in_int
IF ((in_int .LT. 1000) .and. (in_int .GE. 100)  )WRITE(out_char,'(I3)')in_int
IF ((in_int .LT. 10000) .and. (in_int .GE. 1000) )WRITE(out_char,'(I4)')in_int
IF ((in_int .LT. 100000) .and. (in_int .GE. 10000) )WRITE(out_char,'(I5)')in_int

end subroutine int2char


subroutine diasym(a,eig,n,D)
 implicit none

 integer,intent(in)  :: n
 complex(kind(0d0)),intent(in)  :: a(n,n)
 complex(kind(0d0))              :: a_tmp(n,n)
 complex(kind(0d0)),intent(out)  :: D(n,n)

 integer  :: l,inf
 complex(kind(0d0))  :: Lwork(2*n-1)
 real(kind(0d0)),intent(out)      :: eig(n)
 real(kind(0d0))                 :: Rwork(3*n-2)

 a_tmp(:,:)=a(:,:)

 l=2*n-1
 call zheev('V','U',n,a_tmp,n,eig,Lwork,l,Rwork,inf)
 D(:,:)=a_tmp(:,:)



end subroutine diasym


end program

