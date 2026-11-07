program spin
        implicit none
        character(len=1024) :: the_input

        integer :: ix,iy,iz,nx,ny,nz,i,is,j,io,js
        integer :: period_x,period_y,period_z,jx,jy,jz
        integer :: ini_x,ini_y,ini_z,fin_x,fin_y,fin_z
        integer :: iRVR,jRVR,H_dim,i_dim,j_dim,iwan,jwan
        integer :: super_x,super_y,super_z,jrept
        integer,allocatable :: cell(:,:),super_ir_vec(:,:)
        real(kind(0d0)),allocatable        :: x_list(:,:,:),coefficient(:)
        real(kind(0d0)),allocatable        :: spin_vec(:,:,:),off_set(:,:)
        real(kind(0d0))         :: A2B(3),r_tmp(3),Mz
        real(kind(0d0))    :: rad,pi,x,y_com,x_com,theta
        real(kind(0d0))    :: a(3,3),a0 ,avB(3),avA(3),rrBdr(3),rrAdr(3)     
        real(kind(0d0)) :: a_next,b_next,s_abs
        real(kind(0d0)) ::vec6(6,3),vec3(3),axial
        integer :: int3(3),super_nx,super_ny,super_nz,n_func,convert_size_x,convert_size_y,convert_size_z
        integer :: irept,nrept,super_nrept,i_locat,j_locat
        integer :: max_x,max_y,max_z,min_x,min_y,min_z,super_i,super_j,d_ext
        integer :: nspin, nwan_sp, nwan, wan_max,wan_min,i_dir,j_dir,k_dir
        integer,allocatable     :: dir_vec(:,:),inv_dir(:)
        character(len=200)      :: useless_string,func_dir,nscf_dir,THE_dir
        character(len=200)      :: RVR_dir
        character(len=20)       :: rept_char
        complex(kind(0d0)) :: sx(3),sy(3),sz(3),cI
        complex(kind(0d0)),allocatable  :: RVR(:,:,:),RVR_POT_ave(:,:),RVR_tmp(:,:)
        complex(kind(0d0)),allocatable  :: H_sk(:,:),read_H(:,:),full_H(:,:)
        complex(kind(0d0)),allocatable  :: read_AAR(:,:,:,:),full_AAR(:,:,:,:)
        complex(kind(0d0)),allocatable  :: full_H_sp(:,:),full_AAR_sp(:,:,:,:)

call get_command_argument(1,the_input)
if (len_trim(the_input)==0) error stop 'Supply explicit THE.in path as argument 1'
open(unit=101,file=trim(the_input),action='read',access='sequential',status='old')
read(101,'(A9,A200)')useless_string,func_dir
read(101,'(A9,A200)')useless_string,nscf_dir
read(101,'(A8,A200)')useless_string,THE_dir
read(101,'(A7,I3)')  useless_string,n_func
read(101,'(A15,I2)') useless_string,convert_size_x
read(101,'(A15,I2)') useless_string,convert_size_y
read(101,'(A15,I2)') useless_string,convert_size_z
read(101,'(A5,I3)') useless_string, nwan
read(101,'(A6,I2)') useless_string, nspin
read(101,'(A3,I3)') useless_string, nx
read(101,'(A3,I3)') useless_string, ny
read(101,'(A3,I3)') useless_string, nz
read(101,*)!'(A6,F8.4)')  useless_string,axial
read(101,'(A9,I3)') useless_string,period_x
read(101,'(A9,I3)') useless_string, period_y
read(101,'(A9,I3)') useless_string, period_z
read(101,*)
read(101,*)
read(101,*)
read(101,*)
read(101,*)
read(101,'(A8,I3)') useless_string, wan_min
read(101,'(A8,I3)') useless_string, wan_max
close(101)

call load_super_size()
call execute_command_line('mkdir -p tot_HR')
if  (nspin .eq. 1) then
        nwan_sp=2*(wan_max-wan_min+1)
else
        nwan_sp=(wan_max-wan_min+1)
endif

H_dim=nwan*2/nspin*period_x*period_y*period_z
allocate(H_sk(H_dim,H_dim))
allocate(full_H(H_dim,H_dim))
allocate(read_H(H_dim*nspin/2,H_dim*nspin/2))
allocate(read_AAR(H_dim*nspin/2,H_dim*nspin/2,super_nrept,3))
allocate(full_AAR(H_dim,H_dim,super_nrept,3))
allocate(full_H_sp(H_dim,H_dim))
allocate(full_AAR_sp(H_dim,H_dim,super_nrept,2))

do irept=1,super_nrept
        print*,irept
        full_H(:,:)=dcmplx(0.0d0,0.0d0)
        call int2char(irept,rept_char)
        open(unit=102,file=trim(THE_dir)//'6_RVR/super_H/hr_HYC_'//trim(rept_char),action='read',access='sequential'&
                &,status='old',form='unformatted')
        read(102)   read_H(:,:)
        close(102)
        open(unit=102,file='./data/Hsk_'//trim(rept_char),action='read',access='sequential'&
                &,status='old',form='unformatted')
        read(102)   H_sk(:,:)
        close(102)
        
        do iwan=1,H_dim*nspin/2
        do jwan=1,H_dim*nspin/2
        do is=1,2/nspin
                i_dim=is+2/nspin*(iwan-1)
                j_dim=is+2/nspin*(jwan-1)
                full_H(i_dim,j_dim)=read_H(iwan,jwan)
        enddo
        enddo
        enddo
        full_H(:,:)=full_H(:,:)+H_sk(:,:)*13.605693122994d0

        full_H_sp(:,:)=full_H(:,:) ! complex128 stream, version 2



        open(unit=102,file='./tot_HR/hr_'//trim(rept_char),action='write',access='stream'&
                &,status='replace',form='unformatted')
        write(102)      full_H_sp(:,:)
        close(102)
        !if (.not. mod(irept,10) .eq. 8) cycle
       ! print*,irept,sum(abs(H_sk(:,:,irept))**2),super_ir_vec(irept,:)
enddo


!open(unit=102,file=trim(THE_dir)//'6_RVR/super_H/AAR',action='read',access='sequential'&
!                &,status='old',form='unformatted')
!        read(102)   read_AAR(:,:,:,:)
!close(102)


!full_AAR(:,:,:,:)=dcmplx(0.0d0,0.0d0)
!do iwan=1,H_dim*nspin/2
!do jwan=1,H_dim*nspin/2
!do is=1,2/nspin
!        i_dim=is+2/nspin*(iwan-1)
!        j_dim=is+2/nspin*(jwan-1)
!        full_AAR(i_dim,j_dim,:,:)=read_AAR(iwan,jwan,:,:)
!enddo
!enddo
!enddo

!full_AAR_sp(:,:,:,1:2)=full_AAR(:,:,:,1:2)


!do irept=1,super_nrept
!        call int2char(irept,rept_char)

!        open(unit=102,file='./tot_HR_2/AAR_'//trim(rept_char),action='write',access='stream'&
!                &,status='replace',form='unformatted')
!                write(102)      full_AAR_sp(:,:,irept,1:2)
!        close(102)
!enddo



contains

subroutine int2char(in_int,out_char)
        integer,intent(in)              :: in_int
        character(len=20),intent(out)    :: out_char

IF (in_int .LT. 10 )   WRITE(out_char,'(I1)')in_int
IF ((in_int .LT. 100)  .and. (in_int .GE. 10)   ) WRITE(out_char,'(I2)')in_int
IF ((in_int .LT. 1000) .and. (in_int .GE. 100)  ) WRITE(out_char,'(I3)')in_int
IF ((in_int .LT. 10000) .and. (in_int .GE. 1000) ) WRITE(out_char,'(I4)')in_int
IF ((in_int .LT. 100000) .and. (in_int .GE. 10000) )WRITE(out_char,'(I5)')in_int

end subroutine int2char




subroutine load_super_size()
        irept=0
        open(unit=102,file=trim(THE_dir)//'6_RVR/super_H/ir_vec_HYC',action='read',access='sequential',status='old')
        do
        read(102,'(3I8)',iostat=io) int3(:)
        if (io/=0) exit
        irept=irept+1
        enddo
        close(102)
        super_nrept=irept
        allocate(super_ir_vec(super_nrept,3))

        irept=1
        open(unit=102,file=trim(THE_dir)//'6_RVR/super_H/ir_vec_HYC',action='read',access='sequential',status='old')
        do
        read(102,'(3I8)',iostat=io) int3(:)
        if (io/=0) exit
                super_ir_vec(irept,:)=int3(:)
                irept=irept+1
        enddo
        close(102)






end subroutine




end program
