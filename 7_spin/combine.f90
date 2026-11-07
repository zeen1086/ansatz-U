program spin
  use parameter_mod
        implicit none
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
        integer :: int3(3),super_nx,super_ny,super_nz,n_func
        integer :: convert_size_x,convert_size_y,convert_size_z
        integer :: irept,nrept,super_nrept,i_locat,j_locat
        integer :: max_x,max_y,max_z,min_x,min_y,min_z,super_i,super_j,d_ext
        integer :: nspin, nwan_sp, nwan, wan_max,wan_min,i_dir,j_dir,k_dir
        integer,allocatable     :: dir_vec(:,:),inv_dir(:)
        character(len=200)      :: useless_string,func_dir,nscf_dir,rvr_root
        character(len=200)      :: RVR_dir
        character(len=20)       :: rept_char
        complex(kind(0d0)) :: sx(3),sy(3),sz(3),cI
        complex(kind(0d0)),allocatable  :: RVR(:,:,:),RVR_POT_ave(:,:),RVR_tmp(:,:)
        complex(kind(0d0)),allocatable  :: H_sk(:,:),read_H(:,:),full_H(:,:)
        complex(kind(0d0)),allocatable  :: full_H_sp(:,:)
        real(kind(0d0)) :: q_vec(3)

func_dir = p_func_dir
nscf_dir = p_nscf_dir
rvr_root = p_rvr_dir
n_func = p_n_func
convert_size_x = p_convert_size_x
convert_size_y = p_convert_size_y
convert_size_z = p_convert_size_z
nwan = p_nwan
nspin = p_nspin
nx = p_nx
ny = p_ny
nz = p_nz
axial = p_axial
period_x = p_period_x
period_y = p_period_y
period_z = p_period_z
wan_min = p_wan_min
wan_max = p_wan_max

call load_super_size()
if  (nspin .eq. 1) then
        nwan_sp=2*(wan_max-wan_min+1)
else
        nwan_sp=(wan_max-wan_min+1)
endif

H_dim=nwan*2/nspin*period_x*period_y*period_z
allocate(H_sk(H_dim,H_dim))
allocate(full_H(H_dim,H_dim))
allocate(read_H(H_dim*nspin/2,H_dim*nspin/2))
allocate(full_H_sp(H_dim,H_dim))

do irept=1,super_nrept
        print*,irept
        full_H(:,:)=dcmplx(0.0d0,0.0d0)
        call int2char(irept,rept_char)
        open(unit=102,file=trim(p_super_h_dir)//'hr_HYC_'//trim(rept_char),action='read',access='sequential'&
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

        full_H_sp(:,:)=full_H(:,:)

        open(unit=102,file='./tot_HR/hr_'//trim(rept_char),action='write',access='stream'&
                &,status='replace',form='unformatted')
        write(102)      full_H_sp(:,:)
        close(102)

        if (all (super_ir_vec(irept,:).eq. (/0,0,0/) )) then
                open(unit=102,file='./R0',action='write',access='append',status='replace')
                do iwan=1,H_dim
                        do jwan=1,H_dim
                                write(102,*) iwan,jwan, abs(H_sk(iwan,jwan)*13.605693122994d0)
                        enddo
                enddo
                close(102)
        endif
enddo

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
        open(unit=102,file=trim(p_super_h_dir)//'ir_vec_HYC',action='read',access='sequential',status='old')
        do
        read(102,'(3I8)',iostat=io) int3(:)
        if (io/=0) exit
        irept=irept+1
        enddo
        close(102)
        super_nrept=irept
        allocate(super_ir_vec(super_nrept,3))

        irept=1
        open(unit=102,file=trim(p_super_h_dir)//'ir_vec_HYC',action='read',access='sequential',status='old')
        do
        read(102,'(3I8)',iostat=io) int3(:)
        if (io/=0) exit
                super_ir_vec(irept,:)=int3(:)
                irept=irept+1
        enddo
        close(102)

end subroutine

end program
