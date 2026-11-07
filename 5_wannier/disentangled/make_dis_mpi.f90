program make_wan_func
  use parameter_mod
        implicit none
        integer         :: ini_bnd,fin_bnd,ibnd,jbnd,nbnd,is
        integer         :: nx,ny,nz,nk,ik,jk,i,j,io,nk1,nk2,nk3,z_cut,ix,iy,iz,ir
        integer         :: wan_bnd,read_int3(3),nwan,nspin
        complex(kind(0d0)),allocatable :: ini_wfn(:,:),fin_wfn(:,:,:),read_wfn(:,:)
        complex(kind(0d0)),allocatable :: psi_in(:,:),psi_out(:,:),psi_mid(:,:),unmk(:,:,:)
        complex(kind(0d0)),allocatable :: unmk_dis(:,:,:),U(:,:)
        complex(kind(0d0)),allocatable :: unmk_dis_tmp(:)
        complex(kind(0d0)),allocatable :: dis_fun(:,:,:),to_write_wan(:,:,:)
        complex(kind(0d0)),allocatable :: phase(:,:)
        real(kind(0d0)),allocatable :: to_write(:)
        real(kind(0d0)),allocatable ::  k_list(:,:)
        real(kind(0d0)) :: read3(3),r(3)
        complex(kind(0d0)) :: cI
        real(kind(0d0)):: pi
        character(len=20) :: k_char,bnd_char,wan_char
        character(len=200):: file_name
        integer,allocatable:: cell_grid(:,:)
        integer,allocatable :: dir_vec(:,:)
        integer :: convert_size,i_dir,j_dir,k_dir,iwan,bnd_cut
        real :: start, finish,start2, finish2
        complex(kind(0d0)),allocatable ::wan_fun_local(:,:,:)
        integer,allocatable     :: pool(:,:)
        integer                 :: muti_pool,imuti_pool,ierr,my_rank,size
        integer                 :: npool ,ipool

        integer,allocatable     :: include_TF(:,:)
        integer                 :: i_count
        integer,allocatable     :: win_min(:),win_max(:)

include 'mpif.h'

cI=dcmplx(0.0d0,1.0d0)
pi=acos(-1d0)

nk1=p_nk1
nk2=p_nk2
nk3=p_nk3
nk=nk1*nk2*nk3
nbnd=p_num_bands
nspin=p_nspin
bnd_cut=p_first_band-1
allocate(win_min(nk))
allocate(win_max(nk))

open(unit=101,file=p_nscf_dir//'u_mat_HYC',action='read',access='sequential',form='unformatted')
read(101) read_int3
if (any(read_int3 /= [p_nwan,p_nwan,nk])) error stop 'u_mat_HYC dimensions disagree with parameter'
allocate(unmk(read_int3(1),read_int3(2),read_int3(3)))
read(101) unmk
close(101)

open(unit=101,file=p_nscf_dir//'u_dis_mat_HYC',action='read',access='sequential',form='unformatted')
read(101) read_int3
if (any(read_int3 /= [p_num_bands,p_nwan,nk])) error stop 'u_dis_mat_HYC dimensions disagree with parameter'
nbnd=read_int3(1)
nwan=read_int3(2)
if (.not. nk .eq. read_int3(3)) then
        print*, 'nk not correct'
        stop
endif
allocate(unmk_dis(read_int3(1),read_int3(2),read_int3(3)))
allocate(U(read_int3(1),read_int3(2)))
allocate(unmk_dis_tmp(read_int3(2)))
read(101) unmk_dis  !  nbnd, nwan, nk
close(101)

allocate(include_TF(nbnd,nk))

include_TF(:,:)=0

open(unit=101,file=p_nscf_dir//'window_TF_HYC',action='read',access='sequential',form='unformatted')
read(101) include_TF
close(101)

win_min(:)=999
win_max(:)=0

do ik=1,nk
do i=1,nbnd-1
        if ((include_TF(i,ik)   .eq. 1).and. ( include_TF(i+1,ik) .eq. 0)) win_max(ik)=i
        if ((include_TF(i+1,ik) .eq. 1).and. ( include_TF(i  ,ik) .eq. 0)) win_min(ik)=i+1
enddo
        if (include_TF(1,ik) .eq. 1) win_min(ik)=1
        if (include_TF(nbnd,ik) .eq. 1) win_max(ik)=nbnd
        if (win_min(ik) .eq. 999) then
                print*, 'win min not found for k=',ik
                stop
        endif
        if (win_max(ik) .eq. 0) then
                print*, 'win max not found for k=',ik
                stop
        endif
enddo

file_name=p_nscf_dir//'V_HYC_dp'
open(unit=101,file=trim(file_name),access='sequential',form='unformatted')
read(101) nx,ny,nz
print*,nx,ny,nz
close(101)

allocate(read_wfn(nx*ny*nz,nspin))
allocate(dis_fun(nwan,nx*ny*nz,nspin))

npool=nk
allocate(pool(npool,1))
ipool=1
do ik=1,nk
        pool(ipool,1)=ik
        ipool=ipool+1
enddo

call mpi_init(ierr)
call mpi_comm_rank(MPI_COMM_WORLD,my_rank,ierr)
call mpi_comm_size(MPI_COMM_WORLD,size,ierr)
if (size < 2) call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)

muti_pool=0
DO
        muti_pool=muti_pool+1
        IF(muti_pool*(size-1) .gt. npool) EXIT
ENDDO
muti_pool=muti_pool-1

if (my_rank.eq.0) call EXECUTE_COMMAND_LINE("mkdir -p WFN")
call MPI_BARRIER(MPI_COMM_WORLD,ierr)

if (.not. my_rank.eq.0)  then
        DO imuti_pool=0,muti_pool
                ipool=imuti_pool*(size-1)+my_rank
                if (ipool .gt. npool) cycle
                ik=pool(ipool,1)
                call int2char(ik,k_char)
                dis_fun(:,:,:)=dcmplx(0.0d0,0.0d0)

                do ibnd=1,nbnd
                do iwan=1,nwan

                        U(ibnd,iwan)=sum((unmk(:,iwan,ik))*unmk_dis(ibnd,:,ik))
                enddo
                enddo

                do ibnd=win_min(ik),win_max(ik)
                        print*,ik,ibnd
                        call int2char(ibnd+bnd_cut,bnd_char)

                        file_name=p_nscf_dir//'WFN/wfn_k'//trim(k_char)//'_b'//trim(bnd_char)//'_vec_dp'
                        open(unit=101,file=trim(file_name)&
                                &,action='read',access='sequential',form='unformatted')
                                read(101) read_wfn
                        close(101)
                        unmk_dis_tmp(:)=U(1+ibnd-win_min(ik),:)

                        Do i=1,nwan
                                dis_fun(i,:,:)=dis_fun(i,:,:)+unmk_dis_tmp(i)*read_wfn(:,:)
                        enddo

                enddo

                do i=1,nwan
                        call int2char(i,bnd_char)
                        file_name='./WFN/wfn_k'//trim(k_char)//'_w'//trim(bnd_char)//'_vec_dp'
                        open(unit=101,file=trim(file_name)&
                                &,action='write',access='sequential',form='unformatted',status='replace')
                        write(101) dis_fun(i,:,:)

                        close(101)

                enddo

        enddo

endif !my_rank
if (my_rank.eq.0)  then
endif

call mpi_finalize(ierr)

contains

subroutine int2char(in_int,out_char)
        integer,intent(in)              :: in_int
        character(len=20),intent(out)    :: out_char

IF (in_int .LT. 10 )   WRITE(out_char,'(I1)')in_int
IF ((in_int .LT. 100)  .and. (in_int .GE. 10)   )    WRITE(out_char,'(I2)')in_int
IF ((in_int .LT. 1000) .and. (in_int .GE. 100)  )    WRITE(out_char,'(I3)')in_int
IF ((in_int .LT. 10000) .and. (in_int .GE. 1000) )   WRITE(out_char,'(I4)')in_int
IF ((in_int .LT. 100000) .and. (in_int .GE. 10000) ) WRITE(out_char,'(I5)')in_int

end subroutine int2char

end program
