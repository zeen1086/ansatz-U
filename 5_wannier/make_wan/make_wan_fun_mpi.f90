program make_wan_func
  use parameter_mod
        implicit none
        integer         :: ini_bnd,fin_bnd,ibnd,jbnd,nbnd,is
        integer         :: nx,ny,nz,nk,ik,jk,i,j,io,nk1,nk2,nk3,z_cut,ix,iy,iz,ir
        integer         :: nspin,i1,i2
        integer         :: wan_bnd,read_int3(3),nwan
        complex(kind(0d0)),allocatable :: ini_wfn(:,:),fin_wfn(:,:,:),read_wfn(:,:)
        complex(kind(0d0)),allocatable :: unmk(:,:,:)
        complex(kind(0d0)),allocatable :: unmk_dis(:,:,:)
        complex(kind(0d0)),allocatable :: wan_fun(:,:,:)!,to_write_wan(:,:,:)
        complex(kind(0d0)),allocatable :: wan_full_grid(:,:)
        complex(kind(0d0)),allocatable :: w_dx(:,:,:),w_dz(:,:,:),w_dy(:,:,:)
        complex(kind(0d0)),allocatable :: phase(:,:)
        real(kind(0d0)),allocatable ::  k_list(:,:)
        real(kind(0d0)) :: read3(3),r(3)
        complex(kind(0d0)) :: cI
        real(kind(0d0)):: pi
        character(len=20) :: k_char,bnd_char,wan_char
        character(len=200):: file_name
        integer,allocatable:: cell_grid(:,:)
        integer,allocatable :: dir_vec(:,:)
        integer,allocatable :: full_grid(:,:)
        integer :: i_dir,j_dir,k_dir,iwan
        integer :: convert_size_x,convert_size_y,convert_size_z
        real :: start, finish,start2, finish2
        complex(kind(0d0)),allocatable ::wan_fun_local(:,:,:)
        complex(kind(0d0))      :: ctmp,d3(3,2)
        integer,allocatable     :: pool(:,:)
        integer                 :: muti_pool,imuti_pool,ierr,my_rank,size
        integer                 :: npool ,ipool

        integer,allocatable     :: include_TF(:,:)
        integer                 :: i_count
        integer,allocatable     :: win_min(:),win_max(:),R_shift(:,:)

include 'mpif.h'
call read_row_index(iwan)
if (iwan > p_nwan) error stop 'Wannier row exceeds nwan'

cI=dcmplx(0.0d0,1.0d0)
pi=acos(-1d0)

nk1=p_nk1
nk2=p_nk2
nk3=p_nk3
nk=nk1*nk2*nk3
convert_size_x=p_convert_size_x
convert_size_y=p_convert_size_y
convert_size_z=p_convert_size_z
nspin=p_nspin
if (2*convert_size_x+1 .gt. nk1) then
print*, "too large x  wannier size "
stop
endif
if (2*convert_size_y+1 .gt. nk2) then
print*, "too large y  wannier size "
stop
endif
if (2*convert_size_z+1 .gt. nk3) then
print*, "too large z  wannier size "
stop
endif

allocate(dir_vec((2*convert_size_x+1)*(2*convert_size_y+1)*(2*convert_size_z+1),3))
i=1
do i_dir=-convert_size_x,convert_size_x
do j_dir=-convert_size_y,convert_size_y
do k_dir=-convert_size_z,convert_size_z
        dir_vec(i,:)=(/i_dir,j_dir,k_dir/)
        i=i+1
enddo
enddo
enddo

allocate(k_list(nk,3))
allocate(cell_grid(nk,3))
allocate(win_min(nk))
allocate(win_max(nk))
do iz=1,nk3
do iy=1,nk2
do ix=1,nk1
        i=1+(ix-1)+nk1*(iy-1)+nk1*nk2*(iz-1)
        cell_grid(i,:)=(/ix,iy,iz/)
enddo
enddo
enddo

ik=0
open(unit=101,file=p_nscf_dir//'k_list',action='read',access='sequential',status='old')
do
        read(101,'(3F12.8)',iostat=io) read3(:)
        if (io /=0) exit
        if (ik >= nk) error stop 'too many k_list rows'
        ik=ik+1
        k_list(ik,:)=read3(:)
enddo
close(101)
if (ik /= nk) error stop 'k_list row count disagrees with parameter'

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
read(101) unmk_dis  !  nbnd, nwan, nk
close(101)

allocate(R_shift(nwan,3))

R_shift(:,:)=0
ir=1

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

allocate(full_grid(nx*ny*nz*(2*convert_size_x+1)*(2*convert_size_y+1)*(2*convert_size_z+1),4))

do j=1,(2*convert_size_x+1)*(2*convert_size_y+1)*(2*convert_size_z+1)
do ix=0,nx-1
do iy=0,ny-1
do iz=0,nz-1
        i=1+ix+(dir_vec(j,1)+convert_size_x)*nx+&
             &(2*convert_size_x+1)*nx*(iy+(dir_vec(j,2)+convert_size_y)*ny)+&
             &(2*convert_size_x+1)*nx*(2*convert_size_y+1)*ny*(iz+(dir_vec(j,3)+convert_size_z)*nz)
        full_grid(i,:)=(/ix,iy,iz,j/)

enddo
enddo
enddo
enddo

allocate(read_wfn(nx*ny*nz,nspin))
allocate(ini_wfn(nx*ny*nz,nspin))!1+nx*ny*z_cut:nx*ny*(nz-z_cut)
allocate(wan_fun((2*convert_size_x+1)*(2*convert_size_y+1)*(2*convert_size_z+1),nx*ny*nz,nspin))

allocate(wan_fun_local((2*convert_size_x+1)*(2*convert_size_y+1)*(2*convert_size_z+1),nx*ny*nz,nspin))
allocate(phase(nk,nx*ny*nz))

do ix=0,nx-1
do iy=0,ny-1
do iz=0,nz-1
        r(:)=(/1.0d0*ix/nx,1.0d0*iy/ny,1.0d0*iz/nz/)
        ir=1+(ix)+nx*(iy)+nx*ny*(iz)
        do ik=1,nk
                phase(ik,ir)=exp(cI*2.0d0*pi*sum(r(:)*k_list(ik,:)))
        enddo
enddo
enddo
enddo

ipool=0
do ik=1,nk
        ipool=ipool+1
enddo

npool=ipool
allocate(pool(npool,1))
ipool=1
do ik=1,nk
        pool(ipool,1)=ik!(/ik,ibnd/)
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

if ( my_rank.eq.1)  then
open(unit=103,file='./report',action='write',access='append',status='replace')
write(103,*) 'start calculation '
write(103,*) 'iwan=', iwan
close(103)
endif

if (.not. my_rank.eq.0)  then
        wan_fun_local(:,:,:)=dcmplx(0.0d0,0.0d0)
        DO imuti_pool=0,muti_pool
                ipool=imuti_pool*(size-1)+my_rank
                if (ipool .gt. npool) cycle
                ik=pool(ipool,1)

                if ( my_rank.eq.1 .and. mod(imuti_pool,10) .eq.0 )  then
                open(unit=103,file='./report',action='write',access='append')
                write(103,*) 'imuti_pool / muti_pool' , imuti_pool,' /',muti_pool
                write(103,*) 'current ik, ibnd=' , ik,ibnd
                close(103)
                endif

                k_char=row_text(ik)
                bnd_char=row_text(iwan)
                file_name=p_nscf_dir//'disentangled/WFN/wfn_k'//trim(k_char)//'_w'//trim(bnd_char)//'_vec_dp'
                open(unit=101,file=trim(file_name)&
                    &,action='read',access='sequential',form='unformatted',status='old')
                read(101) read_wfn
                close(101)

                do is=1,nspin
                        ini_wfn(:,is)=read_wfn(:,is)*phase(ik,:)
                enddo
                ini_wfn(:,:)=exp(-cI*2*pi*sum(R_shift(iwan,:)*k_list(ik,:)))*&
                                        &ini_wfn(:,:)

                do j=1,(2*convert_size_x+1)*(2*convert_size_y+1)*(2*convert_size_z+1)
                ctmp=exp(cI*2*pi*sum(dir_vec(j,:)*k_list(ik,:)))
                wan_fun_local(j,:,:)=wan_fun_local(j,:,:)+&
                        &ctmp*ini_wfn(:,:)  !ini_wfn(:,is)
                enddo
        enddo
        call MPI_SEND(wan_fun_local(:,:,:),(2*convert_size_x+1)*(2*convert_size_y+1)*(2*convert_size_z+1)*nx*ny*nz*nspin&
                                &,MPI_DOUBLE_COMPLEX,0,my_rank,MPI_COMM_WORLD,ierr)

endif !my_rank
if (my_rank.eq.0)  then
        wan_fun(:,:,:)=dcmplx(0.0d0,0.0d0)
        do ipool=1,size-1
                call MPI_RECV(wan_fun_local(:,:,:),(2*convert_size_x+1)*(2*convert_size_y+1)*(2*convert_size_z+1)*nx*ny*nz*nspin,MPI_DOUBLE_COMPLEX,&
                    &MPI_ANY_SOURCE,ipool,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
                wan_fun(:,:,:)=wan_fun(:,:,:)+wan_fun_local(:,:,:)
        enddo

        wan_fun(:,:,:)=wan_fun(:,:,:)/sqrt(dble(nx)*dble(ny)*dble(nz))/dble(nk)

        wan_char=row_text(iwan)
        open(unit=101,file='./wan_'//trim(wan_char),&
                        &action='write',access='append',status='replace',form='unformatted')
                write(101) wan_fun
        close(101)
        print*,sum(abs(wan_fun)**2)

endif

call mpi_finalize(ierr)

end program
