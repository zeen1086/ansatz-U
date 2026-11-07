program R_sym
  use parameter_mod
        ! This code creates the mapping from position index to 3D vector
        ! crystal symmetry operation can be included,
        ! a rotation symmetry making z -> -z in presented here

        ! Note: Rather than using full BZ, we just take a relavent subspace
        !       countted as ni-points

        implicit none
        integer :: ix,iy,iz,jx,jy,jz,nx,ny,nz,n_sym
        integer :: super_x, super_y,super_z,n_spinor
        integer :: n_atom,i_atom
        integer,allocatable :: map(:),sym_atom(:,:)
        integer,allocatable :: sym_BZ(:,:),vec(:,:)
        integer,allocatable :: sym_op(:,:,:)
        real,allocatable    :: frac_trans(:,:),pos(:,:)
        real    :: tmp(3)
        integer :: i,j,ni,i_sym,io
        integer :: orig(3)
        character(len=20) ::useless_string,pos_format
        real(kind(0d0)),allocatable  :: B_general_old(:,:,:,:)
        real(kind(0d0)),allocatable  :: B_general(:,:,:,:)

  n_spinor=p_n_func/3
  open(unit=101,file=p_scf_dir//'V_HYC_dp',access='sequential',action='read'&
                &,form='unformatted')
        read(101) nx,ny,nz
        print*,'nx ny nz=',nx,ny,nz
  close(101)

  ni=nx*ny*nz

  allocate(vec(ni,3))

  do ix=0,nx-1
  do iy=0,ny-1
  do iz=0,nz-1
        i=1+ix+nx*iy+nx*ny*iz
        vec(i,:)=(/ix,iy,iz/)
  enddo
  enddo
  enddo

  open(unit=101,file=p_func_dir//'R_vec',action='write'&
                &,status='replace',form='unformatted')
        write(101) ni
        write(101) vec
  close(101)

  allocate(B_general_old((nx+1)*(ny+1)*(nz+1),3,n_spinor,3))
  allocate(B_general(ni,3,n_spinor,3))

  open(unit=102,file=p_func_dir//'V_func_tmp_tmp',action='read',access='sequential',form='unformatted')
        read(102)
        read(102)
        read(102)
        read(102) B_general_old
        read(102)
  close(102)

  do ix=0,nx-1
  do iy=0,ny-1
  do iz=0,nz-1
        i=1+ix+nx*iy+nx*ny*iz
        j=1+ix+(nx+1)*iy+(nx+1)*(ny+1)*iz
        B_general(i,:,:,:)=B_general_old(j,:,:,:)
  enddo
  enddo
  enddo

open(unit=102,file=p_func_dir//'V_func_tmp',action='write',access='append',status='replace',form='unformatted')
        write(102)
        write(102)
        write(102)
        write(102) B_general
        write(102)
close(102)

contains

subroutine int2char(in_int,out_char)
        integer,intent(in)              :: in_int
        character(len=20),intent(out)    :: out_char

IF (in_int .LT. 10 )            WRITE(out_char,'(I1)')in_int
IF ((in_int .LT. 100)  .and. (in_int .GE. 10)   )   WRITE(out_char,'(I2)')in_int
IF ((in_int .LT. 1000) .and. (in_int .GE. 100)  )   WRITE(out_char,'(I3)')in_int
IF ((in_int .LT. 10000) .and. (in_int .GE. 1000) )  WRITE(out_char,'(I4)')in_int
IF ((in_int .LT. 100000) .and. (in_int .GE. 10000) )WRITE(out_char,'(I5)')in_int

end subroutine int2char

end program
