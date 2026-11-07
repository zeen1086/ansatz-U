program make_R_vec
  use parameter_mod
  use iso_fortran_env, only: int64
  implicit none
  integer :: nx,ny,nz,ni,ne,u,v,i,j,x,y,z,nspinor
  integer,allocatable :: vec(:,:),take(:)
  integer(int64) :: model_bytes
  real(8),allocatable :: ave(:,:),v0(:,:),vn(:,:),b(:,:,:,:),jh(:,:)
  character(1024) :: grid,input,output,rpath
  call get_command_argument(1,grid)
  call get_command_argument(2,input)
  call get_command_argument(3,output)
  call get_command_argument(4,rpath)
  if(min(len_trim(grid),len_trim(input),len_trim(output),len_trim(rpath))==0) &
    error stop 'Usage: make_R_vec base_V_HYC_dp fitted_extended_model output_model output_R_vec'
  open(newunit=u,file=trim(grid),status='old',action='read',form='unformatted')
  read(u) nx,ny,nz
  close(u)
  if(any([nx,ny,nz]<1)) error stop 'Invalid base grid'
  ni=nx*ny*nz;ne=(nx+1)*(ny+1)*(nz+1);nspinor=p_n_func/3
  inquire(file=trim(input),size=model_bytes)
  if(model_bytes/=40_int64+8_int64*ne*(15_int64+10_int64*nspinor)) &
    error stop 'Fitted model grid/channel count or float64 five-record format mismatch'
  allocate(vec(ni,3),take(ni),ave(ne,4),v0(ne,3),vn(ne,8),b(ne,3,nspinor,3),jh(ne,nspinor))
  open(newunit=u,file=trim(input),status='old',action='read',form='unformatted')
  read(u) ave
  read(u) v0
  read(u) vn
  read(u) b
  read(u) jh
  close(u)
  do z=0,nz-1
  do y=0,ny-1
  do x=0,nx-1
    i=1+x+nx*(y+ny*z)
    vec(i,:)=[x,y,z]
    take(i)=1+x+(nx+1)*(y+(ny+1)*z)
  enddo
  enddo
  enddo
  open(newunit=u,file=trim(rpath),status='new',action='write',form='unformatted')
  write(u) ni
  write(u) vec
  close(u)
  open(newunit=u,file=trim(output),status='new',action='write',form='unformatted')
  write(u) ave(take,:)
  write(u) v0(take,:)
  write(u) vn(take,:)
  write(u) b(take,:,:,:)
  write(u) jh(take,:)
  close(u)
end program
