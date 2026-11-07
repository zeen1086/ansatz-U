program extract
  use parameter_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  integer :: nx,ny,nz,sx,sy,sz,natom,ncell,ndep,u,v,ios,i,c,j,ix,iy,iz,x,y,z,local,global
  integer,allocatable :: mapping(:,:)
  real(8),allocatable :: spins(:,:),field(:,:),part(:,:),coeff(:)
  real(8) :: amp,theta,phi,pi
  character(1024) :: map_path,line
  character(32) :: number
  call get_command_argument(1,map_path)
  if(len_trim(map_path)==0) error stop 'Usage: V_split explicit-channel-map.dat'
  sx=p_split_super_x;sy=p_split_super_y;sz=p_split_super_z;ndep=p_split_nspinor
  ncell=sx*sy*sz
  open(newunit=u,file=trim(map_path),status='old',action='read')
  read(u,*,iostat=ios) natom,j,c
  if(ios/=0) error stop 'channel map header: number_of_spins number_of_channels number_of_cells'
  if(natom<1.or.j/=ndep.or.c/=ncell) error stop 'channel map dimensions disagree with parameter'
  allocate(mapping(ndep,ncell),spins(natom,3),coeff(3*ndep))
  do c=1,ncell
    read(u,*,iostat=ios) mapping(:,c)
    if(ios/=0) error stop 'channel map requires one row per split cell, x,y,z order (z fastest)'
  enddo
  if(any(mapping<1).or.any(mapping>natom)) error stop 'channel map spin index out of range'
  read(u,'(a)',iostat=ios) line
  if(ios==0) error stop 'Unexpected extra channel map rows'
  close(u)
  pi=acos(-1d0)
  open(newunit=u,file='polar2',status='old',action='read')
  do i=1,natom
    read(u,*,iostat=ios) amp,theta,phi
    if(ios/=0) error stop 'polar2 requires amplitude theta_deg phi_deg for every mapped spin'
    if(.not.all(ieee_is_finite([amp,theta,phi])).or.amp<0) error stop 'Invalid polar2 values'
    theta=theta*pi/180;phi=phi*pi/180
    spins(i,:)=amp*[sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)]
  enddo
  close(u)
  open(newunit=u,file='V_HYC_dp',status='old',form='unformatted',action='read')
  read(u) nx,ny,nz
  if(any([nx,ny,nz]<=0).or.mod(nx,sx)/=0.or.mod(ny,sy)/=0.or.mod(nz,sz)/=0) &
    error stop 'Invalid or indivisible split grid'
  allocate(field(nx*ny*nz,4),part((nx/sx+1)*(ny/sy+1)*(nz/sz+1),4))
  read(u) field
  close(u)
  if(.not.all(ieee_is_finite(field))) error stop 'Nonfinite potential input'
  open(newunit=v,file='spin_info',status='new',action='write')
  c=0
  do x=0,sx-1
  do y=0,sy-1
  do z=0,sz-1
    c=c+1
    do iz=0,nz/sz
    do iy=0,ny/sy
    do ix=0,nx/sx
      local=1+ix+(nx/sx+1)*(iy+(ny/sy+1)*iz)
      global=1+modulo(x*nx/sx+ix,nx)+nx*(modulo(y*ny/sy+iy,ny)+ny*modulo(z*nz/sz+iz,nz))
      part(local,:)=field(global,:)
    enddo
    enddo
    enddo
    write(number,'(i0)') c
    open(newunit=u,file='V'//trim(number),status='new',form='unformatted',action='write')
    write(u) size(part,1)
    write(u) part
    close(u)
    call execute_command_line('mkdir -p spin'//trim(number),exitstat=ios)
    if(ios/=0) error stop 'Cannot create fit label directory'
    open(newunit=u,file='spin'//trim(number)//'/V_HYC_dp',status='new',form='unformatted',action='write')
    write(u) nx/sx+1,ny/sy+1,nz/sz+1
    write(u) part
    close(u)
    do j=1,ndep
      coeff(3*j-2:3*j)=spins(mapping(j,c),:)
    enddo
    write(v,'(*(es26.16e3,1x))') coeff
  enddo
  enddo
  enddo
  close(v)
end program
