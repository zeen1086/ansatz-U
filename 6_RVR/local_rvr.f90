module local_rvr
 use local_vfunc
 implicit none
 character(4096)::model_path='',wan_dir='',translations_path='',rvr_dir='',kz_path=''
 character(64)::basis_id=''
 integer::nwan=30
 real(dp)::grid_weight=-1
 namelist /local_projection/ model_path,wan_dir,translations_path,rvr_dir,kz_path,basis_id,nwan,grid_weight
 type::projection_data
  type(local_potential)::m
  integer::ng,nc,nr,ns,nsxy,ndiff,nd,nk
  integer,allocatable::source(:,:),xy(:,:),difference(:,:),domain(:,:),source_xy(:,:)
  real(dp),allocatable::weights(:),kz(:),means(:),fields(:,:,:)
  integer,allocatable::mapping(:,:),overlap_map(:,:)
 end type
 type::local_projected
  integer::np,nc,nr,nw,nk,ndiff
  character(64)::model_id,basis_id
  real(dp)::cell(3,3),alpha
  integer,allocatable::xy(:,:),difference(:,:)
  real(dp),allocatable::kz(:),means(:)
  integer(int64)::data_position,block_bytes
 end type
contains
 subroutine control(path,t,loadfields)
 character(*),intent(in)::path
 type(projection_data),intent(out)::t
 logical,intent(in)::loadfields
 integer::u,ios,i,j,x,y,r,lo(2),hi(2),a,b
 character(1024)::msg
 open(newunit=u,file=path,status='old',action='read',iostat=ios)
 call require(ios==0,'Cannot open local_projection namelist')
 read(u,nml=local_projection,iostat=ios,iomsg=msg);close(u)
 call require(ios==0,'Invalid local_projection: '//trim(msg))
 call require(min(len_trim(model_path),len_trim(wan_dir),len_trim(translations_path), &
 len_trim(rvr_dir),len_trim(kz_path))>0.and.len_trim(basis_id)==64,'Explicit paths and 64-character basis_id required')
 call require(nwan>=1.and.nwan<=30.and.grid_weight>0.and.ieee_is_finite(grid_weight),'Invalid nwan/grid_weight')
 call read_vfunc(trim(model_path),t%m,.false.)
 t%ng=product(t%m%grid);t%nc=2*t%m%np+1
 open(newunit=u,file=trim(translations_path),status='old',action='read',iostat=ios)
 call require(ios==0,'Missing explicit weighted Wannier translations')
 read(u,*,iostat=ios) t%ns;call require(ios==0.and.t%ns>0.and.t%ns<=100000,'Invalid Wannier support count')
 allocate(t%source(3,t%ns),t%weights(t%ns))
 do i=1,t%ns
  read(u,*,iostat=ios) t%source(:,i),t%weights(i)
  call require(ios==0,'Every Wannier translation needs x y z weight; implicit edge weights are forbidden')
  call require(t%weights(i)>0.and.ieee_is_finite(t%weights(i)),'Invalid Wannier weight')
  do j=1,i-1
   call require(any(t%source(:,i)/=t%source(:,j)),'Duplicate source Wannier translation')
  enddo
 enddo
 close(u)
 lo=minval(t%source(:2,:),dim=2);hi=maxval(t%source(:2,:),dim=2)
 t%nsxy=product(hi-lo+1);allocate(t%source_xy(2,t%nsxy));r=0
 do x=lo(1),hi(1)
 do y=lo(2),hi(2)
  r=r+1;t%source_xy(:,r)=[x,y]
 enddo
 enddo
 lo=floor((t%m%center(:2)-t%m%cutoff)/[t%m%cell(1,1),t%m%cell(2,2)])
 hi=floor((t%m%center(:2)+t%m%cutoff)/[t%m%cell(1,1),t%m%cell(2,2)])
 t%nd=product(hi-lo+1);allocate(t%domain(2,t%nd));r=0
 do x=lo(1),hi(1)
 do y=lo(2),hi(2)
  r=r+1;t%domain(:,r)=[x,y]
 enddo
 enddo
 lo=minval(t%domain,dim=2)-maxval(t%source_xy,dim=2)
 hi=maxval(t%domain,dim=2)-minval(t%source_xy,dim=2)
 t%nr=product(hi-lo+1);allocate(t%xy(2,t%nr));r=0
 do x=lo(1),hi(1)
 do y=lo(2),hi(2)
  r=r+1;t%xy(:,r)=[x,y]
 enddo
 enddo
 allocate(t%mapping(t%nsxy,t%nd))
 do i=1,t%nd
 do j=1,t%nsxy
  do r=1,t%nr
   if(all(t%xy(:,r)==t%domain(:,i)-t%source_xy(:,j))) exit
  enddo
  call require(r<=t%nr,'Projection translation map is incomplete');t%mapping(j,i)=r
 enddo
 enddo
 hi=maxval(t%source_xy,dim=2)-minval(t%source_xy,dim=2);lo=-hi
 t%ndiff=product(hi-lo+1);allocate(t%difference(2,t%ndiff));r=0
 do x=lo(1),hi(1)
 do y=lo(2),hi(2)
  r=r+1;t%difference(:,r)=[x,y]
 enddo
 enddo
 allocate(t%overlap_map(t%nsxy,t%nsxy))
 do a=1,t%nsxy
 do b=1,t%nsxy
  do r=1,t%ndiff
   if(all(t%difference(:,r)==t%source_xy(:,a)-t%source_xy(:,b))) exit
  enddo
  call require(r<=t%ndiff,'Overlap difference map is incomplete');t%overlap_map(a,b)=r
 enddo
 enddo
 open(newunit=u,file=trim(kz_path),status='old',action='read',iostat=ios)
 call require(ios==0,'Missing kz points');read(u,*,iostat=ios) t%nk
 call require(ios==0.and.t%nk>0.and.t%nk<=10000,'Invalid kz count');allocate(t%kz(t%nk))
 do i=1,t%nk
  read(u,*,iostat=ios) t%kz(i);call require(ios==0.and.ieee_is_finite(t%kz(i)),'Invalid kz point')
  do j=1,i-1
   call require(abs(t%kz(i)-t%kz(j)-anint(t%kz(i)-t%kz(j)))>1e-12_dp,'Duplicate periodic kz')
  enddo
 enddo
 close(u)
 if(loadfields) call read_fields(t)
 end subroutine
 subroutine prepare_fields(path)
 character(*),intent(in)::path
 type(projection_data)::t
 type(local_potential)::m
 integer::d,g,x,y,z,u,ng
 real(dp)::point(3)
 real(dp),allocatable::v(:),f(:,:),means(:)
 call control(path,t,.false.);call read_vfunc(trim(model_path),m,.true.)
 ng=t%ng;allocate(v(t%nc),f(ng,t%nc),means(m%np+1));means=0
 open(newunit=u,file=trim(rvr_dir)//'/spatial_fields',status='new',access='stream', &
 form='unformatted',convert='little_endian')
 write(u) 'LANIO2_LGRID_V2 ',int([1,t%ng,t%nc,t%nd,m%np,0,0,0],int32),m%id,int(t%domain,int32)
 do d=1,t%nd
 do g=1,ng
  x=mod(g-1,m%grid(1));y=mod((g-1)/m%grid(1),m%grid(2));z=(g-1)/(m%grid(1)*m%grid(2))
  point=matmul(transpose(m%cell),real([x,y,z],dp)/m%grid+real([t%domain(:,d),0],dp))-m%center
  call plane_value(m,point,v);f(g,:)=v
 enddo
 means=means+sum(f(:,:m%np+1),dim=1)/ng
 write(u) f
 enddo
 write(u) means;close(u)
 write(*,'(a,3i12)') 'SPATIAL_FIELDS ng channels domains=',ng,t%nc,t%nd
 end subroutine
 subroutine read_fields(t)
 type(projection_data),intent(inout)::t
 integer::u,ios
 integer(int32)::head(8)
 integer,allocatable::domain(:,:)
 integer(int64)::bytes,expected
 character(16)::magic
 character(64)::id
 open(newunit=u,file=trim(rvr_dir)//'/spatial_fields',status='old',access='stream', &
 form='unformatted',convert='little_endian',iostat=ios)
 call require(ios==0,'Run cal_RVR_split projection.nml 0 once before row projection')
 read(u,iostat=ios) magic,head,id
 call require(ios==0.and.magic=='LANIO2_LGRID_V2 '.and.id==t%m%id,'Wrong spatial cache model')
 call require(all(head==[1,t%ng,t%nc,t%nd,t%m%np,0,0,0]),'Wrong spatial cache dimensions')
 inquire(unit=u,size=bytes)
 expected=112_int64+8_int64*t%nd+8_int64*t%ng*t%nc*t%nd+8_int64*(t%m%np+1)
 call require(bytes==expected,'Spatial cache byte count mismatch')
 allocate(domain(2,t%nd));read(u) domain;call require(all(domain==t%domain),'Spatial domain changed')
 allocate(t%fields(t%ng,t%nc,t%nd),t%means(t%m%np+1));read(u,iostat=ios) t%fields,t%means;close(u)
 call require(ios==0.and.all(ieee_is_finite(t%fields)).and.all(ieee_is_finite(t%means)), &
 'Truncated/nonfinite spatial cache')
 end subroutine
 function row_path(i,k) result(path)
 integer,intent(in)::i,k
 character(4096)::path
 character(32)::a,b
 write(a,'(i0)') i;write(b,'(i0)') k
 path=trim(rvr_dir)//'/B_general_'//trim(a)//'_kz'//trim(b)
 end function
 subroutine read_wan(t,orb,raw)
 type(projection_data),intent(in)::t
 integer,intent(in)::orb
 complex(dp),allocatable,intent(out)::raw(:,:)
 integer::u,ios
 integer(int64)::bytes
 character(32)::a
 write(a,'(i0)') orb
 open(newunit=u,file=trim(wan_dir)//'/wan_'//trim(a),status='old',access='stream', &
 form='unformatted',convert='little_endian',iostat=ios)
 call require(ios==0,'Missing primitive complex128 wan_'//trim(a));inquire(unit=u,size=bytes)
 call require(bytes==16_int64*t%ns*t%ng,'Wannier byte count/grid/support mismatch')
 allocate(raw(t%ns,t%ng));read(u,iostat=ios) raw;close(u)
 call require(ios==0.and.all(ieee_is_finite(real(raw,dp))).and.all(ieee_is_finite(aimag(raw))), &
 'Truncated/nonfinite Wannier data')
 end subroutine
 subroutine plane_wan(t,raw,k,w)
 type(projection_data),intent(in)::t
 complex(dp),intent(in)::raw(:,:)
 integer,intent(in)::k
 complex(dp),intent(out)::w(t%ng,t%nsxy)
 integer::s,j
 complex(dp)::phase
 w=0
 do s=1,t%ns
 do j=1,t%nsxy
  if(any(t%source(:2,s)/=t%source_xy(:,j))) cycle
  phase=exp(-2*pi*zi*t%kz(k)*t%source(3,s))*t%weights(s)
  w(:,j)=w(:,j)+phase*raw(s,:);exit
 enddo
 enddo
 end subroutine
 subroutine project_row(path,column,kfirst,klast)
 character(*),intent(in)::path
 integer,intent(in)::column,kfirst,klast
 type(projection_data)::t
 complex(dp),allocatable::rawi(:,:),rawj(:,:),left(:,:),right(:,:),work(:,:),mat(:,:),kernel(:,:,:),overlap(:)
 integer::i,j,k,d,c,a,b,u,first,last
 call control(path,t,.true.);i=column
 call require(i>=1.and.i<=nwan,'Orbital row outside 1..nwan')
 first=kfirst;last=klast;if(first==0) then
 first=1;last=t%nk
 endif
 call require(first>=1.and.last>=first.and.last<=t%nk,'Invalid kz shard')
 allocate(left(t%ng,t%nsxy),right(t%ng,t%nsxy),work(t%ng,t%nsxy),mat(t%nsxy,t%nsxy), &
 kernel(t%nr,t%nr,t%nc),overlap(t%ndiff))
 call read_wan(t,i,rawi)
 do k=first,last
  call plane_wan(t,rawi,k,left)
  open(newunit=u,file=trim(row_path(i,k)),status='new',access='stream',form='unformatted',convert='little_endian')
  write(u) 'LANIO2_LROW_V2  ',int([1,t%m%np,t%nr,nwan,t%nk,t%ndiff,i,k],int32),t%m%id,basis_id
  write(u) int(t%xy,int32),int(t%difference,int32),t%kz(k),grid_weight
  do j=1,nwan
   call read_wan(t,j,rawj);call plane_wan(t,rawj,k,right);deallocate(rawj)
   call gemm(left,right,mat,ta='C',alpha=cmplx(grid_weight,0._dp,dp));overlap=0
   do a=1,t%nsxy
   do b=1,t%nsxy
    overlap(t%overlap_map(a,b))=overlap(t%overlap_map(a,b))+mat(a,b)
   enddo
   enddo
   kernel=0
   do d=1,t%nd
   do c=1,t%nc
    if(maxval(abs(t%fields(:,c,d)))==0) cycle
    do b=1,t%nsxy
     work(:,b)=right(:,b)*(grid_weight*t%fields(:,c,d))
    enddo
    call gemm(left,work,mat,ta='C')
    do a=1,t%nsxy
    do b=1,t%nsxy
     kernel(t%mapping(a,d),t%mapping(b,d),c)=kernel(t%mapping(a,d),t%mapping(b,d),c)+mat(a,b)
    enddo
    enddo
   enddo
   enddo
   write(u) overlap,kernel
  enddo
  close(u)
  write(*,'(a,2i8)') 'RVR_ROW_COMPLETE orbital,kz=',i,k
 enddo
 end subroutine
 subroutine combine_rows(path)
 character(*),intent(in)::path
 type(projection_data)::t
 integer::u,v,i,j,k,ios
 integer(int32)::head(8)
 integer,allocatable::xy(:,:),difference(:,:)
 complex(dp),allocatable::kernel(:,:,:),overlap(:)
 character(16)::magic
 character(64)::id,bid
 real(dp)::kz,weight
 integer(int64)::bytes,expected
 call control(path,t,.true.)
 allocate(kernel(t%nr,t%nr,t%nc),overlap(t%ndiff),xy(2,t%nr),difference(2,t%ndiff))
 open(newunit=v,file=trim(rvr_dir)//'/B_general',status='new',access='stream', &
 form='unformatted',convert='little_endian')
 write(v) 'LANIO2_LRVR_V2  ',int([1,t%m%np,t%nr,nwan,t%nk,t%ndiff,0,0],int32),t%m%id,basis_id
 write(v) t%m%cell,t%m%alpha,t%means,int(t%xy,int32),int(t%difference,int32),t%kz
 do k=1,t%nk
 do i=1,nwan
  open(newunit=u,file=trim(row_path(i,k)),status='old',access='stream',form='unformatted',convert='little_endian',iostat=ios)
  call require(ios==0,'Missing RVR orbital/kz row')
  read(u,iostat=ios) magic,head,id,bid,xy,difference,kz,weight
  call require(ios==0.and.magic=='LANIO2_LROW_V2  '.and.id==t%m%id.and.bid==basis_id,'Wrong RVR row provenance')
  call require(all(head==[1,t%m%np,t%nr,nwan,t%nk,t%ndiff,i,k]).and.all(xy==t%xy).and. &
   all(difference==t%difference).and.kz==t%kz(k).and.weight==grid_weight,'RVR row contract mismatch')
  inquire(unit=u,size=bytes)
  expected=192_int64+8_int64*(t%nr+t%ndiff)+16_int64*nwan*(t%ndiff+int(t%nr,int64)**2*t%nc)
  call require(bytes==expected,'RVR row byte count mismatch')
  do j=1,nwan
   read(u,iostat=ios) overlap,kernel
   call require(ios==0.and.all(ieee_is_finite(real(kernel,dp))).and.all(ieee_is_finite(aimag(kernel))).and. &
    all(ieee_is_finite(real(overlap,dp))).and.all(ieee_is_finite(aimag(overlap))),'Nonfinite/truncated RVR row')
   write(v) overlap,kernel
  enddo
  close(u)
 enddo
 enddo
 close(v)
 write(*,'(a)') 'B_GENERAL_COMPLETE'
 end subroutine
 subroutine open_projected(path,u,t)
 character(*),intent(in)::path
 integer,intent(out)::u
 type(local_projected),intent(out)::t
 character(16)::magic
 integer(int32)::head(8)
 integer::ios
 integer(int64)::bytes,expected
 open(newunit=u,file=path,status='old',access='stream',form='unformatted',convert='little_endian',iostat=ios)
 call require(ios==0,'Cannot open local B_general')
 read(u,iostat=ios) magic,head,t%model_id,t%basis_id
 call require(ios==0.and.magic=='LANIO2_LRVR_V2  '.and.head(1)==1.and.all(head(2:6)>0),'Invalid local B_general header')
 t%np=head(2);t%nc=2*t%np+1;t%nr=head(3);t%nw=head(4);t%nk=head(5);t%ndiff=head(6)
 call require(t%np<=10000.and.t%nr<=100000.and.t%nw<=30.and.t%nk<=10000.and.t%ndiff<=100000,'Implausible RVR dimensions')
 allocate(t%means(t%np+1),t%xy(2,t%nr),t%difference(2,t%ndiff),t%kz(t%nk))
 read(u,iostat=ios) t%cell,t%alpha,t%means,t%xy,t%difference,t%kz
 call require(ios==0.and.all(ieee_is_finite(t%means)).and.all(ieee_is_finite(t%kz)),'Truncated RVR metadata')
 inquire(unit=u,pos=t%data_position,size=bytes)
 t%block_bytes=16_int64*(t%ndiff+int(t%nr,int64)**2*t%nc)
 expected=t%data_position-1+t%block_bytes*t%nw*t%nw*t%nk
 call require(bytes==expected,'B_general byte count mismatch')
 end subroutine
 subroutine projected_block(u,t,k,i,j,s,kernel)
 integer,intent(in)::u,k,i,j
 type(local_projected),intent(in)::t
 complex(dp),intent(out)::s(t%ndiff),kernel(t%nr,t%nr,t%nc)
 integer(int64)::offset
 integer::ios
 offset=t%data_position+t%block_bytes*((int(k,int64)-1)*t%nw*t%nw+(i-1)*t%nw+j-1)
 read(u,pos=offset,iostat=ios) s,kernel
 call require(ios==0,'Truncated projected block')
 end subroutine
end module
