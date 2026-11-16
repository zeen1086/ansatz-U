module local_vfunc
 use numerics
 use joint_npy, only: npy_header,exact_bytes
 implicit none
 type::local_potential
  integer::np,nn,grid(3),half(3)
  real(dp)::cell(3,3),center(3),spacing(3),cutoff,alpha
  character(64)::id
  integer,allocatable::lookup(:,:,:)
  real(dp),allocatable::cv(:,:),cb(:,:)
 end type
contains
 subroutine read_vfunc(path,m,coefficients)
 character(*),intent(in)::path
 type(local_potential),intent(out)::m
 logical,intent(in)::coefficients
 character(16)::magic
 integer(int32)::head(8)
 integer,allocatable::nodes(:,:)
 integer::u,ios,i,ix(3)
 integer(int64)::bytes
 open(newunit=u,file=trim(path)//'/native_geometry.bin',status='old',access='stream', &
  form='unformatted',convert='little_endian',iostat=ios)
 call require(ios==0,'Missing native V_func_tmp geometry; run local_vfunc.py export')
 read(u,iostat=ios) magic,head
 call require(ios==0.and.magic=='LANIO2_VFUNC_V2 ','Invalid V_func_tmp native header')
 call require(head(1)==1.and.head(2)>0.and.head(3)>0.and.all(head(4:6)>0), 'Invalid local dimensions')
 m%np=head(2);m%nn=head(3);m%grid=head(4:6)
 read(u,iostat=ios) m%half,m%cell,m%center,m%spacing,m%cutoff,m%alpha,m%id
 call require(ios==0.and.all(m%half>0).and.all(m%spacing>0).and.m%cutoff>0,'Invalid local geometry')
 call require(all(ieee_is_finite(m%cell)).and.all(ieee_is_finite(m%center)).and. &
 all(ieee_is_finite(m%spacing)).and.ieee_is_finite(m%cutoff).and.ieee_is_finite(m%alpha), &
 'Nonfinite local geometry')
 inquire(unit=u,size=bytes)
 call require(bytes==260_int64+12_int64*m%nn,'V_func_tmp geometry byte count mismatch')
 allocate(nodes(3,m%nn));read(u,iostat=ios) nodes;close(u)
 call require(ios==0,'Truncated local nodes')
 allocate(m%lookup(-m%half(1):m%half(1),-m%half(2):m%half(2),-m%half(3):m%half(3)))
 m%lookup=0
 do i=1,m%nn
  ix=nodes(:,i)
  call require(all(abs(ix)<=m%half),'Local node outside grid')
  call require(m%lookup(ix(1),ix(2),ix(3))==0,'Duplicate local node')
  m%lookup(ix(1),ix(2),ix(3))=i
 enddo
 if(coefficients) then
  call read_kernel(trim(path)//'/C_v_local.npy',m%np+1,m%nn,m%cv)
  call read_kernel(trim(path)//'/C_B_local.npy',m%np,m%nn,m%cb)
 endif
 end subroutine
 subroutine read_kernel(path,nc,nn,a)
 character(*),intent(in)::path
 integer,intent(in)::nc,nn
 real(dp),allocatable,intent(out)::a(:,:)
 integer::u,ios
 integer,allocatable::shape(:)
 character(:),allocatable::dtype
 logical::fort
 call npy_header(path,u,dtype,shape,fort)
 call require(size(shape)==2,'Local kernel must be rank two')
 call require(dtype=='<f8'.and..not.fort.and.all(shape==[nn,nc]),'Local kernel dtype/order/dimensions mismatch')
 call exact_bytes(u,shape,8);allocate(a(nc,nn));read(u,iostat=ios) a;close(u)
 call require(ios==0.and.all(ieee_is_finite(a)),'Truncated/nonfinite local kernel')
 end subroutine
 subroutine local_value(m,point,value)
 type(local_potential),intent(in)::m
 real(dp),intent(in)::point(3)
 real(dp),intent(out)::value(2*m%np+1)
 real(dp)::r,pos(3),fraction(3),weight,win
 integer::low(3),index(3),x,y,z,row,d(3)
 value=0;r=sqrt(sum(point**2));if(r>=m%cutoff) return
 pos=point/m%spacing;low=floor(pos);fraction=pos-low
 win=cos(pi*r/(2*m%cutoff))**2
 do x=0,1
 do y=0,1
 do z=0,1
  d=[x,y,z];index=low+d
  if(any(abs(index)>m%half)) cycle
  row=m%lookup(index(1),index(2),index(3));if(row==0) cycle
  weight=product(merge(fraction,1-fraction,d==1))*win
  value(:m%np+1)=value(:m%np+1)+weight*m%cv(:,row)
  value(m%np+2:)=value(m%np+2:)+weight*m%cb(:,row)
 enddo
 enddo
 enddo
 end subroutine
 subroutine plane_value(m,point,value)
 type(local_potential),intent(in)::m
 real(dp),intent(in)::point(3)
 real(dp),intent(out)::value(2*m%np+1)
 real(dp)::p(3),v(2*m%np+1),c
 integer::z,lo,hi
 call require(maxval(abs(m%cell-reshape([m%cell(1,1),0._dp,0._dp,0._dp,m%cell(2,2), &
  0._dp,0._dp,0._dp,m%cell(3,3)],[3,3])))<1e-12_dp,'Orthorhombic cell required')
 c=m%cell(3,3);call require(c>0,'Positive z period required')
 lo=ceiling((-point(3)-m%cutoff)/c);hi=floor((-point(3)+m%cutoff)/c);value=0
 do z=lo,hi
  p=point;p(3)=p(3)+z*c;call local_value(m,p,v);value=value+v
 enddo
 end subroutine
end module
