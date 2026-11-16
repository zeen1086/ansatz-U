module joint_npy
 use numerics
 use iso_fortran_env, only:int8
 implicit none
contains
 function quoted_value(h,key) result(value)
 character(*),intent(in)::h,key
 character(:),allocatable::value
 integer::p,q,r
 character::quote
 p=index(h,key);call require(p>0,'NPY header missing '//key)
 p=p+len(key);q=index(h(p:),':');call require(q>0,'Invalid NPY dictionary')
 p=p+q
 do while(p<=len(h))
  if(h(p:p)=="'".or.h(p:p)=='"') exit
  p=p+1
 enddo
 call require(p<len(h),'Invalid NPY string');quote=h(p:p)
 r=index(h(p+1:),quote);call require(r>0,'Unterminated NPY string')
 value=h(p+1:p+r-1)
 end function
 subroutine npy_header(path,u,descriptor,shape,fort)
 character(*),intent(in)::path
 integer,intent(out)::u
 character(:),allocatable,intent(out)::descriptor
 integer,allocatable,intent(out)::shape(:)
 logical,intent(out)::fort
 character(6)::magic
 character(:),allocatable::h,body
 integer(int8)::version(2),bytes(4)
 integer(int64)::length
 integer::ios,n,p,q,i,j,dims(4),rank
 open(newunit=u,file=path,status='old',access='stream',form='unformatted',convert='little_endian',iostat=ios)
 call require(ios==0,'Cannot open DFT NPY: '//trim(path))
 read(u,iostat=ios) magic,version
 call require(ios==0.and.magic==achar(147)//'NUMPY','Bad NPY magic: '//trim(path))
 n=2;if(version(1)==2.or.version(1)==3) n=4
 call require(version(1)>=1.and.version(1)<=3,'Unsupported NPY version')
 read(u,iostat=ios) bytes(:n);call require(ios==0,'Truncated NPY length')
 length=0
 do i=1,n
  length=length+int(iand(int(bytes(i)),255),int64)*256_int64**(i-1)
 enddo
 call require(length>0.and.length<=1048576,'Invalid NPY header length')
 allocate(character(int(length))::h);read(u,iostat=ios) h
 call require(ios==0,'Truncated NPY header');descriptor=quoted_value(h,'descr')
 p=index(h,'fortran_order');call require(p>0,'NPY order missing')
 q=index(h(p:),':');call require(q>0,'Invalid NPY order');p=p+q
 h(p:)=adjustl(h(p:));fort=index(h(p:),'True')==1
 call require(fort.or.index(h(p:),'False')==1,'Invalid NPY order value')
 p=index(h,'shape');call require(p>0,'NPY shape missing')
 q=index(h(p:),'(');call require(q>0,'Invalid NPY shape');p=p+q
 q=index(h(p:),')');call require(q>0,'Invalid NPY shape');body=h(p:p+q-2)
 rank=0;i=1
 do while(i<=len(body))
  if(body(i:i)==' '.or.body(i:i)==',') then
   i=i+1;cycle
  endif
  j=i
  do while(j<=len(body))
   if(body(j:j)==','.or.body(j:j)==' ') exit
   j=j+1
  enddo
  rank=rank+1;call require(rank<=4,'NPY rank exceeds4')
  read(body(i:j-1),*,iostat=ios) dims(rank)
  call require(ios==0.and.dims(rank)>0,'Invalid NPY dimension');i=j
 enddo
 call require(rank>0,'Scalar/object NPY is not a DFT array')
 shape=dims(:rank)
 end subroutine
 subroutine exact_bytes(u,shape,itemsize)
 integer,intent(in)::u,shape(:),itemsize
 integer(int64)::position,bytes,expected
 integer::i
 inquire(unit=u,pos=position,size=bytes)
 expected=itemsize
 do i=1,size(shape)
  call require(expected<=huge(1_int32)/int(shape(i),int64),'DFT array too large for this label reader')
  expected=expected*shape(i)
 enddo
 call require(bytes==position-1+expected,'NPY data byte count/shape mismatch')
 end subroutine
 subroutine npy_real(path,a,shape,fort)
 character(*),intent(in)::path
 real(dp),allocatable,intent(out)::a(:)
 integer,allocatable,intent(out)::shape(:)
 logical,intent(out)::fort
 character(:),allocatable::dtype
 integer::u,ios
 call npy_header(path,u,dtype,shape,fort)
 call require(dtype=='<f8','DFT real arrays must be little-endian float64: '//trim(path))
 call exact_bytes(u,shape,8);allocate(a(product(shape)));read(u,iostat=ios) a;close(u)
 call require(ios==0.and.all(ieee_is_finite(a)),'Truncated/nonfinite DFT real array')
 end subroutine
 subroutine npy_complex(path,a,shape,fort)
 character(*),intent(in)::path
 complex(dp),allocatable,intent(out)::a(:)
 integer,allocatable,intent(out)::shape(:)
 logical,intent(out)::fort
 character(:),allocatable::dtype
 integer::u,ios
 call npy_header(path,u,dtype,shape,fort)
 call require(dtype=='<c16','DFT complex arrays must be little-endian complex128: '//trim(path))
 call exact_bytes(u,shape,16);allocate(a(product(shape)));read(u,iostat=ios) a;close(u)
 call require(ios==0.and.all(ieee_is_finite(real(a,dp))).and.all(ieee_is_finite(aimag(a))), &
  'Truncated/nonfinite DFT complex array')
 end subroutine
 subroutine npy_integer(path,a)
 character(*),intent(in)::path
 integer,allocatable,intent(out)::a(:)
 integer,allocatable::shape(:)
 integer(int32),allocatable::small(:)
 integer(int64),allocatable::large(:)
 character(:),allocatable::dtype
 logical::fort
 integer::u,ios
 call npy_header(path,u,dtype,shape,fort);call require(size(shape)==1,'Expected integer vector')
 allocate(a(shape(1)))
 select case(dtype)
 case('<i4')
  call exact_bytes(u,shape,4);allocate(small(size(a)));read(u,iostat=ios) small;a=int(small)
 case('<i8')
  call exact_bytes(u,shape,8);allocate(large(size(a)));read(u,iostat=ios) large
  call require(all(abs(large)<=huge(1_int32)),'Integer label outside int32 range');a=int(large)
 case default
  call require(.false.,'DFT integer arrays require little-endian int32/int64')
 end select
 close(u);call require(ios==0,'Truncated DFT integer vector')
 end subroutine
 subroutine read_occupation(path,n)
 character(*),intent(in)::path
 complex(dp),intent(out)::n(10,10,4)
 complex(dp),allocatable::a(:)
 integer,allocatable::shape(:)
 logical::fort
 integer::i,j,s,index
 call npy_complex(path,a,shape,fort)
 call require(size(shape)==3,'Occupation must have rank3')
 call require(all(shape==[4,10,10]),'Occupation order must be (Ni5..8, up d1..5/down d1..5, same)')
 do s=1,4
 do j=1,10
 do i=1,10
  index=j+10*(i-1)+100*(s-1)
  if(fort) index=s+4*(i-1)+40*(j-1)
  n(i,j,s)=a(index)
 enddo
 enddo
 enddo
 end subroutine
 subroutine read_dft_label(folder,grid,atom_ids,unit,normalize,ueff,tolerance,e,amplitude,field,n)
 character(*),intent(in)::folder,unit
 integer,intent(in)::grid(3),atom_ids(4)
 logical,intent(in)::normalize
 real(dp),intent(in)::ueff,tolerance
 real(dp),intent(out)::e(3,4),amplitude(4),field(product(grid),4)
 complex(dp),intent(out)::n(10,10,4)
 real(dp),allocatable::a(:)
 complex(dp)::vu(10,10,4),expected(10,10)
 integer,allocatable::shape(:),ids(:)
 logical::fort,exists
 integer::i,j,s,ng,index
 real(dp)::scale
 call npy_integer(trim(folder)//'/grid.npy',ids)
 call require(size(ids)==3,'DFT grid must have3 dimensions');call require(all(ids==grid),'DFT grids differ')
 call npy_integer(trim(folder)//'/atom_ids.npy',ids)
 call require(size(ids)==4,'Exactly four DFT Ni atoms required');call require(all(ids==atom_ids),'DFT Ni order mismatch')
 ng=product(grid);scale=1
 call require(unit=='Ry'.or.unit=='eV','Declare each DFT field unit as Ry or eV')
 if(unit=='Ry') scale=ry_ev
 call npy_real(trim(folder)//'/scalar.npy',a,shape,fort)
 call require(size(shape)==1.and.size(a)==ng,'Scalar field must be a flattened x-fast grid')
 field(:,1)=a*scale
 call npy_real(trim(folder)//'/B.npy',a,shape,fort)
 call require(size(shape)==2,'B field must have rank2');call require(all(shape==[ng,3]),'B shape must be (ngrid,3)')
 do i=1,ng
 do j=1,3
  index=j+3*(i-1);if(fort) index=i+ng*(j-1)
  field(i,j+1)=a(index)*scale
 enddo
 enddo
 call npy_real(trim(folder)//'/moments.npy',a,shape,fort)
 call require(size(shape)==2,'Moments must have rank2');call require(all(shape==[4,3]),'Moments order must be (Ni5..8,xyz)')
 do s=1,4
 do i=1,3
  index=i+3*(s-1);if(fort) index=s+4*(i-1)
  e(i,s)=a(index)
 enddo
 enddo
 amplitude=sqrt(sum(e**2,dim=1));call require(all(amplitude>1e-12_dp),'Zero DFT moment direction')
 if(normalize) then
  do s=1,4
   e(:,s)=e(:,s)/amplitude(s)
  enddo
 else
  call require(maxval(abs(amplitude-1))<1e-9_dp,'Expected unit spin directions; set normalize_moments explicitly')
 endif
 call read_occupation(trim(folder)//'/N.npy',n)
 do s=1,4
  call require(herm_error(n(:,:,s))<=tolerance,'Non-Hermitian raw DFT N; no symmetrization applied')
 enddo
 inquire(file=trim(folder)//'/VU.npy',exist=exists)
 if(exists) then
  call read_occupation(trim(folder)//'/VU.npy',vu)
  do s=1,4
   expected=-ueff*n(:,:,s)
   do i=1,10
    expected(i,i)=expected(i,i)+ueff/2
   enddo
   call require(maxval(abs(vu(:,:,s)-expected))<=tolerance,'DFT VU != Ueff*(I/2-N); check training U/unit')
  enddo
 endif
 end subroutine
end module
