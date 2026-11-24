module npy_io
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
end module
