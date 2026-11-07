module frozen_io
 use numerics
 implicit none
 type::frozen_model
 integer::nr=0
 integer(int32),allocatable::rr(:,:)
 complex(dp),allocatable::hc(:,:,:),p(:,:,:),kb(:,:)
 complex(dp)::n0(10,10)
 end type
contains
 subroutine read_model(path,m)
 character(*),intent(in)::path
 type(frozen_model),intent(out)::m
 character(16)::magic
 integer(int32)::head(8)
 integer::u,ios
 integer(int64)::bytes,expected
 character(1024)::msg
 open(newunit=u,file=path,status='old',access='stream',form='unformatted',action='read', &
 convert='little_endian',iostat=ios,iomsg=msg)
 call require(ios==0,'Cannot open explicit model path: '//trim(path)//' '//trim(msg))
 read(u,iostat=ios) magic,head
 call require(ios==0,'Truncated model header')
 call require(magic=='LANIO2_AU_V1    ','Wrong model magic/version')
 call require(all(head(1:7)==[16909060,729,30,5,10,810,1]),'Model dimensions/endian/unit contract mismatch')
 call require(head(8)==0,'Background contains HU or has unknown Hubbard accounting')
 m%nr=head(2)
 expected=48_int64+12_int64*m%nr+16_int64*(900_int64*m%nr+150_int64*m%nr+810_int64**2+100)
 inquire(unit=u,size=bytes)
 call require(bytes==expected,'Model byte count mismatch (complex128/Fortran order required)')
 allocate(m%rr(3,m%nr),m%hc(30,30,m%nr),m%p(5,30,m%nr),m%kb(810,810))
 read(u,iostat=ios) m%rr,m%hc,m%p,m%kb,m%n0
 close(u)
 call require(ios==0,'Truncated model data')
 call require(all(ieee_is_finite(real(m%hc,dp))).and.all(ieee_is_finite(aimag(m%hc))),'Nonfinite hc_R')
 call require(all(ieee_is_finite(real(m%p,dp))).and.all(ieee_is_finite(aimag(m%p))),'Nonfinite p_R')
 call require(all(ieee_is_finite(real(m%kb,dp))).and.all(ieee_is_finite(aimag(m%kb))),'Nonfinite KB')
 call require(all(ieee_is_finite(real(m%n0,dp))).and.all(ieee_is_finite(aimag(m%n0))),'Nonfinite N0')
 call require(herm_error(m%kb)<1e-10_dp.and.herm_error(m%n0)<1e-12_dp,'Non-Hermitian KB/N0')
 end subroutine
 subroutine interpolate(m,q,hc,p)
 type(frozen_model),intent(in)::m
 real(dp),intent(in)::q(3)
 complex(dp),intent(out)::hc(30,30),p(5,30)
 complex(dp)::phase
 integer::r
 hc=0;p=0
 do r=1,m%nr
 phase=exp(2*pi*zi*dot_product(q,real(m%rr(:,r),dp)))
 hc=hc+phase*m%hc(:,:,r);p=p+phase*m%p(:,:,r)
 enddo
 call require(herm_error(hc)<1e-9_dp,'Interpolated H0+Hv is not Hermitian')
 end subroutine
 subroutine plane_kernel(m,kz,ke)
 type(frozen_model),intent(in)::m
 real(dp),intent(in)::kz
 complex(dp),intent(out)::ke(270,270)
 integer::a,b,z,w,i,j
 ke=0
 ! Translation order: x,y,z, each -1..1; orbital fastest (same as NumPy R3).
 do a=0,8
 do b=0,8
 do z=0,2
 do w=0,2
 do j=1,30
 do i=1,30
 ke(a*30+i,b*30+j)=ke(a*30+i,b*30+j)+ &
 exp(2*pi*zi*kz*(w-z))*m%kb((a*3+z)*30+i,(b*3+w)*30+j)
 enddo
 enddo
 enddo
 enddo
 enddo
 enddo
 call require(herm_error(ke)<1e-10_dp,'Contracted B RVR kernel is not Hermitian')
 end subroutine
end module
