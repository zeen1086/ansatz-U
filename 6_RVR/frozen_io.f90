module frozen_io
 use numerics
 implicit none
 type::frozen_model
 integer::nr=0
 integer(int32),allocatable::rr(:,:)
 complex(dp),allocatable::hc(:,:,:),p(:,:,:)
 end type
contains
 subroutine read_model(path,m)
 character(*),intent(in)::path
 type(frozen_model),intent(out)::m
 ! Keep the LANIO2_AU_V1 layout; legacy KB/N0 fields are validated and discarded.
 complex(dp),allocatable::kb(:,:)
 complex(dp)::n0(10,10)
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
 allocate(m%rr(3,m%nr),m%hc(30,30,m%nr),m%p(5,30,m%nr),kb(810,810))
 read(u,iostat=ios) m%rr,m%hc,m%p,kb,n0
 close(u)
 call require(ios==0,'Truncated model data')
 call require(all(ieee_is_finite(real(m%hc,dp))).and.all(ieee_is_finite(aimag(m%hc))),'Nonfinite hc_R')
 call require(all(ieee_is_finite(real(m%p,dp))).and.all(ieee_is_finite(aimag(m%p))),'Nonfinite p_R')
 call require(all(ieee_is_finite(real(kb,dp))).and.all(ieee_is_finite(aimag(kb))),'Nonfinite KB')
 call require(all(ieee_is_finite(real(n0,dp))).and.all(ieee_is_finite(aimag(n0))),'Nonfinite N0')
 call require(herm_error(kb)<1e-10_dp.and.herm_error(n0)<1e-12_dp,'Non-Hermitian KB/N0')
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
 call require(herm_error(hc)<1e-9_dp,'Interpolated H0 is not Hermitian')
 end subroutine
 subroutine project_at(L,k,pg,site,a)
 integer,intent(in)::L,site(2)
 real(dp),intent(in)::k(3)
 complex(dp),intent(in)::pg(:,:,0:,0:)
 complex(dp),intent(out)::a(10,60*L*L)
 integer::x,y,d(2),orb,cell
 complex(dp)::phase
 a=0
 do x=0,L-1
 do y=0,L-1
 d=site-[x,y];cell=x*L+y
 phase=exp(2*pi*zi*dot_product(real(d,dp),k(:2)/L))
 do orb=1,30
 a(:5,cell*60+2*orb-1)=pg(:,orb,modulo(d(1),L),modulo(d(2),L))*phase
 a(6:,cell*60+2*orb)=pg(:,orb,modulo(d(1),L),modulo(d(2),L))*phase
 enddo
 enddo
 enddo
 end subroutine
 subroutine background_projectors(m,L,k,hcv,pg)
 type(frozen_model),intent(in)::m
 integer,intent(in)::L
 real(dp),intent(in)::k(3)
 complex(dp),intent(out)::hcv(:,:),pg(:,:,0:,0:)
 complex(dp),allocatable::hcq(:,:,:),pq(:,:,:)
 complex(dp)::hc(30,30),p(5,30),phase
 real(dp)::q(3)
 integer::x,y,t,xx,yy,s,d(2),i,j,u,iq,nc,nd
 nc=L*L;nd=60*nc
 allocate(hcq(30,30,nc),pq(5,30,nc))
 do x=0,L-1
 do y=0,L-1
 iq=x*L+y+1;q=[(x+k(1))/L,(y+k(2))/L,k(3)]
 call interpolate(m,q,hc,p);hcq(:,:,iq)=hc;pq(:,:,iq)=p
 enddo
 enddo
 pg=0;hcv=0
 ! Finite inverse Bloch transform: the same normalized IFFT and twist as the reference.
 do x=0,L-1
 do y=0,L-1
 do xx=0,L-1
 do yy=0,L-1
 iq=xx*L+yy+1
 phase=exp(2*pi*zi*real(xx*x+yy*y,dp)/L)/nc
 pg(:,:,x,y)=pg(:,:,x,y)+phase*pq(:,:,iq)
 enddo
 enddo
 enddo
 enddo
 do x=0,L-1
 do y=0,L-1
 t=x*L+y
 do xx=0,L-1
 do yy=0,L-1
 s=xx*L+yy;d=[x-xx,y-yy];hc=0
 do iq=1,nc
 q=[( (iq-1)/L+k(1))/L,(modulo(iq-1,L)+k(2))/L,k(3)]
 phase=exp(2*pi*zi*dot_product(q(:2),real(d,dp)))/nc
 hc=hc+phase*hcq(:,:,iq)
 enddo
 do j=1,30
 do i=1,30
 do u=1,2
 hcv(t*60+2*i-2+u,s*60+2*j-2+u)=hc(i,j)
 enddo
 enddo
 enddo
 enddo
 enddo
 enddo
 enddo
 end subroutine
end module
