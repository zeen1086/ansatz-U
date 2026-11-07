module tensor_rvr
 use numerics
 implicit none
 type::tensor_kernel
 integer::range=0
 complex(dp),allocatable::rvr(:,:,:)
 end type
contains
 subroutine read_tensor(path,radius,t,energy_unit)
 character(*),intent(in)::path,energy_unit
 integer,intent(in)::radius
 type(tensor_kernel),intent(out)::t
 integer::u,ios,d
 integer(int32)::leading,trailing
 integer(int64)::bytes,expected
 call require(radius>=0.and.radius<=4,'Invalid explicit tensor RVR range')
 call require(energy_unit=='Ry'.or.energy_unit=='eV','Declare tensor_unit=Ry or eV for B_general')
 d=(2*radius+1)**3;t%range=radius
 expected=16_int64*3*(60_int64*d)**2
 call require(expected<=huge(leading),'Tensor RVR requires one record with a positive int32 byte count')
 open(newunit=u,file=path,status='old',action='read',access='stream',form='unformatted', &
 convert='little_endian',iostat=ios)
 call require(ios==0,'Cannot open original full-component B_general tensor RVR')
 inquire(unit=u,size=bytes)
 call require(bytes==expected+8,'Tensor RVR byte count mismatch: n_func=3, nwan=30, nspin=1 required')
 read(u,iostat=ios) leading
 call require(ios==0,'Truncated tensor RVR record marker')
 call require(int(leading,int64)==expected,'Tensor RVR leading record marker mismatch')
 allocate(t%rvr(3,60*d,60*d))
 read(u,iostat=ios) t%rvr,trailing
 close(u)
 call require(ios==0,'Tensor RVR record mismatch: n_func=3, nwan=30, nspin=1, complex128 required')
 call require(trailing==leading,'Tensor RVR trailing record marker mismatch')
 if(energy_unit=='Ry') t%rvr=t%rvr*ry_ev
 call require(all(ieee_is_finite(real(t%rvr,dp))).and.all(ieee_is_finite(aimag(t%rvr))),'Nonfinite tensor RVR')
 end subroutine
 subroutine tensor_plane(t,kz,ke)
 type(tensor_kernel),intent(in)::t
 real(dp),intent(in)::kz
 complex(dp),intent(out)::ke(:,:,:)
 integer::a,b,z,w,i,j,c,n,d,iz,jz
 n=2*t%range+1;d=n**3;ke=0
 do c=1,3
 do a=0,n*n-1
 do b=0,n*n-1
 do z=0,n-1
 do w=0,n-1
 do j=1,60
 do i=1,60
 iz=a*n+z+1+d*(i-1);jz=b*n+w+1+d*(j-1)
 ke(a*60+i,b*60+j,c)=ke(a*60+i,b*60+j,c)+exp(2*pi*zi*kz*(w-z))*t%rvr(c,iz,jz)
 enddo
 enddo
 enddo
 enddo
 enddo
 enddo
 call require(herm_error(ke(:,:,c))<1e-8_dp,'Non-Hermitian full-component tensor RVR plane')
 enddo
 end subroutine
 subroutine tensor_embedding(L,k,site,radius,v,local)
 integer,intent(in)::L,site(2),radius
 real(dp),intent(in)::k(3)
 complex(dp),intent(in)::v(:,:)
 complex(dp),intent(out)::local(:,:)
 integer::x,y,a,raw(2),wrapped(2),cells(2),n
 complex(dp)::ph
 n=2*radius+1;local=0;a=0
 do x=-radius,radius
 do y=-radius,radius
 raw=site+[x,y];wrapped=modulo(raw,L);cells=(raw-wrapped)/L
 ph=exp(2*pi*zi*dot_product(real(cells,dp),k(:2)))
 local(a*60+1:(a+1)*60,:)=ph*v((wrapped(1)*L+wrapped(2))*60+1:(wrapped(1)*L+wrapped(2)+1)*60,:)
 a=a+1
 enddo
 enddo
 end subroutine
 subroutine tensor_add(h,L,k,site,radius,ke,spin)
 integer,intent(in)::L,site(2),radius
 real(dp),intent(in)::k(3),spin(3)
 complex(dp),intent(in)::ke(:,:,:)
 complex(dp),intent(inout)::h(:,:)
 integer::x,y,xx,yy,a,b,n,raw(2),wrapped(2),cells(2),ids((2*radius+1)**2),i,j,c
 complex(dp)::ph((2*radius+1)**2),phase
 n=2*radius+1;a=0
 do x=-radius,radius
 do y=-radius,radius
 a=a+1;raw=site+[x,y];wrapped=modulo(raw,L);cells=(raw-wrapped)/L
 ids(a)=wrapped(1)*L+wrapped(2);ph(a)=exp(2*pi*zi*dot_product(real(cells,dp),k(:2)))
 enddo
 enddo
 do a=1,n*n
 do b=1,n*n
 phase=conjg(ph(a))*ph(b)
 do c=1,3
 do j=1,60
 do i=1,60
 h(ids(a)*60+i,ids(b)*60+j)=h(ids(a)*60+i,ids(b)*60+j)+ &
 phase*spin(c)*ke((a-1)*60+i,(b-1)*60+j,c)
 enddo
 enddo
 enddo
 enddo
 enddo
 end subroutine
 subroutine tensor_band(L,k,site,radius,ke,dspin,v,der)
 integer,intent(in)::L,site(2),radius
 real(dp),intent(in)::k(3),dspin(3)
 complex(dp),intent(in)::ke(:,:,:),v(:,:)
 real(dp),intent(out)::der(:)
 complex(dp),allocatable::local(:,:),tmp(:,:)
 integer::c
 allocate(local(size(ke,1),size(v,2)),tmp(size(ke,1),size(v,2)))
 call tensor_embedding(L,k,site,radius,v,local);der=0
 do c=1,3
 call gemm(ke(:,:,c),local,tmp)
 der=der+dspin(c)*real(sum(conjg(local)*tmp,dim=1),dp)
 enddo
 end subroutine
end module
