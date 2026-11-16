module joint_model
 use joint_features
 implicit none
 type::joint_potential
  integer::grid(3),atom_ids(4)=[5,6,7,8],training_gate=0
  real(dp)::training_ueff=1,rcond=1e-10_dp
  logical::removed_scalar_means=.false.
  real(dp),allocatable::scalar(:,:),magnetic(:,:)
  real(dp)::nscalar(8,15)=0,nchi(4,10)=0,nvector(12,15)=0,ncross(6,10)=0
 end type
contains
 subroutine new_joint(m,grid)
 type(joint_potential),intent(out)::m
 integer,intent(in)::grid(3)
 type(crystal_action)::g(64)
 call require(all(grid>0).and.all(grid<=1000000),'Invalid joint grid dimensions')
 call crystal_group(grid,g);m%grid=grid
 call require(product(int(grid,int64))<=10000000_int64,'Joint grid exceeds supported10M points')
 allocate(m%scalar(8,product(grid)),m%magnetic(40,product(grid)))
 m%scalar=0;m%magnetic=0
 end subroutine
 subroutine channels_to_n(q,mag,n)
 complex(dp),intent(in)::q(5,5),mag(5,5,3)
 complex(dp),intent(out)::n(10,10)
 n(:5,:5)=(q+mag(:,:,3))/2;n(6:,6:)=(q-mag(:,:,3))/2
 n(:5,6:)=(mag(:,:,1)-zi*mag(:,:,2))/2;n(6:,:5)=(mag(:,:,1)+zi*mag(:,:,2))/2
 end subroutine
 subroutine n_to_channels(n,q,mag)
 complex(dp),intent(in)::n(10,10)
 complex(dp),intent(out)::q(5,5),mag(5,5,3)
 q=n(:5,:5)+n(6:,6:);mag(:,:,1)=n(:5,6:)+n(6:,:5)
 mag(:,:,2)=(n(6:,:5)-n(:5,6:))/zi;mag(:,:,3)=n(:5,:5)-n(6:,6:)
 end subroutine
 subroutine occupation_value(m,s,v,x,chi,n)
 type(joint_potential),intent(in)::m
 real(dp),intent(in)::s(8),v(3,12),x(3,6),chi(4)
 complex(dp),intent(out)::n(10,10)
 real(dp)::sym(5,5,15),anti(5,5,10),qs(15),qc(10),ms(3,15),mc(3,10)
 complex(dp)::q(5,5),mag(5,5,3)
 integer::i,a
 call orbital_bases(sym,anti)
 qs=matmul(s,m%nscalar);qc=matmul(chi,m%nchi)
 ms=matmul(v,m%nvector);mc=matmul(x,m%ncross);q=0;mag=0
 do i=1,15
  q=q+qs(i)*sym(:,:,i)
  do a=1,3
   mag(:,:,a)=mag(:,:,a)+ms(a,i)*sym(:,:,i)
  enddo
 enddo
 do i=1,10
  q=q+zi*qc(i)*anti(:,:,i)
  do a=1,3
   mag(:,:,a)=mag(:,:,a)+zi*mc(a,i)*anti(:,:,i)
  enddo
 enddo
 call channels_to_n(q,mag,n)
 end subroutine
 subroutine joint_value(m,e,de,field,n,dfield,dn)
 type(joint_potential),intent(in)::m
 real(dp),intent(in)::e(3,4),de(3,4)
 real(dp),intent(out)::field(:,:),dfield(:,:)
 complex(dp),intent(out)::n(10,10,4),dn(10,10,4)
 real(dp)::s(8),b(3,40),v(3,12),x(3,6),chi(4),ds(8),db(3,40),dv(3,12),dx(3,6),dc(4)
 call require(size(field,1)==product(m%grid).and.size(field,2)==4,'Joint field output shape mismatch')
 call require(all(ieee_is_finite(e)).and.all(ieee_is_finite(de)),'Nonfinite joint spin input')
 call require(maxval(abs(sum(e**2,dim=1)-1))<1e-9_dp,'Joint spin inputs must be unit directions')
 call spin_features(e,de,s,b,v,x,chi,ds,db,dv,dx,dc)
 field(:,1)=matmul(s,m%scalar);field(:,2:4)=transpose(matmul(b,m%magnetic))
 dfield(:,1)=matmul(ds,m%scalar);dfield(:,2:4)=transpose(matmul(db,m%magnetic))
 call require(all(ieee_is_finite(field)).and.all(ieee_is_finite(dfield)),'Nonfinite joint field/derivative')
 call joint_occupations(m,e,de,n,dn)
 end subroutine
 subroutine joint_occupations(m,e,de,n,dn)
 type(joint_potential),intent(in)::m
 real(dp),intent(in)::e(3,4),de(3,4)
 complex(dp),intent(out)::n(10,10,4),dn(10,10,4)
 type(crystal_action)::group(64)
 real(dp)::s(8),b(3,40),v(3,12),x(3,6),chi(4),ds(8),db(3,40),dv(3,12),dx(3,6),dc(4)
 complex(dp)::base(10,10),der(10,10),d(10,10)
 integer::g,site
 call crystal_group(m%grid,group);n=0;dn=0
 do g=1,64
  call spin_features(e(:,group(g)%perm),de(:,group(g)%perm),s,b,v,x,chi,ds,db,dv,dx,dc)
  call occupation_value(m,s,v,x,chi,base);call occupation_value(m,ds,dv,dx,dc,der)
  d=0;d(:5,:5)=group(g)%orbital;d(6:,6:)=group(g)%orbital
  site=group(g)%perm(1)
  ! Average the16 stabilizer actions; all four sites respond to any changed input spin.
  n(:,:,site)=n(:,:,site)+matmul(d,matmul(base,transpose(d)))/16
  dn(:,:,site)=dn(:,:,site)+matmul(d,matmul(der,transpose(d)))/16
 enddo
 end subroutine
 subroutine hubbard_value(n,dn,ueff,vu,dvu)
 complex(dp),intent(in)::n(10,10,4),dn(10,10,4)
 real(dp),intent(in)::ueff
 complex(dp),intent(out)::vu(10,10,4),dvu(10,10,4)
 integer::i,s
 call require(ieee_is_finite(ueff).and.ueff>=0,'Ueff must be finite and nonnegative in eV')
 vu=-ueff*n;dvu=-ueff*dn
 do s=1,4
 do i=1,10
  vu(i,i,s)=vu(i,i,s)+ueff/2
 enddo
 enddo
 end subroutine
 subroutine occupation_bounds(n,lo,hi,herm)
 complex(dp),intent(in)::n(10,10,4)
 real(dp),intent(out)::lo,hi,herm
 complex(dp)::a(10,10)
 real(dp)::ev(10)
 integer::s
 lo=huge(1._dp);hi=-lo;herm=0
 do s=1,4
  herm=max(herm,herm_error(n(:,:,s)));a=n(:,:,s);call diagonalize(a,ev)
  lo=min(lo,minval(ev));hi=max(hi,maxval(ev))
 enddo
 end subroutine
 subroutine write_joint(path,m)
 character(*),intent(in)::path
 type(joint_potential),intent(in)::m
 integer::u,ios
 open(newunit=u,file=path,status='new',access='stream',form='unformatted',convert='little_endian',iostat=ios)
 call require(ios==0,'Joint model output exists or cannot be created: '//trim(path))
 write(u) 'LANIO2_JOINT_V1 ',int([16909060,1,m%grid,m%atom_ids,8,40,4],int32), &
  m%training_ueff,m%rcond,int([merge(1,0,m%removed_scalar_means),m%training_gate],int32), &
  m%scalar,m%magnetic,m%nscalar,m%nchi,m%nvector,m%ncross
 close(u)
 end subroutine
 subroutine read_joint(path,m)
 character(*),intent(in)::path
 type(joint_potential),intent(out)::m
 character(16)::magic
 integer(int32)::head(12),flags(2)
 integer(int64)::bytes
 real(dp)::ueff,rcond
 integer::u,ios
 open(newunit=u,file=path,status='old',access='stream',form='unformatted',convert='little_endian',iostat=ios)
 call require(ios==0,'Cannot read joint model: '//trim(path))
 read(u,iostat=ios) magic,head,ueff,rcond,flags
 call require(ios==0.and.magic=='LANIO2_JOINT_V1 ','Invalid joint model header')
 call require(all(head(:2)==[16909060,1]).and.all(head(10:)==[8,40,4]),'Joint model endian/schema mismatch')
 call require(all(head(3:5)>0).and.all(head(3:5)<=1000000),'Invalid joint model grid')
 inquire(unit=u,size=bytes)
 call require(bytes==88_int64+8_int64*(48_int64*product(int(head(3:5),int64))+400),'Joint model byte count mismatch')
 call require(all(flags==0.or.flags==1),'Invalid joint model flags')
 call require(ieee_is_finite(ueff).and.ueff>=0.and.ieee_is_finite(rcond).and.rcond>0.and.rcond<1,'Invalid joint model controls')
 call new_joint(m,int(head(3:5)))
 m%atom_ids=head(6:9);m%training_ueff=ueff;m%rcond=rcond
 m%removed_scalar_means=flags(1)==1;m%training_gate=flags(2)
 call require(all(m%atom_ids>0),'Invalid model atom IDs')
 do ios=1,4
  call require(count(m%atom_ids==m%atom_ids(ios))==1,'Duplicate model atom IDs')
 enddo
 read(u,iostat=ios) m%scalar,m%magnetic,m%nscalar,m%nchi,m%nvector,m%ncross;close(u)
 call require(ios==0,'Truncated joint model')
 call require(all(ieee_is_finite(m%scalar)).and.all(ieee_is_finite(m%magnetic)).and. &
  all(ieee_is_finite(m%nscalar)).and.all(ieee_is_finite(m%nchi)).and. &
  all(ieee_is_finite(m%nvector)).and.all(ieee_is_finite(m%ncross)),'Nonfinite joint coefficients')
 call require(m%training_gate==1,'Joint model failed its training occupation gate; inspect fit report')
 end subroutine
end module
