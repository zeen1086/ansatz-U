module joint_rvr
 use joint_model
 use local_ansatz, only: frozen_model,read_model,background_projectors,project_at,pauli
 use iso_fortran_env, only: int8
 implicit none
 integer,parameter::nw_joint=120,nd_joint=240,nkernel=49
 character(16),parameter::rvr_magic='LANIO2_RVR_V1   ',row_magic='LANIO2_ROW_V1   '
 character(4096)::model_path='',frozen_path='',reference_potential_path='',wan_dir='', &
  translations_path='',support_path='',rvr_dir='',primitive_translations_path=''
 character(16)::wannier_layout='four_cell120'
 character(128)::gauge_id='',reference_scalar_id=''
 character(32)::background=''
 real(dp)::grid_weight=1._dp
 namelist /joint_projection/ model_path,frozen_path,reference_potential_path,wan_dir, &
  translations_path,support_path,rvr_dir,gauge_id,reference_scalar_id,background,grid_weight, &
  wannier_layout,primitive_translations_path
 type::projected_joint
  type(joint_potential)::m
  type(frozen_model)::bg
  integer,allocatable::rr(:,:)
  complex(dp),allocatable::kernel(:,:,:,:)
  integer(int64)::sources(4)
  real(dp)::weight
  character(128)::gauge,reference
  character(4096)::paths(4)
 end type
contains
 subroutine projection_control(path)
 character(*),intent(in)::path
 integer::u,ios
 character(1024)::message
 open(newunit=u,file=path,status='old',action='read',iostat=ios)
 call require(ios==0,'Cannot open joint projection namelist')
 read(u,nml=joint_projection,iostat=ios,iomsg=message);close(u)
 call require(ios==0,'Invalid joint_projection: '//trim(message))
 call require(min(len_trim(model_path),len_trim(frozen_path),len_trim(reference_potential_path), &
  len_trim(wan_dir),len_trim(translations_path),len_trim(rvr_dir),len_trim(gauge_id),len_trim(reference_scalar_id))>0, &
  'Declare model, frozen background, actual fixed scalar reference, Wannier/translation/output paths and source IDs')
 call require(background=='hc_plus_fixed_Hv_no_HU', &
  'Declare background=hc_plus_fixed_Hv_no_HU; supply its actual scalar reference (not a new fitted reference)')
 call require(ieee_is_finite(grid_weight).and.grid_weight>0,'Invalid explicit grid_weight')
 call require(wannier_layout=='four_cell120'.or.wannier_layout=='primitive30','Unknown Wannier layout')
 if(wannier_layout=='primitive30') call require(len_trim(primitive_translations_path)>0, &
  'primitive30 requires the explicit original primitive translation list')
 end subroutine
 function fingerprint_file(path) result(fingerprint)
 character(*),intent(in)::path
 integer(int64)::fingerprint,h1,h2,word,bytes,offset
 integer(int8)::block(65536)
 integer::u,ios,i,n
 open(newunit=u,file=path,status='old',action='read',access='stream',form='unformatted',iostat=ios)
 call require(ios==0,'Cannot identify input: '//trim(path))
 inquire(unit=u,size=bytes);h1=0;h2=0;offset=0
 do while(offset<bytes)
  n=int(min(65536_int64,bytes-offset));read(u,iostat=ios) block(:n)
  call require(ios==0,'Input changed/truncated while identifying: '//trim(path))
  do i=1,n
   word=iand(int(block(i),int64),255_int64)
   h1=modulo(h1*257+word,2147483647_int64)
   h2=modulo(h2*263+word,2147483629_int64)
  enddo
  offset=offset+n
 enddo
 close(u)
 fingerprint=h1*2147483648_int64+h2
 end function
 subroutine projection_sources(control,t,coeff,wr,delta,map)
 character(*),intent(in)::control
 type(projected_joint),intent(out)::t
 real(dp),allocatable,intent(out)::coeff(:,:)
 integer,allocatable,intent(out)::wr(:,:),delta(:,:),map(:,:)
 integer::u,ios,nr,i,j,n,k,grid(3)
 integer(int64)::reference_bytes
 integer,allocatable::work(:,:),inverse(:)
 real(dp),allocatable::ref(:,:)
 call projection_control(control)
 t%sources=[fingerprint_file(trim(model_path)),fingerprint_file(trim(frozen_path)), &
  fingerprint_file(trim(reference_potential_path)),fingerprint_file(control)]
 call read_joint(trim(model_path),t%m);call read_model(trim(frozen_path),t%bg)
 call require(all(t%m%atom_ids==[5,6,7,8]),'Frozen projector bridge requires original DFT Ni order 5,6,7,8')
 open(newunit=u,file=trim(translations_path),status='old',action='read',iostat=ios)
 call require(ios==0,'Cannot open explicit Wannier translation list')
 read(u,*,iostat=ios) nr
 call require(ios==0.and.nr>=1.and.nr<=125,'Wannier translations require 1..125 entries')
 allocate(wr(3,nr),inverse(nr),map(nr,nr),work(3,nr*nr))
 do i=1,nr
  read(u,*,iostat=ios) wr(:,i);call require(ios==0,'Truncated Wannier translations')
 enddo
 close(u)
 do i=1,nr
  inverse(i)=0
  do j=1,nr
   if(j/=i) call require(any(wr(:,i)/=wr(:,j)),'Duplicate Wannier translation')
   if(all(wr(:,j)==-wr(:,i))) inverse(i)=j
  enddo
  call require(inverse(i)>0,'Wannier support must include inverse translations')
 enddo
 n=0
 do j=1,nr
 do i=1,nr
  k=1
  do while(k<=n)
   if(all(work(:,k)==wr(:,j)-wr(:,i))) exit
   k=k+1
  enddo
  if(k>n) then
   n=n+1;work(:,n)=wr(:,j)-wr(:,i)
  endif
  map(i,j)=k
 enddo
 enddo
 call require(n<=4096,'Too many RVR difference translations')
 delta=work(:,:n);t%rr=delta
 allocate(coeff(nkernel,product(t%m%grid)),ref(product(t%m%grid),4))
 coeff(:8,:)=t%m%scalar;coeff(9:48,:)=t%m%magnetic
 open(newunit=u,file=trim(reference_potential_path),status='old',action='read',form='unformatted', &
  convert='little_endian',iostat=ios)
 call require(ios==0,'Cannot read actual fixed scalar reference V_HYC_dp (Ry)')
 inquire(unit=u,size=reference_bytes)
 call require(reference_bytes==28_int64+32_int64*product(t%m%grid),'Reference requires two float64/int32 records, four fields in Ry')
 read(u,iostat=ios) grid
 call require(ios==0.and.all(grid==t%m%grid),'Reference scalar grid/order mismatch')
 read(u,iostat=ios) ref;close(u)
 call require(ios==0.and.all(ieee_is_finite(ref)),'Invalid float64 reference V_HYC_dp')
 coeff(49,:)=ref(:,1)*ry_ev
 t%weight=grid_weight;t%gauge=gauge_id;t%reference=reference_scalar_id
 t%paths=[model_path,frozen_path,reference_potential_path,wan_dir]
 end subroutine
 subroutine read_raw_wannier(orb,nr,ng,raw,support)
 integer,intent(in)::orb,nr,ng
 complex(dp),allocatable,intent(out)::raw(:,:)
 integer,allocatable,intent(out)::support(:)
 integer::u,ios
 integer(int64)::bytes
 character(32)::number
 character(4096)::path
 allocate(raw(nr,ng),support(nr));support=1
 write(number,'(i0)') orb;path=trim(wan_dir)//'/wan_'//trim(number)
 open(newunit=u,file=trim(path),status='old',action='read',access='stream',form='unformatted', &
  convert='little_endian',iostat=ios)
 call require(ios==0,'Missing common-gauge spatial Wannier: '//trim(path))
 inquire(unit=u,size=bytes)
 call require(bytes==16_int64*nr*ng,'Wannier byte count mismatch; require complex128 (nR,ngrid,1)')
 read(u,iostat=ios) raw;close(u)
 call require(ios==0.and.all(ieee_is_finite(real(raw,dp))).and.all(ieee_is_finite(aimag(raw))), &
  'Truncated or nonfinite Wannier data')
 if(len_trim(support_path)>0) then
  open(newunit=u,file=trim(support_path)//'/cal_TF_'//trim(number),status='old',action='read',iostat=ios)
  call require(ios==0,'Missing explicit Wannier support mask')
  read(u,*,iostat=ios) support;close(u)
  call require(ios==0.and.all(support==0.or.support==1),'Support mask must contain one 0/1 per translation')
 endif
 end subroutine
 subroutine read_wannier(orb,wr,grid,w)
 integer,intent(in)::orb,wr(:,:),grid(3)
 complex(dp),allocatable,intent(out)::w(:,:)
 complex(dp),allocatable::raw(:,:)
 integer,allocatable::support(:),primitive_rr(:,:)
 integer::nr,ng,i,j,u,ios,np,g,point(3),local(3),cell(3),centre(3),p(3),target(3),candidate(3),ix,iy
 ng=product(grid);nr=size(wr,2);allocate(w(ng,nr));w=0
 if(wannier_layout=='four_cell120') then
  call read_raw_wannier(orb,nr,ng,raw,support)
  do i=1,nr
  do j=1,nr
   if(any(wr(:,j)/=-wr(:,i))) cycle
   w(:,i)=raw(j,:)*support(j)
   exit
  enddo
  enddo
  return
 endif
 p=[grid(1)/2,grid(2)/2,grid(3)]
 open(newunit=u,file=trim(primitive_translations_path),status='old',action='read',iostat=ios)
 call require(ios==0,'Cannot read original primitive translations')
 read(u,*,iostat=ios) np;call require(ios==0.and.np>=1.and.np<=729,'Invalid primitive translation count')
 allocate(primitive_rr(3,np))
 do i=1,np
  read(u,*,iostat=ios) primitive_rr(:,i);call require(ios==0,'Truncated primitive translations')
 enddo
 close(u)
 do i=1,np
 do j=1,i-1
  call require(any(primitive_rr(:,i)/=primitive_rr(:,j)),'Duplicate primitive translation')
 enddo
 enddo
 call read_raw_wannier(mod(orb-1,30)+1,np,product(p),raw,support)
 centre=[mod((orb-1)/30,2),(orb-1)/60,0]
 ! Relabel original Wannier translations; no interpolation, fitting or extra normalization.
 do i=1,np
 if(support(i)==0) cycle
 do ix=0,1
 do iy=0,1
  candidate=[ix,iy,0]-centre-primitive_rr(:,i)
  if(any(modulo(candidate,[2,2,1])/=0)) cycle
  candidate=candidate/[2,2,1]
  call require(any([(all(wr(:,j)==candidate),j=1,nr)]), &
   'Output translations truncate existing primitive Wannier support; enlarge the explicit list')
 enddo
 enddo
 enddo
 do g=1,ng
  point=[mod(g-1,grid(1)),mod((g-1)/grid(1),grid(2)),(g-1)/(grid(1)*grid(2))]
  cell=point/p;local=modulo(point,p);j=1+local(1)+p(1)*local(2)+p(1)*p(2)*local(3)
  do i=1,nr
   target=cell-centre-[2,2,1]*wr(:,i)
   do u=1,np
    if(any(primitive_rr(:,u)/=target)) cycle
    w(g,i)=raw(u,j)*support(u);exit
   enddo
  enddo
 enddo
 end subroutine
 subroutine project_joint_row(control,column)
 character(*),intent(in)::control
 integer,intent(in)::column
 type(projected_joint)::t
 real(dp),allocatable::coeff(:,:)
 integer,allocatable::wr(:,:),rr(:,:),map(:,:)
 complex(dp),allocatable::wi(:,:),wj(:,:),weighted(:,:),mat(:,:),row(:,:,:)
 integer::j,c,a,b,u,ios,nr,ng,first,last,iwan
 character(32)::number
 character(4096)::path
 call require(column>=0.and.column<=nw_joint,'Joint RVR column must be 1..120 spatial orbitals (0 selects all)')
 call projection_sources(control,t,coeff,wr,rr,map)
 first=column;last=column
 if(column==0) then
  first=1;last=nw_joint
 endif
 do iwan=first,last
  write(number,'(i0)') iwan;call require_new(trim(rvr_dir)//'/tmp/B_general_'//trim(number))
 enddo
 nr=size(wr,2);ng=size(coeff,2)
 allocate(row(nkernel,size(rr,2),nw_joint),weighted(ng,nr),mat(nr,nr))
 do iwan=first,last
 row=0;call read_wannier(iwan,wr,t%m%grid,wi)
 do j=1,nw_joint
  call read_wannier(j,wr,t%m%grid,wj)
  do c=1,nkernel
   do a=1,nr
    weighted(:,a)=wj(:,a)*coeff(c,:)*grid_weight
   enddo
   call gemm(wi,weighted,mat,ta='C')
   do b=1,nr
   do a=1,nr
    row(c,map(a,b),j)=row(c,map(a,b),j)+mat(a,b)
   enddo
   enddo
  enddo
 enddo
 write(number,'(i0)') iwan;path=trim(rvr_dir)//'/tmp/B_general_'//trim(number)
 open(newunit=u,file=trim(path),status='new',access='stream',form='unformatted', &
  convert='little_endian',iostat=ios)
 call require(ios==0,'Cannot create RVR row (create external tmp directory first)')
 write(u) row_magic,int([16909060,1,iwan,nw_joint,size(rr,2),nkernel],int32),t%sources,int(rr,int32),row
 close(u)
 print '(a,i0,a)', 'Projected spatial orbital ',iwan,': scalar8, B40 and actual fixed scalar reference (eV)'
 enddo
 end subroutine
 subroutine require_new(path)
 character(*),intent(in)::path
 logical::exists
 inquire(file=path,exist=exists);call require(.not.exists,'Output exists; use a new run: '//path)
 end subroutine
 subroutine combine_joint(control)
 character(*),intent(in)::control
 type(projected_joint)::t
 real(dp),allocatable::coeff(:,:)
 integer,allocatable::wr(:,:),rr(:,:),map(:,:)
 integer(int32),allocatable::check_rr(:,:)
 complex(dp),allocatable::row(:,:,:)
 integer(int32)::header(6)
 integer(int64)::sources(4),bytes,expected
 integer::i,j,u,ios,c,r,inverse
 character(16)::magic
 character(32)::number
 character(4096)::path
 call projection_sources(control,t,coeff,wr,rr,map)
 call require_new(trim(rvr_dir)//'/B_general')
 allocate(t%kernel(nw_joint,nw_joint,size(rr,2),nkernel),row(nkernel,size(rr,2),nw_joint),check_rr(3,size(rr,2)))
 do i=1,nw_joint
  write(number,'(i0)') i;path=trim(rvr_dir)//'/tmp/B_general_'//trim(number)
  open(newunit=u,file=trim(path),status='old',action='read',access='stream',form='unformatted',convert='little_endian',iostat=ios)
  call require(ios==0,'Missing projected row: '//trim(path))
  inquire(unit=u,size=bytes);expected=72_int64+12_int64*size(rr,2)+16_int64*nkernel*size(rr,2)*nw_joint
  call require(bytes==expected,'Projected row byte count mismatch')
  read(u,iostat=ios) magic,header,sources,check_rr,row;close(u)
  call require(ios==0.and.magic==row_magic,'Invalid projected row header')
  call require(all(header==[16909060,1,i,nw_joint,size(rr,2),nkernel]),'Projected row dimensions/version/order mismatch')
  call require(all(sources==t%sources).and.all(check_rr==rr),'Mixed model/background/reference/control in projected rows')
  do j=1,nw_joint
  do r=1,size(rr,2)
   t%kernel(i,j,r,:)=row(:,r,j)
  enddo
  enddo
 enddo
 do r=1,size(rr,2)
  inverse=0
  do j=1,size(rr,2)
   if(all(rr(:,j)==-rr(:,r))) inverse=j
  enddo
  call require(inverse>0,'RVR missing inverse difference translation')
  do c=1,nkernel
   call require(maxval(abs(t%kernel(:,:,r,c)-transpose(conjg(t%kernel(:,:,inverse,c)))))<1e-9_dp, &
    'RVR violates Hermiticity; check immutable Wannier inputs and row completeness')
  enddo
 enddo
 call write_projected(trim(rvr_dir)//'/B_general',t)
 print '(a)', 'Combined B_general: scalar8 + B40 + scalar reference, original background/projectors and full N coefficients'
 end subroutine
 subroutine write_projected(path,t)
 character(*),intent(in)::path
 type(projected_joint),intent(in)::t
 integer::u,ios
 open(newunit=u,file=path,status='new',access='stream',form='unformatted',convert='little_endian',iostat=ios)
 call require(ios==0,'Cannot create new B_general')
 write(u) rvr_magic,int([16909060,1,nw_joint,size(t%rr,2),t%bg%nr,t%m%grid,t%m%atom_ids,8,40,1,0],int32), &
  t%sources,t%weight,t%gauge,t%reference,t%paths,t%m%training_ueff,t%m%rcond, &
  int([merge(1,0,t%m%removed_scalar_means),t%m%training_gate],int32),int(t%rr,int32), &
  t%bg%rr,t%bg%hc,t%bg%p,t%m%nscalar,t%m%nchi,t%m%nvector,t%m%ncross,t%kernel
 close(u)
 end subroutine
 subroutine read_projected(path,t)
 character(*),intent(in)::path
 type(projected_joint),intent(out)::t
 character(16)::magic
 integer(int32)::head(16),flags(2)
 integer(int64)::bytes,position,expected
 integer::u,ios,nr
 type(crystal_action)::group(64)
 open(newunit=u,file=path,status='old',action='read',access='stream',form='unformatted',convert='little_endian',iostat=ios)
 call require(ios==0,'Cannot read explicit joint B_general')
 read(u,iostat=ios) magic,head,t%sources,t%weight,t%gauge,t%reference,t%paths, &
  t%m%training_ueff,t%m%rcond,flags
 call require(ios==0.and.magic==rvr_magic,'Wrong B_general magic/version; legacy and joint formats differ')
 call require(all(head(:3)==[16909060,1,nw_joint]).and.head(5)==729.and. &
  all(head(13:)==[8,40,1,0]),'Joint RVR dimensions/endian/units/background/HU contract mismatch')
 nr=head(4);call require(nr>=1.and.nr<=4096,'Invalid RVR translation count')
 t%m%grid=head(6:8);t%m%atom_ids=head(9:12);t%bg%nr=head(5)
 call crystal_group(t%m%grid,group)
 call require(all(t%m%atom_ids==[5,6,7,8]).and.all(flags==0.or.flags==1).and.flags(2)==1, &
  'Invalid Ni order/occupation gate in joint RVR')
 call require(ieee_is_finite(t%weight).and.t%weight>0.and.len_trim(t%gauge)>0.and.len_trim(t%reference)>0, &
  'Missing joint RVR normalization/gauge/reference declaration')
 call require(ieee_is_finite(t%m%training_ueff).and.t%m%training_ueff>=0.and. &
  ieee_is_finite(t%m%rcond).and.t%m%rcond>0.and.t%m%rcond<1,'Invalid RVR model controls')
 t%m%removed_scalar_means=flags(1)==1;t%m%training_gate=flags(2)
 inquire(unit=u,pos=position,size=bytes)
 expected=position-1+12_int64*(nr+t%bg%nr)+16_int64*(1050_int64*t%bg%nr+ &
  int(nw_joint,int64)**2*nr*nkernel)+3200_int64
 call require(bytes==expected,'Joint B_general byte count mismatch')
 allocate(t%rr(3,nr),t%bg%rr(3,t%bg%nr),t%bg%hc(30,30,t%bg%nr),t%bg%p(5,30,t%bg%nr), &
  t%kernel(nw_joint,nw_joint,nr,nkernel))
 read(u,iostat=ios) t%rr,t%bg%rr,t%bg%hc,t%bg%p,t%m%nscalar,t%m%nchi,t%m%nvector,t%m%ncross,t%kernel
 close(u);call require(ios==0,'Truncated joint RVR data')
 call require(all(ieee_is_finite(real(t%kernel,dp))).and.all(ieee_is_finite(aimag(t%kernel))).and. &
  all(ieee_is_finite(real(t%bg%hc,dp))).and.all(ieee_is_finite(aimag(t%bg%hc))).and. &
  all(ieee_is_finite(real(t%bg%p,dp))).and.all(ieee_is_finite(aimag(t%bg%p))).and. &
  all(ieee_is_finite(t%m%nscalar)).and.all(ieee_is_finite(t%m%nchi)).and. &
  all(ieee_is_finite(t%m%nvector)).and.all(ieee_is_finite(t%m%ncross)),'Nonfinite joint RVR contents')
 end subroutine
 subroutine rvr_at(t,k,kernel,hbg,a)
 type(projected_joint),intent(in)::t
 real(dp),intent(in)::k(3)
 complex(dp),intent(out)::kernel(nw_joint,nw_joint,nkernel),hbg(nd_joint,nd_joint),a(10,nd_joint,4)
 complex(dp),allocatable::old(:,:)
 complex(dp)::pg(5,30,0:1,0:1),ap(10,nd_joint),phase
 integer::r,i,s,permutation(nd_joint),cell,oldcell
 allocate(old(nd_joint,nd_joint))
 kernel=0
 do r=1,size(t%rr,2)
  phase=exp(2*pi*zi*dot_product(k,real(t%rr(:,r),dp)));kernel=kernel+phase*t%kernel(:,:,r,:)
 enddo
 do i=1,nkernel
  call require(herm_error(kernel(:,:,i))<1e-8_dp,'Non-Hermitian interpolated scalar/B/reference RVR')
 enddo
 call background_projectors(t%bg,2,k,old,pg)
 do cell=0,3
  oldcell=mod(cell,2)*2+cell/2
  do i=1,60
   permutation(cell*60+i)=oldcell*60+i
  enddo
 enddo
 hbg=old(permutation,permutation)
 do s=1,4
  call project_at(2,k,pg,[mod(s-1,2),(s-1)/2],ap)
  a(:,:,s)=ap(:,permutation)
 enddo
 end subroutine
 subroutine projected_value(t,kernel,hbg,a,e,de,ueff,hv,hb,hu,dv,db,du,n,dn)
 type(projected_joint),intent(in)::t
 complex(dp),intent(in)::kernel(nw_joint,nw_joint,nkernel),hbg(nd_joint,nd_joint),a(10,nd_joint,4)
 real(dp),intent(in)::e(3,4),de(3,4),ueff
 complex(dp),intent(out)::hv(nd_joint,nd_joint),hb(nd_joint,nd_joint),hu(nd_joint,nd_joint), &
  dv(nd_joint,nd_joint),db(nd_joint,nd_joint),du(nd_joint,nd_joint),n(10,10,4),dn(10,10,4)
 real(dp)::s(8),b(3,40),v(3,12),x(3,6),chi(4),ds(8),dbf(3,40),dfv(3,12),dx(3,6),dc(4),lo,hi,herm
 complex(dp)::vu(10,10,4),dvu(10,10,4),temp(10,nd_joint),cs,cb,vs(2,2),bs(2,2)
 integer::i,j,c,alpha,beta,ii,jj,site
 call require(all(ieee_is_finite(e)).and.all(ieee_is_finite(de)).and. &
  maxval(abs(sum(e**2,dim=1)-1))<1e-12_dp,'Joint RVR requires four finite unit spin directions')
 call spin_features(e,de,s,b,v,x,chi,ds,dbf,dfv,dx,dc)
 hv=0;hb=0;hu=0;dv=0;db=0;du=0
 do j=1,nw_joint
 do i=1,nw_joint
  cs=sum(s*kernel(i,j,1:8))-kernel(i,j,49)
  cb=sum(ds*kernel(i,j,1:8));vs=0;bs=0
  do c=1,40
   vs=vs+kernel(i,j,8+c)*pauli(b(:,c));bs=bs+kernel(i,j,8+c)*pauli(dbf(:,c))
  enddo
  do beta=1,2
  do alpha=1,2
   ii=2*(i-1)+alpha;jj=2*(j-1)+beta
   if(alpha==beta) then
    hv(ii,jj)=cs;dv(ii,jj)=cb
   endif
   hb(ii,jj)=vs(alpha,beta);db(ii,jj)=bs(alpha,beta)
  enddo
  enddo
 enddo
 enddo
 call joint_occupations(t%m,e,de,n,dn);call occupation_bounds(n,lo,hi,herm)
 call require(lo>=-1e-8_dp.and.hi<=1+1e-8_dp.and.herm<1e-10_dp,'Full-model N fails occupation/Hermiticity gate; no clipping')
 call hubbard_value(n,dn,ueff,vu,dvu)
 ! Any input spin can change N at all four Ni; every original A enters both sums.
 do site=1,4
  call gemm(vu(:,:,site),a(:,:,site),temp)
  call gemm(a(:,:,site),temp,hu,ta='C',beta=cmplx(1._dp,0._dp,dp))
  call gemm(dvu(:,:,site),a(:,:,site),temp)
  call gemm(a(:,:,site),temp,du,ta='C',beta=cmplx(1._dp,0._dp,dp))
 enddo
 call require(max(herm_error(hbg),herm_error(hv),herm_error(hb),herm_error(hu), &
  herm_error(dv),herm_error(db),herm_error(du))<1e-8_dp,'Non-Hermitian full-model component/derivative')
 end subroutine
end module
