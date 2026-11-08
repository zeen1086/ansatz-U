module joint_spin
 use joint_rvr
 implicit none
contains
 subroutine run_joint_spin(rvr_path,texture_path,kpoints_path,readouts_path,prefix,ueff,sigma, &
  electrons_per_ni,mode,delta,perturb_site,dump_matrices)
 character(*),intent(in)::rvr_path,texture_path,kpoints_path,readouts_path,prefix,mode
 real(dp),intent(in)::ueff,sigma,electrons_per_ni,delta
 integer,intent(in)::perturb_site(2)
 logical,intent(in)::dump_matrices
 type(projected_joint)::t
 real(dp)::e(3,4),axis(3,4),spin(3,4),de(3,4),angle,mu,count,band_sum,start,finish
 real(dp),allocatable::q(:,:),w(:),ev(:,:),f(:,:),der(:,:,:,:),tau(:,:,:)
 integer,allocatable::sites(:,:)
 complex(dp),allocatable::kernels(:,:,:),hbg(:,:),a(:,:,:),hv(:,:),hb(:,:),hu(:,:),h(:,:), &
  dv(:,:),db(:,:),du(:,:),vec(:,:),tmp(:,:)
 complex(dp)::n(10,10,4),dn(10,10,4)
 integer::u,ios,nk,nread,i,j,k,pass,passes,site,band
 character(64)::suffix
 character(4),parameter::components(5)=[character(4)::'hbg','hv','hb','hu','h']
 character(2),parameter::derivatives(3)=['dv','db','du']
 call cpu_time(start)
 call require(min(len_trim(rvr_path),len_trim(texture_path),len_trim(kpoints_path), &
  len_trim(readouts_path),len_trim(prefix))>0,'Joint run requires explicit RVR/texture/kpoints/readouts/output paths')
 call require(mode=='torque'.or.mode=='table3','mode must be torque or table3')
 call require(ieee_is_finite(ueff).and.ueff>=0.and.ieee_is_finite(sigma).and.sigma>0, &
  'Invalid Ueff or FD sigma in eV')
 call require(ieee_is_finite(electrons_per_ni).and.electrons_per_ni>0.and.electrons_per_ni<60, &
  'Invalid active electrons/Ni; fixed-number target is 4*electrons_per_ni')
 call require(all(perturb_site>=0).and.all(perturb_site<2).and.ieee_is_finite(delta),'Invalid perturbation site/angle')
 if(mode=='table3') call require(delta/=0,'Table3 requires nonzero delta in radians')
 call read_projected(rvr_path,t)
 open(newunit=u,file=texture_path,status='old',action='read',iostat=ios);call require(ios==0,'Cannot read joint texture')
 do i=1,4
  read(u,*,iostat=ios) e(:,i),axis(:,i);call require(ios==0,'Texture needs 4 rows in Ni5,6,7,8 order (x fastest)')
 enddo
 read(u,*,iostat=ios) angle;call require(ios<0,'Unexpected extra joint texture rows; exactly four independent Ni required')
 close(u)
 call require(all(ieee_is_finite(e)).and.all(ieee_is_finite(axis)),'Nonfinite joint texture')
 call require(maxval(abs(sum(e**2,dim=1)-1))<1e-12_dp.and.maxval(abs(sum(axis**2,dim=1)-1))<1e-12_dp, &
  'Joint texture directions and torque axes must be unit vectors; no extra amplitude factor')
 open(newunit=u,file=kpoints_path,status='old',action='read',iostat=ios);call require(ios==0,'Cannot read joint kpoints')
 read(u,*,iostat=ios) nk;call require(ios==0.and.nk>0.and.nk<=1000000,'Invalid joint kpoint count')
 allocate(q(3,nk),w(nk),ev(nd_joint,nk),f(nd_joint,nk))
 do k=1,nk
  read(u,*,iostat=ios) q(:,k),w(k);call require(ios==0,'Truncated joint kpoints: Kx Ky Kz weight')
 enddo
 close(u)
 call require(all(ieee_is_finite(q)).and.all(ieee_is_finite(w)).and.all(w>=0).and. &
  abs(sum(w)-1)<1e-12_dp,'Kpoints/weights must be finite, nonnegative and normalized')
 open(newunit=u,file=readouts_path,status='old',action='read',iostat=ios);call require(ios==0,'Cannot read joint readouts')
 read(u,*,iostat=ios) nread;call require(ios==0.and.nread>=1.and.nread<=4,'Use 1..4 independent joint readout sites')
 allocate(sites(2,nread),der(nd_joint,3,nread,nk),tau(3,nread,2))
 do i=1,nread
  read(u,*,iostat=ios) sites(:,i);call require(ios==0,'Readouts require zero-based x y')
 enddo
 close(u);call require(all(sites>=0).and.all(sites<2),'Joint readouts are independent sites of the original 2x2 cell')
 passes=1;if(mode=='table3') passes=2
 ! Refuse all collisions before creating the first band, matrix or torque output.
 do pass=1,passes
  write(suffix,'(a,i0)') '.p',pass
  call require_new(prefix//trim(suffix)//'.bands');call require_new(prefix//trim(suffix)//'.result')
  if(.not.dump_matrices) cycle
  do k=1,nk
   write(suffix,'(a,i0,a,i0)') '.p',pass,'.k',k
   do j=1,5
    call require_new(prefix//trim(suffix)//'.'//trim(components(j)))
   enddo
   do i=1,nread
    write(suffix,'(a,i0,a,i0,a,i0)') '.p',pass,'.k',k,'.r',i
    do j=1,3
     call require_new(prefix//trim(suffix)//'.'//derivatives(j))
    enddo
   enddo
  enddo
 enddo
 if(mode=='table3') call require_new(prefix//'.phi')
 allocate(kernels(nw_joint,nw_joint,nkernel),hbg(nd_joint,nd_joint),a(10,nd_joint,4), &
  hv(nd_joint,nd_joint),hb(nd_joint,nd_joint),hu(nd_joint,nd_joint),h(nd_joint,nd_joint), &
  dv(nd_joint,nd_joint),db(nd_joint,nd_joint),du(nd_joint,nd_joint),vec(nd_joint,nd_joint),tmp(nd_joint,nd_joint))
 do pass=1,passes
  angle=0;if(pass==2.or.mode=='torque') angle=delta
  spin=e;site=perturb_site(1)+2*perturb_site(2)+1
  spin(:,site)=cos(angle)*e(:,site)+sin(angle)*cross3(axis(:,site),e(:,site))+ &
   (1-cos(angle))*axis(:,site)*dot_product(axis(:,site),e(:,site))
  do k=1,nk
   call rvr_at(t,q(:,k),kernels,hbg,a);de=0
   call projected_value(t,kernels,hbg,a,spin,de,ueff,hv,hb,hu,dv,db,du,n,dn)
   h=hbg+hv+hb+hu;vec=h;call diagonalize(vec,ev(:,k))
   if(dump_matrices) then
    write(suffix,'(a,i0,a,i0)') '.p',pass,'.k',k
    call dump_joint(prefix//trim(suffix)//'.hbg',hbg);call dump_joint(prefix//trim(suffix)//'.hv',hv)
    call dump_joint(prefix//trim(suffix)//'.hb',hb);call dump_joint(prefix//trim(suffix)//'.hu',hu)
    call dump_joint(prefix//trim(suffix)//'.h',h)
   endif
   do i=1,nread
    site=sites(1,i)+2*sites(2,i)+1;de=0;de(:,site)=cross3(axis(:,site),spin(:,site))
    call projected_value(t,kernels,hbg,a,spin,de,ueff,hv,hb,hu,dv,db,du,n,dn)
    call gemm(dv,vec,tmp);der(:,1,i,k)=real(sum(conjg(vec)*tmp,dim=1),dp)
    call gemm(db,vec,tmp);der(:,2,i,k)=real(sum(conjg(vec)*tmp,dim=1),dp)
    call gemm(du,vec,tmp);der(:,3,i,k)=real(sum(conjg(vec)*tmp,dim=1),dp)
    if(dump_matrices) then
     write(suffix,'(a,i0,a,i0,a,i0)') '.p',pass,'.k',k,'.r',i
     call dump_joint(prefix//trim(suffix)//'.dv',dv);call dump_joint(prefix//trim(suffix)//'.db',db)
     call dump_joint(prefix//trim(suffix)//'.du',du)
    endif
   enddo
  enddo
  call filling(ev,w,4*electrons_per_ni,sigma,mu,f,count)
  tau(:,:,pass)=0;band_sum=0
  do k=1,nk
   band_sum=band_sum+w(k)*sum(f(:,k)*ev(:,k))
   do i=1,nread
   do j=1,3
    tau(j,i,pass)=tau(j,i,pass)-w(k)*sum(f(:,k)*der(:,j,i,k))
   enddo
   enddo
  enddo
  write(suffix,'(a,i0)') '.p',pass
  open(newunit=u,file=prefix//trim(suffix)//'.bands',status='new',iostat=ios)
  call require(ios==0,'Cannot create joint band output')
  write(u,'(a)') '# k band Kx Ky Kz E_eV FD_occupation weight; same full H for all three torques'
  do k=1,nk
  do band=1,nd_joint
   write(u,'(2i8,7es26.16e3)') k,band,q(:,k),ev(band,k),f(band,k),w(k)
  enddo
  enddo
  close(u)
  open(newunit=u,file=prefix//trim(suffix)//'.result',status='new',iostat=ios)
  call require(ios==0,'Cannot create joint result')
  write(u,'(a)') '# Ueff_eV sigma_eV active_electrons EF_eV band_sum_eV angle_rad'
  write(u,'(6es26.16e3)') ueff,sigma,count,mu,band_sum,angle
  write(u,'(a)') '# x y tau_v tau_B tau_U tau_total; eV/rad, same full H eigenstates and FD occupations'
  do i=1,nread
   write(u,'(2i8,4es26.16e3)') sites(:,i),tau(:,i,pass),sum(tau(:,i,pass))
  enddo
  write(u,'(a)') '# H=hc_R+[Hv(e)-Hv_ref]+HB(e)+sum_all_Ni A^dagger Ueff*(I/2-N(e)) A'
  write(u,'(a)') '# hv/dv dumps are scalar difference/derivative; hbg contains the actual fixed scalar reference, no HU'
  write(u,'(a)') '# active basis: (cell x fastest, spatial orbital, up/down); projector rows up d1..5,down d1..5'
  write(u,'(a)') '# gauge='//trim(t%gauge)//'; fixed_scalar_reference='//trim(t%reference)
  write(u,'(a,4(z16.16,1x))') '# source fingerprints (V_func_tmp,frozen,reference,control)=',t%sources
  write(u,'(a,es24.15)') '# grid_weight=',t%weight
  write(u,'(a)') '# band_sum=sum_k w Tr[f_FD(H)H]; not the complete DFT total energy; no U on/off interpretation'
  close(u)
 enddo
 if(mode=='table3') then
  open(newunit=u,file=prefix//'.phi',status='new',iostat=ios);call require(ios==0,'Cannot create joint Phi')
  write(u,'(a)') '# x y Phi_v Phi_B Phi_U Phi_total; meV/rad^2; [tau(delta)-tau(0)]/delta; no S^2'
  do i=1,nread
   write(u,'(2i8,4es26.16e3)') sites(:,i),(tau(:,i,2)-tau(:,i,1))*1000/delta, &
    sum(tau(:,i,2)-tau(:,i,1))*1000/delta
  enddo
  close(u)
 endif
 call cpu_time(finish);write(*,'(a,es14.6)') 'Completed projected three-potential Fortran readout, CPU seconds=',finish-start
 end subroutine
 subroutine dump_joint(path,a)
 character(*),intent(in)::path
 complex(dp),intent(in)::a(:,:)
 integer::u,ios
 open(newunit=u,file=path,status='new',access='stream',form='unformatted',convert='little_endian',iostat=ios)
 call require(ios==0,'Cannot create new joint matrix dump: '//path)
 write(u) a;close(u)
 end subroutine
end module
