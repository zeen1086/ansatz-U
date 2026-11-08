program run_ansatzU
 use local_ansatz
 use joint_spin, only: run_joint_spin
 implicit none
 type(frozen_model)::m
 type(tensor_kernel)::tk
 integer::tensor_range=1
 character(1024)::tensor_path=''
 character(2)::tensor_unit=''
 real(dp)::tensor_spin_scale=-1._dp
 complex(dp),allocatable::tke(:,:,:)
 integer::L=2,perturb_site(2)=[0,0],u,ios,nk,nread,nc,nd,i,j,x,y,k,pass,passes,t,band
 character(1024)::input='',model_path='',texture_path='',kpoints_path='',readouts_path='',output_prefix=''
 character(1024)::joint_rvr_path=''
 character(16)::mode='torque'
 real(dp)::ueff=1._dp,sigma=.05_dp,electrons_per_ni=21._dp,delta=0._dp,mu,count,band_energy,seconds,start
 logical::dump_matrices=.false.
 namelist /run/ L,model_path,texture_path,kpoints_path,readouts_path,output_prefix, &
 ueff,sigma,electrons_per_ni,mode,delta,perturb_site,dump_matrices,tensor_path,tensor_range,tensor_spin_scale,tensor_unit,joint_rvr_path
 real(dp),allocatable::e(:,:),axis(:,:),spin(:,:),q(:,:),w(:),ev(:,:),f(:,:),der(:,:,:,:),tau(:,:,:)
 integer,allocatable::sites(:,:)
 complex(dp),allocatable::hcv(:,:),hb(:,:),hu(:,:),h(:,:),vec(:,:),pg(:,:,:,:),ns(:,:,:),frames(:,:,:),db(:,:),du(:,:)
 complex(dp)::ke(270,270),n0(10,10),dn(10,10),rot(10,10),s2(2,2),g2(2,2)
 real(dp)::angle,ax(3),norm
 character(32)::suffix
 character(1024)::msg
 call cpu_time(start)
 call get_command_argument(1,input)
 call require(len_trim(input)>0,'Usage: ansatzU explicit-run.nml')
 open(newunit=u,file=trim(input),status='old',action='read',iostat=ios,iomsg=msg)
 call require(ios==0,'Cannot open run namelist: '//trim(msg))
 read(u,nml=run,iostat=ios,iomsg=msg);close(u)
 call require(ios==0,'Invalid run namelist: '//trim(msg))
 if(len_trim(joint_rvr_path)>0) then
  call require(L==2,'The nonrigid four-spin joint model supports its original 2x2 cell only')
  call require(len_trim(model_path)==0.and.len_trim(tensor_path)==0, &
   'Joint B_general already contains the background/projector; omit model_path and tensor_path')
  call run_joint_spin(trim(joint_rvr_path),trim(texture_path),trim(kpoints_path),trim(readouts_path), &
   trim(output_prefix),ueff,sigma,electrons_per_ni,trim(mode),delta,perturb_site,dump_matrices)
 else
 call require(L>0.and.L<=100,'L must be positive and <=100 (dense solver; memory is caller responsibility)')
 call require(len_trim(model_path)>0.and.len_trim(texture_path)>0.and.len_trim(kpoints_path)>0.and. &
 len_trim(readouts_path)>0.and.len_trim(output_prefix)>0,'All external paths and output_prefix must be explicit')
 call require(mode=='torque'.or.mode=='table3','mode must be torque or table3')
 call require(ieee_is_finite(ueff).and.ueff>=0.and.sigma>0.and.ieee_is_finite(sigma),'Invalid Ueff or FD sigma (eV)')
 call require(electrons_per_ni>0.and.electrons_per_ni<60,'Expected active electrons/Ni in (0,60)')
 call require(all(perturb_site>=0).and.all(perturb_site<L),'perturb_site outside cell')
 call require(ieee_is_finite(delta),'delta must be finite radians')
 if(mode=='table3') call require(delta/=0,'Table3 requires nonzero delta radians')
 nc=L*L;nd=60*nc
 allocate(e(3,nc),axis(3,nc),spin(3,nc),frames(10,10,nc))
 open(newunit=u,file=trim(texture_path),status='old',action='read',iostat=ios)
 call require(ios==0,'Cannot open explicit texture')
 do i=1,nc
 read(u,*,iostat=ios) e(:,i),axis(:,i)
 call require(ios==0,'Texture requires L*L rows: ex ey ez axis_x axis_y axis_z; x then y (y fastest)')
 call require(all(ieee_is_finite(e(:,i))).and.all(ieee_is_finite(axis(:,i))),'Nonfinite texture')
 call require(abs(sum(e(:,i)**2)-1._dp)<1e-12_dp,'Local spin direction must have unit norm; no implicit rescaling')
 call require(abs(sum(axis(:,i)**2)-1._dp)<1e-12_dp,'Torque rotation axis must have unit norm')
 enddo
 close(u)
 open(newunit=u,file=trim(kpoints_path),status='old',action='read',iostat=ios)
 call require(ios==0,'Cannot open explicit kpoint list')
 read(u,*,iostat=ios) nk
 call require(ios==0.and.nk>0,'kpoints: first row must be positive nk')
 allocate(q(3,nk),w(nk))
 do k=1,nk
 read(u,*,iostat=ios) q(:,k),w(k)
 call require(ios==0,'kpoints require Kx Ky Kz weight for each point')
 enddo
 close(u)
 call require(all(ieee_is_finite(q)).and.all(ieee_is_finite(w)),'Nonfinite kpoints/weights')
 call require(all(w>=0).and.abs(sum(w)-1._dp)<1e-12_dp,'k weights must be nonnegative and sum to one')
 open(newunit=u,file=trim(readouts_path),status='old',action='read',iostat=ios)
 call require(ios==0,'Cannot open explicit readout sites')
 read(u,*,iostat=ios) nread
 call require(ios==0.and.nread>0,'readouts: first row must be positive count')
 allocate(sites(2,nread))
 do i=1,nread
 read(u,*,iostat=ios) sites(:,i)
 call require(ios==0,'Readout sites require x y')
 enddo
 close(u)
 call require(all(sites>=0).and.all(sites<L),'Readout outside cell; use an independent periodic site')
 call read_model(trim(model_path),m)
 if(len_trim(tensor_path)>0) then
 call require(tensor_spin_scale>0.and.ieee_is_finite(tensor_spin_scale), &
 'Tensor RVR requires explicit positive tensor_spin_scale; local KB uses unit directions without this factor')
 call read_tensor(trim(tensor_path),tensor_range,tk,trim(tensor_unit))
 allocate(tke(60*(2*tensor_range+1)**2,60*(2*tensor_range+1)**2,3))
 endif
 allocate(hcv(nd,nd),hb(nd,nd),hu(nd,nd),h(nd,nd),vec(nd,nd),pg(5,30,0:L-1,0:L-1),ns(10,10,nc))
 allocate(ev(nd,nk),f(nd,nk),der(nd,2,nread,nk),tau(2,nread,2))
 if(dump_matrices) then
 call require(nd<=600,'Matrix dump is restricted to small validation cases')
 allocate(db(nd,nd),du(nd,nd))
 endif
 passes=1;if(mode=='table3') passes=2
 do pass=1,passes
 angle=0
 if(pass==2.or.mode=='torque') angle=delta
 ! In torque mode delta is an optional explicit rotation of perturb_site; set tiny/zero via perturb switch below.
 spin=e
 do i=1,nc
 call rotate_n(m%n0,e(:,i),frames(:,:,i))
 enddo
 t=perturb_site(1)*L+perturb_site(2)+1;ax=axis(:,t)
 ! Apply the very same SU(2) left rotation used by the commutator derivative.
 spin(:,t)=cos(angle)*e(:,t)+sin(angle)*cross(ax,e(:,t))+(1-cos(angle))*ax*dot_product(ax,e(:,t))
 g2=pauli(ax);s2=-zi*sin(angle/2)*g2
 s2(1,1)=s2(1,1)+cos(angle/2);s2(2,2)=s2(2,2)+cos(angle/2)
 rot=0
 do i=1,5
 rot(i,i)=s2(1,1);rot(i,i+5)=s2(1,2);rot(i+5,i)=s2(2,1);rot(i+5,i+5)=s2(2,2)
 enddo
 frames(:,:,t)=matmul(rot,matmul(frames(:,:,t),transpose(conjg(rot))))
 do k=1,nk
 if(allocated(tke)) then
 call tensor_plane(tk,q(3,k),tke)
 call assemble(m,L,q(:,k),spin,ueff,hcv,hb,hu,ke,pg,ns,frames,tk,tke,tensor_spin_scale)
 else
 call assemble(m,L,q(:,k),spin,ueff,hcv,hb,hu,ke,pg,ns,frames)
 endif
 h=hcv+hb+hu
 if(dump_matrices) then
 write(suffix,'(a,i0,a,i0)') '.p',pass,'.k',k
 call dump(trim(output_prefix)//trim(suffix)//'.hcv',hcv)
 call dump(trim(output_prefix)//trim(suffix)//'.hb',hb)
 call dump(trim(output_prefix)//trim(suffix)//'.hu',hu)
 call dump(trim(output_prefix)//trim(suffix)//'.h',h)
 endif
 vec=h;call diagonalize(vec,ev(:,k))
 do i=1,nread
 t=sites(1,i)*L+sites(2,i)+1
 call band_derivatives(L,q(:,k),sites(:,i),spin(:,t),axis(:,t),ueff,ke,pg,ns(:,:,t),vec,der(:,1,i,k),der(:,2,i,k))
 if(allocated(tke)) call tensor_band(L,q(:,k),sites(:,i),tk%range,tke, &
 tensor_spin_scale*cross(axis(:,t),spin(:,t)),vec,der(:,1,i,k))
 if(dump_matrices) then
 call derivative_matrices(L,q(:,k),sites(:,i),spin(:,t),axis(:,t),ueff,ke,pg,ns(:,:,t),db,du)
 if(allocated(tke)) then
 db=0
 call tensor_add(db,L,q(:,k),sites(:,i),tk%range,tke,tensor_spin_scale*cross(axis(:,t),spin(:,t)))
 endif
 write(suffix,'(a,i0,a,i0,a,i0)') '.p',pass,'.k',k,'.r',i
 call dump(trim(output_prefix)//trim(suffix)//'.db',db)
 call dump(trim(output_prefix)//trim(suffix)//'.du',du)
 endif
 enddo
 enddo
 call filling(ev,w,electrons_per_ni*nc,sigma,mu,f,count)
 tau(:,:,pass)=0;band_energy=0
 do k=1,nk
 band_energy=band_energy+w(k)*sum(f(:,k)*ev(:,k))
 do i=1,nread
 do j=1,2
 tau(j,i,pass)=tau(j,i,pass)-w(k)*sum(f(:,k)*der(:,j,i,k))
 enddo
 enddo
 enddo
 write(suffix,'(a,i0)') '.p',pass
 open(newunit=u,file=trim(output_prefix)//trim(suffix)//'.bands',status='new',action='write',iostat=ios)
 call require(ios==0,'Band output exists or output directory missing; choose a new output_prefix')
 write(u,'(a)') '# k band Kx Ky Kz E_eV FD_occupation weight'
 do k=1,nk
 do band=1,nd
 write(u,'(2i8,7es26.16e3)') k,band,q(:,k),ev(band,k),f(band,k),w(k)
 enddo
 enddo
 close(u)
 open(newunit=u,file=trim(output_prefix)//trim(suffix)//'.result',status='new',action='write',iostat=ios)
 call require(ios==0,'Result output exists or output directory missing')
 write(u,'(a)') '# Ueff_eV sigma_eV active_electrons EF_eV band_sum_eV angle_rad'
 write(u,'(6es26.16e3)') ueff,sigma,count,mu,band_energy,angle
 write(u,'(a)') '# x y tau_B_eV/rad tau_U_eV/rad tau_total_eV/rad (same full H states)'
 do i=1,nread
 write(u,'(2i8,3es26.16e3)') sites(:,i),tau(:,i,pass),sum(tau(:,i,pass))
 enddo
 write(u,'(a)') '# band_sum=sum_k w Tr[f_FD(H)H]; band sum, eV'
 close(u)
 enddo
 if(mode=='table3') then
 open(newunit=u,file=trim(output_prefix)//'.phi',status='new',action='write',iostat=ios)
 call require(ios==0,'Phi output exists or output directory missing')
 write(u,'(a)') '# x y Phi_B Phi_U Phi_total meV/rad^2; [tau(delta)-tau(0)]/delta; no S^2'
 do i=1,nread
 write(u,'(2i8,3es26.16e3)') sites(:,i),(tau(:,i,2)-tau(:,i,1))*1000/delta, &
 sum(tau(:,i,2)-tau(:,i,1))*1000/delta
 enddo
 close(u)
 endif
 call cpu_time(seconds)
 write(*,'(a,es14.6)') 'Completed Fortran ansatz+U, CPU seconds=',seconds-start
 endif
contains
 subroutine dump(path,a)
 character(*),intent(in)::path
 complex(dp),intent(in)::a(:,:)
 integer::f,stat
 open(newunit=f,file=path,status='new',access='stream',form='unformatted',convert='little_endian',iostat=stat)
 call require(stat==0,'Matrix dump exists or output directory missing: '//path)
 write(f) a;close(f)
 end subroutine
end program
