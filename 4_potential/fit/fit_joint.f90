program fit_joint
 use joint_npy
 use joint_model
 use joint_svd
 implicit none
 type(joint_potential)::model
 type(crystal_action)::group(64)
 type(svd_fit)::fits,fitb,fitv,fitx,fitc
 character(4096)::control,manifest,output,report,base,line
 character(128),allocatable::names(:)
 character(2),allocatable::units(:)
 character(4096),allocatable::folders(:)
 integer,allocatable::grid(:),anchor_flag(:)
 integer::u,r,ios,nt,ng,rows,i,g,a,k,j,p,first,last,nblock,block_points=128,atom_ids(4)=[5,6,7,8]
 real(dp)::rcond=1e-10_dp,ueff=1._dp,anchor_tolerance_eV=1e-9_dp,anchor_tolerance_N=1e-10_dp, &
  hermiticity_tolerance=1e-10_dp,occupation_tolerance=1e-8_dp,lo,hi,herm,nlo,nhi,scalar_error,B_error,N_error
 logical::normalize_moments=.true.,remove_scalar_mean=.false.,exists,passed
 namelist /joint_fit_control/ rcond,ueff,anchor_tolerance_eV,anchor_tolerance_N,hermiticity_tolerance, &
  occupation_tolerance,normalize_moments,remove_scalar_mean,block_points,atom_ids
 real(dp),allocatable::weight(:),spin(:,:,:),amplitude(:,:),fields(:,:,:),removed_mean(:), &
  xs(:,:),xb(:,:),xv(:,:),xx(:,:),xc(:,:),ws(:),wv(:),ys(:,:),yb(:,:), &
  yq(:,:),ychi(:,:),ym(:,:),yx(:,:),pred(:,:),dfield(:,:)
 logical,allocatable::as(:),av(:)
 complex(dp),allocatable::occupations(:,:,:,:)
 real(dp)::f(8),b(3,40),v(3,12),x(3,6),chi(4),df(8),db(3,40),dv(3,12),dx(3,6),dc(4), &
  sym(5,5,15),anti(5,5,10),zero(3,4)
 complex(dp)::d(10,10),n(10,10),q(5,5),mag(5,5,3),pn(10,10,4),dn(10,10,4)
 call get_command_argument(1,control);call get_command_argument(2,manifest)
 call get_command_argument(3,output);call get_command_argument(4,report)
 call require(command_argument_count()>=2.and.command_argument_count()<=4.and. &
  min(len_trim(control),len_trim(manifest))>0, &
  'Usage: fit_joint controls.nml dataset.manifest [V_func_tmp [fit_report.txt]]')
 if(len_trim(output)==0) output='V_func_tmp'
 if(len_trim(report)==0) report='fit_report.txt'
 inquire(file=trim(output),exist=exists);call require(.not.exists,'Model output exists; select a new run')
 open(newunit=u,file=trim(control),status='old',action='read');read(u,nml=joint_fit_control);close(u)
 call require(all(ieee_is_finite([rcond,ueff,anchor_tolerance_eV,anchor_tolerance_N, &
  hermiticity_tolerance,occupation_tolerance])),'Nonfinite joint fit control')
 call require(rcond>0.and.rcond<1.and.ueff>=0.and.anchor_tolerance_eV>0.and.anchor_tolerance_N>0.and.hermiticity_tolerance>0 &
  .and.occupation_tolerance>=0.and.block_points>=1.and.block_points<=4096,'Invalid joint fit control')
 call require(all(atom_ids>0),'Invalid atom IDs')
 do i=1,4
  call require(count(atom_ids==atom_ids(i))==1,'Duplicate Ni atom ID')
 enddo
 open(newunit=r,file=trim(report),status='new',action='write',iostat=ios)
 call require(ios==0,'Report output exists or cannot be created')
 write(r,'(a)') 'JOINT_SVD_V1: scalar8 + B12/28 + full complex nonrigid N(scalar/chi/vector/cross)'
 write(r,'(a)') 'Objective: state-weighted raw DFT fields and occupation blocks; not the old projected-operator loss'
 write(r,'(a,es24.15)') 'rcond=',rcond;write(r,'(a,es24.15)') 'training_Ueff_eV=',ueff
 write(r,'(a,2es24.15)') 'anchor_tolerances_eV_N=',anchor_tolerance_eV,anchor_tolerance_N
 write(r,'(a,2es24.15)') 'Hermiticity_occupation_tolerances=',hermiticity_tolerance,occupation_tolerance
 write(r,'(a,4i8)') 'DFT_atom_order=',atom_ids;write(r,'(a,i0)') 'field_block_points=',block_points
 write(r,'(a,2l3)') 'normalize_moments,remove_scalar_mean=',normalize_moments,remove_scalar_mean
 open(newunit=u,file=trim(manifest),status='old',action='read');read(u,*) nt
 call require(nt>=1.and.nt<=1024,'Expected1..1024 explicitly listed DFT states')
 allocate(names(nt),units(nt),folders(nt),weight(nt),anchor_flag(nt))
 base='.';p=index(trim(manifest),'/',back=.true.);if(p>0) base=manifest(:p)
 do i=1,nt
  read(u,*,iostat=ios) names(i),weight(i),anchor_flag(i),units(i)
  call require(ios==0,'Invalid manifest state row: name weight anchor unit')
  call require(ieee_is_finite(weight(i)).and.weight(i)>=0.and. &
   (anchor_flag(i)==0.or.anchor_flag(i)==1),'Invalid manifest weight/anchor')
  call require(units(i)=='Ry'.or.units(i)=='eV','Manifest field units must be Ry or eV')
  read(u,'(a)',iostat=ios) folders(i);call require(ios==0.and.len_trim(folders(i))>0,'Missing explicit DFT directory')
  if(folders(i)(1:1)/='/') folders(i)=trim(base)//'/'//trim(folders(i))
 enddo
 read(u,'(a)',iostat=ios) line;close(u);call require(ios<0,'Unexpected trailing manifest content')
 call npy_integer(trim(folders(1))//'/grid.npy',grid);call require(size(grid)==3,'Expected3 DFT grid dimensions')
 call new_joint(model,grid);model%atom_ids=atom_ids;model%training_ueff=ueff;model%rcond=rcond
 model%removed_scalar_means=remove_scalar_mean;ng=product(grid);rows=64*nt
 call require(int(ng,int64)*nt<=100000000_int64,'DFT label batch exceeds supported100M state-grid points')
 allocate(fields(ng,4,nt),spin(3,4,nt),amplitude(4,nt),occupations(10,10,4,nt),removed_mean(nt))
 do i=1,nt
  call read_dft_label(trim(folders(i)),grid,atom_ids,trim(units(i)),normalize_moments,ueff, &
   hermiticity_tolerance,spin(:,:,i),amplitude(:,i),fields(:,:,i),occupations(:,:,:,i))
  removed_mean(i)=0
  if(remove_scalar_mean) then
   removed_mean(i)=sum(fields(:,1,i))/ng;fields(:,1,i)=fields(:,1,i)-removed_mean(i)
  endif
  write(r,'(a,1x,a,1x,es18.9,1x,i0,1x,a)') 'state=',trim(names(i)),weight(i),anchor_flag(i),trim(folders(i))
  write(r,'(a,5es22.13)') '  original_moment_norms,removed_scalar_mean_eV=',amplitude(:,i),removed_mean(i)
 enddo
 print '(a,i0,a,i0)', 'Read ',nt,' DFT states; grid points=',ng
 call crystal_group(grid,group);call orbital_bases(sym,anti);zero=0
 allocate(xs(rows,8),xb(3*rows,40),xv(3*rows,12),xx(3*rows,6),xc(rows,4), &
  ws(rows),wv(3*rows),as(rows),av(3*rows),yq(rows,15),ychi(rows,10),ym(3*rows,15),yx(3*rows,10))
 do i=1,nt
 do g=1,64
  p=(i-1)*64+g
  call spin_features(spin(:,group(g)%perm,i),zero,f,b,v,x,chi,df,db,dv,dx,dc)
  xs(p,:)=f;xc(p,:)=chi;ws(p)=weight(i);as(p)=anchor_flag(i)==1
  do a=1,3
   k=3*(p-1)+a;xb(k,:)=b(a,:);xv(k,:)=v(a,:);xx(k,:)=x(a,:);wv(k)=weight(i);av(k)=as(p)
  enddo
  d=0;d(:5,:5)=group(g)%orbital;d(6:,6:)=group(g)%orbital
  n=matmul(transpose(d),matmul(occupations(:,:,group(g)%perm(1),i),d))
  call n_to_channels(n,q,mag)
  do j=1,15
   yq(p,j)=sum(real(q,dp)*sym(:,:,j))
   do a=1,3
    ym(3*(p-1)+a,j)=sum(real(mag(:,:,a),dp)*sym(:,:,j))
   enddo
  enddo
  do j=1,10
   ychi(p,j)=sum(aimag(q)*anti(:,:,j))
   do a=1,3
    yx(3*(p-1)+a,j)=sum(aimag(mag(:,:,a))*anti(:,:,j))
   enddo
  enddo
 enddo
 enddo
 call prepare_svd(xs,ws,as,rcond,fits);call prepare_svd(xb,wv,av,rcond,fitb)
 call prepare_svd(xv,wv,av,rcond,fitv);call prepare_svd(xx,wv,av,rcond,fitx)
 call prepare_svd(xc,ws,as,rcond,fitc)
 call report_svd(r,'scalar8/N_scalar',fits);call report_svd(r,'B12plus28',fitb)
 call report_svd(r,'N_vector',fitv);call report_svd(r,'N_cross',fitx);call report_svd(r,'N_chi',fitc)
 call solve_svd(fits,yq,model%nscalar,anchor_tolerance_N,'N scalar')
 call solve_svd(fitc,ychi,model%nchi,anchor_tolerance_N,'N chirality')
 call solve_svd(fitv,ym,model%nvector,anchor_tolerance_N,'N vector')
 call solve_svd(fitx,yx,model%ncross,anchor_tolerance_N,'N cross')
 allocate(ys(rows,min(block_points,ng)),yb(3*rows,min(block_points,ng)))
 do first=1,ng,block_points
  last=min(ng,first+block_points-1);nblock=last-first+1
  do i=1,nt
  do g=1,64
   p=(i-1)*64+g
   do j=1,nblock
    k=mapped_point(first+j-1,grid,group(g));ys(p,j)=fields(k,1,i)
    do a=1,3
     yb(3*(p-1)+a,j)=fields(k,a+1,i)
    enddo
   enddo
  enddo
  enddo
  call solve_svd(fits,ys(:,:nblock),model%scalar(:,first:last),anchor_tolerance_eV,'scalar field')
  call solve_svd(fitb,yb(:,:nblock),model%magnetic(:,first:last),anchor_tolerance_eV,'B field')
 enddo
 allocate(pred(ng,4),dfield(ng,4));passed=.true.;nlo=huge(1._dp);nhi=-nlo
 write(r,'(a)') 'Residuals: state scalar_RMS_eV B_component_RMS_eV N_complex_RMS eig_min eig_max'
 do i=1,nt
  call joint_value(model,spin(:,:,i),zero,pred,pn,dfield,dn)
  scalar_error=sqrt(sum((pred(:,1)-fields(:,1,i))**2)/ng)
  B_error=sqrt(sum((pred(:,2:4)-fields(:,2:4,i))**2)/(3*ng))
  N_error=sqrt(sum(abs(pn-occupations(:,:,:,i))**2)/400)
  call occupation_bounds(pn,lo,hi,herm);nlo=min(nlo,lo);nhi=max(nhi,hi)
  passed=passed.and.lo>=-occupation_tolerance.and.hi<=1+occupation_tolerance.and.herm<=hermiticity_tolerance
  write(r,'(a,1x,5es24.15)') trim(names(i)),scalar_error,B_error,N_error,lo,hi
  if(anchor_flag(i)==1) then
   call require(maxval(abs(pred-fields(:,:,i)))<=anchor_tolerance_eV,'Predicted field hard anchor failed')
   call require(maxval(abs(pn-occupations(:,:,:,i)))<=anchor_tolerance_N,'Predicted complete-N hard anchor failed')
  endif
 enddo
 model%training_gate=merge(1,0,passed)
 write(r,'(a,2es24.15)') 'predicted_N_eigenvalue_range=',nlo,nhi
 write(r,'(a,l2)') 'training_occupation_gate=',passed
 write(r,'(a)') 'No eigenvalue clipping. Rank-deficient directions use the stated SVD minimum norm.'
 call write_joint(trim(output),model);close(r)
 call require(passed,'Training N occupation gate failed; candidate retained with gate=0 and cannot be evaluated')
 print '(a)', 'Joint scalar/B/nonrigid Hubbard model written: '//trim(output)
end program
