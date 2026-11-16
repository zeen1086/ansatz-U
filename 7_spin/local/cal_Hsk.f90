program cal_Hsk
 use local_rvr, only: local_projected,open_projected,projected_block
 use frozen_io, only: frozen_model,read_model,background_projectors,project_at
 use numerics
 !$ use omp_lib, only: omp_get_max_threads
 implicit none
 type(local_projected)::rvr
 type(frozen_model)::bg
 character(4096)::control,rvr_path='',background_path='',coefficient_path='',kpoints_path='',output_prefix=''
 character(64)::background_kind='',model_id
 character(16)::magic
 character(1024)::message
 logical::dump_components=.false.,derivative
 namelist /local_hamiltonian/ rvr_path,background_path,background_kind,coefficient_path, &
 kpoints_path,output_prefix,dump_components
 integer,parameter::max_k_batch=16
 integer(int32)::head(8)
 integer::u,v,w,ios,L,nsite,nd,nk,ik,kz,s,c,first,last
 integer,allocatable::sectors(:),members(:)
 real(dp)::ueff,offset,mean,bias
 real(dp),allocatable::p(:,:,:),points(:,:),weights(:),errors(:),magnitudes(:)
 complex(dp),allocatable::n(:,:,:),h0(:,:),hu(:,:),h(:,:),pg(:,:,:,:),projector(:,:),vp(:,:),local_u(:,:)
 integer(int64)::bytes,expected
 call require(command_argument_count()==1,'Usage: cal_Hsk local_run.nml (coefficients from local_vfunc.py evaluate)')
 call get_command_argument(1,control)
 open(newunit=u,file=trim(control),status='old',action='read',iostat=ios)
 call require(ios==0,'Cannot open local Hamiltonian control')
 read(u,nml=local_hamiltonian,iostat=ios,iomsg=message);close(u)
 call require(ios==0,'Invalid local_hamiltonian: '//trim(message))
 call require(min(len_trim(rvr_path),len_trim(background_path),len_trim(coefficient_path), &
 len_trim(kpoints_path),len_trim(output_prefix))>0,'Explicit local Hamiltonian paths required')
 call require(background_kind=='pure_H0_no_local_fields_no_HU','Explicit pure-H0 background contract required')
 call open_projected(trim(rvr_path),v,rvr)
 call require(rvr%nw==30,'The LaNiO2 background/projector interface requires 30 spatial orbitals')
 open(newunit=u,file=trim(coefficient_path),status='old',access='stream',form='unformatted',convert='little_endian',iostat=ios)
 call require(ios==0,'Missing complete V1 prediction coefficients')
 read(u,iostat=ios) magic,head,model_id,ueff,offset
 call require(ios==0.and.magic=='LANIO2_COEFF_V2 '.and.head(1)==1,'Invalid local coefficient header')
 call require(model_id==rvr%model_id.and.head(2)==rvr%np,'V_func_tmp / B_general model mismatch')
 L=head(3);nsite=head(4);nd=60*L*L
 call require(L>=1.and.L<=100.and.nsite==L*L.and.(head(5)==0.or.head(5)==1),'Invalid local coefficient cell/mode')
 call require(ueff>=0.and.ieee_is_finite(ueff).and.ieee_is_finite(offset),'Invalid Ueff or gauge offset')
 derivative=head(5)==1;call require(.not.derivative.or.offset==0,'A derivative cannot contain a constant energy offset')
 inquire(unit=u,size=bytes)
 expected=128_int64+32_int64*rvr%np*nsite+1600_int64*nsite
 call require(bytes==expected,'Local prediction byte count mismatch')
 allocate(p(rvr%np,4,nsite),n(10,10,nsite));read(u,iostat=ios) p,n;close(u)
 call require(ios==0.and.all(ieee_is_finite(p)).and.all(ieee_is_finite(real(n,dp))).and. &
 all(ieee_is_finite(aimag(n))),'Nonfinite/truncated local prediction')
 do s=1,nsite
  call require(herm_error(n(:,:,s))<1e-10_dp,'Non-Hermitian full N/dN')
 enddo
 bias=rvr%alpha;if(derivative) bias=0
 mean=(sum(p(:,1,:)*spread(rvr%means(:rvr%np),2,nsite))+nsite*bias*rvr%means(rvr%np+1))/nsite
 call read_model(trim(background_path),bg)
 allocate(h0(nd,nd),hu(nd,nd),h(nd,nd),pg(5,30,0:L-1,0:L-1),projector(10,nd),vp(10,nd),local_u(10,10))
 open(newunit=u,file=trim(kpoints_path),status='old',action='read',iostat=ios)
 call require(ios==0,'Missing k points');read(u,*,iostat=ios) nk
 call require(ios==0.and.nk>0,'Invalid k point count')
 allocate(points(3,nk),weights(nk),sectors(nk),errors(nk),magnitudes(nk))
 do ik=1,nk
  read(u,*,iostat=ios) points(:,ik),weights(ik)
  call require(ios==0.and.all(ieee_is_finite(points(:,ik))).and.ieee_is_finite(weights(ik)).and. &
   weights(ik)>=0,'Invalid weighted k point')
  sectors(ik)=0
  do c=1,rvr%nk
   if(abs(points(3,ik)-rvr%kz(c)-anint(points(3,ik)-rvr%kz(c)))<1e-10_dp) sectors(ik)=c
  enddo
  call require(sectors(ik)>0,'This kz has not been projected; add it to projection input')
 enddo
 close(u)
 call require(abs(sum(weights)-1)<1e-10_dp,'k weights must sum to one')
 open(newunit=w,file=trim(output_prefix)//'.index',status='new',action='write')
 write(w,'(a)') 'LANIO2_LOCAL_H_V2 complex128 little_endian Fortran_order eV'
 write(w,*) L,nd,nk,derivative
 write(w,'(a)') rvr%model_id
 write(w,'(a)') rvr%basis_id
 do kz=1,rvr%nk
  members=pack([(ik,ik=1,nk)],sectors==kz)
  do first=1,size(members),max_k_batch
   last=min(first+max_k_batch-1,size(members))
   call assemble_sector(kz,members(first:last))
  enddo
 enddo
 do ik=1,nk
  write(w,'(i8,7es25.16)') ik,points(:,ik),weights(ik),mean,errors(ik),magnitudes(ik)
 enddo
 close(v);close(w)
contains
 subroutine assemble_sector(sector,indices)
 integer,intent(in)::sector,indices(:)
 integer::nr,np,nk_batch,pairs,ik,slot,i,j,s,x,y,a,b,pair,ii,jj,raw(2),wrapped(2),shift(2)
 integer::block_reads,ns,threads
 integer,allocatable::site(:,:),mean_site(:,:)
 complex(dp),allocatable,target::kernel(:,:,:)
 complex(dp),pointer::matrix(:,:)
 complex(dp),allocatable::hv_batch(:,:,:),hb_batch(:,:,:),phase(:,:,:),mean_phase(:,:,:),overlap(:)
 complex(dp),allocatable::charge(:,:),spin(:,:),v_pair(:,:),b_pair(:,:)
 complex(dp)::ph,values(4)
 nr=rvr%nr;np=rvr%np;ns=rvr%ndiff;nk_batch=size(indices);pairs=nr*nr
 threads=1
 !$ threads=min(nk_batch,omp_get_max_threads())
 allocate(hv_batch(nd,nd,nk_batch),hb_batch(nd,nd,nk_batch),kernel(nr,nr,rvr%nc),overlap(ns), &
  phase(nr,nsite,nk_batch),mean_phase(ns,nsite,nk_batch),site(nr,nsite),mean_site(ns,nsite), &
  charge(np+1,nsite),spin(np,3*nsite),v_pair(pairs,nsite),b_pair(pairs,3*nsite))
 matrix(1:pairs,1:rvr%nc)=>kernel
 do s=1,nsite
  charge(:np,s)=p(:,1,s);charge(np+1,s)=bias
  spin(:,3*s-2:3*s)=p(:,2:4,s)
  x=(s-1)/L;y=mod(s-1,L)
  do a=1,nr
   raw=[x,y]+rvr%xy(:,a);wrapped=modulo(raw,L);shift=(raw-wrapped)/L
   site(a,s)=wrapped(1)*L+wrapped(2)
   do slot=1,nk_batch
    phase(a,s,slot)=exp(2*pi*zi*dot_product(points(:2,indices(slot)),real(shift,dp)))
   enddo
  enddo
  do a=1,ns
   raw=[x,y]+rvr%difference(:,a);wrapped=modulo(raw,L);shift=(raw-wrapped)/L
   mean_site(a,s)=wrapped(1)*L+wrapped(2)
   do slot=1,nk_batch
    mean_phase(a,s,slot)=exp(2*pi*zi*dot_product(points(:2,indices(slot)),real(shift,dp)))
   enddo
  enddo
 enddo
 hv_batch=0;hb_batch=0;block_reads=0
 do i=1,rvr%nw
 do j=1,rvr%nw
  call projected_block(v,rvr,sector,i,j,overlap,kernel)
  block_reads=block_reads+1
  ! The two contractions depend on the coefficient file, orbital pair and kz only.
  call gemm(matrix(:,:np+1),charge,v_pair)
  call gemm(matrix(:,np+2:),spin,b_pair)
  !$omp parallel do if(threads>1) num_threads(threads) private(s,a,b,pair,ii,jj,ph,values)
  do slot=1,nk_batch
   do s=1,nsite
    do a=1,nr
     ii=site(a,s)*60+2*i-1
     do b=1,nr
      pair=a+(b-1)*nr
      ph=conjg(phase(a,s,slot))*phase(b,s,slot)
      values=[v_pair(pair,s),b_pair(pair,3*s-2:3*s)]*ph
      jj=site(b,s)*60+2*j-1
      hv_batch(ii,jj,slot)=hv_batch(ii,jj,slot)+values(1)
      hv_batch(ii+1,jj+1,slot)=hv_batch(ii+1,jj+1,slot)+values(1)
      hb_batch(ii,jj,slot)=hb_batch(ii,jj,slot)+values(4)
      hb_batch(ii+1,jj+1,slot)=hb_batch(ii+1,jj+1,slot)-values(4)
      hb_batch(ii,jj+1,slot)=hb_batch(ii,jj+1,slot)+values(2)-zi*values(3)
      hb_batch(ii+1,jj,slot)=hb_batch(ii+1,jj,slot)+values(2)+zi*values(3)
     enddo
    enddo
    ! Subtract the scalar mean in the same Wannier basis.
    ii=(s-1)*60+2*i-1
    do a=1,ns
     jj=mean_site(a,s)*60+2*j-1
     ph=mean*overlap(a)*mean_phase(a,s,slot)
     hv_batch(ii,jj,slot)=hv_batch(ii,jj,slot)-ph
     hv_batch(ii+1,jj+1,slot)=hv_batch(ii+1,jj+1,slot)-ph
    enddo
   enddo
  enddo
  !$omp end parallel do
 enddo
 enddo
 write(*,'(a,3i10)') 'RVR_KZ_BATCH sector, points, blocks=',sector,nk_batch,block_reads
 do slot=1,nk_batch
  ik=indices(slot)
  call finish_point(ik,hv_batch(:,:,slot),hb_batch(:,:,slot))
 enddo
 end subroutine

 subroutine finish_point(ik,hv,hb)
 integer,intent(in)::ik
 complex(dp),intent(inout)::hv(:,:)
 complex(dp),intent(in)::hb(:,:)
 integer::i,x,y,s
 real(dp)::k(3)
 character(32)::tag
 k=points(:,ik)
 call background_projectors(bg,L,k,h0,pg)
 if(derivative) h0=0
 hu=0
 do i=1,nd
  hv(i,i)=hv(i,i)+offset
 enddo
 do x=0,L-1
 do y=0,L-1
  s=x*L+y+1;call project_at(L,k,pg,[x,y],projector)
  local_u=-ueff*n(:,:,s)
  if(.not.derivative) then
   do i=1,10
    local_u(i,i)=local_u(i,i)+ueff/2
   enddo
  endif
  call gemm(local_u,projector,vp);call gemm(projector,vp,hu,ta='C',beta=cmplx(1._dp,0._dp,dp))
 enddo
 enddo
 h=h0+hv+hb+hu
 errors(ik)=herm_error(h);magnitudes(ik)=maxval(abs(h))
 call require(max(errors(ik),herm_error(hv),herm_error(hb),herm_error(hu))<1e-8_dp, &
  'Non-Hermitian V1 Hamiltonian/component; refusing to symmetrize')
 write(tag,'(a,i6.6)') '.k',ik
 call dump(trim(output_prefix)//trim(tag)//'.h',h)
 if(dump_components) then
  call dump(trim(output_prefix)//trim(tag)//'.h0',h0);call dump(trim(output_prefix)//trim(tag)//'.hv',hv)
  call dump(trim(output_prefix)//trim(tag)//'.hb',hb);call dump(trim(output_prefix)//trim(tag)//'.hu',hu)
 endif
 write(*,'(a,i8,a,i8)') 'HAMILTONIAN_K_COMPLETE ',ik,' dimension ',nd
 end subroutine

 subroutine dump(path,matrix)
 character(*),intent(in)::path
 complex(dp),intent(in)::matrix(:,:)
 integer::fd
 open(newunit=fd,file=path,status='new',access='stream',form='unformatted',convert='little_endian')
 write(fd) matrix;close(fd)
 end subroutine
end program
