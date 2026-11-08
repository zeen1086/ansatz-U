module joint_svd
 use numerics
 implicit none
 type::svd_fit
  real(dp),allocatable::map(:,:),design(:,:),spectrum(:),anchor_spectrum(:)
  integer,allocatable::anchor_rows(:)
  integer::rank=0,anchor_rank=0,columns=0
  real(dp)::anchor_relative_cutoff=0,training_relative_cutoff=0
 end type
contains
 subroutine real_product(a,b,c)
 real(dp),intent(in)::a(:,:),b(:,:)
 real(dp),intent(out)::c(:,:)
 call require(size(a,2)==size(b,1).and.size(c,1)==size(a,1).and.size(c,2)==size(b,2),'DGEMM shape mismatch')
 call dgemm('N','N',size(c,1),size(c,2),size(a,2),1._dp,a,size(a,1),b,size(b,1),0._dp,c,size(c,1))
 end subroutine
 subroutine inverse_and_null(x,rcond,inverse,null,singular,rank)
 real(dp),intent(in)::x(:,:),rcond
 real(dp),allocatable,intent(out)::inverse(:,:),null(:,:),singular(:)
 integer,intent(out)::rank
 real(dp),allocatable::a(:,:),u(:,:),vt(:,:),work(:),scaled(:,:)
 real(dp)::query(1),cutoff
 integer::m,n,k,info,lwork,i
 m=size(x,1);n=size(x,2);k=min(m,n)
 allocate(inverse(n,m),singular(k));inverse=0;rank=0
 if(m==0) then
  allocate(null(n,n));null=0
  do i=1,n
   null(i,i)=1
  enddo
  return
 endif
 a=x;allocate(u(m,k),vt(n,n))
 call dgesvd('S','A',m,n,a,m,singular,u,m,vt,n,query,-1,info)
 call require(info==0,'DGESVD workspace query failed')
 lwork=max(1,int(query(1)));allocate(work(lwork));a=x
 call dgesvd('S','A',m,n,a,m,singular,u,m,vt,n,work,lwork,info)
 call require(info==0,'DGESVD failed')
 cutoff=max(rcond,epsilon(1._dp)*max(m,n))*singular(1)
 rank=count(singular>cutoff);allocate(null(n,n-rank));null=transpose(vt(rank+1:,:))
 if(rank>0) then
  scaled=transpose(vt(:rank,:))
  do i=1,rank
   scaled(:,i)=scaled(:,i)/singular(i)
  enddo
  inverse=matmul(scaled,transpose(u(:,:rank)))
 endif
 end subroutine
 subroutine prepare_svd(x,weight,anchor,rcond,fit)
 real(dp),intent(in)::x(:,:),weight(:),rcond
 logical,intent(in)::anchor(:)
 type(svd_fit),intent(out)::fit
 integer,allocatable::ia(:),it(:)
 real(dp),allocatable::ca(:,:),null(:,:),xa(:,:),xt(:,:),reduced(:,:),pinv(:,:),unused(:,:),left(:,:),right(:,:)
 integer::i,n,r
 call require(size(weight)==size(x,1).and.size(anchor)==size(weight),'SVD row metadata mismatch')
 call require(all(ieee_is_finite(x)).and.all(ieee_is_finite(weight)).and.all(weight>=0),'Invalid SVD data/weights')
 call require(rcond>0.and.rcond<1.and.ieee_is_finite(rcond),'rcond must lie strictly between0 and1')
 n=size(x,2);ia=pack([(i,i=1,size(weight))],anchor);it=pack([(i,i=1,size(weight))],.not.anchor.and.weight>0)
 call require(size(ia)+size(it)>0,'No positive-weight training rows or anchors')
 xa=x(ia,:);xt=x(it,:)
 do i=1,size(it)
  xt(i,:)=sqrt(weight(it(i)))*xt(i,:)
 enddo
 call inverse_and_null(xa,rcond,ca,null,fit%anchor_spectrum,fit%anchor_rank)
 allocate(fit%map(n,size(x,1)));fit%map=0
 if(size(null,2)>0.and.size(it)>0) then
  reduced=matmul(xt,null)
  call inverse_and_null(reduced,rcond,pinv,unused,fit%spectrum,r)
  left=matmul(null,pinv);right=ca-matmul(left,matmul(xt,ca))
  fit%map(:,ia)=right
  do i=1,size(it)
   fit%map(:,it(i))=sqrt(weight(it(i)))*left(:,i)
  enddo
  fit%rank=fit%anchor_rank+r
 else
  fit%map(:,ia)=ca;fit%rank=fit%anchor_rank;allocate(fit%spectrum(0))
 endif
 fit%design=x;fit%anchor_rows=ia;fit%columns=n
 fit%anchor_relative_cutoff=max(rcond,epsilon(1._dp)*max(size(ia),n))
 fit%training_relative_cutoff=max(rcond,epsilon(1._dp)*max(size(it),size(null,2)))
 end subroutine
 subroutine solve_svd(fit,y,coefficient,tolerance,label)
 type(svd_fit),intent(in)::fit
 real(dp),intent(in)::y(:,:),tolerance
 real(dp),intent(out)::coefficient(:,:)
 character(*),intent(in)::label
 real(dp),allocatable::difference(:,:)
 real(dp)::maximum
 call real_product(fit%map,y,coefficient)
 if(size(fit%anchor_rows)>0) then
  difference=matmul(fit%design(fit%anchor_rows,:),coefficient)-y(fit%anchor_rows,:)
  maximum=maxval(abs(difference))
  if(maximum>tolerance) write(error_unit,'(a,1x,a,2es24.15)') 'Hard-anchor residual/tolerance:',trim(label),maximum,tolerance
  call require(maximum<=tolerance, &
   'Incompatible hard anchors for '//label//'; inspect raw labels/symmetry/cutoff, no retuning applied')
 endif
 end subroutine
 subroutine report_svd(u,label,fit)
 integer,intent(in)::u
 character(*),intent(in)::label
 type(svd_fit),intent(in)::fit
 write(u,'(a,3(a,i0))') trim(label),' columns=',fit%columns,' identified_rank=',fit%rank,' anchor_rank=',fit%anchor_rank
 write(u,'(a,2es22.13)') '  effective_anchor_training_relative_cutoffs=',fit%anchor_relative_cutoff,fit%training_relative_cutoff
 write(u,'(a,*(es22.13,1x))') '  anchor_singular_values=',fit%anchor_spectrum
 write(u,'(a,*(es22.13,1x))') '  weighted_nullspace_singular_values=',fit%spectrum
 end subroutine
end module
