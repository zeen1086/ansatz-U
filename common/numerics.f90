module numerics
 use iso_fortran_env, only: real64,int32,int64,error_unit
 use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
 implicit none
 integer,parameter::dp=real64
 real(dp),parameter::pi=acos(-1._dp),ry_ev=13.605693122994_dp
 complex(dp),parameter::zi=(0._dp,1._dp)
contains
 subroutine require(ok,message)
 logical,intent(in)::ok
 character(*),intent(in)::message
 if(.not.ok) then
 write(error_unit,'(a)') trim(message)
 error stop 1
 endif
 end subroutine
 real(dp) function herm_error(a)
 complex(dp),intent(in)::a(:,:)
 herm_error=maxval(abs(a-transpose(conjg(a))))
 end function
 subroutine gemm(a,b,c,ta,tb,alpha,beta)
 complex(dp),intent(in)::a(:,:),b(:,:)
 complex(dp),intent(inout)::c(:,:)
 character,intent(in),optional::ta,tb
 complex(dp),intent(in),optional::alpha,beta
 character::aa,bb
 integer::k
 complex(dp)::al,be
 aa='N';bb='N';al=1;be=0
 if(present(ta)) aa=ta
 if(present(tb)) bb=tb
 if(present(alpha)) al=alpha
 if(present(beta)) be=beta
 k=size(a,2);if(aa/='N') k=size(a,1)
 call zgemm(aa,bb,size(c,1),size(c,2),k,al,a,size(a,1),b,size(b,1),be,c,size(c,1))
 end subroutine
 subroutine diagonalize(a,e)
 complex(dp),intent(inout)::a(:,:)
 real(dp),intent(out)::e(:)
 complex(dp),allocatable::work(:)
 real(dp),allocatable::rw(:)
 complex(dp)::query(1)
 integer::n,info,lw
 n=size(a,1)
 call require(herm_error(a)<1e-8_dp,'Non-Hermitian Hamiltonian; refusing to symmetrize')
 allocate(rw(max(1,3*n-2)))
 call zheev('V','U',n,a,n,e,query,-1,rw,info)
 call require(info==0,'ZHEEV workspace query failed')
 lw=max(2*n-1,int(real(query(1),dp)));allocate(work(lw))
 call zheev('V','U',n,a,n,e,work,lw,rw,info)
 call require(info==0,'ZHEEV diagonalization failed')
 end subroutine
 pure real(dp) function fd(e,mu,sigma)
 real(dp),intent(in)::e,mu,sigma
 real(dp)::x,t
 x=(e-mu)/sigma;t=exp(-abs(x))
 if(x>=0) then
 fd=t/(1+t)
 else
 fd=1/(1+t)
 endif
 end function
 subroutine kahan_add(value,total,compensation)
 real(dp),intent(in)::value
 real(dp),intent(inout)::total,compensation
 real(dp)::y,t
 y=value-compensation;t=total+y
 compensation=(t-total)-y;total=t
 end subroutine
 subroutine filling(e,w,number,sigma,mu,f,count)
 real(dp),intent(in)::e(:,:),w(:),number,sigma
 real(dp),intent(out)::mu,f(:,:),count
 real(dp)::lo,hi,base,corr,cb,cc,x,t,offset
 integer::iter,k,n,filled
 call require(sigma>0.and.number>0.and.number<size(e,1),'Invalid FD width or active electron number')
 call require(size(w)==size(e,2).and.all(w>=0).and.abs(sum(w)-1._dp)<1e-12_dp,'k weights must sum to one')
 lo=minval(e)-max(5._dp,80*sigma);hi=maxval(e)+max(5._dp,80*sigma)
 do iter=1,100
 mu=(lo+hi)/2;base=0;corr=0;cb=0;cc=0
 do k=1,size(e,2)
 filled=0
 do n=1,size(e,1)
 x=(e(n,k)-mu)/sigma;t=exp(-abs(x))
 f(n,k)=fd(e(n,k),mu,sigma)
 if(x<0) then
 filled=filled+1
 call kahan_add(-w(k)*t/(1+t),corr,cc)
 else
 call kahan_add(w(k)*t/(1+t),corr,cc)
 endif
 enddo
 call kahan_add(w(k)*filled,base,cb)
 enddo
 ! Sum thermal particle/hole corrections separately from the integer filled sea.
 offset=(base-number)+corr;count=number+offset
 if(offset>0) then
 hi=mu
 else
 lo=mu
 endif
 enddo
 call require(abs(count-number)<1e-8_dp,'Fixed-number FD solve did not converge')
 end subroutine
end module
