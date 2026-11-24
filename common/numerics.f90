module numerics
 use iso_fortran_env, only: real64,int32,int64,error_unit
 use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
 implicit none
 integer,parameter::dp=real64
 real(dp),parameter::pi=acos(-1._dp)
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
end module
