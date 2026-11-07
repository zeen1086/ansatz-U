module local_ansatz
 use frozen_io
 use tensor_rvr
 implicit none
contains
 function pauli(e) result(s)
 real(dp),intent(in)::e(3)
 complex(dp)::s(2,2)
 s(1,1)=e(3);s(2,2)=-e(3)
 s(1,2)=e(1)-zi*e(2);s(2,1)=e(1)+zi*e(2)
 end function
 subroutine rotate_n(n0,e,n)
 complex(dp),intent(in)::n0(10,10)
 real(dp),intent(in)::e(3)
 complex(dp),intent(out)::n(10,10)
 complex(dp)::d(10,10),r(2,2)
 real(dp)::theta,phi,c,s
 integer::i
 call require(abs(sum(e**2)-1._dp)<1e-12_dp,'Local model requires unit directions, B amplitude is already in KB')
 theta=atan2(sqrt(e(1)**2+e(2)**2),e(3));phi=atan2(e(2),e(1))
 if(abs(e(1))+abs(e(2))<1e-14_dp) phi=0
 c=cos(theta/2);s=sin(theta/2)
 ! D=Rz(phi) Ry(theta), spin-major projector order [up d1..5,down d1..5].
 r(1,1)=exp(-zi*phi/2)*c;r(1,2)=-exp(-zi*phi/2)*s
 r(2,1)=exp(zi*phi/2)*s;r(2,2)=exp(zi*phi/2)*c
 d=0
 do i=1,5
 d(i,i)=r(1,1);d(i,i+5)=r(1,2);d(i+5,i)=r(2,1);d(i+5,i+5)=r(2,2)
 enddo
 n=matmul(d,matmul(n0,transpose(conjg(d))))
 call require(herm_error(n)<1e-12_dp,'Rotated N is not Hermitian')
 end subroutine
 subroutine n_derivative(n,axis,dn)
 complex(dp),intent(in)::n(10,10)
 real(dp),intent(in)::axis(3)
 complex(dp),intent(out)::dn(10,10)
 complex(dp)::g(10,10),s(2,2)
 integer::i
 s=-zi*.5_dp*pauli(axis);g=0
 do i=1,5
 g(i,i)=s(1,1);g(i,i+5)=s(1,2);g(i+5,i)=s(2,1);g(i+5,i+5)=s(2,2)
 enddo
 dn=matmul(g,n)-matmul(n,g)
 end subroutine
 function cross(a,b) result(c)
 real(dp),intent(in)::a(3),b(3)
 real(dp)::c(3)
 c=[a(2)*b(3)-a(3)*b(2),a(3)*b(1)-a(1)*b(3),a(1)*b(2)-a(2)*b(1)]
 end function
 subroutine project_at(L,k,pg,site,a)
 integer,intent(in)::L,site(2)
 real(dp),intent(in)::k(3)
 complex(dp),intent(in)::pg(:,:,0:,0:)
 complex(dp),intent(out)::a(10,60*L*L)
 integer::x,y,d(2),orb,cell
 complex(dp)::phase
 a=0
 do x=0,L-1
 do y=0,L-1
 d=site-[x,y];cell=x*L+y
 phase=exp(2*pi*zi*dot_product(real(d,dp),k(:2)/L))
 do orb=1,30
 a(:5,cell*60+2*orb-1)=pg(:,orb,modulo(d(1),L),modulo(d(2),L))*phase
 a(6:,cell*60+2*orb)=pg(:,orb,modulo(d(1),L),modulo(d(2),L))*phase
 enddo
 enddo
 enddo
 end subroutine
 subroutine local_indices(L,k,site,ids,ph)
 integer,intent(in)::L,site(2)
 real(dp),intent(in)::k(3)
 integer,intent(out)::ids(9)
 complex(dp),intent(out)::ph(9)
 integer::x,y,t,raw(2),wrapped(2),cells(2)
 t=0
 do x=-1,1
 do y=-1,1
 t=t+1;raw=site+[x,y];wrapped=modulo(raw,L);cells=(raw-wrapped)/L
 ids(t)=wrapped(1)*L+wrapped(2)
 ph(t)=exp(2*pi*zi*dot_product(real(cells,dp),k(:2)))
 enddo
 enddo
 end subroutine
 subroutine add_b(h,L,k,site,ke,e)
 integer,intent(in)::L,site(2)
 real(dp),intent(in)::k(3),e(3)
 complex(dp),intent(in)::ke(270,270)
 complex(dp),intent(inout)::h(:,:)
 complex(dp)::ph(9),s(2,2),v
 integer::ids(9),a,b,i,j,u,w,ii,jj
 call local_indices(L,k,site,ids,ph);s=pauli(e)
 do a=1,9
 do b=1,9
 do j=1,30
 do i=1,30
 v=conjg(ph(a))*ph(b)*ke((a-1)*30+i,(b-1)*30+j)
 do u=1,2
 do w=1,2
 ii=ids(a)*60+2*(i-1)+u;jj=ids(b)*60+2*(j-1)+w
 h(ii,jj)=h(ii,jj)+v*s(u,w)
 enddo
 enddo
 enddo
 enddo
 enddo
 enddo
 end subroutine
 subroutine assemble(m,L,k,e,ueff,hcv,hb,hu,ke,pg,nsite,n_override,tk,tke,spin_scale)
 type(frozen_model),intent(in)::m
 integer,intent(in)::L
 real(dp),intent(in)::k(3),e(:,:),ueff
 complex(dp),intent(out)::hcv(:,:),hb(:,:),hu(:,:),ke(270,270),pg(:,:,0:,0:),nsite(:,:,:)
 complex(dp),intent(in),optional::n_override(:,:,:),tke(:,:,:)
 type(tensor_kernel),intent(in),optional::tk
 real(dp),intent(in),optional::spin_scale
 complex(dp),allocatable::hcq(:,:,:),pq(:,:,:),a(:,:),va(:,:)
 complex(dp)::hc(30,30),p(5,30),v(10,10),phase
 real(dp)::q(3)
 integer::x,y,t,xx,yy,s,d(2),i,j,u,iq,nc,nd
 nc=L*L;nd=60*nc
 call require(size(e,2)==nc,'Spin texture size mismatch')
 allocate(hcq(30,30,nc),pq(5,30,nc),a(10,nd),va(10,nd))
 do x=0,L-1
 do y=0,L-1
 iq=x*L+y+1;q=[(x+k(1))/L,(y+k(2))/L,k(3)]
 call interpolate(m,q,hc,p);hcq(:,:,iq)=hc;pq(:,:,iq)=p
 enddo
 enddo
 pg=0;hcv=0;hb=0;hu=0
 ! Finite inverse Bloch transform: the same normalized IFFT and twist as the reference.
 do x=0,L-1
 do y=0,L-1
 do xx=0,L-1
 do yy=0,L-1
 iq=xx*L+yy+1
 phase=exp(2*pi*zi*real(xx*x+yy*y,dp)/L)/nc
 pg(:,:,x,y)=pg(:,:,x,y)+phase*pq(:,:,iq)
 enddo
 enddo
 enddo
 enddo
 do x=0,L-1
 do y=0,L-1
 t=x*L+y
 do xx=0,L-1
 do yy=0,L-1
 s=xx*L+yy;d=[x-xx,y-yy];hc=0
 do iq=1,nc
 q=[( (iq-1)/L+k(1))/L,(modulo(iq-1,L)+k(2))/L,k(3)]
 phase=exp(2*pi*zi*dot_product(q(:2),real(d,dp)))/nc
 hc=hc+phase*hcq(:,:,iq)
 enddo
 do j=1,30
 do i=1,30
 do u=1,2
 hcv(t*60+2*i-2+u,s*60+2*j-2+u)=hc(i,j)
 enddo
 enddo
 enddo
 enddo
 enddo
 enddo
 enddo
 call plane_kernel(m,k(3),ke)
 do x=0,L-1
 do y=0,L-1
 t=x*L+y+1
 if(present(n_override)) then
 nsite(:,:,t)=n_override(:,:,t)
 else
 call rotate_n(m%n0,e(:,t),nsite(:,:,t))
 endif
 if(present(tk)) then
 call tensor_add(hb,L,k,[x,y],tk%range,tke,spin_scale*e(:,t))
 else
 call add_b(hb,L,k,[x,y],ke,e(:,t))
 endif
 if(ueff==0) cycle
 call project_at(L,k,pg,[x,y],a)
 v=-ueff*nsite(:,:,t)
 do i=1,10
 v(i,i)=v(i,i)+ueff/2
 enddo
 call gemm(v,a,va)
 call gemm(a,va,hu,ta='C',beta=cmplx(1._dp,0._dp,dp))
 enddo
 enddo
 call require(herm_error(hcv)<1e-8_dp.and.herm_error(hb)<1e-8_dp.and.herm_error(hu)<1e-8_dp, &
 'H0+Hv/HB/HU assembly is not Hermitian')
 end subroutine
 subroutine derivative_matrices(L,k,site,e,axis,ueff,ke,pg,n,db,du)
 integer,intent(in)::L,site(2)
 real(dp),intent(in)::k(3),e(3),axis(3),ueff
 complex(dp),intent(in)::ke(270,270),pg(:,:,0:,0:),n(10,10)
 complex(dp),intent(out)::db(:,:),du(:,:)
 complex(dp)::dn(10,10),a(10,60*L*L),v(10,60*L*L)
 db=0;du=0
 call add_b(db,L,k,site,ke,cross(axis,e))
 call n_derivative(n,axis,dn)
 call project_at(L,k,pg,site,a)
 call gemm(dn,a,v,alpha=cmplx(-ueff,0._dp,dp))
 call gemm(a,v,du,ta='C')
 end subroutine
 subroutine band_derivatives(L,k,site,e,axis,ueff,ke,pg,n,vec,db,du)
 integer,intent(in)::L,site(2)
 real(dp),intent(in)::k(3),e(3),axis(3),ueff
 complex(dp),intent(in)::ke(270,270),pg(:,:,0:,0:),n(10,10),vec(:,:)
 real(dp),intent(out)::db(:),du(:)
 complex(dp),allocatable::up(:,:),down(:,:),tmp(:,:),av(:,:),nv(:,:),a(:,:)
 complex(dp)::ph(9),s(2,2),dn(10,10)
 integer::ids(9),j,i,row,nb
 nb=size(vec,2)
 allocate(up(270,nb),down(270,nb),tmp(270,nb),a(10,60*L*L),av(10,nb),nv(10,nb))
 call local_indices(L,k,site,ids,ph)
 do j=1,9
 do i=1,30
 row=ids(j)*60+2*i-1
 up((j-1)*30+i,:)=ph(j)*vec(row,:)
 down((j-1)*30+i,:)=ph(j)*vec(row+1,:)
 enddo
 enddo
 s=pauli(cross(axis,e));db=0
 call gemm(ke,up,tmp)
 db=real(sum(conjg(up)*tmp*s(1,1)+conjg(down)*tmp*s(2,1),dim=1),dp)
 call gemm(ke,down,tmp)
 db=db+real(sum(conjg(up)*tmp*s(1,2)+conjg(down)*tmp*s(2,2),dim=1),dp)
 call project_at(L,k,pg,site,a);call gemm(a,vec,av)
 call n_derivative(n,axis,dn);call gemm(dn,av,nv)
 du=-ueff*real(sum(conjg(av)*nv,dim=1),dp)
 end subroutine
end module
