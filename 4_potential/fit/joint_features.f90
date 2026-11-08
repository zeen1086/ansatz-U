module joint_features
 use numerics
 implicit none
 integer,parameter::ns=8,nb=40,nv=12,ncross=6,nchi=4
 integer,parameter::pairs(2,6)=reshape([1,2,1,3,1,4,2,3,2,4,3,4],[2,6])
 integer,parameter::triples(3,4)=reshape([1,2,3,1,2,4,1,3,4,2,3,4],[3,4])
 type::crystal_action
  integer::a(3,3),shift(3),perm(4)
  real(dp)::orbital(5,5)
 end type
contains
 pure function cross3(a,b) result(c)
 real(dp),intent(in)::a(3),b(3)
 real(dp)::c(3)
 c=[a(2)*b(3)-a(3)*b(2),a(3)*b(1)-a(1)*b(3),a(1)*b(2)-a(2)*b(1)]
 end function
 subroutine spin_features(e,de,s,b,v,x,chi,ds,db,dv,dx,dchi)
 real(dp),intent(in)::e(3,4),de(3,4)
 real(dp),intent(out)::s(ns),b(3,nb),v(3,nv),x(3,ncross),chi(nchi), &
                         ds(ns),db(3,nb),dv(3,nv),dx(3,ncross),dchi(nchi)
 real(dp)::q(6),dq(6),alpha,beta,da,dbeta,phi,dphi
 integer::p,i,j,k,t
 do p=1,6
  i=pairs(1,p);j=pairs(2,p)
  q(p)=dot_product(e(:,i),e(:,j))-1
  dq(p)=dot_product(de(:,i),e(:,j))+dot_product(e(:,i),de(:,j))
  x(:,p)=cross3(e(:,i),e(:,j));dx(:,p)=cross3(de(:,i),e(:,j))+cross3(e(:,i),de(:,j))
 enddo
 alpha=-sum(q([1,2,5,6]))/8;beta=-sum(q([3,4]))/4
 da=-sum(dq([1,2,5,6]))/8;dbeta=-sum(dq([3,4]))/4
 phi=alpha*(1-alpha);dphi=da*(1-2*alpha)
 s=[1._dp,q,phi];ds=[0._dp,dq,dphi]
 v(:,1:4)=e;v(:,5:8)=alpha*e;v(:,9:12)=beta*e
 dv(:,1:4)=de;dv(:,5:8)=da*e+alpha*de;dv(:,9:12)=dbeta*e+beta*de
 b(:,1:12)=v;db(:,1:12)=dv
 b(:,13:16)=phi*e;db(:,13:16)=dphi*e+phi*de
 do p=1,6
 do i=1,4
  t=16+(p-1)*4+i
  b(:,t)=phi*q(p)*e(:,i)
  db(:,t)=(dphi*q(p)+phi*dq(p))*e(:,i)+phi*q(p)*de(:,i)
 enddo
 enddo
 do p=1,4
  i=triples(1,p);j=triples(2,p);k=triples(3,p)
  chi(p)=dot_product(e(:,i),cross3(e(:,j),e(:,k)))
  dchi(p)=dot_product(de(:,i),cross3(e(:,j),e(:,k)))+ &
   dot_product(e(:,i),cross3(de(:,j),e(:,k))+cross3(e(:,j),de(:,k)))
 enddo
 end subroutine
 subroutine crystal_group(grid,group)
 integer,intent(in)::grid(3)
 type(crystal_action),intent(out)::group(64)
 integer::sites(3,4),sw,sx,sy,sz,tx,ty,g,i,j,pos(3)
 call require(grid(1)==grid(2).and.mod(grid(1),4)==0.and.mod(grid(3),2)==0, &
  'Joint model requires the original square four-Ni cell: nx=ny divisible by4, nz even')
 sites(1,:)=[1,3,1,3]*grid(1)/4;sites(2,:)=[1,1,3,3]*grid(2)/4;sites(3,:)=grid(3)/2
 g=0
 do sw=0,1
 do sx=1,-1,-2
 do sy=1,-1,-2
 do sz=1,-1,-2
 do tx=0,1
 do ty=0,1
  g=g+1;group(g)%a=0;group(g)%orbital=0
  if(sw==0) then
   group(g)%a(1,1)=sx;group(g)%a(2,2)=sy
   group(g)%orbital(2,2)=sx*sz;group(g)%orbital(3,3)=sy*sz;group(g)%orbital(4,4)=1
  else
   group(g)%a(1,2)=sx;group(g)%a(2,1)=sy
   group(g)%orbital(2,3)=sx*sz;group(g)%orbital(3,2)=sy*sz;group(g)%orbital(4,4)=-1
  endif
  group(g)%a(3,3)=sz;group(g)%shift=[tx*grid(1)/2,ty*grid(2)/2,0]
  group(g)%orbital(1,1)=1;group(g)%orbital(5,5)=sx*sy
  do i=1,4
   pos=modulo(matmul(group(g)%a,sites(:,i))+group(g)%shift,grid)
   group(g)%perm(i)=0
   do j=1,4
    if(all(pos==sites(:,j))) group(g)%perm(i)=j
   enddo
  enddo
  call require(all(group(g)%perm>0),'Invalid four-Ni crystal permutation')
 enddo
 enddo
 enddo
 enddo
 enddo
 enddo
 end subroutine
 integer function mapped_point(point,grid,g) result(mapped)
 integer,intent(in)::point,grid(3)
 type(crystal_action),intent(in)::g
 integer::xyz(3)
 xyz=[mod(point-1,grid(1)),mod((point-1)/grid(1),grid(2)),(point-1)/(grid(1)*grid(2))]
 xyz=modulo(matmul(g%a,xyz)+g%shift,grid)
 mapped=1+xyz(1)+grid(1)*xyz(2)+grid(1)*grid(2)*xyz(3)
 end function
 subroutine orbital_bases(sym,anti)
 real(dp),intent(out)::sym(5,5,15),anti(5,5,10)
 integer::i,j,k
 sym=0;anti=0
 do i=1,5
  sym(i,i,i)=1
 enddo
 k=0
 do i=1,5
 do j=i+1,5
  k=k+1;sym(i,j,5+k)=1/sqrt(2._dp);sym(j,i,5+k)=sym(i,j,5+k)
  anti(i,j,k)=1/sqrt(2._dp);anti(j,i,k)=-anti(i,j,k)
 enddo
 enddo
 end subroutine
end module
