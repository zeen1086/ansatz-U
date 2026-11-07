program test_spin_channels
  use spin_channels_mod
  implicit none
  real(8) :: spin(3,6), ds(3,6), chi(6), dc(6), h=1d-6, lambda, err
  real(8) :: cp(36), cm(36), c(36), d(36), x(3,6), a(6)
  integer :: nf, mode, j, mask(6), source(6)
  do j=1,6
    spin(:,j) = [sin(0.3d0*j),cos(0.7d0*j),0.2d0*j]
    ds(:,j) = [0.1d0*j,-0.03d0*j,0.04d0*j]
  enddo
  do nf=18,36,18
    do mode=1,4
      mask=1; source=[1,2,3,4,5,6]; lambda=0d0
      if (mode==2) mask=[0,1,-1,0,1,1]
      if (mode==3) source=[6,2,4,1,3,3]
      if (mode==4) lambda=0.7d0
      call spin_chirality(spin,chi,ds,dc)
      call spin_channels(spin,chi,mask,source,lambda,c(:nf),ds,dc,d(:nf))
      x=spin+h*ds
      call spin_chirality(x,a)
      call spin_channels(x,a,mask,source,lambda,cp(:nf))
      x=spin-h*ds
      call spin_chirality(x,a)
      call spin_channels(x,a,mask,source,lambda,cm(:nf))
      err=maxval(abs((cp(:nf)-cm(:nf))/(2*h)-d(:nf)))
      if (err>2d-9) error stop 'spin channel finite difference failed'
      if (mode==1) then
        do j=1,6
          if (maxval(abs(c(3*j-2:3*j)-spin(:,j)))>1d-14) error stop 'default coefficients changed'
        enddo
      endif
    enddo
  enddo
  spin=0d0; ds=0d0
  call spin_chirality(spin,chi,ds,dc)
  if (any(chi/=0d0).or.any(dc/=0d0)) error stop 'zero spin handling failed'
  print *, 'PASS: 18/36 channels, masks, remapping, lambda and finite differences'
end program test_spin_channels
