! Known-answer and equivalence self test for the packed/ZGEMM RVR kernel.
! Covers: complex conjugation, Pauli structure (pot1 = B1 - i B2,
! pot2 = B1 + i B2, pot3 = B3), spin-slot pairing, matrix orientation
! (i = bra/R, j = ket/R'), support masks (dir_active, cal_TF, dir_inv),
! index back-fill into the full d_ext x d_ext block and chunk accumulation.
program test_known_answer
  use rvr_kernel_mod
  implicit none
  integer, parameter :: d_ext = 9, ni = 257, ndep = 2, nsc = 3
  integer, parameter :: nfunc = ndep*nsc
  complex(dp) :: fin_wan(d_ext, ni, 2), ini_wan(d_ext, ni, 2)
  real(dp) :: B4(ni, 3, ndep, nsc)
  integer :: dir_vec(d_ext, 3), dir_inv(d_ext), dir_active(d_ext)
  integer :: fin_TF(d_ext), ini_TF(d_ext)
  integer :: n_i, idx_i(d_ext), n_j, idx_j(d_ext)
  integer :: cs, i, j, s, k, jd, sc, nfail, seed
  complex(dp) :: full_ref(d_ext, d_ext), full_new(d_ext, d_ext)
  complex(dp) :: A(2*ni, d_ext), I1(ni, d_ext), I2(ni, d_ext)
  complex(dp), allocatable :: C(:, :)
  complex(dp) :: rtmp(d_ext, d_ext)
  complex(dp) :: expect, v1, v2, p1, p2, p3, f1, f2, g1, g2
  real(dp) :: err, maxerr, denom

  nfail = 0
  call build_dirs()
  print *, '=== RVR kernel known-answer test ==='

  ! ---------------- test 1: single grid point, analytic Pauli check --------
  block
    integer, parameter :: ni1 = 1
    complex(dp) :: fw1(d_ext, ni1, 2), iw1(d_ext, ni1, 2)
    complex(dp) :: A1(2*ni1, d_ext), J1(ni1, d_ext), J2(ni1, d_ext)
    complex(dp) :: c1(d_ext, d_ext)
    real(dp) :: b1(ni1, 3, ndep, nsc)
    integer :: ii(d_ext), jj(d_ext), n1, n2
    fin_TF = 1; ini_TF = 1; dir_active = 1
    call build_support(d_ext, dir_active, dir_inv, fin_TF, ini_TF, n1, ii, n2, jj)
    fw1 = (0.0d0, 0.0d0); iw1 = (0.0d0, 0.0d0); b1 = 0.0d0
    do k = 1, d_ext
      fw1(k, 1, 1) = dcmplx(0.3d0*k, -0.7d0)
      fw1(k, 1, 2) = dcmplx(1.1d0, 0.2d0*k)
      iw1(k, 1, 1) = dcmplx(-0.4d0, 0.9d0*k)
      iw1(k, 1, 2) = dcmplx(0.8d0*k, 1.3d0)
    end do
    b1(1, 1, 2, 3) = 0.55d0      ! B1  (dependency 2, spin component 3 => channel 6)
    b1(1, 2, 2, 3) = -0.25d0     ! B2
    b1(1, 3, 2, 3) = 0.90d0      ! B3
    A1 = (0.0d0, 0.0d0); J1 = (0.0d0, 0.0d0); J2 = (0.0d0, 0.0d0)
    c1 = (0.0d0, 0.0d0)
    call pack_bra(fw1, d_ext, ni1, dir_inv, n1, ii, A1)
    call pack_ket(iw1, d_ext, ni1, dir_inv, n2, jj, J1, J2)
    call apply_channels(ni1, n1, n2, A1, J1, J2, b1, ndep, nsc, (/6/), 1, &
                        ni1, c1)
    p1 = dcmplx(b1(1, 1, 2, 3), -b1(1, 2, 2, 3))
    p2 = dcmplx(b1(1, 1, 2, 3),  b1(1, 2, 2, 3))
    p3 = dcmplx(b1(1, 3, 2, 3),  0.0d0)
    v1 = p1*iw1(dir_inv(jj(1)), 1, 2) + p3*iw1(dir_inv(jj(1)), 1, 1)
    v2 = p2*iw1(dir_inv(jj(1)), 1, 1) - p3*iw1(dir_inv(jj(1)), 1, 2)
    expect = dconjg(fw1(dir_inv(ii(1)), 1, 1))*v1 &
           + dconjg(fw1(dir_inv(ii(1)), 1, 2))*v2
    call chk('analytic Pauli value at (i=1,j=1)', abs(c1(1, 1) - expect), 1d-14)
    v1 = p1*iw1(dir_inv(jj(1)), 1, 2) + p3*iw1(dir_inv(jj(1)), 1, 1)
    v2 = p2*iw1(dir_inv(jj(1)), 1, 1) - p3*iw1(dir_inv(jj(1)), 1, 2)
    expect = dconjg(fw1(dir_inv(ii(3)), 1, 1))*v1 &
           + dconjg(fw1(dir_inv(ii(3)), 1, 2))*v2
    call chk('analytic value at (i=3,j=1)', abs(c1(3, 1) - expect), 1d-14)
    ! orientation: the (i=1,j=3) entry must use jj(3) on the ket side only
    v1 = p1*iw1(dir_inv(jj(3)), 1, 2) + p3*iw1(dir_inv(jj(3)), 1, 1)
    v2 = p2*iw1(dir_inv(jj(3)), 1, 1) - p3*iw1(dir_inv(jj(3)), 1, 2)
    expect = dconjg(fw1(dir_inv(ii(1)), 1, 1))*v1 &
           + dconjg(fw1(dir_inv(ii(1)), 1, 2))*v2
    call chk('analytic value at (i=1,j=3) orientation', &
             abs(c1(1, 3) - expect), 1d-14)
  end block

  ! ---------------- test 2: random data, packed vs production loops --------
  call random_seed(size=seed)
  call random_seed(put=(/(k*7919 + 13, k=1, seed)/))
  call fill_random()
  fin_TF = 1; ini_TF = 1; dir_active = 1
  call run_compare('full support (all 9 R, R'')', 0)

  ! ---------------- test 3: partial support -------------------------------
  dir_active = 1
  dir_active(3) = 0; dir_active(dir_inv(3)) = 0
  dir_active(7) = 0; dir_active(dir_inv(7)) = 0
  fin_TF = 1
  fin_TF(2) = 0; fin_TF(5) = 0
  ini_TF = 1
  ini_TF(6) = 0; ini_TF(8) = 0
  call run_compare('masked support (inactive dirs + cal_TF)', 0)

  ! ---------------- test 4: chunk invariance ------------------------------
  fin_TF = 1; ini_TF = 1; dir_active = 1
  call run_compare('full support, chunk=1', 1)
  call run_compare('full support, chunk=7', 7)
  call run_compare('full support, chunk=64', 64)
  call run_compare('full support, chunk=257', 257)

  print *, '----------------------------------------'
  if (nfail == 0) then
    print *, 'KAT RESULT: PASS'
  else
    print *, 'KAT RESULT: FAIL nfail=', nfail
    stop 1
  end if

contains

  subroutine chk(label, e, tol)
    character(len=*), intent(in) :: label
    real(dp), intent(in) :: e, tol
    if (e <= tol) then
      write (*, '(A,A,A,ES10.3)') '  PASS ', label, ' err=', e
    else
      write (*, '(A,A,A,ES10.3)') '  FAIL ', label, ' err=', e
      nfail = nfail + 1
    end if
  end subroutine chk

  subroutine build_dirs()
    integer :: a, b, c, n
    n = 1
    do a = -1, 1
      do b = -1, 1
        do c = 0, 0
          dir_vec(n, :) = (/a, b, c/)
          n = n + 1
        end do
      end do
    end do
    do a = 1, d_ext
      do b = 1, d_ext
        if (all(dir_vec(a, :) .eq. -dir_vec(b, :))) dir_inv(a) = b
      end do
    end do
  end subroutine build_dirs

  subroutine fill_random()
    real(dp) :: r(4)
    do k = 1, d_ext
      do i = 1, ni
        call random_number(r)
        fin_wan(k, i, 1) = dcmplx(r(1) - 0.5d0, r(2) - 0.5d0)
        fin_wan(k, i, 2) = dcmplx(r(3) - 0.5d0, r(4) - 0.5d0)
        call random_number(r)
        ini_wan(k, i, 1) = dcmplx(r(1) - 0.5d0, r(2) - 0.5d0)
        ini_wan(k, i, 2) = dcmplx(r(3) - 0.5d0, r(4) - 0.5d0)
        do jd = 1, ndep
          do sc = 1, nsc
            call random_number(r)
            B4(i, 1, jd, sc) = 2.0d0*(r(1) - 0.5d0)
            B4(i, 2, jd, sc) = 2.0d0*(r(2) - 0.5d0)
            B4(i, 3, jd, sc) = 2.0d0*(r(3) - 0.5d0)
          end do
        end do
      end do
    end do
  end subroutine fill_random

  subroutine run_compare(label, chunk)
    character(len=*), intent(in) :: label
    integer, intent(in) :: chunk
    integer :: ch, cbase, ncb, blk
    integer :: ich(8)
    maxerr = 0.0d0
    denom = 0.0d0
    call build_support(d_ext, dir_active, dir_inv, fin_TF, ini_TF, n_i, idx_i, n_j, idx_j)
    if (allocated(C)) deallocate (C)
    allocate (C(max(1, n_i), max(1, n_j)*nfunc))
    C = (0.0d0, 0.0d0)
    A = (0.0d0, 0.0d0); I1 = (0.0d0, 0.0d0); I2 = (0.0d0, 0.0d0)
    if (n_i > 0) call pack_bra(fin_wan, d_ext, ni, dir_inv, n_i, idx_i, A)
    if (n_j > 0) call pack_ket(ini_wan, d_ext, ni, dir_inv, n_j, idx_j, I1, I2)
    do ch = 1, nfunc
      full_ref = (0.0d0, 0.0d0)
      do j = 1, d_ext
        if (dir_active(j) == 0) cycle
        if (dir_active(dir_inv(j)) == 0) cycle
        if (ini_TF(dir_inv(j)) /= 1) cycle
        do i = 1, d_ext
          if (dir_active(i) == 0) cycle
          if (dir_active(dir_inv(i)) == 0) cycle
          if (.not. (fin_TF(dir_inv(i)) == 1 .and. ini_TF(dir_inv(j)) == 1)) cycle
          call prod_element(ch, i, j, full_ref(i, j))
        end do
      end do
      ! packed path, in blocks of up to 8 channels
      do cbase = 1, nfunc, 8
        ncb = min(8, nfunc - cbase + 1)
        do blk = 1, ncb
          ich(blk) = cbase + blk - 1
        end do
        call apply_channels(ni, n_i, n_j, A, I1, I2, B4, ndep, nsc, ich, &
                            ncb, max(1, chunk), C)
        do blk = 1, ncb
          if (cbase + blk - 1 /= ch) cycle
          rtmp = (0.0d0, 0.0d0)
          do i = 1, n_i
            do j = 1, n_j
              rtmp(idx_i(i), idx_j(j)) = C(i, (blk - 1)*n_j + j)
            end do
          end do
          do j = 1, d_ext
            do i = 1, d_ext
              err = abs(rtmp(i, j) - full_ref(i, j))
              maxerr = max(maxerr, err)
              denom = max(denom, abs(full_ref(i, j)))
              if (rtmp(i, j) /= (0.0d0, 0.0d0) .and. full_ref(i, j) == (0.0d0, 0.0d0)) then
                print *, '  FAIL support leak at', ch, i, j
                nfail = nfail + 1
              end if
            end do
          end do
        end do
      end do
    end do
    call chk(label, maxerr/max(denom, 1d-30), 1d-13)
  end subroutine run_compare

  subroutine prod_element(ch, i, j, val)
    integer, intent(in) :: ch, i, j
    complex(dp), intent(out) :: val
    integer :: jd2, sc2, ig, iv, jv
    complex(dp) :: vv1(ni), vv2(ni), q1, q2
    jd2 = (ch - 1)/nsc + 1
    sc2 = mod(ch - 1, nsc) + 1
    jv = dir_inv(j)
    iv = dir_inv(i)
    do ig = 1, ni
      q1 = dcmplx(B4(ig, 1, jd2, sc2), -B4(ig, 2, jd2, sc2))
      q2 = dcmplx(B4(ig, 1, jd2, sc2), B4(ig, 2, jd2, sc2))
      vv1(ig) = q1*ini_wan(jv, ig, 2) + dcmplx(B4(ig, 3, jd2, sc2), 0d0)*ini_wan(jv, ig, 1)
      vv2(ig) = q2*ini_wan(jv, ig, 1) - dcmplx(B4(ig, 3, jd2, sc2), 0d0)*ini_wan(jv, ig, 2)
    end do
    val = sum(dconjg(fin_wan(iv, :, 1))*vv1) + sum(dconjg(fin_wan(iv, :, 2))*vv2)
  end subroutine prod_element

end program test_known_answer
