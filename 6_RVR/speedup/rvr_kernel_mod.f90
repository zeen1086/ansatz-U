! RVR contraction kernel: packed contiguous columns + ZGEMM over grid chunks.
! Definition preserved from cal_RVR_split.f90::cal_each_func:
!
!              + sum_ig conjg(fin_wan(iv,ig,2))*Vpsi_2(ig)
!
! Matrix form: C = A^T Y with the grid/spin index k = 2*(ig-1)+s packed
! Only the k order of the summation changes (BLAS blocking), not the algebra.
module rvr_kernel_mod
  use parameter_mod
  implicit none
  integer, parameter :: dp = kind(0d0)
  complex(dp), parameter :: cone = (1.0d0, 0.0d0)
  complex(dp), parameter :: czero = (0.0d0, 0.0d0)
  external zgemm
  real(8) :: t_build = 0.0d0, t_gemm = 0.0d0
contains
  real(8) function wall_clock()
    implicit none
    integer(8) :: c, r
    call system_clock(c, r)
    wall_clock = real(c, 8)/real(r, 8)
  end function wall_clock

  ! public helpers
  real(8) function rvr_kernel_build_time()
    rvr_kernel_build_time = t_build
  end function rvr_kernel_build_time
  real(8) function rvr_kernel_gemm_time()
    rvr_kernel_gemm_time = t_gemm
  end function rvr_kernel_gemm_time
  subroutine rvr_kernel_reset_timers()
    t_build = 0.0d0
    t_gemm = 0.0d0
  end subroutine rvr_kernel_reset_timers

  ! Support lists exactly as in cal_each_func: index d is used only if
  subroutine build_support(d_ext, dir_active, dir_inv, fin_TF, ini_TF, &
                           n_i, idx_i, n_j, idx_j)
    integer, intent(in)  :: d_ext
    integer, intent(in)  :: dir_active(d_ext), dir_inv(d_ext)
    integer, intent(in)  :: fin_TF(d_ext), ini_TF(d_ext)
    integer, intent(out) :: n_i, n_j
    integer, intent(out) :: idx_i(d_ext), idx_j(d_ext)
    integer :: d, dv

    n_i = 0
    do d = 1, d_ext
      if (dir_active(d) == 0) cycle
      dv = dir_inv(d)
      if (dir_active(dv) == 0) cycle
      if (fin_TF(dv) /= 1) cycle
      n_i = n_i + 1
      idx_i(n_i) = d
    end do

    n_j = 0
    do d = 1, d_ext
      if (dir_active(d) == 0) cycle
      dv = dir_inv(d)
      if (dir_active(dv) == 0) cycle
      if (ini_TF(dv) /= 1) cycle
      n_j = n_j + 1
      idx_j(n_j) = d
    end do
  end subroutine build_support

  ! Pack conjg(fin_wan) columns once per iwan row.
  subroutine pack_bra(fin_wan, d_ext, ni, dir_inv, n_i, idx_i, A)
    integer, intent(in) :: d_ext, ni, n_i, idx_i(n_i), dir_inv(d_ext)
    complex(dp), intent(in) :: fin_wan(d_ext, ni, 2)
    complex(dp), intent(out) :: A(2*ni, n_i)
    integer :: m, ig, iv

    do m = 1, n_i
      iv = dir_inv(idx_i(m))
      do ig = 1, ni
        A(2*ig-1, m) = dconjg(fin_wan(iv, ig, 1))
        A(2*ig,   m) = dconjg(fin_wan(iv, ig, 2))
      end do
    end do
  end subroutine pack_bra

  ! Pack the two spin components of ini_wan(jv) once per jwan column.
  subroutine pack_ket(ini_wan, d_ext, ni, dir_inv, n_j, idx_j, I1, I2)
    integer, intent(in) :: d_ext, ni, n_j, idx_j(n_j), dir_inv(d_ext)
    complex(dp), intent(in) :: ini_wan(d_ext, ni, 2)
    complex(dp), intent(out) :: I1(ni, n_j), I2(ni, n_j)
    integer :: n, ig, jv

    do n = 1, n_j
      jv = dir_inv(idx_j(n))
      do ig = 1, ni
        I1(ig, n) = ini_wan(jv, ig, 1)
        I2(ig, n) = ini_wan(jv, ig, 2)
      end do
    end do
  end subroutine pack_ket

  ! A block of ansatz channels: C(:, (b-1)*n_j+1 : b*n_j) holds the
  ! (n_i x n_j) block of channel ichan(b).  Chunk loop outside the channel
  ! loop so that the packed bra and ket chunks are reused from cache.
  subroutine apply_channels(ni, n_i, n_j, A, I1, I2, B4, ndep, nsc, ichan, &
                            ncb, chunk, C)
    integer, intent(in) :: ni, n_i, n_j, ncb, chunk, ndep, nsc, ichan(ncb)
    complex(dp), intent(in)  :: A(2*ni, n_i), I1(ni, n_j), I2(ni, n_j)
    real(dp),    intent(in)  :: B4(ni, 3, ndep, nsc)
    complex(dp), intent(out) :: C(n_i, n_j*ncb)
    complex(dp), allocatable :: Y(:, :), p1(:), p2(:), p3(:)
    complex(dp) :: beta
    integer :: k0, k1, nk, ig, j, b, jd, sc, col
    real(8) :: tw

    C = czero
    if (n_i == 0 .or. n_j == 0 .or. ni == 0) return

    allocate (Y(2*chunk, n_j*ncb), p1(chunk), p2(chunk), p3(chunk))
    beta = czero
    do k0 = 1, ni, chunk
      k1 = min(k0 + chunk - 1, ni)
      nk = k1 - k0 + 1
      tw = wall_clock()
      do b = 1, ncb
        jd = (ichan(b) - 1)/nsc + 1
        sc = mod(ichan(b) - 1, nsc) + 1
        do ig = k0, k1
          p1(ig - k0 + 1) = dcmplx(B4(ig, 1, jd, sc), -B4(ig, 2, jd, sc))
          p2(ig - k0 + 1) = dcmplx(B4(ig, 1, jd, sc),  B4(ig, 2, jd, sc))
          p3(ig - k0 + 1) = dcmplx(B4(ig, 3, jd, sc),  0.0d0)
        end do
        do j = 1, n_j
          col = (b - 1)*n_j + j
          do ig = 1, nk
            Y(2*ig - 1, col) = p1(ig)*I2(k0 + ig - 1, j) + p3(ig)*I1(k0 + ig - 1, j)
            Y(2*ig,     col) = p2(ig)*I1(k0 + ig - 1, j) - p3(ig)*I2(k0 + ig - 1, j)
          end do
        end do
      end do
      t_build = t_build + wall_clock() - tw
      tw = wall_clock()
      call zgemm('T', 'N', n_i, n_j*ncb, 2*nk, cone, A(2*k0 - 1, 1), 2*ni, &
                 Y, 2*chunk, beta, C, n_i)
      t_gemm = t_gemm + wall_clock() - tw
      beta = cone
    end do
    deallocate (Y, p1, p2, p3)
  end subroutine apply_channels

  ! Scalar reference kernel with the exact cal_each_func loop order and
  ! summation order (V_psi built per j column, then two sums over the grid).
  subroutine apply_channel_scalar(ni, n_i, n_j, A, I1, I2, B4, jc, sc, ndep, &
                                  nsc, C)
    integer, intent(in) :: ni, n_i, n_j, jc, sc, ndep, nsc
    complex(dp), intent(in)  :: A(2*ni, n_i), I1(ni, n_j), I2(ni, n_j)
    real(dp),    intent(in)  :: B4(ni, 3, ndep, nsc)
    complex(dp), intent(out) :: C(n_i, n_j)
    complex(dp), allocatable :: v1(:), v2(:)
    complex(dp) :: p1, p2, p3, s1, s2
    integer :: ig, i, j

    C = czero
    allocate (v1(ni), v2(ni))
    do j = 1, n_j
      do ig = 1, ni
        p1 = dcmplx(B4(ig, 1, jc, sc), -B4(ig, 2, jc, sc))
        p2 = dcmplx(B4(ig, 1, jc, sc),  B4(ig, 2, jc, sc))
        p3 = dcmplx(B4(ig, 3, jc, sc),  0.0d0)
        v1(ig) = p1*I2(ig, j) + p3*I1(ig, j)
        v2(ig) = p2*I1(ig, j) - p3*I2(ig, j)
      end do
      do i = 1, n_i
        s1 = czero
        s2 = czero
        do ig = 1, ni
          s1 = s1 + A(2*ig - 1, i)*v1(ig)
          s2 = s2 + A(2*ig,     i)*v2(ig)
        end do
        C(i, j) = s1 + s2
      end do
    end do
    deallocate (v1, v2)
  end subroutine apply_channel_scalar

end module rvr_kernel_mod
