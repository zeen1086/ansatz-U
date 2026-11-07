module spin_channels_mod
  implicit none
contains
  subroutine spin_chirality(spin, chi, dspin, dchi)
    real(8), intent(in) :: spin(3,6)
    real(8), intent(out) :: chi(6)
    real(8), optional, intent(in) :: dspin(3,6)
    real(8), optional, intent(out) :: dchi(6)
    real(8) :: nj, nk, cosine
    integer :: j, k, count
    if (present(dspin) .neqv. present(dchi)) error stop 'chirality derivative arguments are incomplete'
    chi = 0d0
    if (present(dchi)) dchi = 0d0
    do j = 1, 6
      nj = sqrt(sum(spin(:,j)**2))
      if (nj < 1d-12) cycle
      count = 0
      do k = 1, 6
        if (k == j) cycle
        nk = sqrt(sum(spin(:,k)**2))
        if (nk < 1d-12) cycle
        cosine = sum(spin(:,j)*spin(:,k))/(nj*nk)
        chi(j) = chi(j) + cosine - 1d0
        if (present(dchi)) dchi(j) = dchi(j) &
          + (sum(dspin(:,j)*spin(:,k)) + sum(spin(:,j)*dspin(:,k)))/(nj*nk) &
          - cosine*(sum(spin(:,j)*dspin(:,j))/nj**2 + sum(spin(:,k)*dspin(:,k))/nk**2)
        count = count + 1
      enddo
      if (count > 0) then
        chi(j) = chi(j)/count
        if (present(dchi)) dchi(j) = dchi(j)/count
      endif
    enddo
  end subroutine spin_chirality

  subroutine spin_channels(spin, chi, mask, source, lambda, coefficient, dspin, dchi, derivative)
    real(8), intent(in) :: spin(3,6), chi(6), lambda
    integer, intent(in) :: mask(6), source(6)
    real(8), intent(out) :: coefficient(:)
    real(8), optional, intent(in) :: dspin(3,6), dchi(6)
    real(8), optional, intent(out) :: derivative(:)
    integer :: j, src, base, nf
    real(8) :: factor
    logical :: differentiate

    nf = size(coefficient)
    if (nf /= 18 .and. nf /= 36) error stop 'spin channels require 18 or 36 coefficients'
    if (any(abs(mask) > 1) .or. any(source < 1) .or. any(source > 6)) &
      error stop 'invalid spin channel mask or source'
    differentiate = present(derivative)
    if ((present(dspin) .neqv. differentiate) .or. (present(dchi) .neqv. differentiate)) &
      error stop 'spin derivatives require dspin, dchi and derivative together'
    coefficient = 0d0
    if (differentiate) then
      if (size(derivative) /= nf) error stop 'spin derivative dimension mismatch'
      derivative = 0d0
    endif
    do j = 1, 6
      src = source(j)
      base = 1 + 3*(j-1)
      factor = 1d0
      if (nf == 18) factor = factor + lambda*chi(src)
      coefficient(base:base+2) = mask(j)*factor*spin(:,src)
      if (differentiate) then
        derivative(base:base+2) = mask(j)*factor*dspin(:,src)
        if (nf == 18) derivative(base:base+2) = derivative(base:base+2) &
          + mask(j)*lambda*dchi(src)*spin(:,src)
      endif
      if (nf == 36) then
        coefficient(base+18:base+20) = mask(j)*chi(src)*spin(:,src)
        if (differentiate) derivative(base+18:base+20) = &
          mask(j)*(dchi(src)*spin(:,src) + chi(src)*dspin(:,src))
      endif
    enddo
  end subroutine spin_channels
end module spin_channels_mod
