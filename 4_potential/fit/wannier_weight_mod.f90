module wannier_weight_mod
  use parameter_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private

  integer, parameter :: dp = kind(0d0)
  character(len=32), parameter :: weight_magic = 'CRI3_WANNIER_WEIGHT_V1'

  public :: load_wannier_weight

contains

  subroutine load_wannier_weight(path, nx, ny, nz, ni, weight)
    character(len=*), intent(in) :: path
    integer, intent(in) :: nx, ny, nz, ni
    real(dp), intent(out) :: weight(ni)

    character(len=32) :: magic
    integer :: file_nx, file_ny, file_nz, file_ni
    integer :: wan_min, wan_max, nspin, n_translation
    real(dp) :: weight_floor, mean_weight
    logical :: exists

    inquire(file=trim(path), exist=exists)
    if (.not. exists) then
      write(*,'(A,A)') 'ERROR: missing Wannier weight file: ', trim(path)
      write(*,'(A)') 'Run 5_wannier/make_wan/make_w_density.f90 first.'
      stop 1
    end if

    open(unit=97, file=trim(path), access='sequential', action='read',          &
      status='old', form='unformatted')
    read(97) magic
    read(97) file_nx, file_ny, file_nz, file_ni
    read(97) wan_min, wan_max, nspin, n_translation, weight_floor
    if (trim(magic) /= trim(weight_magic)) then
      write(*,'(A,A)') 'ERROR: unsupported weight magic: ', trim(magic)
      stop 1
    end if
    if (file_nx /= nx .or. file_ny /= ny .or. file_nz /= nz .or.             &
        file_ni /= ni) then
      write(*,'(A,4I10)') 'ERROR: fit grid nx ny nz ni = ', nx, ny, nz, ni
      write(*,'(A,4I10)') '       weight grid nx ny nz ni = ',                 &
        file_nx, file_ny, file_nz, file_ni
      stop 1
    end if
    read(97) weight
    close(97)

    if (minval(weight) <= 0d0 .or. any(.not. ieee_is_finite(weight))) then
      write(*,'(A)') 'ERROR: Wannier weights must be finite and positive.'
      stop 1
    end if

    mean_weight = sum(weight) / real(ni,dp)
    if (mean_weight <= 0d0) stop 'invalid mean Wannier weight'
    weight = weight / mean_weight

    write(*,'(A,A)') 'Wannier weight file = ', trim(path)
    write(*,'(A,2I6,A,I4,A,I4)') 'Wannier bands = ', wan_min, wan_max,          &
      '  nspin = ', nspin, '  translations = ', n_translation
    write(*,'(A,3ES14.6)') 'Wannier weight min/mean/max = ',                   &
      minval(weight), sum(weight)/real(ni,dp), maxval(weight)
    write(*,'(A,ES14.6)') 'Wannier density floor fraction = ', weight_floor
  end subroutine load_wannier_weight

end module wannier_weight_mod
