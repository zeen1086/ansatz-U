program make_w_density
  use parameter_mod
  implicit none

  integer, parameter :: dp = kind(0d0)
  character(len=32), parameter :: weight_magic = 'ANSATZ_WANNIER_WEIGHT_V1'

  integer :: nx, ny, nz, fit_nx, fit_ny, fit_nz, fit_ni
  integer :: wan_min, wan_max, nspin
  integer :: convert_x, convert_y, convert_z, n_translation
  integer :: iband, ix, iy, iz, i_base, i_fit, stat
  real(dp) :: weight_floor, mean_density
  real(dp), allocatable :: density(:), weight(:)
  complex(dp), allocatable :: wan(:,:,:)
  character(len=32) :: band_text
  character(len=256) :: argument, output_file, wan_file
  logical :: exists

  output_file = p_func_dir//'w_density'
  wan_min = p_wan_min
  wan_max = p_wan_max
  weight_floor = p_weight_floor
  nspin = p_nspin
  convert_x = p_convert_size_x
  convert_y = p_convert_size_y
  convert_z = p_convert_size_z

  call read_char_arg(1, output_file)
  call read_int_arg(2, wan_min)
  call read_int_arg(3, wan_max)
  call read_real_arg(4, weight_floor)
  call read_int_arg(5, nspin)
  call read_int_arg(6, convert_x)
  call read_int_arg(7, convert_y)
  call read_int_arg(8, convert_z)

  if (wan_min < 1 .or. wan_max < wan_min) stop 'invalid Wannier band range'
  if (weight_floor <= 0d0 .or. weight_floor >= 1d0) then
    stop 'weight floor must be between zero and one'
  end if
  if (nspin < 1) stop 'nspin must be positive'

  open(unit=101, file=p_scf_dir//'V_HYC_dp', access='sequential',             &
    action='read', status='old', form='unformatted')
  read(101) nx, ny, nz
  close(101)

  open(unit=102, file=p_scf_dir//'V_QE', access='sequential', action='read', &
    status='old', form='unformatted')
  read(102) fit_nx, fit_ny, fit_nz
  close(102)

  if (fit_nx /= nx .or. fit_ny /= ny .or. fit_nz /= nz) then
    write(*,'(A,3I8)') 'V_HYC_dp grid = ', nx, ny, nz
    write(*,'(A,3I8)') 'V_QE grid     = ', fit_nx, fit_ny, fit_nz
    stop 'Wannier and fitting base grids differ'
  end if

  n_translation = (2*convert_x+1)*(2*convert_y+1)*(2*convert_z+1)
  fit_ni = (fit_nx+1)*(fit_ny+1)*(fit_nz+1)
  allocate(wan(n_translation,nx*ny*nz,nspin))
  allocate(density(nx*ny*nz))
  allocate(weight(fit_ni))
  density = 0d0

  do iband = wan_min, wan_max
    write(band_text,'(I0)') iband
    wan_file = 'wan_'//trim(band_text)
    inquire(file=trim(wan_file), exist=exists)
    if (.not. exists) then
      write(*,'(A,A)') 'ERROR: missing Wannier function: ', trim(wan_file)
      stop 1
    end if
    open(unit=103, file=trim(wan_file), access='sequential', action='read',    &
      status='old', form='unformatted')
    read(103) wan
    close(103)
    density = density + sum(sum(abs(wan)**2, dim=1), dim=2)
  end do

  mean_density = sum(density) / real(size(density),dp)
  if (mean_density <= 0d0) stop 'zero Wannier density'
  density = density / mean_density

  do iz = 0, fit_nz
    do iy = 0, fit_ny
      do ix = 0, fit_nx
        i_base = 1 + modulo(ix,nx) + nx*modulo(iy,ny) +                      &
          nx*ny*modulo(iz,nz)
        i_fit = 1 + ix + (fit_nx+1)*iy + (fit_nx+1)*(fit_ny+1)*iz
        weight(i_fit) = weight_floor + (1d0-weight_floor)*density(i_base)
      end do
    end do
  end do
  weight = weight / (sum(weight)/real(fit_ni,dp))

  open(unit=104, file=trim(output_file), access='sequential', action='write', &
    status='replace', form='unformatted')
  write(104) weight_magic
  write(104) fit_nx, fit_ny, fit_nz, fit_ni
  write(104) wan_min, wan_max, nspin, n_translation, weight_floor
  write(104) weight
  close(104)

  write(*,'(A,A)') 'Wrote Wannier weight file: ', trim(output_file)
  write(*,'(A,3I8,A,I12)') 'Fit grid = ', fit_nx, fit_ny, fit_nz,             &
    '  points = ', fit_ni
  write(*,'(A,2I6)') 'Wannier band range = ', wan_min, wan_max
  write(*,'(A,3ES14.6)') 'Weight min/mean/max = ', minval(weight),           &
    sum(weight)/real(fit_ni,dp), maxval(weight)

contains

  subroutine read_char_arg(index, value)
    integer, intent(in) :: index
    character(len=*), intent(inout) :: value
    call get_command_argument(index, argument, status=stat)
    if (stat == 0 .and. len_trim(argument) > 0) value = trim(argument)
  end subroutine read_char_arg

  subroutine read_int_arg(index, value)
    integer, intent(in) :: index
    integer, intent(inout) :: value
    call get_command_argument(index, argument, status=stat)
    if (stat == 0 .and. len_trim(argument) > 0) then
      read(argument,*,iostat=stat) value
      if (stat /= 0) stop 'invalid integer argument'
    end if
  end subroutine read_int_arg

  subroutine read_real_arg(index, value)
    integer, intent(in) :: index
    real(dp), intent(inout) :: value
    call get_command_argument(index, argument, status=stat)
    if (stat == 0 .and. len_trim(argument) > 0) then
      read(argument,*,iostat=stat) value
      if (stat /= 0) stop 'invalid real argument'
    end if
  end subroutine read_real_arg

end program make_w_density
