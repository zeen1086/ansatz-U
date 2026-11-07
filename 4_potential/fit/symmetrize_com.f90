program symmetrize_com
  use parameter_mod
  use cri3_block_mask_mod, only: block_mask_mode, make_feature_block_mask,      &
    apply_com_block_mask, count_blocked_features
  implicit none

  integer, parameter :: dp = kind(0d0)
  integer, parameter :: n_spinor = p_fit_nspinor
  integer, parameter :: n_features_expected = 2 * n_spinor
  integer, parameter :: max_symmetry_passes = p_max_symmetry_passes
  real(dp), parameter :: symmetry_tolerance = p_symmetry_tolerance
  character(len=32), parameter :: expected_magic = 'CRI3_COM_FIT_CHANNEL_V1'
  character(len=*), parameter :: default_input_file = 'V_func_com_ortho_channel_v1'
  character(len=*), parameter :: default_output_file = 'V_func_com_ortho_channel_v1_sym'
  character(len=*), parameter :: report_file = 'symmetrize_com_report.dat'

  integer :: nx, ny, nz, ni
  integer :: n_spinor_file, n_features_file, n_file, n_kind
  integer :: max_iter, nonconverged_count
  integer :: i, ix, iy, iz
  integer :: pass, passes_used
  integer :: linear_nz, mod_nz
  integer :: blocked_features
  real(dp) :: l1_relative, l2_relative, alpha1, alpha2, tolerance
  real(dp) :: total_linear_before, total_mod_before
  real(dp) :: total_linear_change, total_mod_change
  real(dp) :: pass_linear_change, pass_mod_change, pass_change
  real(dp) :: zero_threshold
  real(dp) :: pt_pos(n_spinor,3)
  real(dp) :: unit_vec(3,3), origin(3)
  real(dp), allocatable :: k_linear(:,:), k_mod(:,:)
  real(dp), allocatable :: original_linear(:,:), original_mod(:,:)
  real(dp), allocatable :: previous_linear(:,:), previous_mod(:,:)
  real(dp), allocatable :: block_mask(:,:)
  integer, allocatable :: vec(:,:)
  character(len=32) :: magic
  character(len=256) :: input_file, output_file
  logical :: same_output

  input_file = default_input_file
  output_file = default_output_file
  call read_com_file(input_file, magic, nx, ny, nz, ni, n_spinor_file,          &
    n_features_file, n_file, n_kind, max_iter, nonconverged_count,              &
    l1_relative, l2_relative, alpha1, alpha2, tolerance, k_linear, k_mod)
  call validate_header(magic, n_spinor_file, n_features_file)
  allocate(block_mask(ni,n_features_expected))
  call make_feature_block_mask(nx, ny, nz, ni, block_mask)
  blocked_features = count_blocked_features(block_mask)
  call apply_com_block_mask(k_linear, k_mod, block_mask)

  allocate(original_linear(ni,n_spinor))
  allocate(original_mod(ni,n_spinor))
  allocate(previous_linear(ni,n_spinor))
  allocate(previous_mod(ni,n_spinor))
  original_linear = k_linear
  original_mod = k_mod

  allocate(vec(ni,3))
  do ix = 0, nx
    do iy = 0, ny
      do iz = 0, nz
        i = index_from_grid(ix, iy, iz, nx, ny)
        vec(i,:) = (/ix, iy, iz/)
      end do
    end do
  end do

  call build_channel_positions(nx, ny, nz, pt_pos)

  passes_used = max_symmetry_passes
  do pass = 1, max_symmetry_passes
    previous_linear = k_linear
    previous_mod = k_mod

    origin = (/real(nx,dp)/3d0, 2d0*real(ny,dp)/3d0, real(nz,dp)/2d0/)
    unit_vec(1,:) = (/0d0, 1d0, 0d0/)
    unit_vec(2,:) = (/-1d0, -1d0, 0d0/)
    unit_vec(3,:) = (/0d0, 0d0, 1d0/)
    call apply_operation('symmetrize2_local_c3_A', nx, ny, nz, ni, vec, pt_pos, &
      origin, unit_vec, block_mask, k_linear, k_mod, pass == 1)

    origin = (/2d0*real(nx,dp)/3d0, real(ny,dp)/3d0, real(nz,dp)/2d0/)
    unit_vec(1,:) = (/0d0, 1d0, 0d0/)
    unit_vec(2,:) = (/-1d0, -1d0, 0d0/)
    unit_vec(3,:) = (/0d0, 0d0, 1d0/)
    call apply_operation('symmetrize3_local_c3_B', nx, ny, nz, ni, vec, pt_pos, &
      origin, unit_vec, block_mask, k_linear, k_mod, pass == 1)

    origin = (/0d0, 0d0, 0d0/)
    unit_vec(1,:) = (/0d0, 1d0, 0d0/)
    unit_vec(2,:) = (/1d0, 0d0, 0d0/)
    unit_vec(3,:) = (/0d0, 0d0, 1d0/)
    call apply_operation('symmetrize4_xy_mirror', nx, ny, nz, ni, vec, pt_pos,  &
      origin, unit_vec, block_mask, k_linear, k_mod, pass == 1)
    call apply_com_block_mask(k_linear, k_mod, block_mask)

    pass_linear_change = relative_error(sum((k_linear - previous_linear)**2),  &
      sum(previous_linear**2))
    pass_mod_change = relative_error(sum((k_mod - previous_mod)**2),           &
      sum(previous_mod**2))
    pass_change = max(pass_linear_change, pass_mod_change)

    if (pass == 1 .or. mod(pass, 5) == 0 .or.                                &
        pass_change < symmetry_tolerance) then
      write(*,'(A,I0,A,ES12.4,A,ES12.4)') 'symmetry pass ', pass,              &
        ': linear_change = ', pass_linear_change, ', mod_change = ',           &
        pass_mod_change
    end if

    if (pass_change < symmetry_tolerance) then
      passes_used = pass
      exit
    end if
  end do

  total_linear_before = sum(original_linear**2)
  total_mod_before = sum(original_mod**2)
  total_linear_change = sum((k_linear - original_linear)**2)
  total_mod_change = sum((k_mod - original_mod)**2)

  zero_threshold = 1d-12 * max(1d-300, max(maxval(abs(k_linear)),              &
    maxval(abs(k_mod))))
  linear_nz = count(abs(k_linear) > zero_threshold)
  mod_nz = count(abs(k_mod) > zero_threshold)

  same_output = trim(input_file) == trim(output_file)
  call apply_com_block_mask(k_linear, k_mod, block_mask)
  call write_com_file(output_file, magic, nx, ny, nz, ni, n_spinor_file,        &
    n_features_file, n_file, n_kind, max_iter, nonconverged_count,              &
    l1_relative, l2_relative, alpha1, alpha2, tolerance, k_linear, k_mod)
  call write_report(input_file, output_file, same_output, nx, ny, nz, ni,       &
    total_linear_before, total_mod_before, total_linear_change,                 &
    total_mod_change, zero_threshold, linear_nz, mod_nz, passes_used,           &
    pass_change, k_linear, k_mod, blocked_features)

  write(*,'(A)') 'Compressive-kernel symmetrization completed.'
  write(*,'(A,A)') 'Input:  ', trim(input_file)
  write(*,'(A,A)') 'Output: ', trim(output_file)
  write(*,'(A,ES16.8)') 'linear_relative_change = ',                            &
    relative_error(total_linear_change, total_linear_before)
  write(*,'(A,ES16.8)') 'mod_relative_change = ',                               &
    relative_error(total_mod_change, total_mod_before)
  write(*,'(A,I0)') 'symmetry_passes_used = ', passes_used
  write(*,'(A,ES16.8)') 'final_pass_change = ', pass_change
  write(*,'(A,A)') 'block mask mode = ', block_mask_mode
  write(*,'(A,I0)') 'blocked feature coefficients = ', blocked_features
  write(*,'(A,A)') 'Report: ', report_file

contains

  subroutine validate_header(file_magic, spinor_count, feature_count)
    implicit none
    character(len=*), intent(in) :: file_magic
    integer, intent(in) :: spinor_count, feature_count

    if (trim(file_magic) /= trim(expected_magic)) then
      write(*,'(A,A)') 'ERROR: unsupported V_func_com magic: ', trim(file_magic)
      stop 1
    end if
    if (spinor_count /= n_spinor) then
      write(*,'(A,I0)') 'ERROR: expected six spin channels, got ', spinor_count
      stop 1
    end if
    if (feature_count /= n_features_expected) then
      write(*,'(A,I0)') 'ERROR: expected twelve vector features, got ',         &
        feature_count
      stop 1
    end if
  end subroutine validate_header

  subroutine read_com_file(path, file_magic, grid_nx, grid_ny, grid_nz,          &
      grid_size, spinor_count, feature_count, file_count, kind_count, n_iter,   &
      nonconv, l1_rel, l2_rel, l1_abs, l2_abs, conv_tol, linear_part, mod_part)
    implicit none
    character(len=*), intent(in) :: path
    character(len=32), intent(out) :: file_magic
    integer, intent(out) :: grid_nx, grid_ny, grid_nz, grid_size
    integer, intent(out) :: spinor_count, feature_count, file_count, kind_count
    integer, intent(out) :: n_iter, nonconv
    real(dp), intent(out) :: l1_rel, l2_rel, l1_abs, l2_abs, conv_tol
    real(dp), allocatable, intent(out) :: linear_part(:,:), mod_part(:,:)

    open(unit=101, file=trim(path), access='sequential', action='read',         &
      status='old', form='unformatted')
    read(101) file_magic
    read(101) grid_nx, grid_ny, grid_nz
    read(101) grid_size, spinor_count, feature_count, file_count, kind_count
    read(101) n_iter, nonconv
    read(101) l1_rel, l2_rel, l1_abs, l2_abs, conv_tol
    allocate(linear_part(grid_size, n_spinor))
    allocate(mod_part(grid_size, n_spinor))
    read(101) linear_part
    read(101) mod_part
    close(101)
  end subroutine read_com_file

  subroutine write_com_file(path, file_magic, grid_nx, grid_ny, grid_nz,         &
      grid_size, spinor_count, feature_count, file_count, kind_count, n_iter,   &
      nonconv, l1_rel, l2_rel, l1_abs, l2_abs, conv_tol, linear_part, mod_part)
    implicit none
    character(len=*), intent(in) :: path
    character(len=32), intent(in) :: file_magic
    integer, intent(in) :: grid_nx, grid_ny, grid_nz, grid_size
    integer, intent(in) :: spinor_count, feature_count, file_count, kind_count
    integer, intent(in) :: n_iter, nonconv
    real(dp), intent(in) :: l1_rel, l2_rel, l1_abs, l2_abs, conv_tol
    real(dp), intent(in) :: linear_part(grid_size, n_spinor)
    real(dp), intent(in) :: mod_part(grid_size, n_spinor)

    open(unit=102, file=trim(path), access='sequential', action='write',        &
      status='replace', form='unformatted')
    write(102) file_magic
    write(102) grid_nx, grid_ny, grid_nz
    write(102) grid_size, spinor_count, feature_count, file_count, kind_count
    write(102) n_iter, nonconv
    write(102) l1_rel, l2_rel, l1_abs, l2_abs, conv_tol
    write(102) linear_part
    write(102) mod_part
    close(102)
  end subroutine write_com_file

  subroutine build_channel_positions(grid_nx, grid_ny, grid_nz, positions)
    implicit none
    integer, intent(in) :: grid_nx, grid_ny, grid_nz
    real(dp), intent(out) :: positions(n_spinor,3)

    positions(1,:) = (/2d0/3d0, 1d0/3d0, 0.5d0/)
    positions(2,:) = (/1d0/3d0, 2d0/3d0, 0.5d0/)
    positions(3,:) = (/1d0/3d0, -1d0/3d0, 0.5d0/)
    positions(4,:) = (/4d0/3d0, 2d0/3d0, 0.5d0/)
    positions(5,:) = (/-1d0/3d0, 1d0/3d0, 0.5d0/)
    positions(6,:) = (/2d0/3d0, 4d0/3d0, 0.5d0/)

    positions(:,1) = positions(:,1) * real(grid_nx,dp)
    positions(:,2) = positions(:,2) * real(grid_ny,dp)
    positions(:,3) = positions(:,3) * real(grid_nz,dp)
  end subroutine build_channel_positions

  subroutine apply_operation(operation_name, grid_nx, grid_ny, grid_nz,          &
      grid_size, grid_vectors, positions, center, transform, mask, linear_part, &
      mod_part, print_status)
    implicit none
    character(len=*), intent(in) :: operation_name
    integer, intent(in) :: grid_nx, grid_ny, grid_nz, grid_size
    integer, intent(in) :: grid_vectors(grid_size,3)
    real(dp), intent(in) :: positions(n_spinor,3), center(3), transform(3,3)
    real(dp), intent(in) :: mask(grid_size,n_features_expected)
    real(dp), intent(inout) :: linear_part(grid_size,n_spinor)
    real(dp), intent(inout) :: mod_part(grid_size,n_spinor)
    logical, intent(in) :: print_status
    integer :: i, j, channel, mapped_channel
    integer :: ix_map, iy_map, iz_map
    integer :: pt_map(n_spinor), mapped_points, mapped_channels
    integer :: mapped_linear_pairs, mapped_mod_pairs
    integer :: linear_source_active(n_spinor), linear_target_active(n_spinor)
    integer :: linear_active_pairs(n_spinor)
    integer :: mod_source_active(n_spinor), mod_target_active(n_spinor)
    integer :: mod_active_pairs(n_spinor)
    integer, allocatable :: degree_linear(:,:), degree_mod(:,:)
    real(dp) :: mapped_position(3), mapped_grid(3)
    real(dp) :: mapped_positions(n_spinor,3)
    real(dp), allocatable :: old_linear(:,:), old_mod(:,:)
    real(dp), allocatable :: tmp_linear(:,:), tmp_mod(:,:)

    allocate(old_linear(grid_size,n_spinor))
    allocate(old_mod(grid_size,n_spinor))
    allocate(tmp_linear(grid_size,n_spinor))
    allocate(tmp_mod(grid_size,n_spinor))
    allocate(degree_linear(grid_size,n_spinor))
    allocate(degree_mod(grid_size,n_spinor))

    old_linear = linear_part
    old_mod = mod_part
    tmp_linear = old_linear
    tmp_mod = old_mod
    degree_linear = 1
    degree_mod = 1
    pt_map = n_spinor + 1
    mapped_positions = 0d0

    do channel = 1, n_spinor
      mapped_position = center + matmul(transpose(transform),                   &
        positions(channel,:) - center)
      mapped_positions(channel,:) = mapped_position
      do mapped_channel = 1, n_spinor
        if (sum((mapped_position - positions(mapped_channel,:))**2) < 1d-3) then
          pt_map(channel) = mapped_channel
        end if
      end do
    end do

    mapped_points = 0
    mapped_channels = count(pt_map <= n_spinor)
    mapped_linear_pairs = 0
    mapped_mod_pairs = 0
    linear_source_active = 0
    linear_target_active = 0
    linear_active_pairs = 0
    mod_source_active = 0
    mod_target_active = 0
    mod_active_pairs = 0
    do i = 1, grid_size
      mapped_grid = center + matmul(transpose(transform),                       &
        real(grid_vectors(i,:),dp) - center)
      ix_map = nint(mapped_grid(1))
      iy_map = nint(mapped_grid(2))
      iz_map = nint(mapped_grid(3))

      if (abs(mapped_grid(1) - real(ix_map,dp)) > 1d-8) cycle
      if (abs(mapped_grid(2) - real(iy_map,dp)) > 1d-8) cycle
      if (abs(mapped_grid(3) - real(iz_map,dp)) > 1d-8) cycle
      if (ix_map < 0 .or. ix_map > grid_nx) cycle
      if (iy_map < 0 .or. iy_map > grid_ny) cycle
      if (iz_map < 0 .or. iz_map > grid_nz) cycle

      mapped_points = mapped_points + 1
      j = index_from_grid(ix_map, iy_map, iz_map, grid_nx, grid_ny)
      do channel = 1, n_spinor
        mapped_channel = pt_map(channel)
        if (mapped_channel > n_spinor) cycle
        if (mask(i,channel) > 0.5d0) then
          linear_source_active(channel) = linear_source_active(channel) + 1
        end if
        if (mask(j,mapped_channel) > 0.5d0) then
          linear_target_active(channel) = linear_target_active(channel) + 1
        end if
        if (mask(i,n_spinor+channel) > 0.5d0) then
          mod_source_active(channel) = mod_source_active(channel) + 1
        end if
        if (mask(j,n_spinor+mapped_channel) > 0.5d0) then
          mod_target_active(channel) = mod_target_active(channel) + 1
        end if
        if (mask(i,channel) > 0.5d0 .and. mask(j,mapped_channel) > 0.5d0) then
          tmp_linear(j,mapped_channel) = tmp_linear(j,mapped_channel)           &
            + old_linear(i,channel)
          degree_linear(j,mapped_channel) = degree_linear(j,mapped_channel) + 1
          mapped_linear_pairs = mapped_linear_pairs + 1
          linear_active_pairs(channel) = linear_active_pairs(channel) + 1
        end if
        if (mask(i,n_spinor+channel) > 0.5d0 .and.                            &
            mask(j,n_spinor+mapped_channel) > 0.5d0) then
          tmp_mod(j,mapped_channel) = tmp_mod(j,mapped_channel)                 &
            + old_mod(i,channel)
          degree_mod(j,mapped_channel) = degree_mod(j,mapped_channel) + 1
          mapped_mod_pairs = mapped_mod_pairs + 1
          mod_active_pairs(channel) = mod_active_pairs(channel) + 1
        end if
      end do
    end do

    do channel = 1, n_spinor
      tmp_linear(:,channel) = tmp_linear(:,channel)                             &
        / real(degree_linear(:,channel),dp)
      tmp_mod(:,channel) = tmp_mod(:,channel) / real(degree_mod(:,channel),dp)
    end do

    linear_part = tmp_linear
    mod_part = tmp_mod
    if (print_status) then
      write(*,'(A,A,A,I0,A,I0,A,I0,A,I0)') 'Applied ', trim(operation_name),    &
        ': mapped_points = ', mapped_points, ', mapped_channels = ',            &
        mapped_channels, ', active_linear_pairs = ', mapped_linear_pairs,       &
        ', active_mod_pairs = ', mapped_mod_pairs
      call print_operation_diagnostics(operation_name, pt_map, mapped_positions,&
        mapped_points, linear_source_active, linear_target_active,              &
        linear_active_pairs, mod_source_active, mod_target_active,              &
        mod_active_pairs)
    end if

    deallocate(old_linear, old_mod, tmp_linear, tmp_mod, degree_linear, degree_mod)
  end subroutine apply_operation

  subroutine print_operation_diagnostics(operation_name, pt_map,                &
      mapped_positions, mapped_points, linear_source_active,                    &
      linear_target_active, linear_active_pairs, mod_source_active,             &
      mod_target_active, mod_active_pairs)
    implicit none
    character(len=*), intent(in) :: operation_name
    integer, intent(in) :: pt_map(n_spinor), mapped_points
    real(dp), intent(in) :: mapped_positions(n_spinor,3)
    integer, intent(in) :: linear_source_active(n_spinor)
    integer, intent(in) :: linear_target_active(n_spinor)
    integer, intent(in) :: linear_active_pairs(n_spinor)
    integer, intent(in) :: mod_source_active(n_spinor)
    integer, intent(in) :: mod_target_active(n_spinor)
    integer, intent(in) :: mod_active_pairs(n_spinor)
    integer :: channel

    write(*,'(A,A)') '  channel-map diagnostics for ', trim(operation_name)
    write(*,'(A)') '  ch mapped_ch status mapped_x mapped_y mapped_z'
    do channel = 1, n_spinor
      if (pt_map(channel) <= n_spinor) then
        write(*,'(2X,I2,1X,I9,1X,A,3(1X,F12.6))') channel, pt_map(channel),     &
          'mapped ', mapped_positions(channel,1), mapped_positions(channel,2),  &
          mapped_positions(channel,3)
      else
        write(*,'(2X,I2,1X,I9,1X,A,3(1X,F12.6))') channel, 0,                  &
          'missing', mapped_positions(channel,1), mapped_positions(channel,2),  &
          mapped_positions(channel,3)
      end if
    end do

    write(*,'(A,I0)') '  mapped grid points per valid channel = ', mapped_points
    write(*,'(A)') '  ch linear_src linear_tgt linear_both mod_src mod_tgt mod_both'
    do channel = 1, n_spinor
      write(*,'(2X,I2,6(1X,I10))') channel, linear_source_active(channel),       &
        linear_target_active(channel), linear_active_pairs(channel),            &
        mod_source_active(channel), mod_target_active(channel),                 &
        mod_active_pairs(channel)
    end do
  end subroutine print_operation_diagnostics

  integer function index_from_grid(ix, iy, iz, grid_nx, grid_ny)
    implicit none
    integer, intent(in) :: ix, iy, iz, grid_nx, grid_ny

    index_from_grid = 1 + ix + (grid_nx + 1) * iy                               &
      + (grid_nx + 1) * (grid_ny + 1) * iz
  end function index_from_grid

  subroutine write_report(input_path, output_path, is_inplace, grid_nx, grid_ny,&
      grid_nz, grid_size, linear_before, mod_before, linear_change, mod_change, &
      sparse_threshold, linear_nz, mod_nz, passes_used, final_pass_change,      &
      linear_part, mod_part, n_blocked)
    implicit none
    character(len=*), intent(in) :: input_path, output_path
    logical, intent(in) :: is_inplace
    integer, intent(in) :: grid_nx, grid_ny, grid_nz, grid_size
    integer, intent(in) :: linear_nz, mod_nz, passes_used, n_blocked
    real(dp), intent(in) :: linear_before, mod_before, linear_change, mod_change
    real(dp), intent(in) :: sparse_threshold, final_pass_change
    real(dp), intent(in) :: linear_part(grid_size,n_spinor)
    real(dp), intent(in) :: mod_part(grid_size,n_spinor)
    integer :: channel

    open(unit=103, file=report_file, action='write', status='replace')
    write(103,'(A)') 'CrI3 compressive-kernel symmetrization'
    write(103,'(A,A)') 'input_file = ', trim(input_path)
    write(103,'(A,A)') 'output_file = ', trim(output_path)
    write(103,'(A,L1)') 'inplace = ', is_inplace
    write(103,'(A)') 'Applied operations: symmetrize2, symmetrize3, symmetrize4'
    write(103,'(A)') 'Projection mode: mask-aware iterative pair averaging'
    write(103,'(A,I0)') 'max_symmetry_passes = ', max_symmetry_passes
    write(103,'(A,I0)') 'symmetry_passes_used = ', passes_used
    write(103,'(A,ES24.16)') 'symmetry_tolerance = ', symmetry_tolerance
    write(103,'(A,ES24.16)') 'final_pass_change = ', final_pass_change
    write(103,'(A)') 'Action: full per-channel mapping for K_linear and K_mod'
    write(103,'(A)') 'Spin-tensor rotation/projectors are intentionally omitted'
    write(103,'(A,A)') 'block mask mode = ', block_mask_mode
    write(103,'(A,I0)') 'blocked feature coefficients = ', n_blocked
    write(103,'(A,3(I0,1X))') 'grid_intervals = ', grid_nx, grid_ny, grid_nz
    write(103,'(A,I0)') 'grid_points = ', grid_size
    write(103,'(A,ES24.16)') 'linear_relative_change = ',                     &
      relative_error(linear_change, linear_before)
    write(103,'(A,ES24.16)') 'mod_relative_change = ',                         &
      relative_error(mod_change, mod_before)
    write(103,'(A,ES24.16)') 'coefficient_zero_threshold = ', sparse_threshold
    write(103,'(A,I0)') 'linear_nonzero_coefficients = ', linear_nz
    write(103,'(A,I0)') 'mod_nonzero_coefficients = ', mod_nz
    write(103,'(/,A)') 'channel linear_rms mod_rms linear_max_abs mod_max_abs'
    do channel = 1, n_spinor
      write(103,'(I4,4(1X,ES24.16))') channel,                                 &
        sqrt(sum(linear_part(:,channel)**2) / real(grid_size,dp)),              &
        sqrt(sum(mod_part(:,channel)**2) / real(grid_size,dp)),                 &
        maxval(abs(linear_part(:,channel))), maxval(abs(mod_part(:,channel)))
    end do
    close(103)
  end subroutine write_report

  real(dp) function relative_error(squared_error, squared_target)
    implicit none
    real(dp), intent(in) :: squared_error, squared_target

    if (squared_target > 0d0) then
      relative_error = sqrt(max(0d0, squared_error) / squared_target)
    else if (squared_error == 0d0) then
      relative_error = 0d0
    else
      relative_error = huge(1d0)
    end if
  end function relative_error

end program symmetrize_com
