program check_v_func_tmp_c3
  use parameter_mod
  implicit none

  integer, parameter :: dp = kind(0d0)
  integer, parameter :: n_spinor = p_fit_nspinor
  integer, parameter :: n_component = 3
  real(dp), parameter :: strict_tolerance = p_strict_tolerance
  real(dp), parameter :: offdiag_tolerance = p_offdiag_tolerance

  integer :: ni, grid_x, grid_y, grid_z, n_dep
  integer :: sym_x, sym_y, sym_z
  integer :: i, dep, component
  integer, allocatable :: vec(:,:), index_map(:,:,:)
  real(dp), allocatable :: b_general(:,:,:,:), kernel(:,:)
  real(dp), allocatable :: dep_norm(:)
  real(dp) :: offdiag_norm
  real(dp) :: linear_residual(3), mod_residual(3), edge_residual(3), max_residual
  real(dp) :: positions(n_spinor,3), center(3), transform(3,3)
  character(len=256) :: v_file, arg
  logical :: strict_pass

  v_file = 'V_func_tmp'
  n_dep = 12
  call get_command_argument(1, v_file)
  if (len_trim(v_file) == 0) v_file = 'V_func_tmp'
  sym_x = -1
  sym_y = -1
  sym_z = -1
  call get_command_argument(2, arg)
  if (len_trim(arg) > 0) read(arg,*) sym_x
  call get_command_argument(3, arg)
  if (len_trim(arg) > 0) read(arg,*) sym_y
  call get_command_argument(4, arg)
  if (len_trim(arg) > 0) read(arg,*) sym_z
  call get_command_argument(5, arg)
  if (len_trim(arg) > 0) read(arg,*) n_dep
  if (n_dep /= 6 .and. n_dep /= 12 .and. n_dep /= 22) &
    stop 'n_dep argument must be 6, 12, or 22'

  call read_r_vec('R_vec', ni, vec, grid_x, grid_y, grid_z, index_map)
  call set_default_symmetry_intervals(grid_x, grid_y, grid_z, sym_x, sym_y,   &
    sym_z)
  allocate(b_general(ni,n_component,n_dep,n_component))
  allocate(dep_norm(n_dep))
  call read_v_func_tmp(trim(v_file), ni, b_general)

  allocate(kernel(ni,n_dep))
  kernel = 0d0
  offdiag_norm = 0d0
  do dep = 1, n_dep
    do component = 1, n_component
      kernel(:,dep) = kernel(:,dep) + b_general(:,component,dep,component)
    enddo
    kernel(:,dep) = kernel(:,dep) / real(n_component, dp)
    do component = 1, n_component
      do i = 1, n_component
        if (i == component) cycle
        offdiag_norm = offdiag_norm + sum(b_general(:,component,dep,i)**2)
      enddo
    enddo
  enddo

  do dep = 1, n_dep
    dep_norm(dep) = sum(kernel(:,dep)**2)
  enddo

  open(unit=201, file='check_V_func_tmp_c3_report.dat', status='replace', &
    action='write')
  write(201,'(A)') 'V_func_tmp scalar-kernel diagnostics'
  write(201,'(A,A)') 'input_file = ', trim(v_file)
  write(201,'(A,I0)') 'ni = ', ni
  write(201,'(A,3I8)') 'grid max ix iy iz = ', grid_x, grid_y, grid_z
  write(201,'(A,3I8)') 'symmetry grid intervals = ', sym_x, sym_y, sym_z
  write(201,'(A,ES24.16)') 'offdiag_norm2 = ', offdiag_norm
  write(201,'(A)') ''
  write(201,'(A)') 'dep total_scalar_kernel_norm2'
  do dep = 1, n_dep
    write(201,'(I4,ES24.16)') dep, dep_norm(dep)
  enddo
  write(201,'(A)') ''
  call report_channel_balance(201, dep_norm)

  call build_channel_positions(sym_x, sym_y, sym_z, positions)

  center = (/real(sym_x,dp)/3d0, 2d0*real(sym_y,dp)/3d0, real(sym_z,dp)/2d0/)
  transform(1,:) = (/0d0, 1d0, 0d0/)
  transform(2,:) = (/-1d0, -1d0, 0d0/)
  transform(3,:) = (/0d0, 0d0, 1d0/)
  call report_symmetry_residual(201, 'symmetrize2_local_c3_A', grid_x, grid_y, &
    grid_z, index_map, positions, center, transform, kernel,                  &
    linear_residual(1), mod_residual(1), edge_residual(1))

  center = (/2d0*real(sym_x,dp)/3d0, real(sym_y,dp)/3d0, real(sym_z,dp)/2d0/)
  transform(1,:) = (/0d0, 1d0, 0d0/)
  transform(2,:) = (/-1d0, -1d0, 0d0/)
  transform(3,:) = (/0d0, 0d0, 1d0/)
  call report_symmetry_residual(201, 'symmetrize3_local_c3_B', grid_x, grid_y, &
    grid_z, index_map, positions, center, transform, kernel,                  &
    linear_residual(2), mod_residual(2), edge_residual(2))

  center = (/0d0, 0d0, 0d0/)
  transform(1,:) = (/0d0, 1d0, 0d0/)
  transform(2,:) = (/1d0, 0d0, 0d0/)
  transform(3,:) = (/0d0, 0d0, 1d0/)
  call report_symmetry_residual(201, 'symmetrize4_xy_mirror', grid_x, grid_y, &
    grid_z, index_map, positions, center, transform, kernel,                  &
    linear_residual(3), mod_residual(3), edge_residual(3))

  max_residual = max(maxval(linear_residual), maxval(mod_residual))
  if (n_dep == 22) max_residual = max(max_residual, maxval(edge_residual))
  strict_pass = max_residual <= strict_tolerance .and.                        &
    offdiag_norm <= offdiag_tolerance
  write(201,'(A,ES24.16)') 'strict_C3v_tolerance = ', strict_tolerance
  write(201,'(A,ES24.16)') 'strict_C3v_max_scalar_kernel_residual = ',        &
    max_residual
  write(201,'(A,ES24.16)') 'strict_C3v_offdiag_norm2_tolerance = ',           &
    offdiag_tolerance
  if (strict_pass) then
    write(201,'(A)') 'strict_C3v_status = PASS'
  else
    write(201,'(A)') 'strict_C3v_status = FAIL'
  end if

  close(201)
  write(*,'(A)') 'Wrote check_V_func_tmp_c3_report.dat'
  if (.not. strict_pass) stop 1

contains

  subroutine set_default_symmetry_intervals(gx, gy, gz, sx, sy, sz)
    implicit none
    integer, intent(in) :: gx, gy, gz
    integer, intent(inout) :: sx, sy, sz
    integer :: default_x, default_y, default_z
    logical :: exists

    default_x = gx + 1
    default_y = gy + 1
    default_z = gz + 1
    inquire(file=p_scf_dir//'V_HYC_dp', exist=exists)
    if (exists) then
      open(unit=301, file=p_scf_dir//'V_HYC_dp', access='sequential',           &
        action='read', status='old', form='unformatted')
      read(301) default_x, default_y, default_z
      close(301)
    end if
    if (sx < 0) sx = default_x
    if (sy < 0) sy = default_y
    if (sz < 0) sz = default_z
  end subroutine set_default_symmetry_intervals

  subroutine read_r_vec(path, n, coords, gx, gy, gz, idx)
    implicit none
    character(len=*), intent(in) :: path
    integer, intent(out) :: n, gx, gy, gz
    integer, allocatable, intent(out) :: coords(:,:), idx(:,:,:)
    integer :: row

    open(unit=101, file=trim(path), action='read', status='old', &
      form='unformatted')
    read(101) n
    allocate(coords(n,3))
    read(101) coords
    close(101)

    gx = maxval(coords(:,1))
    gy = maxval(coords(:,2))
    gz = maxval(coords(:,3))
    if (minval(coords(:,1)) /= 0 .or. minval(coords(:,2)) /= 0 .or. &
        minval(coords(:,3)) /= 0) stop 'R_vec indices must start at zero'

    allocate(idx(0:gx,0:gy,0:gz))
    idx = 0
    do row = 1, n
      idx(coords(row,1), coords(row,2), coords(row,3)) = row
    enddo
  end subroutine read_r_vec

  subroutine read_v_func_tmp(path, n, b)
    implicit none
    character(len=*), intent(in) :: path
    integer, intent(in) :: n
    real(dp), intent(out) :: b(:,:,:,:)

    open(unit=102, file=trim(path), action='read', status='old', &
      form='unformatted')
    read(102)
    read(102)
    read(102)
    read(102) b
    close(102)
  end subroutine read_v_func_tmp

  subroutine report_channel_balance(unit_id, norm)
    implicit none
    integer, intent(in) :: unit_id
    real(dp), intent(in) :: norm(:)

    write(unit_id,'(A)') 'channel norm balance'
    write(unit_id,'(A,ES24.16)') 'linear dep1_vs_dep2_rel = ', &
      rel_pair(norm(1), norm(2))
    write(unit_id,'(A,ES24.16)') 'linear dep3_to_dep6_spread_over_avg = ', &
      spread_fraction(norm(3:6))
    write(unit_id,'(A,ES24.16)') 'linear dep3_vs_dep5_rel = ', &
      rel_pair(norm(3), norm(5))
    write(unit_id,'(A,ES24.16)') 'linear dep4_vs_dep6_rel = ', &
      rel_pair(norm(4), norm(6))
    if (size(norm) >= 12) then
      write(unit_id,'(A,ES24.16)') 'mod dep7_vs_dep8_rel = ', &
        rel_pair(norm(7), norm(8))
      write(unit_id,'(A,ES24.16)') 'mod dep9_to_dep12_spread_over_avg = ', &
        spread_fraction(norm(9:12))
      write(unit_id,'(A,ES24.16)') 'mod dep9_vs_dep11_rel = ', &
        rel_pair(norm(9), norm(11))
      write(unit_id,'(A,ES24.16)') 'mod dep10_vs_dep12_rel = ', &
        rel_pair(norm(10), norm(12))
    endif
    if (size(norm) == 22) then
      write(unit_id,'(A,ES24.16)') 'edge dep13_to_dep22_spread_over_avg = ', &
        spread_fraction(norm(13:22))
    endif
    write(unit_id,'(A)') ''
  end subroutine report_channel_balance

  subroutine build_channel_positions(gx, gy, gz, pos)
    implicit none
    integer, intent(in) :: gx, gy, gz
    real(dp), intent(out) :: pos(n_spinor,3)

    pos(1,:) = (/2d0/3d0, 1d0/3d0, 0.5d0/)
    pos(2,:) = (/1d0/3d0, 2d0/3d0, 0.5d0/)
    pos(3,:) = (/1d0/3d0, -1d0/3d0, 0.5d0/)
    pos(4,:) = (/4d0/3d0, 2d0/3d0, 0.5d0/)
    pos(5,:) = (/-1d0/3d0, 1d0/3d0, 0.5d0/)
    pos(6,:) = (/2d0/3d0, 4d0/3d0, 0.5d0/)

    pos(:,1) = pos(:,1) * real(gx,dp)
    pos(:,2) = pos(:,2) * real(gy,dp)
    pos(:,3) = pos(:,3) * real(gz,dp)
  end subroutine build_channel_positions

  subroutine report_symmetry_residual(unit_id, label, gx, gy, gz, idx, pos, &
      center, transform, kernel, linear_residual, mod_residual, edge_residual)
    implicit none
    integer, intent(in) :: unit_id, gx, gy, gz, idx(0:gx,0:gy,0:gz)
    character(len=*), intent(in) :: label
    real(dp), intent(in) :: pos(n_spinor,3), center(3), transform(3,3)
    real(dp), intent(in) :: kernel(:,:)
    real(dp), intent(out) :: linear_residual, mod_residual, edge_residual
    integer :: channel, mapped_channel, ix, iy, iz, i_src, i_dst
    integer :: ix_map, iy_map, iz_map, pt_map(n_spinor)
    integer :: mapped_points, mapped_channels, linear_pairs, mod_pairs
    integer :: edge_pairs, edge, mapped_edge
    integer :: edge_src(10), edge_dst(10), edge_map(10)
    real(dp) :: mapped_position(3), mapped_grid(3)
    real(dp) :: linear_diff2, linear_norm2, mod_diff2, mod_norm2
    real(dp) :: edge_diff2, edge_norm2
    real(dp) :: a, b, scale(size(kernel,2))

    pt_map = n_spinor + 1
    do channel = 1, n_spinor
      mapped_position = center + matmul(transpose(transform), &
        pos(channel,:) - center)
      do mapped_channel = 1, n_spinor
        if (sum((mapped_position - pos(mapped_channel,:))**2) < 1d-3) then
          pt_map(channel) = mapped_channel
        endif
      enddo
    enddo
    mapped_channels = count(pt_map <= n_spinor)
    edge_src = (/1,2,1,3,1,4,2,5,2,6/)
    edge_dst = (/2,1,3,1,4,1,5,2,6,2/)
    edge_map = 0
    do edge = 1, 10
      if (pt_map(edge_src(edge)) > n_spinor .or. &
          pt_map(edge_dst(edge)) > n_spinor) cycle
      do mapped_edge = 1, 10
        if (edge_src(mapped_edge) == pt_map(edge_src(edge)) .and. &
            edge_dst(mapped_edge) == pt_map(edge_dst(edge))) edge_map(edge) = mapped_edge
      enddo
    enddo
    do channel = 1, size(kernel,2)
      scale(channel) = maxval(abs(kernel(:,channel)))
    enddo

    linear_diff2 = 0d0
    linear_norm2 = 0d0
    mod_diff2 = 0d0
    mod_norm2 = 0d0
    edge_diff2 = 0d0
    edge_norm2 = 0d0
    mapped_points = 0
    linear_pairs = 0
    mod_pairs = 0
    edge_pairs = 0

    do ix = 0, gx
      do iy = 0, gy
        do iz = 0, gz
          i_src = idx(ix,iy,iz)
          if (i_src == 0) cycle
          mapped_grid = center + matmul(transpose(transform), &
            (/real(ix,dp), real(iy,dp), real(iz,dp)/) - center)
          ix_map = nint(mapped_grid(1))
          iy_map = nint(mapped_grid(2))
          iz_map = nint(mapped_grid(3))

          if (abs(mapped_grid(1) - real(ix_map,dp)) > 1d-8) cycle
          if (abs(mapped_grid(2) - real(iy_map,dp)) > 1d-8) cycle
          if (abs(mapped_grid(3) - real(iz_map,dp)) > 1d-8) cycle
          if (ix_map < 0 .or. ix_map > gx) cycle
          if (iy_map < 0 .or. iy_map > gy) cycle
          if (iz_map < 0 .or. iz_map > gz) cycle

          i_dst = idx(ix_map,iy_map,iz_map)
          if (i_dst == 0) cycle
          mapped_points = mapped_points + 1

          do channel = 1, n_spinor
            mapped_channel = pt_map(channel)
            if (mapped_channel > n_spinor) cycle

            a = kernel(i_src,channel)
            b = kernel(i_dst,mapped_channel)
            if (abs(a) + abs(b) > 1d-30 * &
                max(max(scale(channel), scale(mapped_channel)), 1d-300)) then
              linear_diff2 = linear_diff2 + (b - a)**2
              linear_norm2 = linear_norm2 + b*b + a*a
              linear_pairs = linear_pairs + 1
            endif

            if (size(kernel,2) >= 12) then
              a = kernel(i_src,n_spinor+channel)
              b = kernel(i_dst,n_spinor+mapped_channel)
              if (abs(a) + abs(b) > 1d-30 * &
                  max(max(scale(n_spinor+channel), &
                  scale(n_spinor+mapped_channel)), 1d-300)) then
                mod_diff2 = mod_diff2 + (b - a)**2
                mod_norm2 = mod_norm2 + b*b + a*a
                mod_pairs = mod_pairs + 1
              endif
            endif
          enddo
          if (size(kernel,2) == 22) then
            do edge = 1, 10
              mapped_edge = edge_map(edge)
              if (mapped_edge == 0) cycle
              a = kernel(i_src,12+edge)
              b = kernel(i_dst,12+mapped_edge)
              if (abs(a) + abs(b) > 1d-30 * &
                  max(max(scale(12+edge),scale(12+mapped_edge)),1d-300)) then
                edge_diff2 = edge_diff2 + (b-a)**2
                edge_norm2 = edge_norm2 + b*b + a*a
                edge_pairs = edge_pairs + 1
              endif
            enddo
          endif
        enddo
      enddo
    enddo

    write(unit_id,'(A,A)') 'symmetry residual: ', trim(label)
    write(unit_id,'(A,6I4)') 'source -> target channels = ', pt_map(:)
    write(unit_id,'(A,I0)') 'mapped grid points = ', mapped_points
    write(unit_id,'(A,I0)') 'mapped channels = ', mapped_channels
    write(unit_id,'(A,I0)') 'active linear pairs = ', linear_pairs
    write(unit_id,'(A,I0)') 'active mod pairs = ', mod_pairs
    write(unit_id,'(A,I0)') 'active directed-edge pairs = ', edge_pairs
    linear_residual = sym_rel(linear_diff2, linear_norm2)
    mod_residual = sym_rel(mod_diff2, mod_norm2)
    edge_residual = sym_rel(edge_diff2, edge_norm2)
    write(unit_id,'(A,ES24.16)') 'K_linear scalar-kernel residual = ',        &
      linear_residual
    write(unit_id,'(A,ES24.16)') 'K_mod scalar-kernel residual = ',           &
      mod_residual
    if (size(kernel,2) == 22) then
      write(unit_id,'(A,ES24.16)') 'K_edge scalar-kernel residual = ', edge_residual
    endif
    write(unit_id,'(A)') ''
  end subroutine report_symmetry_residual

  real(dp) function rel_pair(a, b)
    implicit none
    real(dp), intent(in) :: a, b

    rel_pair = abs(a - b) / max(max(abs(a), abs(b)), 1d-300)
  end function rel_pair

  real(dp) function spread_fraction(values)
    implicit none
    real(dp), intent(in) :: values(:)
    real(dp) :: avg

    avg = sum(values) / real(size(values), dp)
    spread_fraction = (maxval(values) - minval(values)) / max(abs(avg), 1d-300)
  end function spread_fraction

  real(dp) function sym_rel(diff2, norm2)
    implicit none
    real(dp), intent(in) :: diff2, norm2

    sym_rel = sqrt(2d0 * diff2 / max(norm2, 1d-300))
  end function sym_rel

end program check_v_func_tmp_c3
