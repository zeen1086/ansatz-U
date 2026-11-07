module cri3_c3_basis_mod
  use parameter_mod
  implicit none
  private

  integer, parameter :: dp = kind(0d0)
  integer, parameter :: c3_n_channel = p_fit_nspinor
  character(len=*), parameter, public :: c3_basis_mode =                     &
    'finite_window_C3v_A_B_mirror_reduced_basis_v1'

  type, public :: c3_reduced_basis
    integer :: n_grid = 0
    integer :: n_channel = 0
    integer :: n_active = 0
    integer :: n_reduced = 0
    integer :: n_relations = 0
    integer, allocatable :: orbit_of(:)
    integer, allocatable :: offset(:)
    integer, allocatable :: member(:)
  end type c3_reduced_basis

  public :: build_c3_reduced_basis
  public :: build_c3_window_basis
  public :: close_c3_orbit_mask
  public :: project_c3_orbit_values
  public :: clear_c3_reduced_basis
  public :: solve_c3_reduced_cd
  public :: c3_basis_residual

contains

  subroutine build_c3_reduced_basis(nx, ny, nz, ni, active_mask, basis)
    implicit none
    integer, intent(in) :: nx, ny, nz, ni
    real(dp), intent(in) :: active_mask(:,:)
    type(c3_reduced_basis), intent(inout) :: basis
    integer :: ncoef, idx, channel, i, root, q
    integer :: relation_count, mismatch_count
    integer, allocatable :: parent(:), root_id(:), count_per_orbit(:), cursor(:)

    if (size(active_mask,1) /= ni) stop 'C3v basis grid-size mismatch'
    if (size(active_mask,2) /= c3_n_channel) then
      stop 'C3v basis requires six ansatz channels'
    end if
    if (ni /= (nx + 1) * (ny + 1) * (nz + 1)) then
      stop 'C3v basis inconsistent finite-window dimensions'
    end if
    if (nx /= ny) stop 'C3v xy mirror requires nx = ny'

    call clear_basis(basis)
    ncoef = ni * c3_n_channel
    allocate(parent(ncoef))
    parent = 0
    do channel = 1, c3_n_channel
      do i = 1, ni
        idx = coefficient_index(i, channel, ni)
        if (active_mask(i,channel) > 0.5d0) parent(idx) = idx
      end do
    end do

    relation_count = 0
    mismatch_count = 0
    call add_c3v_operation(nx, ny, nz, ni, active_mask, 1, parent,            &
      relation_count, mismatch_count)
    call add_c3v_operation(nx, ny, nz, ni, active_mask, 2, parent,            &
      relation_count, mismatch_count)
    call add_c3v_operation(nx, ny, nz, ni, active_mask, 3, parent,            &
      relation_count, mismatch_count)
    if (mismatch_count /= 0) then
      write(*,'(A,I0)')                                                        &
        'ERROR: block mask is not closed under C3v; mismatches = ',            &
        mismatch_count
      stop 1
    end if

    do idx = 1, ncoef
      if (parent(idx) > 0) parent(idx) = find_root(parent, idx)
    end do

    allocate(root_id(ncoef))
    root_id = 0
    basis%n_reduced = 0
    do idx = 1, ncoef
      if (parent(idx) <= 0) cycle
      root = parent(idx)
      if (root_id(root) == 0) then
        basis%n_reduced = basis%n_reduced + 1
        root_id(root) = basis%n_reduced
      end if
      parent(idx) = root_id(root)
    end do
    deallocate(root_id)

    basis%n_grid = ni
    basis%n_channel = c3_n_channel
    basis%n_active = count(parent > 0)
    basis%n_relations = relation_count
    allocate(count_per_orbit(basis%n_reduced))
    count_per_orbit = 0
    do idx = 1, ncoef
      q = parent(idx)
      if (q > 0) count_per_orbit(q) = count_per_orbit(q) + 1
    end do

    allocate(basis%offset(basis%n_reduced + 1))
    basis%offset(1) = 1
    do q = 1, basis%n_reduced
      basis%offset(q+1) = basis%offset(q) + count_per_orbit(q)
    end do
    allocate(basis%member(basis%n_active))
    allocate(cursor(basis%n_reduced))
    cursor = basis%offset(1:basis%n_reduced)
    do idx = 1, ncoef
      q = parent(idx)
      if (q <= 0) cycle
      basis%member(cursor(q)) = idx
      cursor(q) = cursor(q) + 1
    end do
    deallocate(cursor, count_per_orbit)
    call move_alloc(parent, basis%orbit_of)

    write(*,'(A,A)') 'C3v reduced basis mode = ', c3_basis_mode
    write(*,'(A,I0)') 'C3v active coefficients = ', basis%n_active
    write(*,'(A,I0)') 'C3v irreducible coefficients = ', basis%n_reduced
    write(*,'(A,I0)') 'C3v active mapping relations = ', basis%n_relations
    write(*,'(A,F12.6)') 'C3v reduction ratio active/irreducible = ',         &
      real(basis%n_active,dp) / real(max(basis%n_reduced,1),dp)
  end subroutine build_c3_reduced_basis

  subroutine build_c3_window_basis(grid_x, grid_y, grid_z, sym_x, sym_y,     &
      sym_z, ni, basis)
    implicit none
    integer, intent(in) :: grid_x, grid_y, grid_z
    integer, intent(in) :: sym_x, sym_y, sym_z, ni
    type(c3_reduced_basis), intent(inout) :: basis
    integer :: ncoef, idx, channel, i, root, q
    integer :: relation_count
    integer, allocatable :: parent(:), root_id(:), count_per_orbit(:)
    integer, allocatable :: cursor(:)

    if (ni /= (grid_x + 1) * (grid_y + 1) * (grid_z + 1)) then
      stop 'C3v output basis inconsistent grid dimensions'
    end if
    if (sym_x /= sym_y) stop 'C3v output basis requires sym_x = sym_y'
    if (grid_x < 0 .or. grid_y < 0 .or. grid_z < 0) then
      stop 'C3v output basis requires nonnegative grid maxima'
    end if

    call clear_basis(basis)
    ncoef = ni * c3_n_channel
    allocate(parent(ncoef))
    do channel = 1, c3_n_channel
      do i = 1, ni
        idx = coefficient_index(i, channel, ni)
        parent(idx) = idx
      end do
    end do

    relation_count = 0
    call add_c3v_window_operation(grid_x, grid_y, grid_z, sym_x, sym_y,       &
      sym_z, ni, 1, parent, relation_count)
    call add_c3v_window_operation(grid_x, grid_y, grid_z, sym_x, sym_y,       &
      sym_z, ni, 2, parent, relation_count)
    call add_c3v_window_operation(grid_x, grid_y, grid_z, sym_x, sym_y,       &
      sym_z, ni, 3, parent, relation_count)

    do idx = 1, ncoef
      parent(idx) = find_root(parent, idx)
    end do

    allocate(root_id(ncoef))
    root_id = 0
    basis%n_reduced = 0
    do idx = 1, ncoef
      root = parent(idx)
      if (root_id(root) == 0) then
        basis%n_reduced = basis%n_reduced + 1
        root_id(root) = basis%n_reduced
      end if
      parent(idx) = root_id(root)
    end do
    deallocate(root_id)

    basis%n_grid = ni
    basis%n_channel = c3_n_channel
    basis%n_active = ncoef
    basis%n_relations = relation_count
    allocate(count_per_orbit(basis%n_reduced))
    count_per_orbit = 0
    do idx = 1, ncoef
      q = parent(idx)
      count_per_orbit(q) = count_per_orbit(q) + 1
    end do

    allocate(basis%offset(basis%n_reduced + 1))
    basis%offset(1) = 1
    do q = 1, basis%n_reduced
      basis%offset(q+1) = basis%offset(q) + count_per_orbit(q)
    end do
    allocate(basis%member(ncoef))
    allocate(cursor(basis%n_reduced))
    cursor = basis%offset(1:basis%n_reduced)
    do idx = 1, ncoef
      q = parent(idx)
      basis%member(cursor(q)) = idx
      cursor(q) = cursor(q) + 1
    end do
    deallocate(cursor, count_per_orbit)
    call move_alloc(parent, basis%orbit_of)

    write(*,'(A)') 'C3v final-output orbit basis'
    write(*,'(A,3I8)') 'C3v output grid maxima = ', grid_x, grid_y, grid_z
    write(*,'(A,3I8)') 'C3v physical symmetry intervals = ',                &
      sym_x, sym_y, sym_z
    write(*,'(A,I0)') 'C3v output coefficients = ', basis%n_active
    write(*,'(A,I0)') 'C3v output irreducible orbits = ', basis%n_reduced
    write(*,'(A,I0)') 'C3v output mapping relations = ', basis%n_relations
  end subroutine build_c3_window_basis

  subroutine close_c3_orbit_mask(basis, active_mask, removed_count,          &
      active_count)
    implicit none
    type(c3_reduced_basis), intent(in) :: basis
    real(dp), intent(inout) :: active_mask(:,:)
    integer, intent(out) :: removed_count, active_count
    integer :: q, pos, idx, i, channel
    logical :: orbit_active

    if (size(active_mask,1) /= basis%n_grid .or.                             &
        size(active_mask,2) /= basis%n_channel) then
      stop 'C3v output mask shape mismatch'
    end if

    removed_count = 0
    active_count = 0
    do q = 1, basis%n_reduced
      orbit_active = .true.
      do pos = basis%offset(q), basis%offset(q+1) - 1
        idx = basis%member(pos)
        call decode_coefficient_index(idx, basis%n_grid, i, channel)
        if (active_mask(i,channel) <= 0.5d0) orbit_active = .false.
      end do

      do pos = basis%offset(q), basis%offset(q+1) - 1
        idx = basis%member(pos)
        call decode_coefficient_index(idx, basis%n_grid, i, channel)
        if (.not. orbit_active .and. active_mask(i,channel) > 0.5d0) then
          removed_count = removed_count + 1
        end if
        if (orbit_active) then
          active_mask(i,channel) = 1d0
          active_count = active_count + 1
        else
          active_mask(i,channel) = 0d0
        end if
      end do
    end do
  end subroutine close_c3_orbit_mask

  subroutine project_c3_orbit_values(basis, active_mask, values, rel_change, &
      max_change)
    implicit none
    type(c3_reduced_basis), intent(in) :: basis
    real(dp), intent(in) :: active_mask(:,:)
    real(dp), intent(inout) :: values(:,:)
    real(dp), intent(out) :: rel_change, max_change
    integer :: q, pos, idx, i, channel, nmember
    real(dp) :: average, delta, diff2, norm2
    logical :: orbit_active

    if (size(values,1) /= basis%n_grid .or.                                  &
        size(values,2) /= basis%n_channel) then
      stop 'C3v output value shape mismatch'
    end if
    if (size(active_mask,1) /= basis%n_grid .or.                             &
        size(active_mask,2) /= basis%n_channel) then
      stop 'C3v output projection mask shape mismatch'
    end if

    diff2 = 0d0
    norm2 = 0d0
    max_change = 0d0
    do q = 1, basis%n_reduced
      orbit_active = .true.
      average = 0d0
      nmember = basis%offset(q+1) - basis%offset(q)
      do pos = basis%offset(q), basis%offset(q+1) - 1
        idx = basis%member(pos)
        call decode_coefficient_index(idx, basis%n_grid, i, channel)
        norm2 = norm2 + values(i,channel)**2
        average = average + values(i,channel)
        if (active_mask(i,channel) <= 0.5d0) orbit_active = .false.
      end do

      if (orbit_active) then
        average = average / real(nmember,dp)
      else
        average = 0d0
      end if

      do pos = basis%offset(q), basis%offset(q+1) - 1
        idx = basis%member(pos)
        call decode_coefficient_index(idx, basis%n_grid, i, channel)
        delta = average - values(i,channel)
        diff2 = diff2 + delta**2
        max_change = max(max_change, abs(delta))
        values(i,channel) = average
      end do
    end do
    rel_change = sqrt(diff2 / max(norm2, 1d-300))
  end subroutine project_c3_orbit_values

  subroutine solve_c3_reduced_cd(basis, gram, rhs, penalty, l1, l2, n_iter,    &
      tolerance, scale, grid_weight, coef, iter_used, converged, anchor,       &
      prox_mu)
    implicit none
    type(c3_reduced_basis), intent(in) :: basis
    real(dp), intent(in) :: gram(:,:), rhs(:,:), penalty(:)
    real(dp), intent(in) :: l1, l2, tolerance, scale(:), grid_weight(:)
    integer, intent(in) :: n_iter
    real(dp), intent(out) :: coef(:,:)
    integer, intent(out) :: iter_used
    logical, intent(out) :: converged
    real(dp), intent(in), optional :: anchor(:,:), prox_mu
    integer :: q, pos, idx, i, channel, other, other_orbit, iter
    real(dp) :: rho, denom, old_value, new_value, max_change
    real(dp) :: local_diag, mu
    real(dp), allocatable :: value(:), rhs_reduced(:), data_diag(:)
    real(dp), allocatable :: ridge_diag(:), l1_reduced(:)
    real(dp), allocatable :: prox_diag(:), prox_rhs(:)

    call validate_solver_shapes(basis, gram, rhs, penalty, scale, grid_weight, &
      coef)
    if (present(anchor) .neqv. present(prox_mu)) then
      stop 'C3v solver requires both anchor and prox_mu'
    end if
    mu = 0d0
    if (present(prox_mu)) mu = prox_mu
    if (mu < 0d0) stop 'C3v solver prox_mu must be nonnegative'
    if (present(anchor)) then
      if (size(anchor,1) /= basis%n_grid .or.                                &
          size(anchor,2) /= basis%n_channel) then
        stop 'C3v solver anchor shape mismatch'
      end if
    end if

    allocate(value(basis%n_reduced))
    allocate(rhs_reduced(basis%n_reduced))
    allocate(data_diag(basis%n_reduced))
    allocate(ridge_diag(basis%n_reduced))
    allocate(l1_reduced(basis%n_reduced))
    allocate(prox_diag(basis%n_reduced))
    allocate(prox_rhs(basis%n_reduced))
    value = 0d0
    rhs_reduced = 0d0
    data_diag = 0d0
    ridge_diag = 0d0
    l1_reduced = 0d0
    prox_diag = 0d0
    prox_rhs = 0d0

    do q = 1, basis%n_reduced
      do pos = basis%offset(q), basis%offset(q+1) - 1
        idx = basis%member(pos)
        call decode_coefficient_index(idx, basis%n_grid, i, channel)
        rhs_reduced(q) = rhs_reduced(q) + grid_weight(i) * scale(channel)     &
          * rhs(channel,i)
        ridge_diag(q) = ridge_diag(q) + l2 * penalty(channel)                &
          * scale(channel)**2
        l1_reduced(q) = l1_reduced(q) + l1 * penalty(channel)                &
          * abs(scale(channel))

        do other = 1, basis%n_channel
          other_orbit = basis%orbit_of(coefficient_index(i, other,            &
            basis%n_grid))
          if (other_orbit == q) then
            data_diag(q) = data_diag(q) + grid_weight(i) * scale(channel)    &
              * gram(channel,other) * scale(other)
          end if
        end do

        if (present(anchor)) then
          local_diag = (grid_weight(i) * gram(channel,channel)               &
            + l2 * penalty(channel)) * scale(channel)**2
          prox_diag(q) = prox_diag(q) + mu * local_diag
          prox_rhs(q) = prox_rhs(q) + mu * local_diag * anchor(i,channel)
          value(q) = value(q) + anchor(i,channel)
        end if
      end do
      if (present(anchor)) then
        value(q) = value(q) / real(basis%offset(q+1)-basis%offset(q),dp)
      end if
      if (data_diag(q) + ridge_diag(q) + prox_diag(q) <= 0d0) then
        stop 'nonpositive C3v reduced normal-equation diagonal'
      end if
    end do

    converged = .false.
    do iter = 1, n_iter
      max_change = 0d0
      do q = 1, basis%n_reduced
        rho = rhs_reduced(q) + prox_rhs(q)
        do pos = basis%offset(q), basis%offset(q+1) - 1
          idx = basis%member(pos)
          call decode_coefficient_index(idx, basis%n_grid, i, channel)
          do other = 1, basis%n_channel
            other_orbit = basis%orbit_of(coefficient_index(i, other,          &
              basis%n_grid))
            if (other_orbit <= 0 .or. other_orbit == q) cycle
            rho = rho - grid_weight(i) * scale(channel)                      &
              * gram(channel,other) * scale(other) * value(other_orbit)
          end do
        end do
        denom = data_diag(q) + ridge_diag(q) + prox_diag(q)
        old_value = value(q)
        new_value = soft_threshold(rho, l1_reduced(q)) / denom
        value(q) = new_value
        max_change = max(max_change, abs(new_value - old_value))
      end do
      if (max_change < tolerance) then
        converged = .true.
        exit
      end if
    end do
    iter_used = min(iter, n_iter)

    coef = 0d0
    do idx = 1, size(basis%orbit_of)
      q = basis%orbit_of(idx)
      if (q <= 0) cycle
      call decode_coefficient_index(idx, basis%n_grid, i, channel)
      coef(i,channel) = value(q)
    end do

    deallocate(value, rhs_reduced, data_diag, ridge_diag, l1_reduced)
    deallocate(prox_diag, prox_rhs)
  end subroutine solve_c3_reduced_cd

  subroutine c3_basis_residual(basis, coef, max_abs, rms)
    implicit none
    type(c3_reduced_basis), intent(in) :: basis
    real(dp), intent(in) :: coef(:,:)
    real(dp), intent(out) :: max_abs, rms
    integer :: q, pos, idx, i, channel, nterm
    real(dp) :: average, delta, sumsq

    if (size(coef,1) /= basis%n_grid .or.                                     &
        size(coef,2) /= basis%n_channel) then
      stop 'C3v residual coefficient shape mismatch'
    end if
    max_abs = 0d0
    sumsq = 0d0
    nterm = 0
    do q = 1, basis%n_reduced
      average = 0d0
      do pos = basis%offset(q), basis%offset(q+1) - 1
        idx = basis%member(pos)
        call decode_coefficient_index(idx, basis%n_grid, i, channel)
        average = average + coef(i,channel)
      end do
      average = average / real(basis%offset(q+1)-basis%offset(q),dp)
      do pos = basis%offset(q), basis%offset(q+1) - 1
        idx = basis%member(pos)
        call decode_coefficient_index(idx, basis%n_grid, i, channel)
        delta = coef(i,channel) - average
        max_abs = max(max_abs, abs(delta))
        sumsq = sumsq + delta**2
        nterm = nterm + 1
      end do
    end do
    rms = sqrt(sumsq / real(max(nterm,1),dp))
  end subroutine c3_basis_residual

  subroutine add_c3v_operation(nx, ny, nz, ni, active_mask, operation_kind,   &
      parent,                                                                  &
      relation_count, mismatch_count)
    implicit none
    integer, intent(in) :: nx, ny, nz, ni, operation_kind
    real(dp), intent(in) :: active_mask(ni,c3_n_channel)
    integer, intent(inout) :: parent(:), relation_count, mismatch_count
    integer :: ix, iy, iz, i, j, channel, mapped_channel
    integer :: ix_map, iy_map, iz_map
    integer :: channel_map(c3_n_channel)
    real(dp) :: center(3), transform(3,3), mapped_grid(3)
    real(dp) :: positions(c3_n_channel,3)
    logical :: source_active, target_active

    call build_channel_positions(nx, ny, nz, positions)
    if (operation_kind == 1) then
      center = (/real(nx,dp)/3d0, 2d0*real(ny,dp)/3d0, real(nz,dp)/2d0/)
      call c3_plus_transform(transform)
    else if (operation_kind == 2) then
      center = (/2d0*real(nx,dp)/3d0, real(ny,dp)/3d0, real(nz,dp)/2d0/)
      call c3_plus_transform(transform)
    else if (operation_kind == 3) then
      center = 0d0
      call mirror_transform(transform)
    else
      stop 'invalid C3v operation kind'
    end if
    call build_channel_map(center, transform, positions, channel_map)

    do iz = 0, nz
      do iy = 0, ny
        do ix = 0, nx
          i = index_from_grid(ix, iy, iz, nx, ny)
          mapped_grid = center + matmul(transpose(transform),                   &
            (/real(ix,dp), real(iy,dp), real(iz,dp)/) - center)
          ix_map = nint(mapped_grid(1))
          iy_map = nint(mapped_grid(2))
          iz_map = nint(mapped_grid(3))
          if (abs(mapped_grid(1)-real(ix_map,dp)) > 1d-8) cycle
          if (abs(mapped_grid(2)-real(iy_map,dp)) > 1d-8) cycle
          if (abs(mapped_grid(3)-real(iz_map,dp)) > 1d-8) cycle
          if (ix_map < 0 .or. ix_map > nx) cycle
          if (iy_map < 0 .or. iy_map > ny) cycle
          if (iz_map < 0 .or. iz_map > nz) cycle
          j = index_from_grid(ix_map, iy_map, iz_map, nx, ny)

          do channel = 1, c3_n_channel
            mapped_channel = channel_map(channel)
            if (mapped_channel > c3_n_channel) cycle
            source_active = active_mask(i,channel) > 0.5d0
            target_active = active_mask(j,mapped_channel) > 0.5d0
            if (source_active .neqv. target_active) then
              mismatch_count = mismatch_count + 1
              cycle
            end if
            if (.not. source_active) cycle
            call union_sets(parent, coefficient_index(i,channel,ni),           &
              coefficient_index(j,mapped_channel,ni))
            relation_count = relation_count + 1
          end do
        end do
      end do
    end do
  end subroutine add_c3v_operation

  subroutine add_c3v_window_operation(grid_x, grid_y, grid_z, sym_x, sym_y, &
      sym_z, ni, operation_kind, parent, relation_count)
    implicit none
    integer, intent(in) :: grid_x, grid_y, grid_z
    integer, intent(in) :: sym_x, sym_y, sym_z, ni, operation_kind
    integer, intent(inout) :: parent(:), relation_count
    integer :: ix, iy, iz, i, j, channel, mapped_channel
    integer :: ix_map, iy_map, iz_map
    integer :: channel_map(c3_n_channel)
    real(dp) :: center(3), transform(3,3), mapped_grid(3)
    real(dp) :: positions(c3_n_channel,3)

    call build_channel_positions(sym_x, sym_y, sym_z, positions)
    if (operation_kind == 1) then
      center = (/real(sym_x,dp)/3d0, 2d0*real(sym_y,dp)/3d0,               &
        real(sym_z,dp)/2d0/)
      call c3_plus_transform(transform)
    else if (operation_kind == 2) then
      center = (/2d0*real(sym_x,dp)/3d0, real(sym_y,dp)/3d0,               &
        real(sym_z,dp)/2d0/)
      call c3_plus_transform(transform)
    else if (operation_kind == 3) then
      center = 0d0
      call mirror_transform(transform)
    else
      stop 'invalid C3v output operation kind'
    end if
    call build_channel_map(center, transform, positions, channel_map)

    do iz = 0, grid_z
      do iy = 0, grid_y
        do ix = 0, grid_x
          i = index_from_grid(ix, iy, iz, grid_x, grid_y)
          mapped_grid = center + matmul(transpose(transform),                 &
            (/real(ix,dp), real(iy,dp), real(iz,dp)/) - center)
          ix_map = nint(mapped_grid(1))
          iy_map = nint(mapped_grid(2))
          iz_map = nint(mapped_grid(3))
          if (abs(mapped_grid(1)-real(ix_map,dp)) > 1d-8) cycle
          if (abs(mapped_grid(2)-real(iy_map,dp)) > 1d-8) cycle
          if (abs(mapped_grid(3)-real(iz_map,dp)) > 1d-8) cycle
          if (ix_map < 0 .or. ix_map > grid_x) cycle
          if (iy_map < 0 .or. iy_map > grid_y) cycle
          if (iz_map < 0 .or. iz_map > grid_z) cycle
          j = index_from_grid(ix_map, iy_map, iz_map, grid_x, grid_y)

          do channel = 1, c3_n_channel
            mapped_channel = channel_map(channel)
            if (mapped_channel > c3_n_channel) cycle
            call union_sets(parent, coefficient_index(i,channel,ni),         &
              coefficient_index(j,mapped_channel,ni))
            relation_count = relation_count + 1
          end do
        end do
      end do
    end do
  end subroutine add_c3v_window_operation

  subroutine build_channel_positions(nx, ny, nz, positions)
    implicit none
    integer, intent(in) :: nx, ny, nz
    real(dp), intent(out) :: positions(c3_n_channel,3)

    positions(1,:) = (/2d0/3d0, 1d0/3d0, 0.5d0/)
    positions(2,:) = (/1d0/3d0, 2d0/3d0, 0.5d0/)
    positions(3,:) = (/1d0/3d0, -1d0/3d0, 0.5d0/)
    positions(4,:) = (/4d0/3d0, 2d0/3d0, 0.5d0/)
    positions(5,:) = (/-1d0/3d0, 1d0/3d0, 0.5d0/)
    positions(6,:) = (/2d0/3d0, 4d0/3d0, 0.5d0/)
    positions(:,1) = positions(:,1) * real(nx,dp)
    positions(:,2) = positions(:,2) * real(ny,dp)
    positions(:,3) = positions(:,3) * real(nz,dp)
  end subroutine build_channel_positions

  subroutine build_channel_map(center, transform, positions, channel_map)
    implicit none
    real(dp), intent(in) :: center(3), transform(3,3)
    real(dp), intent(in) :: positions(c3_n_channel,3)
    integer, intent(out) :: channel_map(c3_n_channel)
    integer :: channel, mapped_channel
    real(dp) :: mapped_position(3)

    channel_map = c3_n_channel + 1
    do channel = 1, c3_n_channel
      mapped_position = center + matmul(transpose(transform),                  &
        positions(channel,:) - center)
      do mapped_channel = 1, c3_n_channel
        if (sum((mapped_position-positions(mapped_channel,:))**2) < 1d-3) then
          channel_map(channel) = mapped_channel
        end if
      end do
    end do
  end subroutine build_channel_map

  subroutine c3_plus_transform(transform)
    implicit none
    real(dp), intent(out) :: transform(3,3)

    transform(1,:) = (/0d0, 1d0, 0d0/)
    transform(2,:) = (/-1d0, -1d0, 0d0/)
    transform(3,:) = (/0d0, 0d0, 1d0/)
  end subroutine c3_plus_transform

  subroutine mirror_transform(transform)
    implicit none
    real(dp), intent(out) :: transform(3,3)

    transform(1,:) = (/0d0, 1d0, 0d0/)
    transform(2,:) = (/1d0, 0d0, 0d0/)
    transform(3,:) = (/0d0, 0d0, 1d0/)
  end subroutine mirror_transform

  subroutine union_sets(parent, a, b)
    implicit none
    integer, intent(inout) :: parent(:)
    integer, intent(in) :: a, b
    integer :: root_a, root_b

    if (parent(a) <= 0 .or. parent(b) <= 0) then
      stop 'union on blocked C3v coefficient'
    end if
    root_a = find_root(parent, a)
    root_b = find_root(parent, b)
    if (root_a == root_b) return
    if (root_a < root_b) then
      parent(root_b) = root_a
    else
      parent(root_a) = root_b
    end if
  end subroutine union_sets

  integer function find_root(parent, item)
    implicit none
    integer, intent(inout) :: parent(:)
    integer, intent(in) :: item
    integer :: current, next_item, root

    if (parent(item) <= 0) stop 'find_root on inactive coefficient'
    root = item
    do while (parent(root) /= root)
      root = parent(root)
    end do
    current = item
    do while (parent(current) /= current)
      next_item = parent(current)
      parent(current) = root
      current = next_item
    end do
    find_root = root
  end function find_root

  subroutine validate_solver_shapes(basis, gram, rhs, penalty, scale, weight,  &
      coef)
    implicit none
    type(c3_reduced_basis), intent(in) :: basis
    real(dp), intent(in) :: gram(:,:), rhs(:,:), penalty(:), scale(:), weight(:)
    real(dp), intent(in) :: coef(:,:)

    if (size(gram,1) /= basis%n_channel .or.                                 &
        size(gram,2) /= basis%n_channel) stop 'C3v solver Gram shape mismatch'
    if (size(rhs,1) /= basis%n_channel .or.                                  &
        size(rhs,2) /= basis%n_grid) stop 'C3v solver RHS shape mismatch'
    if (size(penalty) /= basis%n_channel) stop 'C3v solver penalty mismatch'
    if (size(scale) /= basis%n_channel) stop 'C3v solver scale mismatch'
    if (size(weight) /= basis%n_grid) stop 'C3v solver weight mismatch'
    if (size(coef,1) /= basis%n_grid .or.                                    &
        size(coef,2) /= basis%n_channel) stop 'C3v solver output shape mismatch'
  end subroutine validate_solver_shapes

  subroutine clear_basis(basis)
    implicit none
    type(c3_reduced_basis), intent(inout) :: basis

    if (allocated(basis%orbit_of)) deallocate(basis%orbit_of)
    if (allocated(basis%offset)) deallocate(basis%offset)
    if (allocated(basis%member)) deallocate(basis%member)
    basis%n_grid = 0
    basis%n_channel = 0
    basis%n_active = 0
    basis%n_reduced = 0
    basis%n_relations = 0
  end subroutine clear_basis

  subroutine clear_c3_reduced_basis(basis)
    implicit none
    type(c3_reduced_basis), intent(inout) :: basis

    call clear_basis(basis)
  end subroutine clear_c3_reduced_basis

  integer function coefficient_index(i, channel, ni)
    implicit none
    integer, intent(in) :: i, channel, ni

    coefficient_index = i + ni * (channel - 1)
  end function coefficient_index

  subroutine decode_coefficient_index(idx, ni, i, channel)
    implicit none
    integer, intent(in) :: idx, ni
    integer, intent(out) :: i, channel

    channel = (idx - 1) / ni + 1
    i = mod(idx - 1, ni) + 1
  end subroutine decode_coefficient_index

  integer function index_from_grid(ix, iy, iz, nx, ny)
    implicit none
    integer, intent(in) :: ix, iy, iz, nx, ny

    index_from_grid = 1 + ix + (nx + 1) * iy + (nx + 1) * (ny + 1) * iz
  end function index_from_grid

  real(dp) function soft_threshold(value, threshold)
    implicit none
    real(dp), intent(in) :: value, threshold

    if (value > threshold) then
      soft_threshold = value - threshold
    else if (value < -threshold) then
      soft_threshold = value + threshold
    else
      soft_threshold = 0d0
    end if
  end function soft_threshold

end module cri3_c3_basis_mod
