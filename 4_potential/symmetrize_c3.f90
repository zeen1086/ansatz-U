program symmetrize_c3
  use parameter_mod
        implicit none

        integer, parameter :: dp = kind(0d0)
        integer, parameter :: n_spinor = p_fit_nspinor
        integer, parameter :: n_ops = 3
        real(dp), parameter :: pi = 3.1415926535897932384626433832795_dp
        real(dp), parameter :: map_tol = p_map_tolerance

        ! Change this if the C3 axis used by the fitted local channels is not
        ! the primitive-cell origin. The channel-map check below will stop if
        ! this center is inconsistent with the six fitted channels.
        real(dp), parameter :: c3_center(3) = (/ 0.0_dp, 0.0_dp, 0.5_dp /)

        integer :: nx, ny, nz, ni
        character(len=256) :: input_file, output_file
        real(dp), allocatable :: VA_ave(:,:), V0(:,:), Vn(:,:), JH(:,:)
        real(dp), allocatable :: B_old(:,:,:,:), B_c3(:,:,:,:)
        real(dp) :: diff_norm, old_norm, residual

        input_file = 'V_func_tmp_tmp'
        output_file = 'V_func_tmp_tmp_c3'
        if (command_argument_count() >= 1) call get_command_argument(1, input_file)
        if (command_argument_count() >= 2) call get_command_argument(2, output_file)

        call read_grid_shape(nx, ny, nz)
        if (nx /= ny) then
                print *, 'ERROR: C3 projection assumes nx == ny on the hexagonal grid.'
                print *, 'nx, ny = ', nx, ny
                stop 1
        endif
        ni = (nx + 1) * (ny + 1) * (nz + 1)

        allocate(VA_ave(ni,4))
        allocate(V0(ni,3))
        allocate(Vn(ni,8))
        allocate(JH(ni,n_spinor))
        allocate(B_old(ni,3,n_spinor,3))
        allocate(B_c3(ni,3,n_spinor,3))

        call read_vfunc(trim(input_file), VA_ave, V0, Vn, B_old, JH)
        call project_b_general(nx, ny, nz, B_old, B_c3)

        old_norm = sqrt(sum(B_old * B_old))
        diff_norm = sqrt(sum((B_c3 - B_old) * (B_c3 - B_old)))
        if (old_norm > 0.0_dp) then
                print '(A,ES14.6)', 'relative change in B_general = ', diff_norm / old_norm
        else
                print '(A)', 'relative change in B_general = NaN (zero input norm)'
        endif

        call c3_residual(nx, ny, nz, B_c3, residual)
        print '(A,ES14.6)', 'max relative C3 residual after projection = ', residual

        call write_vfunc(trim(output_file), VA_ave, V0, Vn, B_c3, JH)
        print *, 'wrote ', trim(output_file)

contains

subroutine read_grid_shape(nx, ny, nz)
        integer, intent(out) :: nx, ny, nz

        open(unit=101, file=p_scf_dir//'V_QE', action='read', access='sequential', &
                form='unformatted', status='old')
        read(101) nx, ny, nz
        close(101)
        print *, 'nx ny nz = ', nx, ny, nz
end subroutine read_grid_shape

subroutine read_vfunc(path, VA_ave, V0, Vn, B_general, JH)
        character(len=*), intent(in) :: path
        real(dp), intent(out) :: VA_ave(:,:), V0(:,:), Vn(:,:), B_general(:,:,:,:), JH(:,:)

        open(unit=102, file=path, action='read', access='sequential', &
                form='unformatted', status='old')
        read(102) VA_ave
        read(102) V0
        read(102) Vn
        read(102) B_general
        read(102) JH
        close(102)
        print *, 'read ', trim(path)
end subroutine read_vfunc

subroutine write_vfunc(path, VA_ave, V0, Vn, B_general, JH)
        character(len=*), intent(in) :: path
        real(dp), intent(in) :: VA_ave(:,:), V0(:,:), Vn(:,:), B_general(:,:,:,:), JH(:,:)

        open(unit=102, file=path, action='write', access='sequential', &
                form='unformatted', status='replace')
        write(102) VA_ave
        write(102) V0
        write(102) Vn
        write(102) B_general
        write(102) JH
        close(102)
end subroutine write_vfunc

subroutine project_b_general(nx, ny, nz, B_old, B_out)
        integer, intent(in) :: nx, ny, nz
        real(dp), intent(in) :: B_old(:,:,:,:)
        real(dp), intent(out) :: B_out(:,:,:,:)

        integer :: op, ix, iy, iz, j_src, j_dst, i_src, i_dst, ix_dst, iy_dst
        integer :: channel_map(n_spinor)
        integer, allocatable :: count(:,:)
        real(dp) :: rot(3,3), block_in(3,3), block_rot(3,3)

        allocate(count(size(B_out,1), n_spinor))
        B_out = 0.0_dp
        count = 0

        do op = 0, 2
                call rotation_matrix(op, rot)
                call build_channel_map(op, channel_map)
                print '(A,I1,A,6I4)', 'op ', op, ' channel map = ', channel_map(:)

                do iz = 0, nz - 1
                do iy = 0, ny - 1
                do ix = 0, nx - 1
                        i_src = index_of(ix, iy, iz, nx, ny)
                        call map_grid_point(op, ix, iy, nx, ny, ix_dst, iy_dst)
                        i_dst = index_of(ix_dst, iy_dst, iz, nx, ny)

                        do j_src = 1, n_spinor
                                j_dst = channel_map(j_src)
                                block_in(:,:) = B_old(i_src,:,j_src,:)
                                call rotate_block(rot, block_in, block_rot)
                                B_out(i_dst,:,j_dst,:) = B_out(i_dst,:,j_dst,:) + block_rot(:,:)
                                count(i_dst,j_dst) = count(i_dst,j_dst) + 1
                        enddo
                enddo
                enddo
                enddo
        enddo

        do iz = 0, nz - 1
        do iy = 0, ny - 1
        do ix = 0, nx - 1
                i_src = index_of(ix, iy, iz, nx, ny)
                do j_src = 1, n_spinor
                        if (count(i_src,j_src) /= n_ops) then
                                print *, 'ERROR: incomplete C3 averaging at grid/channel ', &
                                        ix, iy, iz, j_src, ' count=', count(i_src,j_src)
                                stop 1
                        endif
                        B_out(i_src,:,j_src,:) = B_out(i_src,:,j_src,:) / dble(count(i_src,j_src))
                enddo
        enddo
        enddo
        enddo

        call fill_periodic_boundaries(nx, ny, nz, B_out)
        deallocate(count)
end subroutine project_b_general

subroutine rotation_matrix(op, rot)
        integer, intent(in) :: op
        real(dp), intent(out) :: rot(3,3)
        real(dp) :: angle

        rot = 0.0_dp
        select case (op)
        case (0)
                rot(1,1) = 1.0_dp
                rot(2,2) = 1.0_dp
                rot(3,3) = 1.0_dp
        case (1)
                angle = 2.0_dp * pi / 3.0_dp
                rot(1,1) = cos(angle)
                rot(1,2) = -sin(angle)
                rot(2,1) = sin(angle)
                rot(2,2) = cos(angle)
                rot(3,3) = 1.0_dp
        case (2)
                angle = 4.0_dp * pi / 3.0_dp
                rot(1,1) = cos(angle)
                rot(1,2) = -sin(angle)
                rot(2,1) = sin(angle)
                rot(2,2) = cos(angle)
                rot(3,3) = 1.0_dp
        end select
end subroutine rotation_matrix

subroutine rotate_block(rot, block_in, block_out)
        real(dp), intent(in) :: rot(3,3), block_in(3,3)
        real(dp), intent(out) :: block_out(3,3)
        integer :: a, b, c, d

        block_out = 0.0_dp
        do a = 1, 3
        do b = 1, 3
        do c = 1, 3
        do d = 1, 3
                block_out(a,b) = block_out(a,b) + rot(a,c) * block_in(c,d) * rot(b,d)
        enddo
        enddo
        enddo
        enddo
end subroutine rotate_block

subroutine build_channel_map(op, channel_map)
        integer, intent(in) :: op
        integer, intent(out) :: channel_map(n_spinor)

        integer :: j, k, tx, ty, best_tx, best_ty, best_k
        integer :: used(n_spinor), trial_map(n_spinor), best_map(n_spinor)
        real(dp) :: pos(3), rot_pos(3), shifted(3), dist, best_dist
        real(dp) :: score, best_score, max_err, best_max_err
        real(dp) :: pt_pos(n_spinor,3)

        call init_channel_positions(pt_pos)

        if (op == 0) then
                do j = 1, n_spinor
                        channel_map(j) = j
                enddo
                return
        endif

        best_score = huge(1.0_dp)
        best_max_err = huge(1.0_dp)
        best_map = 0
        best_tx = 0
        best_ty = 0

        do tx = -3, 3
        do ty = -3, 3
                used = 0
                trial_map = 0
                score = 0.0_dp
                max_err = 0.0_dp

                do j = 1, n_spinor
                        pos(:) = pt_pos(j,:)
                        call apply_c3_frac(op, pos, c3_center, rot_pos)
                        shifted(:) = rot_pos(:)
                        shifted(1) = shifted(1) + dble(tx)
                        shifted(2) = shifted(2) + dble(ty)

                        best_dist = huge(1.0_dp)
                        best_k = 0
                        do k = 1, n_spinor
                                if (used(k) /= 0) cycle
                                dist = frac_metric_dist2(shifted(:) - pt_pos(k,:))
                                if (dist < best_dist) then
                                        best_dist = dist
                                        best_k = k
                                endif
                        enddo

                        trial_map(j) = best_k
                        used(best_k) = 1
                        score = score + best_dist
                        max_err = max(max_err, best_dist)
                enddo

                if (score < best_score) then
                        best_score = score
                        best_max_err = max_err
                        best_map(:) = trial_map(:)
                        best_tx = tx
                        best_ty = ty
                endif
        enddo
        enddo

        if (best_max_err > map_tol) then
                print *, 'ERROR: six fitted channels are not closed under this C3 operation.'
                print *, 'op = ', op, ' center = ', c3_center(:)
                print *, 'best lattice shift = ', best_tx, best_ty
                print *, 'best max channel-map error = ', best_max_err
                print *, 'best attempted channel map = ', best_map(:)
                print *, 'Adjust c3_center or extend/remap the fitted channel basis before symmetrizing.'
                stop 1
        endif

        channel_map(:) = best_map(:)
end subroutine build_channel_map

subroutine init_channel_positions(pt_pos)
        real(dp), intent(out) :: pt_pos(n_spinor,3)

        pt_pos(1,:) = (/  2.0_dp/3.0_dp,  1.0_dp/3.0_dp, 0.5_dp /)
        pt_pos(2,:) = (/  1.0_dp/3.0_dp,  2.0_dp/3.0_dp, 0.5_dp /)
        pt_pos(3,:) = (/  1.0_dp/3.0_dp, -1.0_dp/3.0_dp, 0.5_dp /)
        pt_pos(4,:) = (/  4.0_dp/3.0_dp,  2.0_dp/3.0_dp, 0.5_dp /)
        pt_pos(5,:) = (/ -1.0_dp/3.0_dp,  1.0_dp/3.0_dp, 0.5_dp /)
        pt_pos(6,:) = (/  2.0_dp/3.0_dp,  4.0_dp/3.0_dp, 0.5_dp /)
end subroutine init_channel_positions

subroutine map_grid_point(op, ix, iy, nx, ny, ix_dst, iy_dst)
        integer, intent(in) :: op, ix, iy, nx, ny
        integer, intent(out) :: ix_dst, iy_dst
        real(dp) :: pos(3), rot_pos(3)

        if (op == 0) then
                ix_dst = ix
                iy_dst = iy
                return
        endif

        pos(1) = dble(ix) / dble(nx)
        pos(2) = dble(iy) / dble(ny)
        pos(3) = c3_center(3)
        call apply_c3_frac(op, pos, c3_center, rot_pos)

        ix_dst = wrap_index(nint(rot_pos(1) * dble(nx)), nx)
        iy_dst = wrap_index(nint(rot_pos(2) * dble(ny)), ny)
end subroutine map_grid_point

subroutine apply_c3_frac(op, pos_in, center, pos_out)
        integer, intent(in) :: op
        real(dp), intent(in) :: pos_in(3), center(3)
        real(dp), intent(out) :: pos_out(3)
        real(dp) :: x, y

        x = pos_in(1) - center(1)
        y = pos_in(2) - center(2)

        select case (op)
        case (0)
                pos_out(:) = pos_in(:)
        case (1)
                pos_out(1) = -y + center(1)
                pos_out(2) = x - y + center(2)
                pos_out(3) = pos_in(3)
        case (2)
                pos_out(1) = y - x + center(1)
                pos_out(2) = -x + center(2)
                pos_out(3) = pos_in(3)
        end select
end subroutine apply_c3_frac

real(dp) function frac_metric_dist2(delta)
        real(dp), intent(in) :: delta(3)
        real(dp) :: cx, cy

        cx = delta(1) - 0.5_dp * delta(2)
        cy = sqrt(3.0_dp) * 0.5_dp * delta(2)
        frac_metric_dist2 = cx * cx + cy * cy + delta(3) * delta(3)
end function frac_metric_dist2

integer function wrap_index(i, n)
        integer, intent(in) :: i, n

        wrap_index = modulo(i, n)
end function wrap_index

integer function index_of(ix, iy, iz, nx, ny)
        integer, intent(in) :: ix, iy, iz, nx, ny

        index_of = 1 + ix + (nx + 1) * iy + (nx + 1) * (ny + 1) * iz
end function index_of

subroutine fill_periodic_boundaries(nx, ny, nz, B)
        integer, intent(in) :: nx, ny, nz
        real(dp), intent(inout) :: B(:,:,:,:)
        integer :: ix, iy, iz, i_full, i_base

        do iz = 0, nz
        do iy = 0, ny
        do ix = 0, nx
                i_full = index_of(ix, iy, iz, nx, ny)
                i_base = index_of(wrap_index(ix, nx), wrap_index(iy, ny), &
                        wrap_index(iz, nz), nx, ny)
                B(i_full,:,:,:) = B(i_base,:,:,:)
        enddo
        enddo
        enddo
end subroutine fill_periodic_boundaries

subroutine c3_residual(nx, ny, nz, B, residual)
        integer, intent(in) :: nx, ny, nz
        real(dp), intent(in) :: B(:,:,:,:)
        real(dp), intent(out) :: residual

        integer :: op, ix, iy, iz, j_src, j_dst, i_src, i_dst, ix_dst, iy_dst
        integer :: channel_map(n_spinor)
        real(dp) :: rot(3,3), block_in(3,3), block_rot(3,3)
        real(dp) :: denom, err

        residual = 0.0_dp
        denom = maxval(abs(B))
        if (denom == 0.0_dp) denom = 1.0_dp

        do op = 1, 2
                call rotation_matrix(op, rot)
                call build_channel_map(op, channel_map)
                do iz = 0, nz - 1
                do iy = 0, ny - 1
                do ix = 0, nx - 1
                        i_src = index_of(ix, iy, iz, nx, ny)
                        call map_grid_point(op, ix, iy, nx, ny, ix_dst, iy_dst)
                        i_dst = index_of(ix_dst, iy_dst, iz, nx, ny)
                        do j_src = 1, n_spinor
                                j_dst = channel_map(j_src)
                                block_in(:,:) = B(i_src,:,j_src,:)
                                call rotate_block(rot, block_in, block_rot)
                                err = maxval(abs(block_rot(:,:) - B(i_dst,:,j_dst,:))) / denom
                                residual = max(residual, err)
                        enddo
                enddo
                enddo
                enddo
        enddo
end subroutine c3_residual

end program symmetrize_c3
