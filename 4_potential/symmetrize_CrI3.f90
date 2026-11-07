program symmetrize_CrI3
  use parameter_mod
        implicit none

        integer, parameter :: dp = kind(0d0)
        integer, parameter :: n_spinor = p_split_nspinor
        integer, parameter :: n_sym = 6
        logical, parameter :: use_nonclosed_symmetry = p_use_nonclosed_symmetry
        real(dp), parameter :: spin_isotropy_mix = p_spin_isotropy_mix

        integer :: nx, ny, nz, ni
        integer :: ix, iy, iz, i, j, isym, j_dir, j_map
        integer :: i_dir, ix2, iy2, iz2
        integer :: pt_map(n_spinor, n_sym)
        logical :: valid_sym(n_sym)
        logical :: backup_exists
        integer, allocatable :: deg_tmp(:,:)
        real(dp) :: fin(3), fout(3)
        real(dp) :: block(3,3), block_rot(3,3), trace_ave
        real(dp) :: pt_pos(n_spinor,3)
        real(dp) :: a_lat(3,3), cart_basis(3,3), inv_cart_basis(3,3)
        real(dp) :: sym_crys(3,3,n_sym), sym_cart(3,3,n_sym)
        character(len=24) :: sym_name(n_sym)

        real(dp), allocatable :: VA_ave(:,:), V0(:,:), Vn(:,:), JH(:,:)
        real(dp), allocatable :: B_general(:,:,:,:), B_general_old(:,:,:,:)

        call init_lattice(a_lat)
        cart_basis = transpose(a_lat)
        call inverse3(cart_basis, inv_cart_basis)
        call init_symmetry(sym_crys, sym_name)

        do isym = 1, n_sym
                sym_cart(:,:,isym) = matmul(cart_basis, matmul(sym_crys(:,:,isym), inv_cart_basis))
        enddo

        open(unit=101, file=p_scf_dir//'V_HYC_dp', access='sequential', action='read', &
             form='unformatted', status='old')
                read(101) nx, ny, nz
        close(101)
        print *, 'nx ny nz = ', nx, ny, nz

        ni = (nx+1) * (ny+1) * (nz+1)

        allocate(VA_ave(ni,4))
        allocate(V0(ni,3))
        allocate(Vn(ni,8))
        allocate(JH(ni,n_spinor))
        allocate(B_general_old(ni,3,n_spinor,3))
        allocate(B_general(ni,3,n_spinor,3))
        allocate(deg_tmp(ni,n_spinor))

        open(unit=102, file='V_func_tmp_tmp', action='read', access='sequential', &
             form='unformatted', status='old')
                read(102) VA_ave
                read(102) V0
                read(102) Vn
                read(102) B_general_old
                read(102) JH
        close(102)

        inquire(file='V_func_tmp_tmp.before_sym', exist=backup_exists)
        if (backup_exists) then
                print *, 'backup V_func_tmp_tmp.before_sym already exists; it is kept unchanged'
        else
                call write_vfunc('V_func_tmp_tmp.before_sym', VA_ave, V0, Vn, B_general_old, JH)
                print *, 'backup written to V_func_tmp_tmp.before_sym'
        endif

        pt_pos(1,:) = (/ 1.0_dp/3.0_dp, 2.0_dp/3.0_dp, 0.5_dp /)
        pt_pos(2,:) = (/ 2.0_dp/3.0_dp, 1.0_dp/3.0_dp, 0.5_dp /)
        call build_pt_map(pt_pos, sym_crys, sym_name, pt_map, valid_sym)

        B_general = 0.0_dp
        deg_tmp = 0

        do ix = 0, nx-1
        do iy = 0, ny-1
        do iz = 0, nz-1
                i = index3(ix, iy, iz)
                fin = (/ dble(ix)/dble(nx), dble(iy)/dble(ny), dble(iz)/dble(nz) /)

                do isym = 1, n_sym
                        if ((.not. use_nonclosed_symmetry) .and. (.not. valid_sym(isym))) cycle

                        fout = matmul(sym_crys(:,:,isym), fin)
                        call wrap_frac(fout)
                        ix2 = modulo(nint(fout(1) * dble(nx)), nx)
                        iy2 = modulo(nint(fout(2) * dble(ny)), ny)
                        iz2 = modulo(nint(fout(3) * dble(nz)), nz)
                        j = index3(ix2, iy2, iz2)

                        do j_dir = 1, n_spinor
                                j_map = pt_map(j_dir, isym)
                                block(:,:) = B_general_old(i,:,j_dir,:)
                                block_rot = matmul(sym_cart(:,:,isym), &
                                           matmul(block, transpose(sym_cart(:,:,isym))))
                                B_general(j,:,j_map,:) = B_general(j,:,j_map,:) + block_rot
                                deg_tmp(j,j_map) = deg_tmp(j,j_map) + 1
                        enddo
                enddo
        enddo
        enddo
        enddo

        do i = 1, ni
        do j_dir = 1, n_spinor
                if (deg_tmp(i,j_dir) > 0) then
                        B_general(i,:,j_dir,:) = B_general(i,:,j_dir,:) / dble(deg_tmp(i,j_dir))
                else
                        B_general(i,:,j_dir,:) = B_general_old(i,:,j_dir,:)
                endif
        enddo
        enddo

        call copy_periodic_boundaries(B_general)

        if (spin_isotropy_mix > 0.0_dp) then
                print *, 'spin_isotropy_mix = ', spin_isotropy_mix
                print *, 'each 3x3 spin block is softly mixed toward trace(B)/3 * I'
                do i = 1, ni
                do j_dir = 1, n_spinor
                        block(:,:) = B_general(i,:,j_dir,:)
                        trace_ave = (B_general(i,1,j_dir,1) + B_general(i,2,j_dir,2) + &
                                     B_general(i,3,j_dir,3)) / 3.0_dp
                        B_general(i,:,j_dir,:) = 0.0_dp
                        do i_dir = 1, 3
                                B_general(i,i_dir,j_dir,i_dir) = trace_ave
                        enddo
                        B_general(i,:,j_dir,:) = (1.0_dp - spin_isotropy_mix) * block(:,:) + &
                                                  spin_isotropy_mix * B_general(i,:,j_dir,:)
                enddo
                enddo
        else
                print *, 'spin_isotropy_mix = 0: spin tensor is not projected'
        endif

        call write_vfunc('V_func_tmp_tmp', VA_ave, V0, Vn, B_general, JH)
        print *, 'symmetrized V_func_tmp_tmp written'

contains

subroutine init_lattice(a)
        real(dp), intent(out) :: a(3,3)

        a(1,:) = (/  1.000000_dp, 0.000000_dp, 0.000000_dp /)
        a(2,:) = (/ -0.500000_dp, 0.866025_dp, 0.000000_dp /)
        a(3,:) = (/  0.000000_dp, 0.000000_dp, p_celldm3 /)
end subroutine init_lattice

subroutine init_symmetry(sym, names)
        real(dp), intent(out) :: sym(3,3,n_sym)
        character(len=24), intent(out) :: names(n_sym)

        sym = 0.0_dp

        names(1) = 'E'
        sym(:,:,1) = reshape((/ &
             1.0_dp,  0.0_dp, 0.0_dp, &
             0.0_dp,  1.0_dp, 0.0_dp, &
             0.0_dp,  0.0_dp, 1.0_dp /), (/3,3/), order=(/2,1/))

        names(2) = 'C3'
        sym(:,:,2) = reshape((/ &
             0.0_dp,  1.0_dp, 0.0_dp, &
            -1.0_dp, -1.0_dp, 0.0_dp, &
             0.0_dp,  0.0_dp, 1.0_dp /), (/3,3/), order=(/2,1/))

        names(3) = 'C3^2'
        sym(:,:,3) = reshape((/ &
            -1.0_dp, -1.0_dp, 0.0_dp, &
             1.0_dp,  0.0_dp, 0.0_dp, &
             0.0_dp,  0.0_dp, 1.0_dp /), (/3,3/), order=(/2,1/))

        names(4) = 'Mxy'
        sym(:,:,4) = reshape((/ &
             0.0_dp,  1.0_dp, 0.0_dp, &
             1.0_dp,  0.0_dp, 0.0_dp, &
             0.0_dp,  0.0_dp, 1.0_dp /), (/3,3/), order=(/2,1/))

        names(5) = 'Mxy*C3'
        sym(:,:,5) = matmul(sym(:,:,4), sym(:,:,2))

        names(6) = 'Mxy*C3^2'
        sym(:,:,6) = matmul(sym(:,:,4), sym(:,:,3))
end subroutine init_symmetry

subroutine build_pt_map(pt, sym, names, map_out, valid_out)
        real(dp), intent(in) :: pt(n_spinor,3), sym(3,3,n_sym)
        character(len=24), intent(in) :: names(n_sym)
        integer, intent(out) :: map_out(n_spinor,n_sym)
        logical, intent(out) :: valid_out(n_sym)

        integer :: is, ia, ja
        real(dp) :: p2(3), diff(3), best_dist, dist

        do is = 1, n_sym
        valid_out(is) = .true.
        do ia = 1, n_spinor
                p2 = matmul(sym(:,:,is), pt(ia,:))
                call wrap_frac(p2)
                best_dist = huge(1.0_dp)
                map_out(ia,is) = ia
                do ja = 1, n_spinor
                        diff = p2 - pt(ja,:)
                        diff = diff - dnint(diff)
                        dist = sum(diff*diff)
                        if (dist < best_dist) then
                                best_dist = dist
                                map_out(ia,is) = ja
                        endif
                enddo
                if (best_dist > 1.0d-8) then
                        valid_out(is) = .false.
                        print *, 'warning: ', trim(names(is)), ' maps center ', ia, &
                                 ' outside the 2-center basis; nearest center would be ', map_out(ia,is)
                endif
        enddo
        if (valid_out(is)) then
                print '(A,A,2I4,A)', 'pt_map for ', trim(names(is)), map_out(1,is), &
                        map_out(2,is), '  [used]'
        else if (use_nonclosed_symmetry) then
                print '(A,A,2I4,A)', 'pt_map for ', trim(names(is)), map_out(1,is), &
                        map_out(2,is), '  [forced]'
        else
                print '(A,A,2I4,A)', 'pt_map for ', trim(names(is)), map_out(1,is), &
                        map_out(2,is), '  [skipped]'
        endif
        enddo
end subroutine build_pt_map

subroutine wrap_frac(x)
        real(dp), intent(inout) :: x(3)
        integer :: k

        do k = 1, 3
                x(k) = x(k) - floor(x(k))
                if (x(k) >= 1.0_dp - 1.0d-10) x(k) = 0.0_dp
                if (x(k) < 1.0d-10) x(k) = 0.0_dp
        enddo
end subroutine wrap_frac

integer function index3(ix_in, iy_in, iz_in)
        integer, intent(in) :: ix_in, iy_in, iz_in

        index3 = 1 + ix_in + (nx+1)*iy_in + (nx+1)*(ny+1)*iz_in
end function index3

subroutine copy_periodic_boundaries(B)
        real(dp), intent(inout) :: B(:,:,:,:)
        integer :: lx, ly, lz

        do ly = 0, ny
        do lz = 0, nz
                B(index3(nx,ly,lz),:,:,:) = B(index3(0,ly,lz),:,:,:)
        enddo
        enddo

        do lx = 0, nx
        do lz = 0, nz
                B(index3(lx,ny,lz),:,:,:) = B(index3(lx,0,lz),:,:,:)
        enddo
        enddo

        do lx = 0, nx
        do ly = 0, ny
                B(index3(lx,ly,nz),:,:,:) = B(index3(lx,ly,0),:,:,:)
        enddo
        enddo
end subroutine copy_periodic_boundaries

subroutine write_vfunc(filename, VA, V0_in, Vn_in, B_in, JH_in)
        character(len=*), intent(in) :: filename
        real(dp), intent(in) :: VA(:,:), V0_in(:,:), Vn_in(:,:), B_in(:,:,:,:), JH_in(:,:)

        open(unit=202, file=filename, action='write', access='sequential', &
             status='replace', form='unformatted')
                write(202) VA
                write(202) V0_in
                write(202) Vn_in
                write(202) B_in
                write(202) JH_in
        close(202)
end subroutine write_vfunc

subroutine inverse3(a, ainv)
        real(dp), intent(in) :: a(3,3)
        real(dp), intent(out) :: ainv(3,3)
        real(dp) :: det

        det = a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2)) &
            - a(1,2)*(a(2,1)*a(3,3)-a(2,3)*a(3,1)) &
            + a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))

        if (abs(det) < 1.0d-14) then
                print *, 'inverse3 error: singular matrix'
                stop
        endif

        ainv(1,1) =  (a(2,2)*a(3,3)-a(2,3)*a(3,2)) / det
        ainv(1,2) = -(a(1,2)*a(3,3)-a(1,3)*a(3,2)) / det
        ainv(1,3) =  (a(1,2)*a(2,3)-a(1,3)*a(2,2)) / det
        ainv(2,1) = -(a(2,1)*a(3,3)-a(2,3)*a(3,1)) / det
        ainv(2,2) =  (a(1,1)*a(3,3)-a(1,3)*a(3,1)) / det
        ainv(2,3) = -(a(1,1)*a(2,3)-a(1,3)*a(2,1)) / det
        ainv(3,1) =  (a(2,1)*a(3,2)-a(2,2)*a(3,1)) / det
        ainv(3,2) = -(a(1,1)*a(3,2)-a(1,2)*a(3,1)) / det
        ainv(3,3) =  (a(1,1)*a(2,2)-a(1,2)*a(2,1)) / det
end subroutine inverse3

end program symmetrize_CrI3
