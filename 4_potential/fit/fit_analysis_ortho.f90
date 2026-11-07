program fit_analysis_ortho
  use parameter_mod
        use wannier_weight_mod, only: load_wannier_weight
        use cri3_block_mask_mod, only: block_mask_mode, make_feature_block_mask, &
                apply_com_block_mask, count_blocked_features

        ! Orthogonal-basis diagnostic for files written by fit_com_ortho.f90.
        !
        ! The file format stores raw coefficients:
        !
        ! fit_com_ortho.f90 fitted instead:
        !
        ! Therefore:
        !
        ! This program reconstructs P from the training spins, transforms the
        ! stored coefficients back to the orthogonal basis, and reports the
        ! real linear/modulated contributions in that basis.

        implicit none

        integer, parameter :: dp = kind(0d0)
        integer, parameter :: n_spinor = p_fit_nspinor
        integer, parameter :: n_file = p_n_file
        integer, parameter :: n_kind = p_n_kind
        integer, parameter :: n_sample = n_file * n_kind
        integer, parameter :: n_comp = 3
        integer, parameter :: n_feat = 2 * n_spinor
        character(len=32), parameter :: expected_magic = 'CRI3_COM_FIT_CHANNEL_V1'

        integer :: nx, ny, nz, ni
        integer :: file_nx, file_ny, file_nz, file_grid
        integer :: file_spinor, file_feat, file_count, kind_count
        integer :: max_iter_file, nonconv_file
        integer :: i_sample, i_file, i_kind, j, c, n_read
        integer :: i_target, blocked_features
        real(dp) :: l1_rel, l2_rel, alpha1, alpha2, tol
        real(dp) :: target_norm, scalar_norm
        real(dp) :: linear_err, full_err, mod_gain_err, mod_gain_rms
        real(dp) :: comp_target(3), comp_linear_err(3), comp_full_err(3)
        real(dp) :: raw_lambda, raw_lambda_resid
        real(dp) :: ortho_lambda, ortho_lambda_resid
        real(dp) :: corr_before, corr_after, max_mod_mod_corr
        real(dp) :: projection(n_spinor,n_spinor)
        real(dp) :: raw_K0_rms(6), raw_K1_rms(6)
        real(dp) :: ortho_K0_rms(6), ortho_K1_rms(6)
        real(dp) :: raw_K0_max(6), raw_K1_max(6)
        real(dp) :: ortho_K0_max(6), ortho_K1_max(6)
        real(dp) :: lin_eff_norm(6), mod_eff_norm(6), eff_cross(6)
        real(dp) :: raw_lin_eff_norm(6), raw_mod_eff_norm(6), raw_eff_cross(6)
        real(dp), allocatable :: K0_raw(:,:), K1_raw(:,:)
        real(dp), allocatable :: K0_ortho(:,:), K1_ortho(:,:)
        real(dp), allocatable :: V_read(:,:)
        real(dp), allocatable :: spin(:,:,:)
        real(dp), allocatable :: feature_raw(:,:,:), feature_ortho(:,:,:)
        real(dp), allocatable :: lin_pred(:,:), mod_pred(:,:)
        real(dp), allocatable :: raw_lin_pred(:,:), raw_mod_pred(:,:)
        real(dp), allocatable :: grid_weight(:)
        real(dp), allocatable :: block_mask(:,:)
        real(dp) :: lin_piece, mod_piece
        character(len=32) :: magic
        character(len=200) :: input_file, path, weight_file, report_file

        weight_file = 'w_density'

        call read_grid(nx, ny, nz, ni)

        do i_target = 1, 2
        if (i_target == 1) then
                input_file = 'V_func_com_ortho_channel_v1'
                report_file = 'fit_analysis_ortho_raw_report.dat'
        else
                input_file = 'V_func_com_ortho_channel_v1_sym'
                report_file = 'fit_analysis_ortho_sym_report.dat'
        endif

        allocate(K0_raw(ni,n_spinor))
        allocate(K1_raw(ni,n_spinor))
        allocate(K0_ortho(ni,n_spinor))
        allocate(K1_ortho(ni,n_spinor))
        allocate(V_read(ni,4))
        allocate(spin(n_sample,n_spinor,n_comp))
        allocate(feature_raw(n_sample,n_feat,n_comp))
        allocate(feature_ortho(n_sample,n_feat,n_comp))
        allocate(lin_pred(ni,3))
        allocate(mod_pred(ni,3))
        allocate(raw_lin_pred(ni,3))
        allocate(raw_mod_pred(ni,3))
        allocate(grid_weight(ni))
        call load_wannier_weight(trim(weight_file), nx, ny, nz, ni, grid_weight)

        call read_com_file(trim(input_file), magic, file_nx, file_ny, file_nz,    &
                file_grid, file_spinor, file_feat, file_count, kind_count,        &
                max_iter_file, nonconv_file, l1_rel, l2_rel, alpha1, alpha2,      &
                tol, K0_raw, K1_raw)
        call validate_header(magic, file_nx, file_ny, file_nz, file_grid,         &
                file_spinor, file_feat)
        allocate(block_mask(file_grid,n_feat))
        call make_feature_block_mask(file_nx, file_ny, file_nz, file_grid,        &
                block_mask)
        blocked_features = count_blocked_features(block_mask)
        call apply_com_block_mask(K0_raw, K1_raw, block_mask)

        call read_spins(spin)
        do i_sample = 1, n_sample
                call make_features(spin(i_sample,:,:), feature_raw(i_sample,:,:))
        enddo

        call build_projection(feature_raw, projection, corr_before)
        call make_ortho_features(feature_raw, projection, feature_ortho)
        call max_linear_mod_corr(feature_ortho, corr_after)
        call max_mod_mod_corr_sub(feature_ortho, max_mod_mod_corr)
        call raw_to_ortho_coefficients(K0_raw, K1_raw, projection,                &
                K0_ortho, K1_ortho)
        call apply_com_block_mask(K0_ortho, K1_ortho, block_mask)

        call kernel_summary(K0_raw, K1_raw, raw_K0_rms, raw_K1_rms,               &
                raw_K0_max, raw_K1_max)
        call kernel_summary(K0_ortho, K1_ortho, ortho_K0_rms, ortho_K1_rms,       &
                ortho_K0_max, ortho_K1_max)
        call lambda_like_diagnostic(K0_raw, K1_raw, raw_lambda, raw_lambda_resid)
        call lambda_like_diagnostic(K0_ortho, K1_ortho, ortho_lambda,             &
                ortho_lambda_resid)

        target_norm = 0d0
        scalar_norm = 0d0
        linear_err = 0d0
        full_err = 0d0
        comp_target(:) = 0d0
        comp_linear_err(:) = 0d0
        comp_full_err(:) = 0d0
        lin_eff_norm(:) = 0d0
        mod_eff_norm(:) = 0d0
        eff_cross(:) = 0d0
        raw_lin_eff_norm(:) = 0d0
        raw_mod_eff_norm(:) = 0d0
        raw_eff_cross(:) = 0d0

        i_sample = 0
        do i_file = 1, n_file
        do i_kind = 1, n_kind
                i_sample = i_sample + 1
                call make_path(i_file, i_kind, path)
                call read_V(path, V_read, n_read)
                if (n_read .ne. ni) then
                        print *, 'warning: V grid record differs from ni:',      &
                                n_read, ni
                endif

                lin_pred(:,:) = 0d0
                mod_pred(:,:) = 0d0
                raw_lin_pred(:,:) = 0d0
                raw_mod_pred(:,:) = 0d0
                do j = 1, n_spinor
                do c = 1, 3
                        lin_piece = feature_ortho(i_sample,j,c)
                        mod_piece = feature_ortho(i_sample,n_spinor+j,c)
                        lin_pred(:,c) = lin_pred(:,c)                            &
                                + K0_ortho(:,j) * lin_piece
                        mod_pred(:,c) = mod_pred(:,c)                            &
                                + K1_ortho(:,j) * mod_piece
                        lin_eff_norm(j) = lin_eff_norm(j)                        &
                                + sum(grid_weight*(K0_ortho(:,j) * lin_piece)**2)
                        mod_eff_norm(j) = mod_eff_norm(j)                        &
                                + sum(grid_weight*(K1_ortho(:,j) * mod_piece)**2)
                        eff_cross(j) = eff_cross(j)                              &
                                + sum(grid_weight*(K0_ortho(:,j) * lin_piece)   &
                                * (K1_ortho(:,j) * mod_piece))

                        lin_piece = feature_raw(i_sample,j,c)
                        mod_piece = feature_raw(i_sample,n_spinor+j,c)
                        raw_lin_pred(:,c) = raw_lin_pred(:,c)                    &
                                + K0_raw(:,j) * lin_piece
                        raw_mod_pred(:,c) = raw_mod_pred(:,c)                    &
                                + K1_raw(:,j) * mod_piece
                        raw_lin_eff_norm(j) = raw_lin_eff_norm(j)                &
                                + sum(grid_weight*(K0_raw(:,j) * lin_piece)**2)
                        raw_mod_eff_norm(j) = raw_mod_eff_norm(j)                &
                                + sum(grid_weight*(K1_raw(:,j) * mod_piece)**2)
                        raw_eff_cross(j) = raw_eff_cross(j)                      &
                                + sum(grid_weight*(K0_raw(:,j) * lin_piece)     &
                                * (K1_raw(:,j) * mod_piece))
                enddo
                enddo

                scalar_norm = scalar_norm + sum(grid_weight*V_read(:,1)**2)
                do c = 1, 3
                        target_norm = target_norm                              &
                                + sum(grid_weight*V_read(:,c+1)**2)
                        comp_target(c) = comp_target(c)                        &
                                + sum(grid_weight*V_read(:,c+1)**2)
                        linear_err = linear_err                                  &
                                + sum(grid_weight*(V_read(:,c+1)-lin_pred(:,c))**2)
                        full_err = full_err                                      &
                                + sum(grid_weight*(V_read(:,c+1)-lin_pred(:,c)  &
                                - mod_pred(:,c))**2)
                        comp_linear_err(c) = comp_linear_err(c)                  &
                                + sum(grid_weight*(V_read(:,c+1)-lin_pred(:,c))**2)
                        comp_full_err(c) = comp_full_err(c)                      &
                                + sum(grid_weight*(V_read(:,c+1)-lin_pred(:,c)  &
                                - mod_pred(:,c))**2)
                enddo
        enddo
        enddo

        mod_gain_err = 1d0 - full_err / max(linear_err, 1d-300)
        mod_gain_rms = 1d0 - sqrt(full_err / max(linear_err, 1d-300))

        call write_report(trim(report_file), trim(input_file), file_nx, file_ny,  &
                file_nz,                                                         &
                max_iter_file, nonconv_file, l1_rel, l2_rel, alpha1, alpha2,      &
                tol, target_norm, scalar_norm, linear_err, full_err,              &
                mod_gain_err, mod_gain_rms, comp_target, comp_linear_err,         &
                comp_full_err, corr_before, corr_after, max_mod_mod_corr,         &
                projection, raw_lambda, raw_lambda_resid, ortho_lambda,           &
                ortho_lambda_resid, raw_K0_rms, raw_K1_rms, raw_K0_max,           &
                raw_K1_max, ortho_K0_rms, ortho_K1_rms, ortho_K0_max,             &
                ortho_K1_max, lin_eff_norm, mod_eff_norm, eff_cross,              &
                raw_lin_eff_norm, raw_mod_eff_norm, raw_eff_cross, weight_file,   &
                blocked_features)

        print *, 'fit_analysis_ortho done'
        print *, 'input  = ', trim(input_file)
        print *, 'report = ', trim(report_file)
        print *, 'block mask mode = ', block_mask_mode
        print *, 'blocked feature coefficients = ', blocked_features
        print *, 'corr before/after = ', corr_before, corr_after
        print *, 'ortho linear relative error = ', sqrt(linear_err / target_norm)
        print *, 'ortho full relative error   = ', sqrt(full_err / target_norm)
        deallocate(K0_raw, K1_raw, K0_ortho, K1_ortho, V_read, spin,             &
                feature_raw, feature_ortho, lin_pred, mod_pred, raw_lin_pred,    &
                raw_mod_pred, grid_weight, block_mask)
        enddo

contains

subroutine read_grid(grid_nx, grid_ny, grid_nz, grid_size)
        implicit none
        integer, intent(out) :: grid_nx, grid_ny, grid_nz, grid_size

        open(unit=101,file=p_scf_dir//'V_QE',access='sequential',action='read',     &
                status='old',form='unformatted')
        read(101) grid_nx, grid_ny, grid_nz
        close(101)
        grid_size = (grid_nx + 1) * (grid_ny + 1) * (grid_nz + 1)
end subroutine read_grid

subroutine read_com_file(filename, read_magic, gx, gy, gz, grid_size, spin_count, &
                feat_count, file_count_out, kind_count_out, max_iter_out,         &
                nonconv_out, l1_out, l2_out, a1_out, a2_out, tol_out,             &
                linear_part, mod_part)
        implicit none
        character(len=*), intent(in) :: filename
        character(len=32), intent(out) :: read_magic
        integer, intent(out) :: gx, gy, gz, grid_size, spin_count, feat_count
        integer, intent(out) :: file_count_out, kind_count_out
        integer, intent(out) :: max_iter_out, nonconv_out
        real(dp), intent(out) :: l1_out, l2_out, a1_out, a2_out, tol_out
        real(dp), intent(out) :: linear_part(ni,n_spinor), mod_part(ni,n_spinor)

        open(unit=102,file=trim(filename),access='sequential',action='read',      &
                status='old',form='unformatted')
        read(102) read_magic
        read(102) gx, gy, gz
        read(102) grid_size, spin_count, feat_count, file_count_out,              &
                kind_count_out
        read(102) max_iter_out, nonconv_out
        read(102) l1_out, l2_out, a1_out, a2_out, tol_out
        read(102) linear_part
        read(102) mod_part
        close(102)
end subroutine read_com_file

subroutine validate_header(read_magic, gx, gy, gz, grid_size, spin_count,         &
                feat_count)
        implicit none
        character(len=*), intent(in) :: read_magic
        integer, intent(in) :: gx, gy, gz, grid_size, spin_count, feat_count

        if (trim(read_magic) .ne. trim(expected_magic)) stop 'bad magic'
        if (gx .ne. nx .or. gy .ne. ny .or. gz .ne. nz) stop 'bad grid'
        if (grid_size .ne. ni) stop 'bad grid size'
        if (spin_count .ne. n_spinor) stop 'bad spinor count'
        if (feat_count .ne. n_feat) stop 'bad feature count'
end subroutine validate_header

subroutine read_spins(all_spin)
        implicit none
        real(dp), intent(out) :: all_spin(n_sample,n_spinor,n_comp)
        real(dp) :: flat(3*n_spinor)
        integer :: ifile, ikind, isample, ispin
        character(len=20) :: file_char
        character(len=200) :: spin_path

        isample = 0
        do ifile = 1, n_file
                call int2char(ifile, file_char)
                spin_path = './angle_'//trim(file_char)//'/spin_info'
                open(unit=103,file=trim(spin_path),access='sequential',           &
                        action='read',status='old')
                do ikind = 1, n_kind
                        isample = isample + 1
                        read(103,*) flat(:)
                        do ispin = 1, n_spinor
                                all_spin(isample,ispin,:) =                     &
                                        flat(3*(ispin-1)+1:3*ispin)
                        enddo
                enddo
                close(103)
        enddo
end subroutine read_spins

subroutine make_features(spin_in, feat_out)
        implicit none
        real(dp), intent(in) :: spin_in(n_spinor,n_comp)
        real(dp), intent(out) :: feat_out(n_feat,n_comp)
        real(dp) :: snorm(n_spinor), chi_site(n_spinor), cosab
        integer :: ia, ib

        feat_out(:,:) = 0d0
        do ia = 1, n_spinor
                snorm(ia) = sqrt(sum(spin_in(ia,:)**2))
                if (snorm(ia) .lt. 1d-12) stop 'zero spin in spin_info'
                feat_out(ia,:) = spin_in(ia,:)
        enddo

        chi_site(:) = 0d0
        do ia = 1, n_spinor
                do ib = 1, n_spinor
                        if (ib .eq. ia) cycle
                        cosab = dot_product(spin_in(ia,:),spin_in(ib,:))         &
                                / (snorm(ia) * snorm(ib))
                        chi_site(ia) = chi_site(ia) + cosab - 1d0
                enddo
                chi_site(ia) = chi_site(ia) / dble(n_spinor - 1)
        enddo

        do ia = 1, n_spinor
                feat_out(n_spinor+ia,:) = chi_site(ia) * spin_in(ia,:)
        enddo
end subroutine make_features

subroutine build_projection(raw_feature, proj, max_corr_before)
        implicit none
        real(dp), intent(in) :: raw_feature(n_sample,n_feat,n_comp)
        real(dp), intent(out) :: proj(n_spinor,n_spinor)
        real(dp), intent(out) :: max_corr_before
        real(dp) :: GLL(n_spinor,n_spinor), GLM(n_spinor,n_spinor)
        integer :: is, ia, ib

        GLL(:,:) = 0d0
        GLM(:,:) = 0d0
        do is = 1, n_sample
        do ia = 1, n_spinor
        do ib = 1, n_spinor
                GLL(ia,ib) = GLL(ia,ib)                                         &
                        + dot_product(raw_feature(is,ia,:), raw_feature(is,ib,:))
                GLM(ia,ib) = GLM(ia,ib)                                         &
                        + dot_product(raw_feature(is,ia,:),                      &
                        raw_feature(is,n_spinor+ib,:))
        enddo
        enddo
        enddo

        max_corr_before = 0d0
        do ia = 1, n_spinor
        do ib = 1, n_spinor
                max_corr_before = max(max_corr_before,                           &
                        abs(GLM(ia,ib)) / sqrt(max(GLL(ia,ia)                   &
                        * raw_feature_norm2(raw_feature,n_spinor+ib), 1d-300)))
        enddo
        enddo

        call solve_matrix_rhs(GLL, GLM, proj)
end subroutine build_projection

real(dp) function raw_feature_norm2(raw_feature, idx)
        implicit none
        real(dp), intent(in) :: raw_feature(n_sample,n_feat,n_comp)
        integer, intent(in) :: idx

        raw_feature_norm2 = sum(raw_feature(:,idx,:)**2)
end function raw_feature_norm2

subroutine solve_matrix_rhs(A_in, B_in, X_out)
        implicit none
        real(dp), intent(in) :: A_in(n_spinor,n_spinor), B_in(n_spinor,n_spinor)
        real(dp), intent(out) :: X_out(n_spinor,n_spinor)
        real(dp) :: aug(n_spinor,2*n_spinor)
        real(dp) :: pivot, factor
        integer :: row, piv, best

        aug(:,1:n_spinor) = A_in(:,:)
        aug(:,n_spinor+1:2*n_spinor) = B_in(:,:)

        do piv = 1, n_spinor
                best = piv
                do row = piv + 1, n_spinor
                        if (abs(aug(row,piv)) .gt. abs(aug(best,piv))) best = row
                enddo
                if (abs(aug(best,piv)) .lt. 1d-14) stop 'singular GLL'
                if (best .ne. piv) call swap_rows(aug, best, piv)

                pivot = aug(piv,piv)
                aug(piv,:) = aug(piv,:) / pivot
                do row = 1, n_spinor
                        if (row .eq. piv) cycle
                        factor = aug(row,piv)
                        aug(row,:) = aug(row,:) - factor * aug(piv,:)
                enddo
        enddo

        X_out(:,:) = aug(:,n_spinor+1:2*n_spinor)
end subroutine solve_matrix_rhs

subroutine swap_rows(mat, r1, r2)
        implicit none
        real(dp), intent(inout) :: mat(n_spinor,2*n_spinor)
        integer, intent(in) :: r1, r2
        real(dp) :: tmp(2*n_spinor)

        tmp(:) = mat(r1,:)
        mat(r1,:) = mat(r2,:)
        mat(r2,:) = tmp(:)
end subroutine swap_rows

subroutine make_ortho_features(raw_feature, proj, ortho_feature)
        implicit none
        real(dp), intent(in) :: raw_feature(n_sample,n_feat,n_comp)
        real(dp), intent(in) :: proj(n_spinor,n_spinor)
        real(dp), intent(out) :: ortho_feature(n_sample,n_feat,n_comp)
        integer :: is, ia, ib

        ortho_feature(:,:,:) = raw_feature(:,:,:)
        do is = 1, n_sample
        do ia = 1, n_spinor
                ortho_feature(is,n_spinor+ia,:) = raw_feature(is,n_spinor+ia,:)
                do ib = 1, n_spinor
                        ortho_feature(is,n_spinor+ia,:)                         &
                                = ortho_feature(is,n_spinor+ia,:)                &
                                - proj(ib,ia) * raw_feature(is,ib,:)
                enddo
        enddo
        enddo
end subroutine make_ortho_features

subroutine raw_to_ortho_coefficients(raw0, raw1, proj, ortho0, ortho1)
        implicit none
        real(dp), intent(in) :: raw0(ni,n_spinor), raw1(ni,n_spinor)
        real(dp), intent(in) :: proj(n_spinor,n_spinor)
        real(dp), intent(out) :: ortho0(ni,n_spinor), ortho1(ni,n_spinor)
        integer :: ip, ia, ib

        ortho1(:,:) = raw1(:,:)
        ortho0(:,:) = raw0(:,:)
        do ip = 1, ni
        do ib = 1, n_spinor
                do ia = 1, n_spinor
                        ortho0(ip,ib) = ortho0(ip,ib)                            &
                                + raw1(ip,ia) * proj(ib,ia)
                enddo
        enddo
        enddo
end subroutine raw_to_ortho_coefficients

subroutine max_linear_mod_corr(feature_in, corr)
        implicit none
        real(dp), intent(in) :: feature_in(n_sample,n_feat,n_comp)
        real(dp), intent(out) :: corr
        real(dp) :: num, den
        integer :: ia, ib

        corr = 0d0
        do ia = 1, n_spinor
        do ib = 1, n_spinor
                num = sum(feature_in(:,ia,:) * feature_in(:,n_spinor+ib,:))
                den = sqrt(max(sum(feature_in(:,ia,:)**2)                        &
                        * sum(feature_in(:,n_spinor+ib,:)**2), 1d-300))
                corr = max(corr, abs(num / den))
        enddo
        enddo
end subroutine max_linear_mod_corr

subroutine max_mod_mod_corr_sub(feature_in, corr)
        implicit none
        real(dp), intent(in) :: feature_in(n_sample,n_feat,n_comp)
        real(dp), intent(out) :: corr
        real(dp) :: num, den
        integer :: ia, ib

        corr = 0d0
        do ia = 1, n_spinor
        do ib = 1, n_spinor
                if (ia .eq. ib) cycle
                num = sum(feature_in(:,n_spinor+ia,:)                            &
                        * feature_in(:,n_spinor+ib,:))
                den = sqrt(max(sum(feature_in(:,n_spinor+ia,:)**2)               &
                        * sum(feature_in(:,n_spinor+ib,:)**2), 1d-300))
                corr = max(corr, abs(num / den))
        enddo
        enddo
end subroutine max_mod_mod_corr_sub

subroutine kernel_summary(linear_part, mod_part, lin_norm, mod_norm, lin_max,     &
                mod_max)
        implicit none
        real(dp), intent(in) :: linear_part(ni,n_spinor), mod_part(ni,n_spinor)
        real(dp), intent(out) :: lin_norm(6), mod_norm(6), lin_max(6), mod_max(6)
        integer :: jf

        do jf = 1, n_spinor
                lin_norm(jf) = sqrt(sum(linear_part(:,jf)**2) / dble(ni))
                mod_norm(jf) = sqrt(sum(mod_part(:,jf)**2) / dble(ni))
                lin_max(jf) = maxval(abs(linear_part(:,jf)))
                mod_max(jf) = maxval(abs(mod_part(:,jf)))
        enddo
end subroutine kernel_summary

subroutine lambda_like_diagnostic(linear_part, mod_part, best_ratio, residual)
        implicit none
        real(dp), intent(in) :: linear_part(ni,n_spinor), mod_part(ni,n_spinor)
        real(dp), intent(out) :: best_ratio, residual
        real(dp) :: denom, resnorm, modnorm

        denom = sum(linear_part**2)
        if (denom .lt. 1d-300) then
                best_ratio = 0d0
        else
                best_ratio = sum(linear_part * mod_part) / denom
        endif

        resnorm = sum((mod_part - best_ratio * linear_part)**2)
        modnorm = sum(mod_part**2)
        if (modnorm .lt. 1d-300) then
                residual = 0d0
        else
                residual = sqrt(resnorm / modnorm)
        endif
end subroutine lambda_like_diagnostic

subroutine make_path(ifile, ikind, out_path)
        implicit none
        integer, intent(in) :: ifile, ikind
        character(len=*), intent(out) :: out_path
        character(len=20) :: file_char, kind_char

        call int2char(ifile, file_char)
        call int2char(ikind, kind_char)
        out_path = './angle_'//trim(file_char)//'/V'//trim(kind_char)
end subroutine make_path

subroutine read_V(in_path, V_out, n_record)
        implicit none
        character(len=*), intent(in) :: in_path
        real(dp), intent(out) :: V_out(ni,4)
        integer, intent(out) :: n_record

        open(unit=104,file=trim(in_path),access='sequential',action='read',       &
                status='old',form='unformatted')
        read(104) n_record
        read(104) V_out
        close(104)
end subroutine read_V

subroutine write_report(report_path, filename, gx, gy, gz, max_iter_in,           &
                nonconv_in, l1_in, l2_in, a1_in, a2_in, tol_in, target, scalar,   &
                err_lin, err_full,                                                &
                gain_err, gain_rms, comp_target_in, comp_lin_in, comp_full_in,    &
                corr0, corr1, mm_corr, proj, raw_ratio, raw_resid, ortho_ratio,   &
                ortho_resid, raw_lin_norm, raw_mod_norm, raw_lin_max,             &
                raw_mod_max, ortho_lin_norm, ortho_mod_norm, ortho_lin_max,       &
                ortho_mod_max, lin_eff, mod_eff, cross_eff, raw_lin_eff,          &
                raw_mod_eff, raw_cross_eff, weight_path, n_blocked)
        implicit none
        character(len=*), intent(in) :: report_path, filename, weight_path
        integer, intent(in) :: gx, gy, gz, max_iter_in, nonconv_in
        integer, intent(in) :: n_blocked
        real(dp), intent(in) :: l1_in, l2_in, a1_in, a2_in, tol_in
        real(dp), intent(in) :: target, scalar, err_lin, err_full
        real(dp), intent(in) :: gain_err, gain_rms
        real(dp), intent(in) :: comp_target_in(3), comp_lin_in(3), comp_full_in(3)
        real(dp), intent(in) :: corr0, corr1, mm_corr
        real(dp), intent(in) :: proj(n_spinor,n_spinor)
        real(dp), intent(in) :: raw_ratio, raw_resid, ortho_ratio, ortho_resid
        real(dp), intent(in) :: raw_lin_norm(6), raw_mod_norm(6)
        real(dp), intent(in) :: raw_lin_max(6), raw_mod_max(6)
        real(dp), intent(in) :: ortho_lin_norm(6), ortho_mod_norm(6)
        real(dp), intent(in) :: ortho_lin_max(6), ortho_mod_max(6)
        real(dp), intent(in) :: lin_eff(6), mod_eff(6), cross_eff(6)
        real(dp), intent(in) :: raw_lin_eff(6), raw_mod_eff(6), raw_cross_eff(6)
        integer :: ia, ib, ic
        real(dp) :: eff_ratio, eff_corr

        open(unit=105,file=trim(report_path),action='write',                     &
                status='replace')
        write(105,'(A)') 'CrI3 orthogonal-basis K_mod diagnostic report'
        write(105,'(A)') 'Modulation: independent per-channel all-other average'
        write(105,'(A,A)') 'input_file = ', trim(filename)
        write(105,'(A,A)') 'Wannier weight file = ', trim(weight_path)
        write(105,'(A,A)') 'block mask mode = ', block_mask_mode
        write(105,'(A,I0)') 'blocked feature coefficients = ', n_blocked
        write(105,'(A,3I8)') 'grid nx ny nz = ', gx, gy, gz
        write(105,'(A,I0)') 'grid_points = ', ni
        write(105,'(A,I0)') 'file_max_iter = ', max_iter_in
        write(105,'(A,I0)') 'file_nonconverged = ', nonconv_in
        write(105,'(A,ES24.16)') 'file_l1_relative = ', l1_in
        write(105,'(A,ES24.16)') 'file_l2_relative = ', l2_in
        write(105,'(A,ES24.16)') 'file_alpha1 = ', a1_in
        write(105,'(A,ES24.16)') 'file_alpha2 = ', a2_in
        write(105,'(A,ES24.16)') 'file_tolerance = ', tol_in
        write(105,*)
        write(105,'(A)') 'Basis definition: L=S, M_ortho=chi*S-L*P.'
        write(105,'(A)') 'Stored raw K0/K1 are analyzed directly; K0 is refit after K1.'
        write(105,'(A,ES24.16)') 'max_linear_mod_corr_raw = ', corr0
        write(105,'(A,ES24.16)') 'max_linear_mod_corr_ortho = ', corr1
        write(105,'(A,ES24.16)') 'max_offdiag_mod_corr_ortho = ', mm_corr
        write(105,*)
        write(105,'(A,ES24.16)') 'weighted_target_magnetic_norm = ', sqrt(target)
        write(105,'(A,ES24.16)') 'weighted_scalar_norm = ', sqrt(scalar)
        write(105,'(A,ES24.16)') 'weighted_ortho_linear_relative_error = ',      &
                sqrt(err_lin / max(target,1d-300))
        write(105,'(A,ES24.16)') 'weighted_ortho_full_relative_error = ',        &
                sqrt(err_full / max(target,1d-300))
        write(105,'(A,ES24.16)') 'ortho_mod_gain_error_fraction = ', gain_err
        write(105,'(A,ES24.16)') 'ortho_mod_gain_rms_fraction = ', gain_rms
        write(105,*)
        write(105,'(A)') 'component ortho_linear_rel ortho_full_rel gain_error_fraction'
        do ic = 1, 3
                write(105,'(I4,3ES24.16)') ic,                                   &
                        sqrt(comp_lin_in(ic) / max(comp_target_in(ic),1d-300)),  &
                        sqrt(comp_full_in(ic) / max(comp_target_in(ic),1d-300)), &
                        1d0 - comp_full_in(ic) / max(comp_lin_in(ic),1d-300)
        enddo
        write(105,*)
        write(105,'(A,ES24.16)') 'raw_global_lambda_like_ratio = ', raw_ratio
        write(105,'(A,ES24.16)') 'raw_global_lambda_like_residual = ', raw_resid
        write(105,'(A,ES24.16)') 'ortho_global_lambda_like_ratio = ', ortho_ratio
        write(105,'(A,ES24.16)') 'ortho_global_lambda_like_residual = ',          &
                ortho_resid
        write(105,*)
        write(105,'(A)') 'projection P(row=linear channel, col=mod channel)'
        do ia = 1, n_spinor
                write(105,'(6ES18.9)') (proj(ia,ib), ib=1,n_spinor)
        enddo
        write(105,*)
        write(105,'(A)') 'channel raw_K0_rms raw_K1_rms raw_K1_over_K0 raw_K0_max raw_K1_max'
        do ia = 1, n_spinor
                write(105,'(I4,5ES24.16)') ia, raw_lin_norm(ia),                 &
                        raw_mod_norm(ia),                                        &
                        raw_mod_norm(ia) / max(raw_lin_norm(ia),1d-300),         &
                        raw_lin_max(ia), raw_mod_max(ia)
        enddo
        write(105,*)
        write(105,'(A)') 'channel ortho_K0_rms ortho_K1_rms ortho_K1_over_K0 ortho_K0_max ortho_K1_max'
        do ia = 1, n_spinor
                write(105,'(I4,5ES24.16)') ia, ortho_lin_norm(ia),               &
                        ortho_mod_norm(ia),                                      &
                        ortho_mod_norm(ia) / max(ortho_lin_norm(ia),1d-300),     &
                        ortho_lin_max(ia), ortho_mod_max(ia)
        enddo
        write(105,*)
        write(105,'(A)') 'channel ortho_linear_effect_rms ortho_mod_effect_rms mod_over_linear effect_corr'
        do ia = 1, n_spinor
                eff_ratio = sqrt(mod_eff(ia) / max(lin_eff(ia),1d-300))
                eff_corr = cross_eff(ia) / sqrt(max(lin_eff(ia)*mod_eff(ia),     &
                        1d-300))
                write(105,'(I4,4ES24.16)') ia, sqrt(lin_eff(ia)),                &
                        sqrt(mod_eff(ia)), eff_ratio, eff_corr
        enddo
        write(105,*)
        write(105,'(A)') 'channel raw_linear_effect_rms raw_mod_effect_rms mod_over_linear effect_corr'
        do ia = 1, n_spinor
                eff_ratio = sqrt(raw_mod_eff(ia) / max(raw_lin_eff(ia),1d-300))
                eff_corr = raw_cross_eff(ia) / sqrt(max(raw_lin_eff(ia)          &
                        * raw_mod_eff(ia), 1d-300))
                write(105,'(I4,4ES24.16)') ia, sqrt(raw_lin_eff(ia)),            &
                        sqrt(raw_mod_eff(ia)), eff_ratio, eff_corr
        enddo
        write(105,*)
        write(105,'(A)') 'Decision guide:'
        write(105,'(A)') '1. Use ortho_* quantities for fit_com_ortho files.'
        write(105,'(A)') '2. Raw K0/K1 can look large because raw terms cancel.'
        write(105,'(A)') '3. If ortho_mod_gain is large, residual really needs the orthogonal high-order basis.'
        write(105,'(A)') '4. If ortho_K1 is large but ortho_mod_gain is small, suspect residual overfit.'
        close(105)
end subroutine write_report

subroutine int2char(in_int, out_char)
        implicit none
        integer, intent(in) :: in_int
        character(len=*), intent(out) :: out_char

        if (in_int .lt. 10) write(out_char,'(I1)') in_int
        if (in_int .lt. 100 .and. in_int .ge. 10) write(out_char,'(I2)') in_int
        if (in_int .lt. 1000 .and. in_int .ge. 100) write(out_char,'(I3)') in_int
        if (in_int .lt. 10000 .and. in_int .ge. 1000) write(out_char,'(I4)') in_int
        if (in_int .lt. 100000 .and. in_int .ge. 10000) write(out_char,'(I5)') in_int
end subroutine int2char

end program fit_analysis_ortho
