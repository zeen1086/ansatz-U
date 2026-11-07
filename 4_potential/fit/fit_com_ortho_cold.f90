program fitting_compressive_ortho
  use parameter_mod
        use wannier_weight_mod, only: load_wannier_weight
        use cri3_block_mask_mod, only: block_mask_mode, make_feature_block_mask, &
                apply_feature_block_mask, count_blocked_features,               &
                make_feature_active_fraction, masked_absmax_rhs
        use cri3_c3_basis_mod, only: c3_basis_mode, c3_reduced_basis,            &
                build_c3_reduced_basis, solve_c3_reduced_cd, c3_basis_residual

        !
        ! Raw features:
        !
        ! Before fitting K_mod, remove from M its projection onto the linear
        ! subspace:
        !
        ! K_mod is fitted in the orthogonal basis.  Final raw K_linear is then
        ! conditionally refitted.  Both written kernels are expanded from the
        ! finite-window C3v(A/B C3 plus xy mirror) irreducible orbit basis.

        implicit none

        integer, parameter :: dp = kind(0d0)
        integer, parameter :: n_spinor = p_fit_nspinor
        integer, parameter :: n_file = p_n_file
        integer, parameter :: n_kind = p_n_kind
        integer, parameter :: n_sample = n_file * n_kind
        integer, parameter :: n_comp = 3
        integer, parameter :: n_feat = 2 * n_spinor
        character(len=32), parameter :: magic = 'CRI3_COM_FIT_CHANNEL_V1'

        integer :: nx, ny, nz, ni
        integer :: i_sample
        integer :: max_iter, linear_iter, raw_linear_iter, mod_iter
        integer :: linear_nonconv, mod_nonconv, nonconv_total
        integer :: blocked_features
        real(dp) :: l1_relative, l2_relative, tol
        real(dp) :: alpha1_linear, alpha1_mod, alpha2
        real(dp) :: err_comp(3), norm_comp(3), err_total, norm_total
        real(dp) :: lambda_ratio, lambda_residual
        real(dp), allocatable :: spin(:,:,:)
        real(dp), allocatable :: feat_raw(:,:,:)
        real(dp), allocatable :: feat_ortho(:,:,:)
        real(dp), allocatable :: feat_fit(:,:,:)
        real(dp), allocatable :: scale(:)
        real(dp), allocatable :: gram(:,:), rhs(:,:)
        real(dp), allocatable :: coef_ortho(:,:), coef_raw(:,:)
        real(dp), allocatable :: V_read(:,:), V_pred(:,:)
        real(dp), allocatable :: grid_weight(:)
        real(dp), allocatable :: block_mask(:,:)
        real(dp), allocatable :: active_fraction(:)
        real(dp) :: projection(n_spinor,n_spinor)
        real(dp) :: corr_before, corr_after
        real(dp) :: penalty(n_spinor)
        real(dp) :: c3_linear_max, c3_linear_rms, c3_mod_max, c3_mod_rms
        character(len=200) :: output_file, weight_file
        logical :: linear_converged, raw_linear_converged, mod_converged
        type(c3_reduced_basis) :: c3_basis

        l1_relative = p_fit_l1_relative
        l2_relative = p_fit_l2_relative
        max_iter = p_fit_max_iter
        tol = p_fit_tolerance
        output_file = 'V_func_com_ortho_channel_v1'
        weight_file = 'w_density'
        call read_grid(nx, ny, nz, ni)

        allocate(spin(n_sample,n_spinor,n_comp))
        allocate(feat_raw(n_sample,n_feat,n_comp))
        allocate(feat_ortho(n_sample,n_feat,n_comp))
        allocate(feat_fit(n_sample,n_feat,n_comp))
        allocate(scale(n_feat))
        allocate(gram(n_feat,n_feat))
        allocate(rhs(n_feat,ni))
        allocate(coef_ortho(ni,n_feat))
        allocate(coef_raw(ni,n_feat))
        allocate(V_read(ni,4))
        allocate(V_pred(ni,4))
        allocate(grid_weight(ni))
        allocate(block_mask(ni,n_feat))
        allocate(active_fraction(n_feat))
        call load_wannier_weight(trim(weight_file), nx, ny, nz, ni, grid_weight)
        call make_feature_block_mask(nx, ny, nz, ni, block_mask)
        call make_feature_active_fraction(block_mask, active_fraction)
        blocked_features = count_blocked_features(block_mask)
        if (maxval(abs(block_mask(:,1:n_spinor)                              &
                - block_mask(:,n_spinor+1:n_feat))) > 1d-12) then
                stop 'linear/mod block masks differ under shared C3v basis'
        endif
        call build_c3_reduced_basis(nx, ny, nz, ni,                           &
                block_mask(:,1:n_spinor), c3_basis)

        call read_spins(spin)
        do i_sample = 1, n_sample
                call make_features(spin(i_sample,:,:), feat_raw(i_sample,:,:))
        enddo

        call build_projection(feat_raw, projection, corr_before)
        call make_ortho_features(feat_raw, projection, feat_ortho)
        call max_linear_mod_corr(feat_ortho, corr_after)

        feat_fit(:,:,:) = feat_ortho(:,:,:)
        call normalize_features(feat_fit, scale)
        call make_gram(feat_fit, gram)
        call make_rhs(feat_fit, rhs)

        alpha1_linear = l1_relative                                             &
                * masked_absmax_rhs(rhs(1:n_spinor,:), block_mask(:,1:n_spinor))
        alpha2 = l2_relative

        coef_ortho(:,:) = 0d0
        linear_iter = 0
        mod_iter = 0
        linear_nonconv = 0
        mod_nonconv = 0
        penalty(:) = active_fraction(1:n_spinor)

        call solve_c3_reduced_cd(c3_basis, gram(1:n_spinor,1:n_spinor),       &
                rhs(1:n_spinor,:), penalty, alpha1_linear, alpha2, max_iter,  &
                tol, scale(1:n_spinor), grid_weight,                          &
                coef_ortho(:,1:n_spinor), linear_iter, linear_converged)
        if (.not. linear_converged) linear_nonconv = linear_nonconv + 1

        call make_mod_rhs(coef_ortho(:,1:n_spinor), feat_ortho, feat_fit,         &
                rhs(n_spinor+1:n_feat,:))
        alpha1_mod = l1_relative                                                  &
                * masked_absmax_rhs(rhs(n_spinor+1:n_feat,:),                    &
                block_mask(:,n_spinor+1:n_feat))
        penalty(:) = active_fraction(n_spinor+1:n_feat)

        call solve_c3_reduced_cd(c3_basis,                                    &
                gram(n_spinor+1:n_feat,n_spinor+1:n_feat),                    &
                rhs(n_spinor+1:n_feat,:), penalty, alpha1_mod, alpha2,         &
                max_iter, tol, scale(n_spinor+1:n_feat), grid_weight,         &
                coef_ortho(:,n_spinor+1:n_feat), mod_iter, mod_converged)
        if (.not. mod_converged) mod_nonconv = mod_nonconv + 1

        ! K_mod is unchanged by the orthogonal-to-raw transformation.  Refit
        ! final raw K_linear conditionally on K_mod so the written coefficients,
        ! rather than only the intermediate orthogonal coefficients, obey C3v.
        coef_raw(:,:) = 0d0
        coef_raw(:,n_spinor+1:n_feat) = coef_ortho(:,n_spinor+1:n_feat)
        feat_fit(:,:,:) = feat_raw(:,:,:)
        call normalize_features(feat_fit, scale)
        call make_gram(feat_fit, gram)
        call make_raw_linear_rhs(coef_raw(:,n_spinor+1:n_feat), feat_raw,     &
                feat_fit, rhs(1:n_spinor,:))
        alpha1_linear = l1_relative                                           &
                * masked_absmax_rhs(rhs(1:n_spinor,:),                        &
                block_mask(:,1:n_spinor))
        penalty(:) = active_fraction(1:n_spinor)
        call solve_c3_reduced_cd(c3_basis, gram(1:n_spinor,1:n_spinor),       &
                rhs(1:n_spinor,:), penalty, alpha1_linear, alpha2, max_iter,  &
                tol, scale(1:n_spinor), grid_weight, coef_raw(:,1:n_spinor),  &
                raw_linear_iter, raw_linear_converged)
        linear_iter = max(linear_iter, raw_linear_iter)
        if (.not. raw_linear_converged) linear_nonconv = linear_nonconv + 1

        call apply_feature_block_mask(coef_raw, block_mask)
        call raw_to_ortho_coefficients(coef_raw, projection, coef_ortho)
        call apply_feature_block_mask(coef_ortho, block_mask)
        call c3_basis_residual(c3_basis, coef_raw(:,1:n_spinor),              &
                c3_linear_max, c3_linear_rms)
        call c3_basis_residual(c3_basis, coef_raw(:,n_spinor+1:n_feat),       &
                c3_mod_max, c3_mod_rms)
        call training_error(coef_raw, feat_raw, grid_weight, err_comp, norm_comp, &
                err_total, norm_total)
        call lambda_like(coef_raw(:,1:n_spinor), coef_raw(:,n_spinor+1:n_feat),   &
                lambda_ratio, lambda_residual)

        nonconv_total = linear_nonconv + mod_nonconv
        call write_V_func_com(trim(output_file), nx, ny, nz, ni, l1_relative,     &
                l2_relative, alpha1_linear, alpha2, tol, max_iter, nonconv_total, &
                coef_raw)
        call write_fit_sample(coef_raw, feat_raw)
        call write_report(trim(output_file), projection, corr_before, corr_after, &
                l1_relative, l2_relative, alpha1_linear, alpha1_mod, alpha2, tol, &
                linear_iter, mod_iter, linear_nonconv, mod_nonconv, err_comp,     &
                norm_comp, err_total, norm_total, lambda_ratio, lambda_residual,  &
                coef_ortho, coef_raw, weight_file, blocked_features,              &
                c3_basis%n_active, c3_basis%n_reduced, c3_linear_max,             &
                c3_linear_rms, c3_mod_max, c3_mod_rms)

        print *, 'fit_com_ortho done'
        print *, 'output = ', trim(output_file)
        print *, 'report = fit_com_ortho_report.dat'
        print *, 'block mask mode = ', block_mask_mode
        print *, 'C3v basis mode = ', c3_basis_mode
        print *, 'C3v active/reduced = ', c3_basis%n_active, c3_basis%n_reduced
        print *, 'C3v residual max K0/K1 = ', c3_linear_max, c3_mod_max
        print *, 'blocked feature coefficients = ', blocked_features
        print *, 'corr before/after = ', corr_before, corr_after
        print *, 'relative error = ', sqrt(err_total / norm_total)

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
                open(unit=102,file=trim(spin_path),access='sequential',          &
                        action='read',status='old')
                do ikind = 1, n_kind
                        isample = isample + 1
                        read(102,*) flat(:)
                        do ispin = 1, n_spinor
                                all_spin(isample,ispin,:) =                     &
                                        flat(3*(ispin-1)+1:3*ispin)
                        enddo
                enddo
                close(102)
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

subroutine normalize_features(feature_inout, col_scale)
        implicit none
        real(dp), intent(inout) :: feature_inout(n_sample,n_feat,n_comp)
        real(dp), intent(out) :: col_scale(n_feat)
        integer :: jf

        do jf = 1, n_feat
                col_scale(jf) = sqrt(sum(feature_inout(:,jf,:)**2))
                if (col_scale(jf) .lt. 1d-14) stop 'zero feature column'
                feature_inout(:,jf,:) = feature_inout(:,jf,:) / col_scale(jf)
        enddo
end subroutine normalize_features

subroutine make_gram(feature_in, gram_out)
        implicit none
        real(dp), intent(in) :: feature_in(n_sample,n_feat,n_comp)
        real(dp), intent(out) :: gram_out(n_feat,n_feat)
        integer :: is, jf, kf

        gram_out(:,:) = 0d0
        do is = 1, n_sample
        do jf = 1, n_feat
        do kf = 1, n_feat
                gram_out(jf,kf) = gram_out(jf,kf)                                &
                        + dot_product(feature_in(is,jf,:), feature_in(is,kf,:))
        enddo
        enddo
        enddo
end subroutine make_gram

subroutine make_rhs(feature_in, rhs_out)
        implicit none
        real(dp), intent(in) :: feature_in(n_sample,n_feat,n_comp)
        real(dp), intent(out) :: rhs_out(n_feat,ni)
        integer :: ifile, ikind, isample, jf, nrec
        character(len=200) :: v_path

        rhs_out(:,:) = 0d0
        isample = 0
        do ifile = 1, n_file
        do ikind = 1, n_kind
                isample = isample + 1
                call make_path(ifile, ikind, v_path)
                call read_V(v_path, V_read, nrec)
                do jf = 1, n_feat
                        rhs_out(jf,:) = rhs_out(jf,:)                            &
                                + feature_in(isample,jf,1) * V_read(:,2)         &
                                + feature_in(isample,jf,2) * V_read(:,3)         &
                                + feature_in(isample,jf,3) * V_read(:,4)
                enddo
        enddo
        enddo
end subroutine make_rhs

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

        open(unit=103,file=trim(in_path),access='sequential',action='read',       &
                status='old',form='unformatted')
        read(103) n_record
        read(103) V_out
        close(103)
end subroutine read_V

subroutine make_mod_rhs(linear_part, feature_raw, feature_norm, rhs_mod)
        implicit none
        real(dp), intent(in) :: linear_part(ni,n_spinor)
        real(dp), intent(in) :: feature_raw(n_sample,n_feat,n_comp)
        real(dp), intent(in) :: feature_norm(n_sample,n_feat,n_comp)
        real(dp), intent(out) :: rhs_mod(n_spinor,ni)
        real(dp) :: residual(ni,3)
        integer :: ifile, ikind, isample, jf, ic, nrec
        character(len=200) :: v_path

        rhs_mod(:,:) = 0d0
        isample = 0
        do ifile = 1, n_file
        do ikind = 1, n_kind
                isample = isample + 1
                call make_path(ifile, ikind, v_path)
                call read_V(v_path, V_read, nrec)
                do ic = 1, 3
                        residual(:,ic) = V_read(:,ic+1)
                        do jf = 1, n_spinor
                                residual(:,ic) = residual(:,ic)                 &
                                        - linear_part(:,jf)                      &
                                        * feature_raw(isample,jf,ic)
                        enddo
                enddo
                do jf = 1, n_spinor
                do ic = 1, 3
                        rhs_mod(jf,:) = rhs_mod(jf,:)                           &
                                + feature_norm(isample,n_spinor+jf,ic)          &
                                * residual(:,ic)
                enddo
                enddo
        enddo
        enddo
end subroutine make_mod_rhs

subroutine make_raw_linear_rhs(mod_part, raw_feature, normalized_feature,       &
        rhs_linear)
        implicit none
        real(dp), intent(in) :: mod_part(ni,n_spinor)
        real(dp), intent(in) :: raw_feature(n_sample,n_feat,n_comp)
        real(dp), intent(in) :: normalized_feature(n_sample,n_feat,n_comp)
        real(dp), intent(out) :: rhs_linear(n_spinor,ni)
        real(dp) :: residual(ni,3)
        integer :: ifile, ikind, isample, jf, ic, nrec
        character(len=200) :: v_path

        rhs_linear(:,:) = 0d0
        isample = 0
        do ifile = 1, n_file
        do ikind = 1, n_kind
                isample = isample + 1
                call make_path(ifile, ikind, v_path)
                call read_V(v_path, V_read, nrec)
                do ic = 1, n_comp
                        residual(:,ic) = V_read(:,ic+1)
                        do jf = 1, n_spinor
                                residual(:,ic) = residual(:,ic)                 &
                                        - mod_part(:,jf)                        &
                                        * raw_feature(isample,n_spinor+jf,ic)
                        enddo
                enddo
                do jf = 1, n_spinor
                do ic = 1, n_comp
                        rhs_linear(jf,:) = rhs_linear(jf,:)                     &
                                + normalized_feature(isample,jf,ic)             &
                                * residual(:,ic)
                enddo
                enddo
        enddo
        enddo
end subroutine make_raw_linear_rhs

subroutine raw_to_ortho_coefficients(coef_in, proj, coef_out)
        implicit none
        real(dp), intent(in) :: coef_in(ni,n_feat)
        real(dp), intent(in) :: proj(n_spinor,n_spinor)
        real(dp), intent(out) :: coef_out(ni,n_feat)
        integer :: ip, ia, ib

        coef_out(:,:) = coef_in(:,:)
        do ip = 1, ni
        do ib = 1, n_spinor
                coef_out(ip,ib) = coef_in(ip,ib)
                do ia = 1, n_spinor
                        coef_out(ip,ib) = coef_out(ip,ib)                        &
                                + coef_in(ip,n_spinor+ia) * proj(ib,ia)
                enddo
        enddo
        enddo
end subroutine raw_to_ortho_coefficients

subroutine training_error(fit_coef, raw_feature, weight, e_comp, n_comp_out,  &
        e_tot, n_tot)
        implicit none
        real(dp), intent(in) :: fit_coef(ni,n_feat)
        real(dp), intent(in) :: raw_feature(n_sample,n_feat,n_comp)
        real(dp), intent(in) :: weight(ni)
        real(dp), intent(out) :: e_comp(3), n_comp_out(3), e_tot, n_tot
        real(dp) :: pred(ni,3)
        integer :: ifile, ikind, isample, jf, ic, nrec
        character(len=200) :: v_path

        e_comp(:) = 0d0
        n_comp_out(:) = 0d0
        isample = 0
        do ifile = 1, n_file
        do ikind = 1, n_kind
                isample = isample + 1
                call make_path(ifile, ikind, v_path)
                call read_V(v_path, V_read, nrec)
                pred(:,:) = 0d0
                do jf = 1, n_feat
                do ic = 1, 3
                        pred(:,ic) = pred(:,ic)                                  &
                                + fit_coef(:,jf) * raw_feature(isample,jf,ic)
                enddo
                enddo
                do ic = 1, 3
                        e_comp(ic) = e_comp(ic)                                  &
                                + sum(weight*(pred(:,ic)-V_read(:,ic+1))**2)
                        n_comp_out(ic) = n_comp_out(ic)                          &
                                + sum(weight*V_read(:,ic+1)**2)
                enddo
        enddo
        enddo
        e_tot = sum(e_comp)
        n_tot = sum(n_comp_out)
end subroutine training_error

subroutine lambda_like(linear_part, mod_part, best_ratio, residual)
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
        residual = sqrt(resnorm / max(modnorm, 1d-300))
end subroutine lambda_like

subroutine write_V_func_com(file_out, gx, gy, gz, grid_size, l1rel, l2rel,        &
                l1abs, l2abs, tolerance, n_iter, nonconv, fit_coef)
        implicit none
        character(len=*), intent(in) :: file_out
        integer, intent(in) :: gx, gy, gz, grid_size, n_iter, nonconv
        real(dp), intent(in) :: l1rel, l2rel, l1abs, l2abs, tolerance
        real(dp), intent(in) :: fit_coef(grid_size,n_feat)

        open(unit=104,file=trim(file_out),access='sequential',action='write',     &
                status='replace',form='unformatted')
        write(104) magic
        write(104) gx, gy, gz
        write(104) grid_size, n_spinor, n_feat, n_file, n_kind
        write(104) n_iter, nonconv
        write(104) l1rel, l2rel, l1abs, l2abs, tolerance
        write(104) fit_coef(:,1:n_spinor)
        write(104) fit_coef(:,n_spinor+1:n_feat)
        close(104)
end subroutine write_V_func_com

subroutine write_fit_sample(fit_coef, raw_feature)
        implicit none
        real(dp), intent(in) :: fit_coef(ni,n_feat)
        real(dp), intent(in) :: raw_feature(n_sample,n_feat,n_comp)
        integer :: isample, jf, ic

        isample = (min(10,n_file) - 1) * n_kind + 1
        V_pred(:,:) = 0d0
        do jf = 1, n_feat
        do ic = 1, 3
                V_pred(:,ic+1) = V_pred(:,ic+1)                                  &
                        + fit_coef(:,jf) * raw_feature(isample,jf,ic)
        enddo
        enddo
        open(unit=105,file='V_fit_com_ortho',access='sequential',action='write',  &
                status='replace',form='unformatted')
        write(105) ni
        write(105) V_pred
        close(105)
end subroutine write_fit_sample

subroutine write_report(file_out, proj, corr0, corr1, l1rel, l2rel, a1lin,        &
                a1mod, a2, tolerance, lin_iter, high_iter, lin_nonconv,           &
                high_nonconv, e_comp, n_comp_out, e_tot, n_tot, best_ratio,       &
                lambda_resid, coef_o, coef_r, weight_path, n_blocked, c3_active,  &
                c3_reduced, c3_lin_max, c3_lin_rms, c3_high_max, c3_high_rms)
        implicit none
        character(len=*), intent(in) :: file_out
        character(len=*), intent(in) :: weight_path
        real(dp), intent(in) :: proj(n_spinor,n_spinor), corr0, corr1
        real(dp), intent(in) :: l1rel, l2rel, a1lin, a1mod, a2, tolerance
        integer, intent(in) :: lin_iter, high_iter, lin_nonconv, high_nonconv
        real(dp), intent(in) :: e_comp(3), n_comp_out(3), e_tot, n_tot
        real(dp), intent(in) :: best_ratio, lambda_resid
        real(dp), intent(in) :: coef_o(ni,n_feat), coef_r(ni,n_feat)
        integer, intent(in) :: n_blocked, c3_active, c3_reduced
        real(dp), intent(in) :: c3_lin_max, c3_lin_rms, c3_high_max, c3_high_rms
        integer :: ia, ib
        real(dp) :: lin_rms, mod_rms, raw_lin_rms, raw_mod_rms

        open(unit=106,file='fit_com_ortho_report.dat',action='write',             &
                status='replace')
        write(106,'(A)') 'CrI3 orthogonalized SU(2) compressive fit'
        write(106,'(A)') 'Modulation: independent per-channel all-other average'
        write(106,'(A,A)') 'output_file = ', trim(file_out)
        write(106,'(A,A)') 'Wannier weight file = ', trim(weight_path)
        write(106,'(A,A)') 'block mask mode = ', block_mask_mode
        write(106,'(A,I0)') 'blocked feature coefficients = ', n_blocked
        write(106,'(A,A)') 'C3v reduced basis mode = ', c3_basis_mode
        write(106,'(A,I0)') 'C3v active coefficients per K block = ', c3_active
        write(106,'(A,I0)') 'C3v irreducible coefficients per K block = ',      &
                c3_reduced
        write(106,'(A,ES24.16)') 'C3v K0 max_abs_residual = ', c3_lin_max
        write(106,'(A,ES24.16)') 'C3v K0 rms_residual = ', c3_lin_rms
        write(106,'(A,ES24.16)') 'C3v K1 max_abs_residual = ', c3_high_max
        write(106,'(A,ES24.16)') 'C3v K1 rms_residual = ', c3_high_rms
        write(106,'(A)') 'Ridge penalty: active-fraction normalized'
        write(106,'(A)') 'Objective: Wannier-weighted data loss with Ridge-dominant regularization'
        write(106,'(A)') 'Fit basis: L=S, Mortho=Mraw-L*P, Mraw=chi*S'
        write(106,'(A)') 'K0 is conditionally refit in raw L after constrained K1.'
        write(106,'(A)') 'Written raw K0 and K1 both lie in the C3v reduced basis.'
        write(106,'(A,ES24.16)') 'l1_relative = ', l1rel
        write(106,'(A,ES24.16)') 'l2_relative = ', l2rel
        write(106,'(A,ES24.16)') 'alpha1_linear = ', a1lin
        write(106,'(A,ES24.16)') 'alpha1_mod = ', a1mod
        write(106,'(A,ES24.16)') 'alpha2 = ', a2
        write(106,'(A,ES24.16)') 'tolerance = ', tolerance
        write(106,'(A,I0)') 'linear_max_iterations_used = ', lin_iter
        write(106,'(A,I0)') 'mod_max_iterations_used = ', high_iter
        write(106,'(A,I0)') 'linear_nonconverged_reduced_solves = ', lin_nonconv
        write(106,'(A,I0)') 'mod_nonconverged_reduced_solves = ', high_nonconv
        write(106,'(A,ES24.16)') 'weighted_training_relative_fit_error_total = ', &
                sqrt(e_tot / n_tot)
        do ia = 1, 3
                write(106,'(A,I0,A,ES24.16)')                                    &
                        'weighted_training_relative_fit_error_component_', ia,  &
                        sqrt(e_comp(ia) / n_comp_out(ia))
        enddo
        write(106,'(A,ES24.16)') 'max_linear_mod_corr_before = ', corr0
        write(106,'(A,ES24.16)') 'max_linear_mod_corr_after = ', corr1
        write(106,'(A,ES24.16)') 'global_lambda_like_ratio_raw = ', best_ratio
        write(106,'(A,ES24.16)') 'global_lambda_like_residual_raw = ',            &
                lambda_resid
        write(106,*)
        write(106,'(A)') 'projection P(row=linear channel, col=mod channel)'
        do ia = 1, n_spinor
                write(106,'(6ES18.9)') (proj(ia,ib), ib=1,n_spinor)
        enddo
        write(106,*)
        write(106,'(A)') 'channel ortho_K0_rms ortho_K1_rms raw_K0_rms raw_K1_rms raw_K1_over_K0'
        do ia = 1, n_spinor
                lin_rms = sqrt(sum(coef_o(:,ia)**2) / dble(ni))
                mod_rms = sqrt(sum(coef_o(:,n_spinor+ia)**2) / dble(ni))
                raw_lin_rms = sqrt(sum(coef_r(:,ia)**2) / dble(ni))
                raw_mod_rms = sqrt(sum(coef_r(:,n_spinor+ia)**2) / dble(ni))
                write(106,'(I4,5ES24.16)') ia, lin_rms, mod_rms, raw_lin_rms,    &
                        raw_mod_rms, raw_mod_rms / max(raw_lin_rms,1d-300)
        enddo
        close(106)
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

end program fitting_compressive_ortho

