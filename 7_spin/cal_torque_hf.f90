program direct_torque_hf
  use parameter_mod
  use spin_channels_mod, only: spin_channels, spin_chirality
    use omp_lib
    implicit none

    integer, parameter :: dp = kind(0d0)
    real(dp), parameter :: hartree_to_ev = 13.6d0

    complex(dp), allocatable :: H_wan(:,:,:), H_k(:,:)
    complex(dp), allocatable :: dH_wan(:,:,:,:), dH_k(:,:)
    complex(dp), allocatable :: unmk(:,:), H_wan_tmp(:,:), dH_wan_tmp(:,:)
    complex(dp), allocatable :: H_step1(:,:), dH_step1(:,:), dH_psi(:)
    complex(dp), allocatable :: phase0(:,:), phase1(:,:), phase_k(:)
    complex(dp), allocatable :: RVR(:,:,:), RVR_tmp(:,:)
    complex(dp), allocatable :: work(:)

    real(dp), allocatable :: DOS_E_list(:), eta(:), torque(:,:,:), torque_local(:,:,:)
    real(dp), allocatable :: k_list(:,:), eig_tmp(:), tau_eig(:,:)
    real(dp), allocatable :: nb(:,:), R_list(:,:), spin_vec(:,:,:)
    real(dp), allocatable :: dspin_vec(:,:,:,:), physical_spin(:,:), dphysical_spin(:,:)
    real(dp), allocatable :: rwork(:)

    integer, allocatable :: ir_vec(:,:), super_ir_vec(:,:), pool(:,:), cell(:,:)
    integer, allocatable :: dir_vec(:,:), inv_dir(:), iwork(:)
    integer, allocatable :: response_sites(:)

    integer :: ik, jk, nk, nkx, nky, nkz, irept, nrept, io, read3_int(3)
    integer :: i, j, ix, iy, iz, nwan_cell, nwan_rvr, ikxy, H_dim
    integer :: nspin, period_x, period_y, period_z, super_size
    integer :: c_ix, c_iy, c_iz, c_jx, c_jy, c_jz, n1, n2
    integer :: i_dim, j_dim, iE, i1, i2, i3, i4, f1, f2, f3, f4
    integer :: n_DOS, nwan, n_func, convert_size_x, convert_size_y, convert_size_z
    integer :: grid_x, grid_y, grid_z, spin_grid_x, spin_grid_y, spin_grid_z
    integer :: n_eta, i_eta, npool, ipool, n_input, n_response, n_dep, isite, ib
    integer :: response_index, selected_eta
    integer :: wan_min, wan_max, d_ext, i_dir, j_dir, k_dir
    integer :: ini_x, ini_y, ini_z, fin_x, fin_y, fin_z
    integer :: super_x, super_y, super_z, i_locat, j_locat, jrept
    integer :: jx, jy, jz
    integer :: iRVR, jRVR, iwan, jwan, is, js
    integer :: max_x, max_y, max_z, min_x, min_y, min_z
    integer :: lwork, lrwork, liwork, info, env_status

    real(dp) :: a0, a(3,3), Ef, fac, celldm3, axial
    real(dp) :: pi, read1, E1, width, inv_eta, delta_E_val, fermi_func, term_exp
    real(dp) :: t_start, t_end, t_diag_start, t_diag_end, time_diag

    complex(dp) :: ctmp1, ctmp2, ctmp3, cI, I2pi

    character(len=20) :: irept_char, eta_char, site_char
    character(len=200) :: useless_string, func_dir, nscf_dir, rvr_root, RVR_dir
    character(len=200) :: unused_dir_override

    call read_global_parameters()
    call read_broadening()
    call read_torque_controls()
    call setup_dimensions()
    call read_response_sites()
    call load_ir_vectors()
    call setup_k_space()
    call read_total_hamiltonian()
    call read_spin_derivatives()
    call read_RVR()
    call build_dH_wan()
    call compute_torque()
    call write_direct_torque()

contains

    subroutine read_global_parameters()
func_dir = p_func_dir
nscf_dir = p_nscf_dir
rvr_root = p_rvr_dir
n_func = p_n_func
convert_size_x = p_convert_size_x
convert_size_y = p_convert_size_y
convert_size_z = p_convert_size_z
nwan = p_nwan
nspin = p_nspin
grid_x = p_nx
grid_y = p_ny
grid_z = p_nz
axial = p_axial
period_x = p_period_x
period_y = p_period_y
period_z = p_period_z
fac = p_fac
wan_min = p_wan_min
wan_max = p_wan_max
a0 = p_a0
celldm3 = p_celldm3
Ef = p_ef

        if (len_trim(rvr_root) .gt. 0) then
            if (rvr_root(len_trim(rvr_root):len_trim(rvr_root)) .ne. '/') then
                rvr_root=trim(rvr_root)//'/'
            endif
        endif
    end subroutine read_global_parameters

    subroutine read_broadening()
        n_eta = ubound(p_broadening,1)
        allocate(eta(n_eta))
        eta = p_broadening
    end subroutine read_broadening

    subroutine read_torque_controls()
        selected_eta = p_selected_eta
    end subroutine read_torque_controls

    subroutine setup_dimensions()
        pi = 3.1415926535d0
        cI = dcmplx(0.0d0, 1.0d0)
        I2pi = 2.0d0*pi*cI

        nwan_cell = nwan*2/nspin
        if (nspin .eq. 1) then
            nwan_rvr = 2*(wan_max-wan_min+1)
        else
            nwan_rvr = wan_max-wan_min+1
        endif

        H_dim = nwan_cell*period_x*period_y*period_z
        super_size = period_x*period_y*period_z
        d_ext = (2*convert_size_x+1)*(2*convert_size_y+1)*(2*convert_size_z+1)
        n_dep = 6
        n_input = 2*period_x*period_y*period_z

        if (n_func .ne. n_dep*3 .and. n_func .ne. 2*n_dep*3) then
            print *, 'Unexpected n_func for six-channel spin ansatz:', n_func
            stop
        endif
    end subroutine setup_dimensions

    subroutine read_response_sites()
        logical :: file_exists
        integer :: site_zero

        inquire(file='response_sites.in', exist=file_exists)
        if (file_exists) then
            open(unit=101,file='response_sites.in',status='old',action='read')
            read(101,*) n_response
            allocate(response_sites(n_response))
            read(101,*) response_sites(:)
            close(101)

            response_sites(:)=response_sites(:)+1
            do i=1,n_response
                if (response_sites(i) .lt. 1 .or. response_sites(i) .gt. n_input) then
                    site_zero=response_sites(i)-1
                    print *, 'response site is outside spin.in:', site_zero
                    stop
                endif
            enddo
        else
            n_response=n_input
            allocate(response_sites(n_response))
            do i=1,n_response
                response_sites(i)=i
            enddo
        endif
    end subroutine read_response_sites

    subroutine load_ir_vectors()
        irept=0
        open(unit=101,file=trim(p_super_h_dir)//'ir_vec_HYC',action='read',access='sequential',status='old')
        do
            read(101,'(3I8)',iostat=io) read3_int(:)
            if (io/=0) exit
            irept=irept+1
        enddo
        close(101)

        nrept=irept
        allocate(super_ir_vec(nrept,3))
        allocate(ir_vec(nrept,3))

        irept=1
        open(unit=101,file=trim(p_super_h_dir)//'ir_vec_HYC',action='read',access='sequential',status='old')
        do
            read(101,'(3I8)',iostat=io) read3_int(:)
            if (io/=0) exit
            super_ir_vec(irept,:)=read3_int(:)
            ir_vec(irept,:)=read3_int(:)
            irept=irept+1
        enddo
        close(101)

        ir_vec(:,1)=ir_vec(:,1)*period_x
        ir_vec(:,2)=ir_vec(:,2)*period_y
        ir_vec(:,3)=ir_vec(:,3)*period_z
    end subroutine load_ir_vectors

    subroutine setup_k_space()
        if (.not. mod(grid_x,period_x) .eq. 0) stop 'not commensurate in x'
        if (.not. mod(grid_y,period_y) .eq. 0) stop 'not commensurate in y'
        if (.not. mod(grid_z,period_z) .eq. 0) stop 'not commensurate in z'

        nkx=grid_x/period_x
        nky=grid_y/period_y
        nkz=grid_z/period_z
        nk=nkx*nky*nkz

        a(1,:) = (/   1.000000d0,   0.000000d0,   0.000000d0 /)
        a(2,:) = (/  -0.500000d0,   0.866025d0,   0.000000d0 /)
        a(3,:) = (/   0.000000d0,   0.000000d0,   celldm3  /)
        a(:,:)=a(:,:)*a0*0.529177d0

        allocate(k_list(nk,3))
        npool=nk/nkz
        allocate(pool(npool,nkz))

        do ix=1,nkx
        do iy=1,nky
        do iz=1,nkz
            ik=1+(iz-1)+(iy-1)*nkz+(ix-1)*nky*nkz
            k_list(ik,1)=1.0d0*(ix-1)/(nkx*period_x)
            k_list(ik,2)=1.0d0*(iy-1)/(nky*period_y)
            k_list(ik,3)=1.0d0*(iz-1)/(nkz*period_z)
            ipool=1+(iy-1)+(ix-1)*nky
            pool(ipool,iz)=ipool
        enddo
        enddo
        enddo

        allocate(nb(super_size,3))
        allocate(phase1(super_size,super_size))
        allocate(R_list(super_size,3))

        do ix=0,period_x-1
        do iy=0,period_y-1
        do iz=0,period_z-1
            i=1+iz+period_z*iy+period_z*period_y*ix
            nb(i,:)=(/ix*1.0d0/period_x, iy*1.0d0/period_y, iz*1.0d0/period_z/)
            R_list(i,1)=sum(a(:,1)*(/ix,iy,iz/))
            R_list(i,2)=sum(a(:,2)*(/ix,iy,iz/))
            R_list(i,3)=sum(a(:,3)*(/ix,iy,iz/))
        enddo
        enddo
        enddo

        do n1=1,super_size
            do c_ix=0,period_x-1
            do c_iy=0,period_y-1
            do c_iz=0,period_z-1
                i=1+c_iz+period_z*c_iy+period_z*period_y*c_ix
                phase1(i,n1)=exp(-I2pi*dot_product(1.0d0*(/c_ix,c_iy,c_iz/),nb(n1,:)))
            enddo
            enddo
            enddo
        enddo

    end subroutine setup_k_space

    subroutine read_total_hamiltonian()
        allocate(H_wan(H_dim,H_dim,nrept))

        do irept=1,nrept
            call int2char(irept, irept_char)
            open(unit=101,file='./tot_HR/hr_'//trim(irept_char),action='read', &
                access='stream',form='unformatted',status='old')
            read(101) H_wan(:,:,irept)
            close(101)
        enddo
    end subroutine read_total_hamiltonian

    subroutine read_spin_derivatives()
        real(dp), allocatable :: theta_nu(:), phi_nu(:)
        real(dp) :: deg2rad, mag_moment
        integer :: map_idx(6)
        logical :: file_exists

        allocate(theta_nu(n_input))
        allocate(phi_nu(n_input))
        allocate(physical_spin(n_input,3))
        allocate(dphysical_spin(n_input,3))

        inquire(file='spin.in', exist=file_exists)
        if (.not. file_exists) then
            print *, "Error: 'spin.in' not found!"
            stop
        endif

        open(unit=201,file='spin.in',status='old',action='read')
        read(201,*) theta_nu(:)
        read(201,*) phi_nu(:)
        close(201)

        deg2rad = pi/180.0d0
        theta_nu(:)=theta_nu(:)*deg2rad
        phi_nu(:)=phi_nu(:)*deg2rad
        mag_moment = p_cr_magnetization

        do isite=1,n_input
            physical_spin(isite,:) = (/ mag_moment*sin(theta_nu(isite))*cos(phi_nu(isite)), &
                                        mag_moment*sin(theta_nu(isite))*sin(phi_nu(isite)), &
                                        mag_moment*cos(theta_nu(isite)) /)
            dphysical_spin(isite,:) = (/ mag_moment*cos(theta_nu(isite))*cos(phi_nu(isite)), &
                                         mag_moment*cos(theta_nu(isite))*sin(phi_nu(isite)), &
                                        -mag_moment*sin(theta_nu(isite)) /)
        enddo

        max_x=-99; max_y=-99; max_z=-99
        min_x=99;  min_y=99;  min_z=99

        do irept=1,nrept
            if (super_ir_vec(irept,1) .gt. max_x) max_x=super_ir_vec(irept,1)
            if (super_ir_vec(irept,2) .gt. max_y) max_y=super_ir_vec(irept,2)
            if (super_ir_vec(irept,3) .gt. max_z) max_z=super_ir_vec(irept,3)
            if (super_ir_vec(irept,1) .lt. min_x) min_x=super_ir_vec(irept,1)
            if (super_ir_vec(irept,2) .lt. min_y) min_y=super_ir_vec(irept,2)
            if (super_ir_vec(irept,3) .lt. min_z) min_z=super_ir_vec(irept,3)
        enddo

        spin_grid_x=(max_x-min_x+1)*period_x
        spin_grid_y=(max_y-min_y+1)*period_y
        spin_grid_z=(max_z-min_z+1)*period_z

        allocate(cell(spin_grid_x*spin_grid_y*spin_grid_z,3))
        allocate(spin_vec(spin_grid_x*spin_grid_y*spin_grid_z,3,n_dep))
        allocate(dspin_vec(spin_grid_x*spin_grid_y*spin_grid_z,3,n_dep,n_input))
        spin_vec(:,:,:)=0.0d0
        dspin_vec(:,:,:,:)=0.0d0

        do super_x=min_x,max_x
        do super_y=min_y,max_y
        do super_z=min_z,max_z
        do ix=0,period_x-1
        do iy=0,period_y-1
        do iz=0,period_z-1
            i=1+(iz+(super_z-min_z)*period_z)+spin_grid_z*(iy+(super_y-min_y)*period_y) &
              +spin_grid_z*spin_grid_y*(ix+(super_x-min_x)*period_x)

            cell(i,1)=ix+super_x*period_x
            cell(i,2)=iy+super_y*period_y
            cell(i,3)=iz+super_z*period_z

            map_idx(1)=spin_input_index(ix,       iy,       iz,2)
            map_idx(2)=spin_input_index(ix,       iy,       iz,1)
            map_idx(3)=spin_input_index(ix,       iy+1,     iz,1)
            map_idx(4)=spin_input_index(ix+1,     iy,       iz,1)
            map_idx(5)=spin_input_index(ix+1,     iy,       iz,2)
            map_idx(6)=spin_input_index(ix,       iy+1,     iz,2)

            do j=1,n_dep
                spin_vec(i,:,j)=physical_spin(map_idx(j),:)
                do isite=1,n_input
                    if (map_idx(j) .eq. isite) then
                        dspin_vec(i,:,j,isite)=dphysical_spin(isite,:)
                    endif
                enddo
            enddo
        enddo
        enddo
        enddo
        enddo
        enddo
        enddo

        deallocate(theta_nu, phi_nu, physical_spin)
    end subroutine read_spin_derivatives

    subroutine read_RVR()
        allocate(RVR(n_func,nwan_rvr*d_ext,nwan_rvr*d_ext))
        allocate(RVR_tmp(nwan_rvr*d_ext,nwan_rvr*d_ext))

        RVR_dir=p_rvr_dir
        open(unit=103,file=trim(RVR_dir)//'./B_general',action='read', &
            access='sequential',form='unformatted',status='old')
        read(103) RVR
        close(103)

        allocate(dir_vec(d_ext,3))
        allocate(inv_dir(d_ext))
        inv_dir(:)=0
        i=1
        do i_dir=-convert_size_x,convert_size_x
        do j_dir=-convert_size_y,convert_size_y
        do k_dir=-convert_size_z,convert_size_z
            dir_vec(i,:)=(/i_dir,j_dir,k_dir/)
            i=i+1
        enddo
        enddo
        enddo

        do i=1,d_ext
        do j=1,d_ext
            if (all(dir_vec(i,:)+dir_vec(j,:) .eq. (/0,0,0/))) then
                inv_dir(i)=j
                exit
            endif
        enddo
            if (inv_dir(i) .eq. 0) then
                print *, i, '-dir inverse not found'
                stop
            endif
        enddo
    end subroutine read_RVR

    subroutine build_dH_wan()
        real(dp), allocatable :: coefficient(:)
        real(dp) :: channel_chi(6), channel_dchi(6), values(n_func)
        integer :: base_index, mod_index

        allocate(coefficient(n_func))
        allocate(dH_wan(H_dim,H_dim,nrept,n_response))
        dH_wan(:,:,:,:)=dcmplx(0.0d0,0.0d0)

        do response_index=1,n_response
            isite=response_sites(response_index)
            print *, 'Building dH/dtheta for physical site', isite-1

            do ix=0,spin_grid_x-1
            do iy=0,spin_grid_y-1
            do iz=0,spin_grid_z-1
                i=1+iz+spin_grid_z*iy+spin_grid_y*spin_grid_z*ix

                call spin_chirality(spin_vec(i,:,:), channel_chi, &
                     dspin_vec(i,:,:,isite), channel_dchi)
                call spin_channels(spin_vec(i,:,:), channel_chi, p_channel_mask, &
                     p_channel_source, p_lambda_high, values, dspin_vec(i,:,:,isite), &
                     channel_dchi, coefficient)
                if (sum(abs(coefficient(:))) .lt. 1.0d-14) cycle

                do i_dir=1,d_ext
                    ini_x=cell(i,1)+dir_vec(i_dir,1)
                    ini_y=cell(i,2)+dir_vec(i_dir,2)
                    ini_z=cell(i,3)+dir_vec(i_dir,3)

                    if (.not. (ini_x .ge. 0 .and. ini_x/period_x .eq. 0)) cycle
                    if (.not. (ini_y .ge. 0 .and. ini_y/period_y .eq. 0)) cycle
                    if (.not. (ini_z .ge. 0 .and. ini_z/period_z .eq. 0)) cycle

                    i_locat=1+ini_z+period_z*ini_y+period_z*period_y*ini_x

                    do j_dir=1,nwan_rvr*d_ext
                        RVR_tmp(:,j_dir)=matmul(coefficient(:),RVR(:,:,j_dir))
                    enddo

                    do j_dir=1,d_ext
                        fin_x=cell(i,1)+dir_vec(j_dir,1)
                        fin_y=cell(i,2)+dir_vec(j_dir,2)
                        fin_z=cell(i,3)+dir_vec(j_dir,3)

                        if (fin_x .ge. 0) then
                            super_x=fin_x/period_x
                        else
                            super_x=(fin_x+1)/period_x-1
                        endif
                        if (fin_y .ge. 0) then
                            super_y=fin_y/period_y
                        else
                            super_y=(fin_y+1)/period_y-1
                        endif
                        if (fin_z .ge. 0) then
                            super_z=fin_z/period_z
                        else
                            super_z=(fin_z+1)/period_z-1
                        endif

                        jrept=0
                        do irept=1,nrept
                            if (all((/super_x,super_y,super_z/) .eq. super_ir_vec(irept,:))) then
                                jrept=irept
                                exit
                            endif
                        enddo
                        if (jrept .eq. 0) then
                            print *, 'supercell vector not found:', super_x, super_y, super_z
                            stop
                        endif

                        jx=fin_x-super_x*period_x
                        jy=fin_y-super_y*period_y
                        jz=fin_z-super_z*period_z
                        j_locat=1+jz+period_z*jy+period_z*period_y*jx

                        do iwan=wan_min,wan_max
                        do jwan=wan_min,wan_max
                        do is=1,2/nspin
                        do js=1,2/nspin
                            iRVR=1+(i_dir-1)+d_ext*((iwan-wan_min)*2+is-1)
                            jRVR=1+(j_dir-1)+d_ext*((jwan-wan_min)*2+js-1)
                            i_dim=(is+2/nspin*(iwan-1))+nwan_cell*(i_locat-1)
                            j_dim=(js+2/nspin*(jwan-1))+nwan_cell*(j_locat-1)

                            if (abs(RVR_tmp(iRVR,jRVR)) .lt. 1.0d-12) cycle
                            dH_wan(i_dim,j_dim,jrept,response_index)= &
                                dH_wan(i_dim,j_dim,jrept,response_index) &
                                + hartree_to_ev*RVR_tmp(iRVR,jRVR)
                        enddo
                        enddo
                        enddo
                        enddo
                    enddo
                enddo
            enddo
            enddo
            enddo
        enddo

        deallocate(coefficient)
    end subroutine build_dH_wan

    subroutine compute_torque()
        n_DOS=1
        allocate(DOS_E_list(n_DOS))
        DOS_E_list(1)=Ef

        allocate(torque(n_DOS,1,n_response))
        torque(:,:,:)=0.0d0
        time_diag=0.0d0
        t_start=omp_get_wtime()

        !$OMP PARALLEL PRIVATE(ipool,n1,n2,ik,ikxy,iz,c_ix,c_iy,c_iz,c_jx,c_jy,c_jz,i,j, &
        !$OMP& irept,i1,i2,i3,i4,f1,f2,f3,f4,ctmp1,ctmp2,ctmp3,E1,i_eta,i_dim,iE, &
        !$OMP& isite,ib,t_diag_start,t_diag_end,width,inv_eta,delta_E_val,fermi_func,term_exp, &
        !$OMP& response_index, &
        !$OMP& H_wan_tmp,dH_wan_tmp,H_k,dH_k,H_step1,dH_step1,unmk,eig_tmp,dH_psi,tau_eig, &
        !$OMP& torque_local,phase0,phase_k,work,rwork,iwork,lwork,lrwork,liwork,info)

        allocate(H_k(H_dim,H_dim))
        allocate(dH_k(H_dim,H_dim))
        allocate(unmk(H_dim,H_dim))
        allocate(eig_tmp(H_dim))
        allocate(H_wan_tmp(H_dim,H_dim))
        allocate(dH_wan_tmp(H_dim,H_dim))
        allocate(H_step1(H_dim,H_dim))
        allocate(dH_step1(H_dim,H_dim))
        allocate(dH_psi(H_dim))
        allocate(tau_eig(H_dim,n_response))
        allocate(torque_local(n_DOS,1,n_response))
        allocate(phase0(super_size,super_size))
        allocate(phase_k(nrept))

        torque_local(:,:,:)=0.0d0

        allocate(work(1), rwork(1), iwork(1))
        call zheevd('V','U',H_dim,H_k,H_dim,eig_tmp,work,-1,rwork,-1,iwork,-1,info)
        lwork=int(real(work(1)))
        lrwork=int(rwork(1))
        liwork=iwork(1)
        deallocate(work,rwork,iwork)

        if (lwork < 2*H_dim + H_dim*H_dim) lwork = 2*H_dim + H_dim*H_dim
        if (lrwork < 1 + 5*H_dim + 2*H_dim*H_dim) lrwork = 1 + 5*H_dim + 2*H_dim*H_dim
        if (liwork < 3 + 5*H_dim) liwork = 3 + 5*H_dim

        allocate(work(lwork))
        allocate(rwork(lrwork))
        allocate(iwork(liwork))

        !$OMP DO SCHEDULE(dynamic)
        do ipool=1,npool
            ikxy=pool(ipool,1)
            ik=1+(1-1)+(ikxy-1)*nkz

            do c_ix=0,period_x-1
            do c_iy=0,period_y-1
            do c_iz=0,period_z-1
                i=1+c_iz+period_z*c_iy+period_z*period_y*c_ix
            do c_jx=0,period_x-1
            do c_jy=0,period_y-1
            do c_jz=0,period_z-1
                j=1+c_jz+period_z*c_jy+period_z*period_y*c_jx
                phase0(i,j)=exp(I2pi*dot_product(1.0d0*(/c_jx-c_ix,c_jy-c_iy,c_jz-c_iz/),k_list(ik,:)))
            enddo
            enddo
            enddo
            enddo
            enddo
            enddo

            H_wan_tmp(:,:)=dcmplx(0.0d0,0.0d0)
            do irept=1,nrept
                phase_k(irept)=exp(I2pi*sum(ir_vec(irept,:)*k_list(ik,:)))
                H_wan_tmp(:,:)=H_wan_tmp(:,:)+H_wan(:,:,irept)*phase_k(irept)
            enddo

            call transform_supercell_matrix(H_wan_tmp,H_k,H_step1,phase0)

            t_diag_start=omp_get_wtime()
            call diasym_opt(H_k,eig_tmp(:),H_dim,unmk(:,:),work,lwork,rwork,lrwork,iwork,liwork)
            t_diag_end=omp_get_wtime()

            !$OMP ATOMIC
            time_diag=time_diag+(t_diag_end-t_diag_start)

            tau_eig(:,:)=0.0d0
            do response_index=1,n_response
                dH_wan_tmp(:,:)=dcmplx(0.0d0,0.0d0)
                do irept=1,nrept
                    dH_wan_tmp(:,:)=dH_wan_tmp(:,:) &
                        + dH_wan(:,:,irept,response_index)*phase_k(irept)
                enddo

                call transform_supercell_matrix(dH_wan_tmp,dH_k,dH_step1,phase0)

                do ib=1,H_dim
                    dH_psi(:)=matmul(dH_k(:,:),unmk(:,ib))
                    tau_eig(ib,response_index)= &
                        -real(dot_product(unmk(:,ib),dH_psi(:)),dp)
                enddo
            enddo

            do ib=1,H_dim
                E1=eig_tmp(ib)
                i_eta=1
                width=20.0d0*eta(selected_eta)
                inv_eta=1.0d0/eta(selected_eta)

                do iE=1,n_DOS
                    delta_E_val=E1-DOS_E_list(iE)

                    if (delta_E_val < -width) then
                        torque_local(iE,i_eta,:)=torque_local(iE,i_eta,:)+tau_eig(ib,:)
                    elseif (delta_E_val > width) then
                        cycle
                    else
                        term_exp=exp(delta_E_val*inv_eta)
                        fermi_func=1.0d0/(1.0d0+term_exp)
                        torque_local(iE,i_eta,:)=torque_local(iE,i_eta,:) &
                            + fermi_func*tau_eig(ib,:)
                    endif
                enddo
            enddo
        enddo
        !$OMP END DO

        !$OMP CRITICAL
        torque(:,:,:)=torque(:,:,:)+torque_local(:,:,:)
        !$OMP END CRITICAL

        deallocate(H_k,dH_k,unmk,eig_tmp,H_wan_tmp,dH_wan_tmp,H_step1,dH_step1,dH_psi,tau_eig)
        deallocate(torque_local,phase0,phase_k,work,rwork,iwork)

        !$OMP END PARALLEL

        t_end=omp_get_wtime()
        torque(:,:,:)=torque(:,:,:)/nk

        print *, 'Total Wall Time:', t_end-t_start
        print *, 'Total Diag Time:', time_diag
    end subroutine compute_torque

    subroutine transform_supercell_matrix(matrix_in,matrix_out,step1,phase0_in)
        complex(dp), intent(in) :: matrix_in(H_dim,H_dim)
        complex(dp), intent(out) :: matrix_out(H_dim,H_dim)
        complex(dp), intent(inout) :: step1(H_dim,H_dim)
        complex(dp), intent(in) :: phase0_in(super_size,super_size)
        integer :: li, lj, ln1, ln2
        integer :: li1, li2, li3, li4, lf1, lf2, lf3, lf4
        complex(dp) :: lctmp1, lctmp2, lctmp3

        matrix_out(:,:)=dcmplx(0.0d0,0.0d0)
        step1(:,:)=dcmplx(0.0d0,0.0d0)

        do li=1,super_size
            li3=nwan_cell*(li-1)+1
            lf3=nwan_cell*li

            do ln2=1,super_size
                li2=nwan_cell*(ln2-1)+1
                lf2=nwan_cell*ln2

                do ln1=1,super_size
                    lctmp2=phase1(li,ln1)
                    li1=nwan_cell*(ln1-1)+1
                    lf1=nwan_cell*ln1
                    step1(li3:lf3,li2:lf2)=step1(li3:lf3,li2:lf2)+lctmp2*matrix_in(li1:lf1,li2:lf2)
                enddo
            enddo
        enddo

        do li=1,super_size
            li3=nwan_cell*(li-1)+1
            lf3=nwan_cell*li

            do lj=1,super_size
                lctmp1=phase0_in(li,lj)
                li4=nwan_cell*(lj-1)+1
                lf4=nwan_cell*lj

                do ln2=1,super_size
                    lctmp3=lctmp1*conjg(phase1(lj,ln2))
                    li2=nwan_cell*(ln2-1)+1
                    lf2=nwan_cell*ln2
                    matrix_out(li3:lf3,li4:lf4)=matrix_out(li3:lf3,li4:lf4)+lctmp3*step1(li3:lf3,li2:lf2)
                enddo
            enddo
        enddo

        matrix_out(:,:)=matrix_out(:,:)/super_size
    end subroutine transform_supercell_matrix

    subroutine write_direct_torque()
        call int2char(selected_eta,eta_char)

        open(unit=101,file='./direct_torque_'//trim(eta_char), &
            action='write',access='sequential',status='replace')

        write(101,'(A)',advance='no') '# DOS_E'
        do response_index=1,n_response
            call int2char(response_sites(response_index)-1,site_char)
            write(101,'(A)',advance='no') ' tau_site_'//trim(site_char)
        enddo
        write(101,*)

        do iE=1,n_DOS
            write(101,'(E24.12)',advance='no') DOS_E_list(iE)
            do response_index=1,n_response
                write(101,'(1X,E24.12)',advance='no') torque(iE,1,response_index)
            enddo
            write(101,*)
        enddo

        close(101)
    end subroutine write_direct_torque

    integer function spin_input_index(ix_in,iy_in,iz_in,atom)
        integer, intent(in) :: ix_in, iy_in, iz_in, atom
        integer :: ix0, iy0, iz0

        ix0=modulo(ix_in,period_x)
        iy0=modulo(iy_in,period_y)
        iz0=modulo(iz_in,period_z)

        spin_input_index=atom+2*(ix0+period_x*(iy0+period_y*iz0))
    end function spin_input_index

    subroutine int2char(in_int,out_char)
        integer, intent(in) :: in_int
        character(len=20), intent(out) :: out_char

        out_char=' '
        if (in_int .lt. 10) write(out_char,'(I1)') in_int
        if ((in_int .lt. 100) .and. (in_int .ge. 10)) write(out_char,'(I2)') in_int
        if ((in_int .lt. 1000) .and. (in_int .ge. 100)) write(out_char,'(I3)') in_int
        if ((in_int .lt. 10000) .and. (in_int .ge. 1000)) write(out_char,'(I4)') in_int
        if ((in_int .lt. 100000) .and. (in_int .ge. 10000)) write(out_char,'(I5)') in_int
    end subroutine int2char

    subroutine diasym_opt(a_in,eig,n,D,work_in,lwork_in,rwork_in,lrwork_in,iwork_in,liwork_in)
        integer, intent(in) :: n, lwork_in, lrwork_in, liwork_in
        complex(dp), intent(in) :: a_in(n,n)
        complex(dp), intent(out) :: D(n,n)
        real(dp), intent(out) :: eig(n)
        complex(dp) :: work_in(lwork_in)
        real(dp) :: rwork_in(lrwork_in)
        integer :: iwork_in(liwork_in)
        integer :: info_local

        D(:,:)=a_in(:,:)
        call zheevd('V','U',n,D,n,eig,work_in,lwork_in,rwork_in,lrwork_in,iwork_in,liwork_in,info_local)
        if (info_local /= 0) stop 'Error in zheevd'
    end subroutine diasym_opt

end program direct_torque_hf
