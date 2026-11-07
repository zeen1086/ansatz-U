program wan_interpolat
  use parameter_mod
    use mpi
    implicit none

    ! 定义双精度参数
    integer, parameter :: dp = 8

    complex(kind=dp), allocatable :: H_wan(:,:,:), H_k(:,:), H_buff(:,:), eig_vec(:,:,:)
    complex(kind=dp), allocatable :: unmk(:,:), mat_tmp(:,:), unmk_Hermit(:,:)
    complex(kind=dp), allocatable :: delHH(:,:,:), D_h(:,:,:)
    complex(kind=dp), allocatable :: H_wan_tmp(:,:), RH_tmp(:,:,:)
    complex(kind=dp), allocatable :: sig_H(:,:,:), sig_AH(:,:,:)
    complex(kind=dp), allocatable :: sig_H_local(:,:,:), sig_AH_local(:,:,:)
    complex(kind=dp), allocatable :: phase0(:,:), phase1(:,:), phase_k(:)

    real(kind=dp), allocatable    :: E_list(:), eig_list(:,:), spin_list(:,:), spin_tmp(:)
    real(kind=dp), allocatable    :: k_list(:,:), eigen(:,:), eigen_local(:,:)
    real(kind=dp), allocatable    :: cr_vec(:,:), eig_tmp(:)
    real(kind=dp), allocatable    :: nb(:,:), R_list(:,:), pt_k(:)

    integer, allocatable :: ir_vec(:,:), pool(:,:), read3_int(:)

    integer :: ik, jk, nk, nkx, nky, nkz, irept, nrept, io
    integer :: i, j, wan_x, wan_y, wan_z, ix, iy, iz
    integer :: ikxy, H_dim, nspin, period_x, period_y, period_z, super_size
    integer :: c_ix, c_iy, c_iz, c_jx, c_jy, c_jz, n1, n2, i_dim, j_dim
    integer :: nE, iE, i1, i2, i3, i4, f1, f2, f3, f4
    integer :: n_cell, nwan, nline, iline, npt, ipool, muti_pool, imuti_pool
    integer :: ierr, my_rank, size, npool

    real(kind=dp) :: a0, a(3,3), Ef, fac, pt_ini_k(3), pt_fin_k(3)
    real(kind=dp) :: length, vec_tmp(3), b(3,3), celldm3, eta
    real(kind=dp) :: t1, t2, t4, E_max, E_min, pi

    complex(kind=dp) :: cI, I2pi, ctmp1, ctmp2, ctmp3
    character(len=20) :: irept_char, kz_char
    character(len=200) :: useless_string, nscf_dir, rvr_root, title

    ! MPI 初始化
    call mpi_init(ierr)
    call mpi_comm_rank(MPI_COMM_WORLD, my_rank, ierr)
    call mpi_comm_size(MPI_COMM_WORLD, size, ierr)

    ! 读取控制文件
nscf_dir = p_nscf_dir
rvr_root = p_rvr_dir
nwan = p_nwan
nspin = p_nspin
nkx = p_nx
nky = p_ny
nkz = p_nz
period_x = p_period_x
period_y = p_period_y
period_z = p_period_z
eta = p_eta
fac = p_fac
a0 = p_a0
celldm3 = p_celldm3
Ef = p_ef

    ! 晶格常数设置
    a(1,:) = [1.0_dp, 0.0_dp, 0.0_dp]
    a(2,:) = [-0.5_dp, 0.866025_dp, 0.0_dp]
    a(3,:) = [0.0_dp, 0.0_dp, celldm3]

    b(1,:) = [1.0_dp, 0.577350_dp, 0.0_dp]
    b(2,:) = [0.0_dp, 1.154701_dp, 0.0_dp]
    b(3,:) = [0.0_dp, 0.0_dp, 1.0_dp/celldm3]

    a = a * a0 * 0.529177_dp

    ! 能量范围
    E_max = p_band_energy_max
    E_min = p_band_energy_min
    nE = p_band_energy_points
    allocate(E_list(nE))
    allocate(sig_H(2,2,nE), sig_AH(2,2,nE))
    allocate(sig_H_local(2,2,nE), sig_AH_local(2,2,nE))

    do iE = 1, nE
        E_list(iE) = E_min + (iE-1)*(E_max-E_min)/nE
    end do

    cI = cmplx(0.0_dp, 1.0_dp, kind=dp)
    pi = 3.1415926535897932_dp
    I2pi = 2.0_dp * pi * cI
    H_dim = nwan * 2 / nspin * period_x * period_y * period_z

    ! 检查周期性
    if (mod(nkx, period_x) /= 0 .or. mod(nky, period_y) /= 0 .or. mod(nkz, period_z) /= 0) then
        print*, 'not commensurate'
        call mpi_finalize(ierr)
        stop
    end if
    nkx = nkx / period_x
    nky = nky / period_y
    nkz = nkz / period_z

    ! 内存分配
    allocate(H_k(H_dim,H_dim), H_buff(H_dim,H_dim))
    allocate(unmk(H_dim,H_dim), eig_tmp(H_dim))
    allocate(delHH(H_dim,H_dim,3), D_h(H_dim,H_dim,3))
    allocate(mat_tmp(H_dim,H_dim), unmk_Hermit(H_dim,H_dim))
    allocate(spin_tmp(H_dim))

    ! 读取不可约矢量
    irept = 0
    allocate(read3_int(3))
    open(unit=101, file=trim(p_super_h_dir)//'ir_vec_HYC', status='old')
    do
        read(101,'(3I8)', iostat=io) read3_int
        if (io /= 0) exit
        irept = irept + 1
    end do
    close(101)

    nrept = irept
    allocate(ir_vec(nrept,3), cr_vec(nrept,3))
    allocate(H_wan(H_dim,H_dim,nrept), H_wan_tmp(H_dim,H_dim))
    allocate(RH_tmp(H_dim,H_dim,3), phase_k(nrept))

    irept = 1
    open(unit=101, file=trim(p_super_h_dir)//'ir_vec_HYC', status='old')
    do irept=1,nrept
        read(101,'(3I8)') ir_vec(irept,:)
    enddo
    close(101)

    ir_vec(:,1) = ir_vec(:,1) * period_x
    ir_vec(:,2) = ir_vec(:,2) * period_y
    ir_vec(:,3) = ir_vec(:,3) * period_z

    do irept = 1, nrept
        cr_vec(irept,1) = sum(a(:,1)*ir_vec(irept,:))
        cr_vec(irept,2) = sum(a(:,2)*ir_vec(irept,:))
        cr_vec(irept,3) = sum(a(:,3)*ir_vec(irept,:))
    end do

    ! 读取 K 路径
    open(unit=101, file='./plt_path.in', status='old')
    read(101,'(A200)') title
    if (index(title,'! Generated')==1) read(101,'(A200)') title
    read(101,'(A6,I3)') useless_string, nline
    read(101,'(A4,I3)') useless_string, npt
    if (min(nline,npt)<1) error stop 'band path needs positive line and point counts'
    nk = nline * npt

    allocate(pt_k(nk), eig_list(nk,H_dim), spin_list(nk,H_dim))
    allocate(k_list(nk,3), pool(nk,1))

    if (my_rank == 0) open(unit=102, file='./path_check.out', access='append')

    do iline = 1, nline
        read(101,'(6F12.6)') pt_ini_k, pt_fin_k
        if (my_rank == 0) write(102,'(6F12.6)') pt_ini_k, pt_fin_k

        vec_tmp(1) = sum(b(:,1)*(pt_fin_k - pt_ini_k))
        vec_tmp(2) = sum(b(:,2)*(pt_fin_k - pt_ini_k))
        vec_tmp(3) = sum(b(:,3)*(pt_fin_k - pt_ini_k))
        length = sqrt(sum(vec_tmp**2)) / npt

        do ik = npt*(iline-1)+1, npt*iline
            k_list(ik,:) = pt_ini_k + (ik-npt*(iline-1))*(pt_fin_k - pt_ini_k)/npt
            pool(ik,1) = ik
            if (ik == 1) then
                pt_k(ik) = length
            else
                pt_k(ik) = pt_k(ik-1) + length
            end if
        end do
    end do
    close(101)
    if (my_rank == 0) close(102)

    k_list(:,1) = k_list(:,1) / period_x
    k_list(:,2) = k_list(:,2) / period_y
    k_list(:,3) = k_list(:,3) / period_z

    super_size = period_x * period_y * period_z
    allocate(nb(super_size,3), R_list(super_size,3))
    allocate(phase0(super_size,super_size), phase1(super_size,super_size))

    do ix = 0, period_x-1
    do iy = 0, period_y-1
    do iz = 0, period_z-1
        i = 1 + iz + period_z*iy + period_z*period_y*ix
        nb(i,:) = [real(ix,dp)/period_x, real(iy,dp)/period_y, real(iz,dp)/period_z]
        R_list(i,1) = sum(a(:,1)*[ix, iy, iz])
        R_list(i,2) = sum(a(:,2)*[ix, iy, iz])
        R_list(i,3) = sum(a(:,3)*[ix, iy, iz])
    end do
    end do
    end do

    muti_pool = (nk-1)/max(1,size-1)
    npool = nk

    ! 主进程读取 Hamiltonian 并广播
    if (my_rank == 0) then
        open(unit=101, file='report', access='append')
        write(101,*) 'start reading'
        close(101)

        do irept = 1, nrept
            call int2char(irept, irept_char)
            ! 注意：这里读取双精度数据
            open(unit=101, file='./tot_HR/hr_'//trim(irept_char), &
                 access='stream', form='unformatted')
            read(101) H_wan(:,:,irept)
            close(101)

            ! 使用 MPI_DOUBLE_COMPLEX
            call MPI_Bcast(H_wan(:,:,irept), H_dim*H_dim, &
                           MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)
            call MPI_Barrier(MPI_COMM_WORLD, ierr)
        end do

        eig_list = 0.0_dp
        do iz = 1, 1
            if (size==1) then
                do ik=1,nk
                    call calculate_point()
                    eig_list(ik,:)=eig_tmp
                    spin_list(ik,:)=spin_tmp
                enddo
            else
            do ipool = 1, nk
                ! 使用 MPI_DOUBLE_PRECISION
                call MPI_RECV(eig_list(ipool,:), H_dim, MPI_DOUBLE_PRECISION, &
                              MPI_ANY_SOURCE, ipool, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
                call MPI_RECV(spin_list(ipool,:), H_dim, MPI_DOUBLE_PRECISION, &
                              MPI_ANY_SOURCE, ipool+npool, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
            end do
            endif
            call MPI_Barrier(MPI_COMM_WORLD, ierr)
        end do

        open(unit=101, file='./plt_bnd_'//trim(title), access='append')
        do i_dim = 1, H_dim
            do ik = 1, nk
                write(101,'(I8,6F18.6)') ik, k_list(ik,:), pt_k(ik), &
                                         eig_list(ik,i_dim), spin_list(ik,i_dim)
            end do
            write(101,*)
        end do
        close(101)
    end if

    ! 子进程计算
    if (my_rank /= 0) then
        do irept = 1, nrept
            call MPI_Bcast(H_wan(:,:,irept), H_dim*H_dim, &
                           MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)
            call MPI_Barrier(MPI_COMM_WORLD, ierr)
        end do

        do iz = 1, 1
            do imuti_pool = 0, muti_pool
                ipool = imuti_pool*(size-1) + my_rank
                if (ipool > npool) cycle

                ikxy = pool(ipool,iz)
                ik = ikxy
                call calculate_point()

                call MPI_SEND(eig_tmp, H_dim, MPI_DOUBLE_PRECISION, &
                              0, ikxy, MPI_COMM_WORLD, ierr)
                call MPI_SEND(spin_tmp, H_dim, MPI_DOUBLE_PRECISION, &
                              0, ikxy+npool, MPI_COMM_WORLD, ierr)
            end do
            call MPI_Barrier(MPI_COMM_WORLD, ierr)
        end do
    end if

    call mpi_finalize(ierr)

contains

    subroutine calculate_point()
                call cpu_time(t1)

                ! 计算相位 phase0
                do c_ix=0,period_x-1
                do c_iy=0,period_y-1
                do c_iz=0,period_z-1
                    i=1+c_iz+period_z*c_iy+period_z*period_y*c_ix
                    do c_jx=0,period_x-1
                    do c_jy=0,period_y-1
                    do c_jz=0,period_z-1
                        j=1+c_jz+period_z*c_jy+period_z*period_y*c_jx
                        phase0(i,j) = exp(I2pi * dot_product([real(c_jx-c_ix,dp), &
                            real(c_jy-c_iy,dp), real(c_jz-c_iz,dp)], k_list(ik,:)))
                    end do
                    end do
                    end do
                end do
                end do
                end do

                ! 计算相位 phase1
                do n1=1,super_size
                do c_ix=0,period_x-1
                do c_iy=0,period_y-1
                do c_iz=0,period_z-1
                    i=1+c_iz+period_z*c_iy+period_z*period_y*c_ix
                    phase1(i,n1) = exp(-I2pi * dot_product([real(c_ix,dp), &
                        real(c_iy,dp), real(c_iz,dp)], nb(n1,:)))
                end do
                end do
                end do
                end do

                H_wan_tmp = cmplx(0.0_dp, 0.0_dp, kind=dp)
                do irept = 1, nrept
                    phase_k(irept) = exp(I2pi * sum(ir_vec(irept,:)*k_list(ik,:)))
                    H_wan_tmp = H_wan_tmp + H_wan(:,:,irept) * phase_k(irept)
                end do

                H_k = cmplx(0.0_dp, 0.0_dp, kind=dp)

                ! 构建 Hamiltonian
                do i = 1, super_size
                    i3 = 2/nspin*nwan*(i-1) + 1
                    f3 = 2/nspin*nwan*i
                    do j = 1, super_size
                        ctmp1 = phase0(i,j)
                        i4 = 2/nspin*nwan*(j-1) + 1
                        f4 = 2/nspin*nwan*j
                        do n1 = 1, super_size
                            ctmp2 = phase1(i,n1) * ctmp1
                            i1 = 2/nspin*nwan*(n1-1) + 1
                            f1 = 2/nspin*nwan*n1
                            do n2 = 1, super_size
                                ctmp3 = conjg(phase1(j,n2)) * ctmp2
                                i2 = 2/nspin*nwan*(n2-1) + 1
                                f2 = 2/nspin*nwan*n2

                                H_k(i1:f1, i2:f2) = H_k(i1:f1, i2:f2) + &
                                                    ctmp3 * H_wan_tmp(i3:f3, i4:f4)
                            end do
                        end do
                    end do
                end do

                H_k = H_k / real(super_size, dp)
                call cpu_time(t2)

                ! 调用对角化
                call diasym(H_k, eig_tmp, H_dim, unmk)

                spin_tmp = 0.0_dp
                do i_dim = 1, H_dim
                    do j_dim = 1, H_dim
                        spin_tmp(i_dim) = spin_tmp(i_dim) + &
                                          (-1.0_dp)**(j_dim+1) * abs(unmk(j_dim,i_dim))**2
                    end do
                end do

                call cpu_time(t4)

                if (my_rank == 1) then
                    open(unit=101, file='report', access='append')
                    write(101,*) 'build', t2-t1, 'diag', t4-t2
                    close(101)
                end if

    end subroutine calculate_point

    subroutine int2char(in_int, out_char)
        integer, intent(in) :: in_int
        character(len=20), intent(out) :: out_char
        if (in_int < 10) write(out_char,'(I1)') in_int
        if (in_int < 100 .and. in_int >= 10) write(out_char,'(I2)') in_int
        if (in_int < 1000 .and. in_int >= 100) write(out_char,'(I3)') in_int
        if (in_int < 10000 .and. in_int >= 1000) write(out_char,'(I4)') in_int
        if (in_int < 100000 .and. in_int >= 10000) write(out_char,'(I5)') in_int
    end subroutine int2char

    subroutine diasym(a, eig, n, D)
        implicit none
        integer, intent(in) :: n
        complex(kind=dp), intent(in) :: a(n,n)
        complex(kind=dp) :: a_tmp(n,n)
        complex(kind=dp), intent(out) :: D(n,n)

        integer :: l, inf
        ! LAPACK 工作数组大小需足够
        complex(kind=dp), allocatable :: Lwork(:)
        real(kind=dp), intent(out) :: eig(n)
        real(kind=dp), allocatable :: Rwork(:)

        allocate(Lwork(2*n-1))
        allocate(Rwork(3*n-2))

        a_tmp = a
        l = 2*n-1

        ! 使用双精度 LAPACK 程序 ZHEEV
        call zheev('V', 'U', n, a_tmp, n, eig, Lwork, l, Rwork, inf)

        if (inf /= 0) print*, "Error in diagonalization: ", inf

        D = a_tmp
        deallocate(Lwork, Rwork)
    end subroutine diasym

end program wan_interpolat
