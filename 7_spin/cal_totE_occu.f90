program wan_interpolat
  use parameter_mod
    use omp_lib
    implicit none

    complex(kind(0d0)),allocatable  :: H_wan(:,:,:),H_k(:,:)
    integer :: ik,nk,nkx,nky,nkz,irept,nrept,io,read3_int(3)
    integer :: i,j,ix,iy,iz,nwan_sp
    integer :: ikxy,H_dim
    integer :: nspin, period_x, period_y, period_z,super_size
    integer :: c_ix,c_iy,c_iz,c_jx,c_jy,c_jz,n1,n2,i_dim
    integer :: iE,i1,i2,i3,i4,f1,f2,f3,f4
    integer :: n_DOS
    real(kind(0d0)),allocatable     :: DOS_E_list(:),occu(:,:),occu_local(:,:)
    real(kind(0d0)),allocatable     :: tot_E(:,:),tot_E_local(:,:)
    real(kind(0d0)),allocatable     :: k_list(:,:)
    real(kind(0d0)),allocatable     :: cr_vec(:,:),eig_tmp(:)
    real(kind(0d0)),allocatable     :: nb(:,:),R_list(:,:)
    real(kind(0d0)) :: a0,a(3,3),Ef,fac,celldm3
    integer :: nwan
    integer,allocatable     :: ir_vec(:,:)
    character(len=20)       :: irept_char
    complex(kind(0d0))      :: cI,ctmp1,ctmp2,ctmp3
    complex(kind(0d0))      :: I2pi
    complex(kind(0d0)),allocatable:: phase0(:,:),phase1(:,:),phase_k(:)
    complex(kind(0d0)),allocatable:: H_wan_tmp(:,:)
    real(kind(0d0))         :: pi,DOS_max,DOS_min
    character(len=200)      :: useless_string,nscf_dir,rvr_root
    character(len=20)       :: eta_char
    integer,allocatable     :: pool(:,:)
    integer                 :: npool ,ipool,n_eta,i_eta
    real(kind(0d0)), allocatable       :: eta(:)
    real(kind(0d0))    :: read1,E1

    real(kind=8) :: t_start, t_end, t_diag_start, t_diag_end
    real(kind=8) :: time_diag

    complex(kind(0d0)), allocatable :: H_step1(:,:)

    complex(kind(0d0)), allocatable :: work(:)
    real(kind(0d0)),    allocatable :: rwork(:)
    integer,            allocatable :: iwork(:)
    integer :: lwork, lrwork, liwork, info

    real(kind=8) :: width, inv_eta, delta_E_val, fermi_func, term_exp

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
fac = p_fac
a0 = p_a0
celldm3 = p_celldm3
Ef = p_ef

        n_eta = ubound(p_broadening,1)
        allocate(eta(n_eta))
        eta = p_broadening

    a(1,:) = (/   1.000000d0,   0.000000d0,   0.000000d0 /)
    a(2,:) = (/   0.000000d0,   1.000000d0,   0.000000d0 /)
    a(3,:) = (/   0.000000d0,   0.000000d0,   celldm3  /)
    a(:,:)=a(:,:)*a0*0.529177d0

    time_diag = 0.0d0
    t_start = omp_get_wtime()

    DOS_max = p_dos_max
    DOS_min = p_dos_min
    n_DOS = p_dos_points

    allocate(DOS_E_list(n_DOS))
    allocate(occu(n_DOS,n_eta))
    allocate(tot_E(n_DOS,n_eta))

    Do iE=1,n_DOS
        DOS_E_list(iE)=DOS_min+(iE-1)*(DOS_max-DOS_min)/n_DOS
    ENDDO

    cI=dcmplx(0.000d0,1.000d0)
    pi=3.1415926535d0
    I2pi=2*pi*cI
    nwan_sp=nwan*2/nspin
    H_dim=nwan_sp*period_x*period_y*period_z

    if (.not. mod(nkx,period_x) .eq. 0) stop 'not conmenserate'
    if (.not. mod(nky,period_y) .eq. 0) stop 'not conmenserate'
    if (.not. mod(nkz,period_z) .eq. 0) stop 'not conmenserate'
    nkx=nkx/period_x
    nky=nky/period_y
    nkz=nkz/period_z
    nk=nkx*nky*nkz

    irept=0
    open(unit=101,file=trim(p_super_h_dir)//'ir_vec_HYC',action='read',access='sequential')
    do
        read(101,'(3I8)',iostat=io) read3_int(:)
        if(io/=0) exit
        irept=irept+1
    enddo
    close(101)
    nrept=irept
    allocate(ir_vec(nrept,3))
    allocate(cr_vec(nrept,3))
    allocate(H_wan(H_dim,H_dim,nrept))

    irept=1
    open(unit=101,file=trim(p_super_h_dir)//'ir_vec_HYC',action='read',access='sequential')
    do
        read(101,'(3I8)',iostat=io) read3_int(:)
        if(io/=0) exit
        ir_vec(irept,:)=read3_int(:)
        irept=irept+1
    enddo
    close(101)

    ir_vec(:,1)=ir_vec(:,1)*period_x
    ir_vec(:,2)=ir_vec(:,2)*period_y
    ir_vec(:,3)=ir_vec(:,3)*period_z

    do irept=1,nrept
        cr_vec(irept,1)=sum(a(:,1)*ir_vec(irept,:))
        cr_vec(irept,2)=sum(a(:,2)*ir_vec(irept,:))
        cr_vec(irept,3)=sum(a(:,3)*ir_vec(irept,:))
    enddo

    allocate(k_list(nk,3))
    npool=nk/nkz
    allocate(pool(npool,nkz))
    do ix=1,nkx
    do iy=1,nky
    do iz=1,nkz
        ik=1+(iz-1)+(iy-1)*nkz+(ix-1)*nky*nkz
        k_list(ik,1)=1.0*(ix-1)/(nkx*period_x)
        k_list(ik,2)=1.0*(iy-1)/(nky*period_y)
        k_list(ik,3)=1.0*(iz-1)/(nkz*period_z)
        ipool=1+(iy-1)+(ix-1)*nky
        pool(ipool,iz)=ipool
    enddo
    enddo
    enddo

    super_size=period_x*period_y*period_z
    allocate(nb(super_size,3))
    allocate(phase1(super_size,super_size))
    allocate(R_list(super_size,3))
    do ix=0,period_x-1
    do iy=0,period_y-1
    do iz=0,period_z-1
        i=1+iz+period_z*iy+period_z*period_y*ix
        nb(i,:)=(/ix*1.000d0/period_x,iy*1.000d0/period_y,iz*1.000d0/period_z/)

        R_list(i,1)=sum(a(:,1)*(/ix,iy,iz/))
        R_list(i,2)=sum(a(:,2)*(/ix,iy,iz/))
        R_list(i,3)=sum(a(:,3)*(/ix,iy,iz/))
    enddo
    enddo
    enddo
    print*,'super_size',super_size

    do irept=1,nrept
        call int2char(irept, irept_char)
        Open(unit=101,file='./tot_HR/hr_'//trim(irept_char),action='read'&
                &,access='stream',form='unformatted')
        read(101) H_wan(:,:,irept)
        close(101)
    enddo

    do n1=1,super_size
        do c_ix=0,period_x-1
        do c_iy=0,period_y-1
        do c_iz=0,period_z-1
            i=1+c_iz+period_z*c_iy+period_z*period_y*c_ix
            phase1(i,n1)=exp(-I2pi*dot_product(1.000d0*(/c_ix,c_iy,c_iz/),nb(n1,:)))
        enddo
        enddo
        enddo
    enddo

    occu(:,:)=0.00d0
    tot_E(:,:)=0.00d0

    !$OMP PARALLEL PRIVATE(ipool,n1,n2,ik,ikxy,iz,c_ix,c_iy,c_iz,c_jx,c_jy,c_jz,i,j,&
    !$OMP& irept,i1,i2,i3,i4,f1,f2,f3,f4,ctmp1,ctmp2,ctmp3,E1,i_eta,i_dim,iE) &
    !$OMP& PRIVATE(H_wan_tmp, H_k, eig_tmp, occu_local, tot_E_local) &
    !$OMP& PRIVATE(t_diag_start, t_diag_end) &
    !$OMP& PRIVATE(H_step1, work, rwork, iwork, lwork, lrwork, liwork, info) &
    !$OMP& PRIVATE(width, inv_eta, delta_E_val, fermi_func, term_exp) &
    !$OMP& PRIVATE(phase0, phase_k) &
    !$OMP& SHARED(H_wan,k_list,npool,super_size,period_x,period_y,period_z,phase1,I2pi,nb,nkz,ir_vec,&
    !$OMP& nwan_sp,H_dim,n_DOS,n_eta,pool,DOS_E_LIST,eta,occu,tot_E,nrept,time_diag,DOS_min)&
    !$OMP& DEFAULT(NONE)

    allocate(H_k(H_dim,H_dim))
    allocate(eig_tmp(H_dim))
    allocate(H_wan_tmp(H_dim,H_dim))
    allocate(occu_local(n_DOS,n_eta))
    allocate(tot_E_local(n_DOS,n_eta))
    allocate(H_step1(H_dim, H_dim))
    allocate(phase0(super_size, super_size))
    allocate(phase_k(nrept))

    occu_local(:,:)=0.00d0
    tot_E_local(:,:)=0.00d0

    allocate(work(1), rwork(1), iwork(1))
    call zheevd('N', 'U', H_dim, H_k, H_dim, eig_tmp, work, -1, rwork, -1, iwork, -1, info)
    lwork  = int(real(work(1)))
    lrwork = int(rwork(1))
    liwork = iwork(1)
    deallocate(work, rwork, iwork)

    if (lwork < H_dim + 1) lwork = H_dim + 1
    if (lrwork < H_dim) lrwork = H_dim
    if (liwork < 1) liwork = 1

    allocate(work(lwork))
    allocate(rwork(lrwork))
    allocate(iwork(liwork))

    !$OMP DO
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
                    phase0(i,j)=exp(I2pi*&
                            &dot_product(1.000d0*(/c_jx-c_ix,c_jy-c_iy,c_jz-c_iz/),k_list(ik,:)))
            enddo
            enddo
            enddo
            enddo
            enddo
            enddo

            H_wan_tmp(:,:)=dcmplx(0.000d0,0.000d0)
            do irept=1,nrept
                    phase_k(irept)=exp(I2pi*sum(ir_vec(irept,:)*k_list(ik,:)))
                    H_wan_tmp(:,:)=H_wan_tmp(:,:)+H_wan(:,:,irept)*phase_k(irept)
            enddo

            H_k(:,:)=dcmplx(0.000d0,0.000d0)
            H_step1(:,:) = dcmplx(0.000d0, 0.000d0)

            do i = 1, super_size
                i3 = nwan_sp*(i-1) + 1
                f3 = nwan_sp*i

                do n2 = 1, super_size
                    i2 = nwan_sp*(n2-1) + 1
                    f2 = nwan_sp*n2

                    do n1 = 1, super_size
                        ctmp2 = phase1(i, n1)
                        i1 = nwan_sp*(n1-1) + 1
                        f1 = nwan_sp*n1
                        H_step1(i3:f3, i2:f2) = H_step1(i3:f3, i2:f2) + &
                                                ctmp2 * H_wan_tmp(i1:f1, i2:f2)
                    enddo
                enddo
            enddo

            do i = 1, super_size
                i3 = nwan_sp*(i-1) + 1
                f3 = nwan_sp*i

                do j = 1, super_size
                    ctmp1 = phase0(i, j)
                    i4 = nwan_sp*(j-1) + 1
                    f4 = nwan_sp*j

                    do n2 = 1, super_size
                        ctmp3 = ctmp1 * conjg(phase1(j, n2))
                        i2 = nwan_sp*(n2-1) + 1
                        f2 = nwan_sp*n2
                        H_k(i3:f3, i4:f4) = H_k(i3:f3, i4:f4) + &
                                            ctmp3 * H_step1(i3:f3, i2:f2)
                    enddo
                enddo
            enddo

            H_k=H_k/super_size

            t_diag_start = omp_get_wtime()
            call diasym_eigvals(H_k, eig_tmp(:), H_dim, work, lwork, rwork, lrwork, iwork, liwork)
            t_diag_end = omp_get_wtime()

            !$OMP ATOMIC
            time_diag = time_diag + (t_diag_end - t_diag_start)

            do i_dim=1,H_dim
                E1=eig_tmp(i_dim)

                do i_eta=1,n_eta
                    width = 20.0d0 * eta(i_eta)
                    inv_eta = 1.0d0 / eta(i_eta)

                    do iE = 1, n_DOS
                        delta_E_val = E1 - DOS_E_list(iE)

                        if (delta_E_val < -width) then
                            occu_local(iE, i_eta) = occu_local(iE, i_eta) + 1.0d0
                            tot_E_local(iE, i_eta) = tot_E_local(iE, i_eta) + E1

                        elseif (delta_E_val > width) then
                            cycle

                        else
                            term_exp = exp(delta_E_val * inv_eta)
                            fermi_func = 1.0d0 / (1.0d0 + term_exp)

                            occu_local(iE, i_eta) = occu_local(iE, i_eta) + fermi_func
                            tot_E_local(iE, i_eta) = tot_E_local(iE, i_eta) + E1 * fermi_func
                        endif
                    enddo
                enddo
            enddo

    ENDDO
    !$OMP END DO

    !$OMP CRITICAL
    occu(:,:)=occu(:,:)+occu_local(:,:)
    tot_E(:,:)=tot_E(:,:)+tot_E_local(:,:)
    !$OMP END CRITICAL

    deallocate(H_k, H_step1, eig_tmp, H_wan_tmp)
    deallocate(occu_local, tot_E_local)
    deallocate(work, rwork, iwork)
    deallocate(phase0, phase_k)

    !$OMP END PARALLEL

    t_end = omp_get_wtime()

    print *, "Total Wall Time: ", t_end - t_start
    print *, "Total Diag Time (Sum of all threads): ", time_diag
    print *, "Ratio (Diag/Total/Threads): ", time_diag / (t_end - t_start) / omp_get_max_threads()

    occu(:,:)=occu(:,:)/nk
    tot_E(:,:)=tot_E(:,:)/nk

    do i_eta=1,n_eta
            call int2char(i_eta,eta_char)

            open(unit=101,file='./occupation_'//trim(eta_char)&
                                    &,action='write',access='append',status='replace')
                    do iE=1,n_DOS
                            write(101,'(2E24.12)')    DOS_E_list(iE),occu(iE,i_eta)
                    enddo
            close(101)

            open(unit=101,file='./tot_energy_'//trim(eta_char)&
                                    &,action='write',access='append',status='replace')
                    do iE=1,n_DOS
                            write(101,'(2E24.12)')    DOS_E_list(iE),tot_E(iE,i_eta)
                    enddo
            close(101)
    enddo

contains

    subroutine int2char(in_int,out_char)
            integer,intent(in)              :: in_int
            character(len=20),intent(out)    :: out_char
            IF (in_int .LT. 10 )    WRITE(out_char,'(I1)')in_int
            IF ((in_int .LT. 100)  .and. (in_int .GE. 10)   )    WRITE(out_char,'(I2)')in_int
            IF ((in_int .LT. 1000) .and. (in_int .GE. 100)  )    WRITE(out_char,'(I3)')in_int
            IF ((in_int .LT. 10000) .and. (in_int .GE. 1000) )   WRITE(out_char,'(I4)')in_int
            IF ((in_int .LT. 100000) .and. (in_int .GE. 10000) ) WRITE(out_char,'(I5)')in_int
    end subroutine int2char

    subroutine diasym_eigvals(a, eig, n, work, lwork, rwork, lrwork, iwork, liwork)
        implicit none
        integer,intent(in)  :: n, lwork, lrwork, liwork
        complex(kind(0d0)),intent(inout)  :: a(n,n)
        real(kind(0d0)),intent(out)       :: eig(n)
        complex(kind(0d0)) :: work(lwork)
        real(kind(0d0))    :: rwork(lrwork)
        integer            :: iwork(liwork)
        integer :: info

        call zheevd('N', 'U', n, a, n, eig, work, lwork, rwork, lrwork, iwork, liwork, info)
        if (info /= 0) stop "Error in zheevd"
    end subroutine diasym_eigvals

end program
