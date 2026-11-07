program rearrange_HR
    use parameter_mod
    implicit none
    integer, parameter :: dp = kind(0d0)

    complex(dp), allocatable :: H_wan(:,:,:), super_H_wan(:,:,:)
    integer, allocatable :: ir_vec(:,:), super_ir_vec(:,:), deg(:)
    character(len=200) :: nscf_dir
    character(len=20) :: irept_char

    integer :: i, j, k, io
    integer :: nwan
    integer :: period_x, period_y, period_z
    integer :: irept, nrept
    integer :: max_x, max_y, max_z, min_x, min_y, min_z
    integer :: super_nx, super_ny, super_nz, super_nrept, super_size
    integer :: ix, iy, iz, jx, jy, jz
    integer :: super_x, super_y, super_z, super_rept
    integer :: iwan, jwan, H_dim, super_i, super_j

    nscf_dir=p_nscf_dir
    nwan=p_nwan
    period_x=p_period_x
    period_y=p_period_y
    period_z=p_period_z

    ! Keep the primitive translation order used by the exported HR blocks.
    irept = 0
    open(unit=101, file=trim(nscf_dir)//'/ir_vec_HYC', action='read', status='old')
    do
        read(101, *, iostat=io) i, j, k
        if (io /= 0) exit
        irept = irept + 1
    end do
    close(101)
    nrept = irept
    if (nrept < 1) error stop 'empty ir_vec_HYC'

    allocate(ir_vec(nrept, 3))
    allocate(H_wan(nwan, nwan, nrept))
    allocate(deg(nrept))

    open(unit=101, file=trim(nscf_dir)//'/ir_vec_HYC', action='read', status='old')
    do i = 1, nrept
        read(101, *) ir_vec(i, :)
    end do
    close(101)

    deg = 0
    open(unit=101, file=trim(nscf_dir)//'/H_deg_HYC', action='read', status='old')
    do
        read(101, '(2I8)', iostat=io) i, j
        if (io /= 0) exit
        if (i > 0 .and. i <= nrept) then
            deg(i) = j
        endif
    end do
    close(101)
    if (any(deg <= 0)) error stop 'invalid or incomplete H_deg_HYC'

    do irept = 1, nrept
        call int2char(irept, irept_char)
        open(unit=101, file=trim(nscf_dir)//'/hr_HYC_'//trim(irept_char), &
             action='read', status='old', form='unformatted')
        read(101) H_wan(:, :, irept)
        close(101)
    end do

    max_x = -9999; max_y = -9999; max_z = -9999
    min_x =  9999; min_y =  9999; min_z =  9999

    do irept = 1, nrept
        if (ir_vec(irept,1)+period_x-1 > max_x) max_x = ir_vec(irept,1)+period_x-1
        if (ir_vec(irept,2)+period_y-1 > max_y) max_y = ir_vec(irept,2)+period_y-1
        if (ir_vec(irept,3)+period_z-1 > max_z) max_z = ir_vec(irept,3)+period_z-1
        if (ir_vec(irept,1) < min_x) min_x = ir_vec(irept,1)
        if (ir_vec(irept,2) < min_y) min_y = ir_vec(irept,2)
        if (ir_vec(irept,3) < min_z) min_z = ir_vec(irept,3)
    end do

    max_x = max_x / period_x
    max_y = max_y / period_y
    max_z = max_z / period_z
    min_x = (min_x + 1) / period_x - 1
    min_y = (min_y + 1) / period_y - 1
    min_z = (min_z + 1) / period_z - 1

    super_nx = max_x - min_x + 1
    super_ny = max_y - min_y + 1
    super_nz = max_z - min_z + 1
    super_nrept = super_nx * super_ny * super_nz
    super_size = period_x * period_y * period_z
    H_dim = nwan * super_size

    allocate(super_ir_vec(super_nrept, 3))
    allocate(super_H_wan(H_dim, H_dim, super_nrept))

    do ix = min_x, max_x
        do iy = min_y, max_y
            do iz = min_z, max_z
                irept = 1 + (iz - min_z) + super_nz * (iy - min_y) + &
                        super_nz * super_ny * (ix - min_x)
                super_ir_vec(irept, :) = (/ix, iy, iz/)
            end do
        end do
    end do

    ! Divide by primitive Wigner-Seitz degeneracy exactly once during folding.
    super_H_wan = dcmplx(0.0_dp, 0.0_dp)

    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ix,iy,iz,iwan,i,super_i,irept,jx,jy,jz,super_x,super_y,super_z,super_rept,j,jwan,super_j)
    do ix = 0, period_x - 1
        do iy = 0, period_y - 1
            do iz = 0, period_z - 1
                do iwan = 1, nwan
                    i = 1 + iz + period_z * iy + period_z * period_y * ix
                    super_i = 1 + (iwan - 1) + nwan * (i - 1)

                    do irept = 1, nrept
                        jx = ix + ir_vec(irept, 1)
                        jy = iy + ir_vec(irept, 2)
                        jz = iz + ir_vec(irept, 3)

                        if (jx < 0) then
                            super_x = (jx + 1) / period_x - 1
                        else
                            super_x = jx / period_x
                        endif
                        if (jy < 0) then
                            super_y = (jy + 1) / period_y - 1
                        else
                            super_y = jy / period_y
                        endif
                        if (jz < 0) then
                            super_z = (jz + 1) / period_z - 1
                        else
                            super_z = jz / period_z
                        endif

                        super_rept = 1 + (super_z - min_z) + &
                                     super_nz * (super_y - min_y) + &
                                     super_nz * super_ny * (super_x - min_x)

                        jx = jx - super_x * period_x
                        jy = jy - super_y * period_y
                        jz = jz - super_z * period_z

                        j = 1 + jz + period_z * jy + period_z * period_y * jx

                        do jwan = 1, nwan
                            super_j = 1 + (jwan - 1) + nwan * (j - 1)
                            super_H_wan(super_i, super_j, super_rept) = &
                                H_wan(iwan, jwan, irept) / deg(irept)
                        end do
                    end do
                end do
            end do
        end do
    end do
    !$OMP END PARALLEL DO

    do irept = 1, super_nrept
        call int2char(irept, irept_char)
        open(unit=102, file=trim(p_super_h_dir)//'hr_HYC_'//trim(irept_char), &
             action='write', status='replace', form='unformatted')
        write(102) super_H_wan(:, :, irept)
        close(102)
    end do

    open(unit=102, file=trim(p_super_h_dir)//'ir_vec_HYC', action='write', status='replace')
    do irept = 1, super_nrept
        write(102, '(3I8)') super_ir_vec(irept, :)
    end do
    close(102)

    write(*,'(A,I0,A,I0)') 'SUPER_H_TORQUE_ONLY_PASS blocks=', super_nrept, ' H_dim=', H_dim

contains

    subroutine int2char(in_int, out_char)
        integer, intent(in) :: in_int
        character(len=20), intent(out) :: out_char
        if (in_int < 10) write(out_char, '(I1)') in_int
        if (in_int < 100 .and. in_int >= 10) write(out_char, '(I2)') in_int
        if (in_int < 1000 .and. in_int >= 100) write(out_char, '(I3)') in_int
        if (in_int < 10000 .and. in_int >= 1000) write(out_char, '(I4)') in_int
        if (in_int < 100000 .and. in_int >= 10000) write(out_char, '(I5)') in_int
    end subroutine int2char

end program rearrange_HR
