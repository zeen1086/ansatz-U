program spin
  use parameter_mod
  use spin_channels_mod, only: spin_channels, spin_chirality
        implicit none
        integer :: ix,iy,iz,nx,ny,nz,i,is,j,io,js
        integer :: period_x,period_y,period_z,jx,jy,jz
        integer :: ini_x,ini_y,ini_z,fin_x,fin_y,fin_z
        integer :: iRVR,jRVR,H_dim,i_dim,j_dim,iwan,jwan
        integer :: super_x,super_y,super_z,jrept
        integer,allocatable :: cell(:,:),super_ir_vec(:,:)
        real(kind(0d0)),allocatable        :: x_list(:,:,:),coefficient(:)
        real(kind(0d0)),allocatable        :: spin_vec(:,:,:),off_set(:,:)
        real(kind(0d0))         :: A2B(3),r_tmp(3),Mz
        real(kind(0d0))    :: rad,pi,x,y_com,x_com,theta
        real(kind(0d0))    :: a(3,3),a0 ,avB(3),avA(3),rrBdr(3),rrAdr(3)
        real(kind(0d0)) :: a_next,b_next,s_abs,celldm3
        real(kind(0d0)) ::vec6(6,3),vec3(3),axial
        integer :: int3(3),super_nx,super_ny,super_nz,n_func
        integer :: convert_size_x, convert_size_y,convert_size_z
        integer :: irept,nrept,super_nrept,i_locat,j_locat
        integer :: max_x,max_y,max_z,min_x,min_y,min_z,super_i,super_j,d_ext
        integer :: nspin, nwan_sp, nwan, wan_max,wan_min,i_dir,j_dir,k_dir
        integer :: env_status
        integer,allocatable     :: dir_vec(:,:),inv_dir(:)
        character(len=200)      :: useless_string,func_dir,nscf_dir,rvr_root
        character(len=200)      :: unused_dir_override
        character(len=200)      :: RVR_dir
        character(len=20)       :: rept_char
        complex(kind(0d0)),allocatable  :: RVR(:,:,:),RVR_POT_ave(:,:),RVR_tmp(:,:)
        complex(kind(0d0)),allocatable  :: H_sk(:,:,:)
        complex(kind(0d0))      :: cI

        real(kind(0d0)) :: q_vec(3)
        real(kind(0d0)) :: r_pos(3), phase_q

        integer :: channel_mask(6), channel_source(6)
        real(8) :: channel_chi(6)
        integer :: src, base_index, mod_index
        real(kind(0d0)) :: lambda_high, chi_mod, mod_factor
        logical :: channel_file_exists, source_file_exists, lambda_file_exists

 call EXECUTE_COMMAND_LINE("mkdir -p  data")

func_dir = p_func_dir
nscf_dir = p_nscf_dir
rvr_root = p_rvr_dir
n_func = p_n_func
convert_size_x = p_convert_size_x
convert_size_y = p_convert_size_y
convert_size_z = p_convert_size_z
nwan = p_nwan
nspin = p_nspin
nx = p_nx
ny = p_ny
nz = p_nz
axial = p_axial
period_x = p_period_x
period_y = p_period_y
period_z = p_period_z
wan_min = p_wan_min
wan_max = p_wan_max
a0 = p_a0
celldm3 = p_celldm3

if (len_trim(rvr_root) .gt. 0) then
        if (rvr_root(len_trim(rvr_root):len_trim(rvr_root)) .ne. '/') then
                rvr_root=trim(rvr_root)//'/'
        endif
endif

call make_spin()
call load_super_size()
if  (nspin .eq. 1) then
        nwan_sp=2*(wan_max-wan_min+1)
else
        nwan_sp=(wan_max-wan_min+1)
endif

H_dim=nwan*2/nspin*period_x*period_y*period_z
allocate(H_sk(H_dim,H_dim,super_nrept))

d_ext=(2*convert_size_x+1)*(2*convert_size_y+1)*(2*convert_size_z+1)
RVR_dir=p_rvr_dir

call read_RVR()

allocate(coefficient(n_func))

channel_mask = p_channel_mask
channel_source = p_channel_source
lambda_high = p_lambda_high

do j=1,6
        if (channel_mask(j) .lt. -1 .or. channel_mask(j) .gt. 1) then
                stop 'parameter channel_mask entries must be -1, 0, or 1'
        endif
        if (channel_source(j) .lt. 1 .or. channel_source(j) .gt. 6) then
                stop 'parameter channel_source entries must be 1..6'
        endif
enddo
print *, 'channel_mask = ', channel_mask(:)
print *, 'channel_source = ', channel_source(:)
if (n_func .eq. 18) then
        print *, 'lambda_high = ', lambda_high
elseif (n_func .eq. 36) then
        print *, 'n_func = 36: using fitted K_mod channels; lambda_high is ignored.'
else
        print *, 'unsupported n_func = ', n_func
        stop
endif

H_sk(:,:,:)=dcmplx(0.0d0,0.0d0)

do ix=0,nx-1
do iy=0,ny-1
do iz=0,nz-1
        i=1+iz+nz*iy+ny*nz*ix
        coefficient(:)=0.0d0
        do i_dir=1,d_ext
                ini_x=cell(i,1)+dir_vec(i_dir,1)
                ini_y=cell(i,2)+dir_vec(i_dir,2)
                ini_z=cell(i,3)+dir_vec(i_dir,3)
                if (.not. (ini_x .ge.0 .and. ini_x/period_x .eq. 0)) cycle
                if (.not. (ini_y .ge.0 .and. ini_y/period_y .eq. 0)) cycle
                if (.not. (ini_z .ge.0 .and. ini_z/period_z .eq. 0)) cycle
                i_locat=1+ini_z+period_z*ini_y+period_z*period_y*ini_x

                call spin_chirality(spin_vec(i,:,:), channel_chi)
                call spin_channels(spin_vec(i,:,:), channel_chi, p_channel_mask, &
                     p_channel_source, p_lambda_high, coefficient)

                do j_dir=1,nwan_sp*d_ext
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

                do irept=1,super_nrept
                        if (all((/super_x,super_y,super_z/) .eq. (super_ir_vec(irept,:)))) then
                                jrept=irept
                                exit
                        endif
                enddo

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
                i_dim=(is+2/nspin*(iwan-1))+nwan*2/nspin*(i_locat-1)
                j_dim=(js+2/nspin*(jwan-1))+nwan*2/nspin*(j_locat-1)

                if (abs(RVR_tmp(iRVR,jRVR)) .lt. 10.0**(-12)) cycle
                H_sk(i_dim,j_dim,jrept)=H_sk(i_dim,j_dim,jrept)+RVR_tmp(iRVR,jRVR)

        enddo
        enddo
        enddo
        enddo

        enddo!j_dir
        enddo!i_dir
enddo!iz
enddo!iy
enddo!ix

!

do irept=1,super_nrept

                do super_i=0,0

                enddo

        call int2char(irept,rept_char)
        open(unit=102,file='./data/Hsk_'//trim(rept_char),action='write',access='append'&
                &,status='replace',form='unformatted')
        write(102)      H_sk(:,:,irept)
        close(102)
        print*,irept,sum(abs(H_sk(:,:,irept))**2),super_ir_vec(irept,:)
enddo

contains

subroutine int2char(in_int,out_char)
        integer,intent(in)              :: in_int
        character(len=20),intent(out)    :: out_char

IF (in_int .LT. 10 )   WRITE(out_char,'(I1)')in_int
IF ((in_int .LT. 100)  .and. (in_int .GE. 10)   ) WRITE(out_char,'(I2)')in_int
IF ((in_int .LT. 1000) .and. (in_int .GE. 100)  ) WRITE(out_char,'(I3)')in_int
IF ((in_int .LT. 10000) .and. (in_int .GE. 1000) ) WRITE(out_char,'(I4)')in_int
IF ((in_int .LT. 100000) .and. (in_int .GE. 10000) )WRITE(out_char,'(I5)')in_int

end subroutine int2char

subroutine make_spin()
        implicit none
        real(kind(0d0)), allocatable :: theta_nu(:)
        real(kind(0d0)), allocatable :: phi_nu(:)
        real(kind(0d0)), allocatable :: physical_spin(:,:)
        real(kind(0d0)) :: mag_moment
        real(kind(0d0)) :: Mz, s_abs
        real(kind(0d0)) :: deg2rad, pi

        integer :: io, int3(3)
        integer :: i, j, i_spin
        integer :: n_dep, n_input
        integer :: map_idx(6)

        real(kind(0d0)) :: a_super(3,3)
        logical :: file_exists

        n_dep = 6
        n_input = 2 * period_x * period_y * period_z

        allocate(theta_nu(n_input))
        allocate(phi_nu(n_input))
        allocate(physical_spin(n_input,3))

        inquire(file='spin.in', exist=file_exists)
        if (.not. file_exists) then
                print *, "Error: 'spin.in' not found!"
                stop
        endif

        open(unit=201, file='spin.in', status='old', action='read')
        read(201, *) theta_nu(:)
        read(201, *) phi_nu(:)
        close(201)

        pi = 3.1415926535d0
        deg2rad = pi / 180.0d0
        theta_nu(:) = theta_nu(:) * deg2rad
        phi_nu(:)  = phi_nu(:) * deg2rad
        mag_moment = p_cr_magnetization

        do i_spin=1,n_input
                physical_spin(i_spin,:) = (/ mag_moment * sin(theta_nu(i_spin))*cos(phi_nu(i_spin)), &
                                             mag_moment * sin(theta_nu(i_spin))*sin(phi_nu(i_spin)), &
                                             mag_moment * cos(theta_nu(i_spin)) /)

                s_abs = sum(physical_spin(i_spin,:)**2)**0.5d0
                physical_spin(i_spin,:) = physical_spin(i_spin,:) / s_abs * mag_moment
        enddo

        max_x=-99; max_y=-99; max_z=-99
        min_x=99;  min_y=99;  min_z=99

        open(unit=102,file=trim(p_super_h_dir)//'ir_vec_HYC',action='read',access='sequential',status='old')
        do
                read(102,'(3I8)',iostat=io) int3(:)
                if (io/=0) exit
                if (int3(1) .gt. max_x) max_x=int3(1)
                if (int3(2) .gt. max_y) max_y=int3(2)
                if (int3(3) .gt. max_z) max_z=int3(3)
                if (int3(1) .lt. min_x) min_x=int3(1)
                if (int3(2) .lt. min_y) min_y=int3(2)
                if (int3(3) .lt. min_z) min_z=int3(3)
        enddo
        close(102)

        nx=(max_x-min_x+1)*period_x
        ny=(max_y-min_y+1)*period_y
        nz=(max_z-min_z+1)*period_z

        print*, "Grid size:", nx, ny, nz

        allocate(cell(nx*ny*nz,3))
        allocate(spin_vec(nx*ny*nz,3,n_dep))
        allocate(off_set(n_dep,3))

        a(1,:) = (/  1.0d0,       0.0d0,       0.0d0    /)
        a(2,:) = (/ -0.5d0,       0.866025d0,  0.0d0    /)
        a(3,:) = (/  0.0d0,       0.0d0,       celldm3  /)

        a_super = a

        off_set(:,:)=0.0d0
        off_set(1,:) = (/ 1.0d0/3.0d0, 2.0d0/3.0d0 + 1.0d0, 0.5d0 /)
        off_set(2,:) = (/ 1.0d0/3.0d0, 2.0d0/3.0d0,         0.5d0 /)
        off_set(3,:) = (/ 2.0d0/3.0d0, 1.0d0/3.0d0,         0.5d0 /)
        off_set(4,:) = off_set(2,:)
        off_set(5,:) = off_set(1,:)
        off_set(6,:) = (/ 2.0d0/3.0d0, 1.0d0/3.0d0 + 1.0d0, 0.5d0 /)

        Mz = 0.0d0

        do super_x=min_x,max_x
        do super_y=min_y,max_y
        do super_z=min_z,max_z
        do ix=0,period_x-1
        do iy=0,period_y-1
        do iz=0,period_z-1

                i=1+(iz+(super_z-min_z)*period_z)+nz*(iy+(super_y-min_y)*period_y)&
                  &+nz*ny*(ix+(super_x-min_x)*period_x)

                cell(i,1) = dble(ix + super_x * period_x)
                cell(i,2) = dble(iy + super_y * period_y)
                cell(i,3) = dble(iz + super_z * period_z)

                ! Channel order:
                ! Cr2, Cr1, y-neighbor Cr1, x-neighbor Cr1,
                ! x-neighbor Cr2, y-neighbor Cr2.
                map_idx(1)=spin_input_index(ix,       iy,       iz,2)
                map_idx(2)=spin_input_index(ix,       iy,       iz,1)
                map_idx(3)=spin_input_index(ix,       iy+1,     iz,1)
                map_idx(4)=spin_input_index(ix+1,     iy,       iz,1)
                map_idx(5)=spin_input_index(ix+1,     iy,       iz,2)
                map_idx(6)=spin_input_index(ix,       iy+1,     iz,2)

                do j=1, n_dep
                        spin_vec(i,:,j)=physical_spin(map_idx(j),:)
                        Mz = Mz + spin_vec(i,3,j)
                enddo

        enddo
        enddo
        enddo
        enddo
        enddo
        enddo

        Mz = Mz / (nx*ny*nz)
        print *, "Average Mz:", Mz / n_dep

        open(unit=101,file='cmap',action='write',access='append',status='replace')
        do i=1,nx*ny*nz
                if (.not. cell(i,3) .eq. 0) cycle
                do j=1,n_dep
                        write(101,'(6F24.6)') &
                        sum(a(:,1)*(cell(i,:)+off_set(j,:)))*15.53313653*0.529177-spin_vec(i,1,j), spin_vec(i,1,j)*2, &
                        sum(a(:,2)*(cell(i,:)+off_set(j,:)))*15.53313653*0.529177-spin_vec(i,2,j), spin_vec(i,2,j)*2, &
                        sum(a(:,3)*(cell(i,:)+off_set(j,:)))*15.53313653*0.529177-spin_vec(i,3,j), spin_vec(i,3,j)*2
                enddo
        enddo
        close(101)

end subroutine

integer function spin_input_index(ix, iy, iz, atom)
        implicit none
        integer, intent(in) :: ix, iy, iz, atom
        integer :: ix0, iy0, iz0

        ix0 = modulo(ix, period_x)
        iy0 = modulo(iy, period_y)
        iz0 = modulo(iz, period_z)

        spin_input_index = atom + 2 * (ix0 + period_x * (iy0 + period_y * iz0))
end function spin_input_index

subroutine load_super_size()
        irept=0
        open(unit=102,file=trim(p_super_h_dir)//'ir_vec_HYC',action='read',access='sequential',status='old')
        do
        read(102,'(3I8)',iostat=io) int3(:)
        if (io/=0) exit
        irept=irept+1
        enddo
        close(102)
        super_nrept=irept
        allocate(super_ir_vec(super_nrept,3))

        irept=1
        open(unit=102,file=trim(p_super_h_dir)//'ir_vec_HYC',action='read',access='sequential',status='old')
        do
        read(102,'(3I8)',iostat=io) int3(:)
        if (io/=0) exit
                super_ir_vec(irept,:)=int3(:)
                irept=irept+1
        enddo
        close(102)

end subroutine

subroutine read_RVR()

        allocate(RVR(n_func,nwan_sp*d_ext,nwan_sp*d_ext))
        allocate(RVR_POT_ave(nwan_sp*d_ext,nwan_sp*d_ext))
        allocate(RVR_tmp(nwan_sp*d_ext,nwan_sp*d_ext))

        open(unit=103,file=trim(RVR_dir)//'./B_general',action='read',access='sequential',form='unformatted')
                read(103) RVR
        close(103)

        RVR_POT_ave(:,:)=0.0d0

        allocate(dir_vec((2*convert_size_x+1)*(2*convert_size_y+1)*(2*convert_size_z+1),3))
        allocate(inv_dir((2*convert_size_x+1)*(2*convert_size_y+1)*(2*convert_size_z+1)))
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

        do i=1,(2*convert_size_x+1)*(2*convert_size_y+1)*(2*convert_size_z+1)
        do j=1,(2*convert_size_x+1)*(2*convert_size_y+1)*(2*convert_size_z+1)
                if (all( dir_vec(i,:)+dir_vec(j,:) .eq. (/0,0,0/))) then
                        inv_dir(i)=j
                        exit
                endif
        enddo
                if (inv_dir(i) .eq. 0) then
                        print*, i, "-dir inverse not found"
                        stop
                endif
        enddo

end subroutine read_RVR

end program
