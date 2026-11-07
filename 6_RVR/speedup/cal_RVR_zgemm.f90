! Candidate RVR generator: identical inputs, support rules, index layout and
! output format as 6_RVR/cal_RVR_split.f90, with
!   (1) the iwan row selected at run time (RVR_IWAN) instead of by sed+compile,
!   (2) one packing of conjg(fin_wan) per row and of ini_wan per column,
!   (3) the per-channel grid contraction evaluated as ZGEMM over grid chunks.
! All physical/numerical definitions are unchanged; only the summation order of
! the grid contraction differs (BLAS blocking / chunk accumulation).
program cal_RVR_zgemm
  use parameter_mod
  use rvr_kernel_mod
  implicit none
  include 'mpif.h'

  integer :: iwan, jwan, nwan, nwan_sp
  integer :: nx, ny, nz, j
  integer :: ni
  complex(dp), allocatable :: read_wan(:, :, :), ini_wan(:, :, :), fin_wan(:, :, :)
  integer, allocatable :: ini_TF(:), fin_TF(:)
  integer :: convert_size_x, convert_size_y, convert_size_z
  integer :: d_ext, i_dir, j_dir, k_dir, i, i2_dir
  integer, allocatable :: dir_vec(:, :), dir_inv(:), dir_active(:)
  character(len=200) :: wan_dir, func_dir, nscf_dir, rvr_root, parameter_file
  character(len=200) :: useless_string
  character(len=20) :: iwan_char, jwan_char, ndir_char
  complex(dp) :: cI
  integer :: n_func, i_func, wan_min, wan_max, n_dep, s_comp, B_comp
  complex(dp), allocatable :: RVR(:, :, :), RVR_local(:, :, :), RVR_tmp(:, :)
  real(dp), allocatable :: Bgen(:, :, :, :)
  integer :: nspin, env_status

  ! packed workspace
  complex(dp), allocatable :: A_pack(:, :), I1(:, :), I2(:, :), C(:, :)
  complex(dp), allocatable :: Cbig(:, :)
  integer, allocatable :: idx_i(:), idx_j(:), pool(:), ichan(:)
  integer :: n_i, n_j, chunk, kernel_mode, nchan, npairs, cb, ncb, c0, b

  ! mpi
  integer :: ierr, my_rank, size, nworker, worker_id, imuti, ipool
  integer :: npool, col0

  real(8) :: t_all0, tw, t_read_pot, t_read_wan, t_pack, t_kernel, t_comm, &
             t_write, t_rank0

  cI = dcmplx(0.0d0, 1.0d0)
  n_i = 0
  n_j = 0

  ! ---------------- inputs (identical parsing to cal_RVR_split.f90) --------
func_dir = p_func_dir
nscf_dir = p_nscf_dir
rvr_root = p_rvr_dir
n_func = p_n_func
convert_size_x = p_convert_size_x
convert_size_y = p_convert_size_y
convert_size_z = p_convert_size_z
nwan = p_nwan
nspin = p_nspin
wan_min = p_wan_min
wan_max = p_wan_max

  s_comp = 3
  B_comp = 3
  if (mod(n_func, s_comp) .ne. 0) then
    print *, "n_function must be a multiple of spin components"
    stop
  endif
  n_dep = n_func/s_comp
  print *, 'Configuration: root parameter (build-time)'
  print *, 'n_func = ', n_func, ' n_dep = ', n_dep, ' nspin = ', nspin, &
           ' wan_min/max = ', wan_min, wan_max
  if (.not. (n_dep .eq. 6 .or. n_dep .eq. 12 .or. n_dep .eq. 22)) then
    print *, "unsupported dependency count n_dep=", n_dep
    stop
  endif

  ! run-time iwan row (compile once)
  call read_row_index(iwan)
  chunk = p_rvr_chunk
  kernel_mode = 1
  cb = p_rvr_channel_block

  if (nspin .eq. 1) then
    if (iwan .lt. wan_min*2 - 1 .or. iwan .gt. wan_max*2) then
      print *, 'iwan out of min,max'
      stop
    endif
  else
    if (iwan .lt. wan_min .or. iwan .gt. wan_max) then
      print *, 'iwan out of min,max'
      stop
    endif
  endif

  wan_dir = trim(nscf_dir)//'./make_wan/'
  d_ext = (2*convert_size_x + 1)*(2*convert_size_y + 1)*(2*convert_size_z + 1)
  call int2char(d_ext, ndir_char)

  allocate (dir_vec(d_ext, 3))
  allocate (dir_inv(d_ext))
  allocate (dir_active(d_ext))
  i = 1
  do i_dir = -convert_size_x, convert_size_x
    do j_dir = -convert_size_y, convert_size_y
      do k_dir = -convert_size_z, convert_size_z
        dir_vec(i, :) = (/i_dir, j_dir, k_dir/)
        i = i + 1
      end do
    end do
  end do
  do i_dir = 1, d_ext
    do j_dir = 1, d_ext
      if (all(dir_vec(i_dir, :) .eq. -dir_vec(j_dir, :))) then
        dir_inv(i_dir) = j_dir
        exit
      endif
    end do
  end do
  call load_dir_active()

  open (unit=101, file=trim(wan_dir)//'../V_HYC_dp', access='sequential', &
        form='unformatted', action='read')
  read (101) nx, ny, nz
  close (101)

  open (unit=101, file=trim(func_dir)//'./R_vec', action='read', &
        access='sequential', form='unformatted')
  read (101) ni
  close (101)

  allocate (Bgen(ni, B_comp, n_dep, s_comp))

  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
  t_all0 = MPI_WTIME()

  ! ---------------- potential (stage 1a) -----------------------------------
  tw = MPI_WTIME()
  open (unit=101, file=trim(func_dir)//'./V_func_tmp', &
        access='sequential', action='read', form='unformatted')
  read (101)   ! ave_read (not used by RVR, as in cal_RVR_split.f90)
  read (101)
  read (101)
  read (101) Bgen
  close (101)
  t_read_pot = MPI_WTIME() - tw

  if (size .gt. 1) then
    nworker = size - 1
    worker_id = my_rank
  else
    nworker = 1
    worker_id = 1
  endif

  if (nspin .eq. 1) then
    nwan_sp = 2*(wan_max - wan_min + 1)
    npool = nwan_sp
    allocate (pool(npool))
    ipool = 1
    do jwan = 2*wan_min - 1, 2*wan_max
      pool(ipool) = jwan
      ipool = ipool + 1
    end do
  else
    nwan_sp = (wan_max - wan_min + 1)
    npool = nwan_sp
    allocate (pool(npool))
    ipool = 1
    do jwan = wan_min, wan_max
      pool(ipool) = jwan
      ipool = ipool + 1
    end do
  endif

  allocate (RVR(n_func, d_ext, nwan_sp*d_ext))
  allocate (RVR_local(n_func, d_ext, nwan_sp*d_ext))
  allocate (RVR_tmp(d_ext, d_ext))
  allocate (idx_i(d_ext), idx_j(d_ext))

  t_read_wan = 0.0d0
  t_pack     = 0.0d0
  t_kernel   = 0.0d0
  t_comm     = 0.0d0
  t_write    = 0.0d0
  nchan      = 0
  npairs     = 0

  RVR_local = dcmplx(0.0d0, 0.0d0)

  if (size .eq. 1 .or. my_rank .ne. 0) then
    allocate (read_wan(d_ext, ni, nspin))
    allocate (ini_wan(d_ext, ni, 2))
    allocate (fin_wan(d_ext, ni, 2))
    allocate (fin_TF(d_ext), ini_TF(d_ext))

    if (nspin .eq. 1) then
      if (mod(iwan, 2) .eq. 0) then
        call int2char(iwan/2, iwan_char)
      else
        call int2char((iwan + 1)/2, iwan_char)
      endif
    else
      call int2char(iwan, iwan_char)
    endif

    tw = MPI_WTIME()
    open (unit=102, file=trim(wan_dir)//'./wan_data/wan_'//trim(iwan_char), &
          access='stream', form='unformatted', action='read')
    read (102) read_wan
    close (102)
    t_read_wan = t_read_wan + MPI_WTIME() - tw

    fin_wan = dcmplx(0.0d0, 0.0d0)
    if (nspin .eq. 1) then
      if (mod(iwan, 2) .eq. 0) then
        fin_wan(:, :, 2) = read_wan(:, :, 1)
      else
        fin_wan(:, :, 1) = read_wan(:, :, 1)
      endif
    else
      fin_wan = read_wan
    endif

    tw = MPI_WTIME()
    open (unit=102, file=trim(wan_dir)//'./cal_TF_'//trim(iwan_char), &
          access='sequential', action='read')
    read (102, '('//trim(ndir_char)//'I8)') fin_TF
    close (102)
    call build_support(d_ext, dir_active, dir_inv, fin_TF, fin_TF, &
                       n_i, idx_i, n_j, idx_j)
    ! leading dimension of every kernel array is fixed by n_i (row of this job)
    allocate (A_pack(2*ni, max(1, n_i)), I1(ni, d_ext), I2(ni, d_ext))
    allocate (C(max(1, n_i), d_ext), Cbig(max(1, n_i), d_ext*cb))
    allocate (ichan(cb))
    A_pack = dcmplx(0.0d0, 0.0d0)
    I1 = dcmplx(0.0d0, 0.0d0)
    I2 = dcmplx(0.0d0, 0.0d0)
    C = dcmplx(0.0d0, 0.0d0)
    Cbig = dcmplx(0.0d0, 0.0d0)
    if (n_i .gt. 0) call pack_bra(fin_wan, d_ext, ni, dir_inv, n_i, idx_i, A_pack)
    t_pack = t_pack + MPI_WTIME() - tw

    do imuti = 0, (npool - 1)/nworker
      ipool = imuti*nworker + worker_id
      if (ipool .gt. npool) cycle
      jwan = pool(ipool)
      npairs = npairs + 1

      if (nspin .eq. 1) then
        if (mod(jwan, 2) .eq. 0) then
          call int2char(jwan/2, jwan_char)
        else
          call int2char((jwan + 1)/2, jwan_char)
        endif
      else
        call int2char(jwan, jwan_char)
      endif

      tw = MPI_WTIME()
      open (unit=104, file=trim(wan_dir)//'./wan_data/wan_'//trim(jwan_char), &
            access='stream', form='unformatted', action='read')
      read (104) read_wan
      close (104)
      t_read_wan = t_read_wan + MPI_WTIME() - tw

      ini_wan = dcmplx(0.0d0, 0.0d0)
      if (nspin .eq. 1) then
        if (mod(jwan, 2) .eq. 0) then
          ini_wan(:, :, 2) = read_wan(:, :, 1)
        else
          ini_wan(:, :, 1) = read_wan(:, :, 1)
        endif
      else
        ini_wan = read_wan
      endif

      tw = MPI_WTIME()
      open (unit=102, file=trim(wan_dir)//'./cal_TF_'//trim(jwan_char), &
            access='sequential', action='read')
      read (102, '('//trim(ndir_char)//'I8)') ini_TF
      close (102)
      call build_support(d_ext, dir_active, dir_inv, fin_TF, ini_TF, &
                         n_i, idx_i, n_j, idx_j)
      if (n_j .gt. 0) call pack_ket(ini_wan, d_ext, ni, dir_inv, n_j, idx_j, I1, I2)
      t_pack = t_pack + MPI_WTIME() - tw

      if (nspin .eq. 1) then
        col0 = (jwan - (2*wan_min - 1))*d_ext + 1
      else
        col0 = (jwan - wan_min)*d_ext + 1
      endif

      tw = MPI_WTIME()
      if (kernel_mode .eq. 1) then
        do c0 = 1, n_func, cb
          ncb = min(cb, n_func - c0 + 1)
          do b = 1, ncb
            ichan(b) = c0 + b - 1
          end do
          call apply_channels(ni, n_i, n_j, A_pack, I1, I2, Bgen, n_dep, &
                              s_comp, ichan, ncb, chunk, Cbig)
          do b = 1, ncb
            i_func = ichan(b)
            RVR_tmp = dcmplx(0.0d0, 0.0d0)
            do i = 1, n_i
              do j = 1, n_j
                RVR_tmp(idx_i(i), idx_j(j)) = Cbig(i, (b - 1)*n_j + j)
              end do
            end do
            RVR_local(i_func, :, col0:col0 + d_ext - 1) = RVR_tmp
            nchan = nchan + 1
          end do
        end do
      else
        do j_dir = 1, n_dep
          do i2_dir = 1, s_comp
            i_func = 1 + (i2_dir - 1) + s_comp*(j_dir - 1)
            call apply_channel_scalar(ni, n_i, n_j, A_pack, I1, I2, Bgen, &
                                      j_dir, i2_dir, n_dep, s_comp, C)
            RVR_tmp = dcmplx(0.0d0, 0.0d0)
            do i = 1, n_i
              do j = 1, n_j
                RVR_tmp(idx_i(i), idx_j(j)) = C(i, j)
              end do
            end do
            RVR_local(i_func, :, col0:col0 + d_ext - 1) = RVR_tmp
            nchan = nchan + 1
          end do
        end do
      endif
      t_kernel = t_kernel + MPI_WTIME() - tw
    end do

    if (size .gt. 1) then
      tw = MPI_WTIME()
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call MPI_SEND(RVR_local, n_func*(nwan_sp*d_ext)*d_ext, &
                    MPI_DOUBLE_COMPLEX, 0, my_rank, MPI_COMM_WORLD, ierr)
      t_comm = MPI_WTIME() - tw
    endif
  endif

  if (my_rank .eq. 0) then
    ! Output names use spin rows; input Wannier files use orbital indices.
    call int2char(iwan, iwan_char)
    tw = MPI_WTIME()
    if (size .eq. 1) then
      RVR = RVR_local
    else
      RVR = dcmplx(0.0d0, 0.0d0)
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      do ipool = 1, size - 1
        call MPI_RECV(RVR_local, n_func*(nwan_sp*d_ext)*d_ext, &
                      MPI_DOUBLE_COMPLEX, MPI_ANY_SOURCE, ipool, &
                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        RVR = RVR + RVR_local
      end do
    endif
    t_comm = t_comm + MPI_WTIME() - tw

    tw = MPI_WTIME()
    open (unit=103, file='./tmp/B_general_'//trim(iwan_char), action='write', &
          access='append', status='replace', form='unformatted')
    write (103) RVR
    close (103)
    t_write = t_write + MPI_WTIME() - tw
  endif

  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  t_rank0 = MPI_WTIME() - t_all0
  call report_peak_rss()

  write (*, '(A,I0,A,I0,A,I0,A,ES16.8,A,ES16.8,A,ES16.8,A,ES16.8,A,ES16.8,A,ES16.8,A,ES16.8)') &
    'RVR_TIMING rank=', my_rank, ' iwan=', iwan, ' npairs=', npairs, &
    ' total=', t_rank0, &
    ' read_pot=', t_read_pot, ' read_wan=', t_read_wan, &
    ' pack=', t_pack, ' kernel=', t_kernel, ' comm=', t_comm, &
    ' write=', t_write
  if (kernel_mode .eq. 1) write (*, '(A,ES16.8,A,ES16.8)') &
    'RVR_KSPLIT build=', rvr_kernel_build_time(), &
    ' gemm=', rvr_kernel_gemm_time()
  write (*, '(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0)') 'RVR_SHAPE rank=', my_rank, &
    ' cb=', cb, ' chunk=', chunk, &
    ' n_i=', n_i, ' n_j=', n_j, ' nchan=', nchan

  call MPI_FINALIZE(ierr)

contains

  subroutine report_peak_rss()
    ! Linux-only diagnostic: peak RSS of this rank from /proc/self/status.
    logical :: ex
    integer :: u, ios
    character(len=256) :: line
    integer :: peak_kb
    peak_kb = 0
    inquire (file='/proc/self/status', exist=ex)
    if (.not. ex) then
      write (*, '(A,I0,A,I0)') 'RVR_MEM rank=', my_rank, ' peak_kb=', 0
      return
    endif
    open (newunit=u, file='/proc/self/status', status='old', action='read', iostat=ios)
    if (ios /= 0) return
    do
      read (u, '(A)', iostat=ios) line
      if (ios /= 0) exit
      if (index(line, 'VmHWM:') == 1) then
        read (line(7:), *, iostat=ios) peak_kb
        exit
      endif
    end do
    close (u)
    write (*, '(A,I0,A,I0)') 'RVR_MEM rank=', my_rank, ' peak_kb=', peak_kb
  end subroutine report_peak_rss

  subroutine read_int_env(name, val)
    character(len=*), intent(in)    :: name
    integer, intent(inout)          :: val
    character(len=64) :: buf
    integer :: l, st, v
    buf = ' '
    call get_environment_variable(name, buf, length=l, status=st)
    if (st .ne. 0 .or. l .eq. 0) return
    read (buf, *, iostat=st) v
    if (st .eq. 0) val = v
  end subroutine read_int_env

  subroutine int2char(in_int, out_char)
    integer, intent(in)             :: in_int
    character(len=20), intent(out)  :: out_char
    IF (in_int .LT. 10) WRITE (out_char, '(I1)') in_int
    IF ((in_int .LT. 100) .and. (in_int .GE. 10)) WRITE (out_char, '(I2)') in_int
    IF ((in_int .LT. 1000) .and. (in_int .GE. 100)) WRITE (out_char, '(I3)') in_int
    IF ((in_int .LT. 10000) .and. (in_int .GE. 1000)) WRITE (out_char, '(I4)') in_int
    IF ((in_int .LT. 100000) .and. (in_int .GE. 10000)) WRITE (out_char, '(I5)') in_int
  end subroutine int2char

  subroutine load_dir_active()
    implicit none
    logical :: file_exists
    integer :: i_dir
    dir_active(:) = 1
    inquire (file='dir_active.in', exist=file_exists)
    if (file_exists) then
      open (unit=778, file='dir_active.in', status='old', action='read')
      read (778, *) dir_active(:)
      close (778)
    endif
    do i_dir = 1, d_ext
      if (.not. (dir_active(i_dir) .eq. 0 .or. dir_active(i_dir) .eq. 1)) then
        stop 'dir_active.in entries must be 0 or 1'
      endif
      if (dir_active(i_dir) .ne. dir_active(dir_inv(i_dir))) then
        stop 'dir_active.in must be inversion symmetric'
      endif
    end do
    print *, 'dir_active active/inactive = ', sum(dir_active(:)), d_ext - sum(dir_active(:))
  end subroutine load_dir_active

end program cal_RVR_zgemm
