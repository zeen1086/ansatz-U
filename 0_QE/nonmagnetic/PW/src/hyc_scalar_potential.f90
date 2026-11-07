! HYC fixed scalar Hamiltonian: canonical V_HYC_dp, no empirical shifts.
MODULE hyc_scalar_potential
  USE kinds, ONLY : DP
  USE, INTRINSIC :: iso_fortran_env, ONLY : int64
  USE, INTRINSIC :: ieee_arithmetic, ONLY : ieee_is_finite
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: hyc_apply_scalar, hyc_clear_scalar
  REAL(DP), ALLOCATABLE, SAVE :: scalar(:)
  INTEGER, SAVE :: saved_grid(3)=0
CONTAINS
  SUBROUTINE hyc_apply_scalar(vrs,nrxx,nspin)
    USE fft_base, ONLY : dffts
    USE io_global, ONLY : ionode, stdout
    USE mp_images, ONLY : nimage
    USE mp_pools, ONLY : npool
    USE mp_bands, ONLY : nbgrp
    USE scatter_mod, ONLY : scatter_grid
    USE ldaU, ONLY : lda_plus_u
    USE uspp, ONLY : okvan
    USE paw_variables, ONLY : okpaw
    USE xc_lib, ONLY : xclib_dft_is
    USE noncollin_module, ONLY : i_cons
    USE control_flags, ONLY : lbfgs, lmd, use_gpu
    USE cellmd, ONLY : lmovecell
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nrxx,nspin
    REAL(DP), INTENT(OUT) :: vrs(nrxx,nspin)
    REAL(DP), ALLOCATABLE :: packed(:),field(:,:)
    INTEGER :: dims(3),ix,iy,iz,ir,ip
    dims=[dffts%nr1,dffts%nr2,dffts%nr3]
    IF (nspin/=1 .AND. nspin/=2 .AND. nspin/=4) CALL errore('hyc_no_mag','unsupported nspin',1)
    IF (nrxx<dffts%nnr) CALL errore('hyc_no_mag','local grid buffer too small',1)
    IF (.NOT.ALLOCATED(scalar)) THEN
       IF (nimage/=1 .OR. npool/=1 .OR. nbgrp/=1) &
          CALL errore('hyc_no_mag','use -ni 1 -nk 1 -nb 1 (MPI FFT remains enabled)',1)
       IF (use_gpu) CALL errore('hyc_no_mag','CPU MPI build required',1)
       IF (lda_plus_u .OR. xclib_dft_is('meta') .OR. xclib_dft_is('hybrid')) &
          CALL errore('hyc_no_mag','local field alone does not define Hubbard/meta/hybrid operators',1)
       IF (okvan .OR. okpaw) CALL errore('hyc_no_mag','norm-conserving pseudopotentials required',1)
       IF (i_cons/=0) CALL errore('hyc_no_mag','magnetic constraints/external B are incompatible with scalar mode',1)
       IF (lbfgs .OR. lmd .OR. lmovecell) &
          CALL errore('hyc_no_mag','fixed geometry only; use the magnetic build for relaxation',1)
       saved_grid=dims
       ! QE scatter_grid checks the full input buffer size on every rank.
       ALLOCATE(packed(dffts%nr1x*dffts%nr2x*dffts%nr3x),scalar(dffts%nnr))
       packed=0.0_DP
       IF (ionode) THEN
          CALL read_field(dims,field)
          ir=0
          DO iz=1,dims(3)
             DO iy=1,dims(2)
                DO ix=1,dims(1)
                   ir=ir+1
                   ip=ix+(iy-1)*dffts%nr1x+(iz-1)*dffts%nr1x*dffts%nr2x
                   packed(ip)=field(ir,1)
                ENDDO
             ENDDO
          ENDDO
          DEALLOCATE(field)
       ENDIF
#if defined(__MPI)
       CALL scatter_grid(dffts,packed,scalar)
#else
       scalar=packed(1:dffts%nnr)
#endif
       DEALLOCATE(packed)
       IF (ionode) THEN
          WRITE(stdout,'(5x,a,3i8,a,i2)') 'HYC_SCALAR_LOADED V_HYC_dp smooth grid=',dims,' nspin=',nspin
          WRITE(stdout,'(5x,a)') 'HYC fixed-potential interface: total energies/forces are not physical SCF observables.'
       ENDIF
    ELSE
       IF (ANY(dims/=saved_grid) .OR. SIZE(scalar)/=dffts%nnr) &
          CALL errore('hyc_no_mag','FFT layout changed while scalar potential is cached',1)
    ENDIF
    vrs=0.0_DP
    vrs(1:dffts%nnr,1)=scalar
    IF (nspin==2) vrs(1:dffts%nnr,2)=scalar
    ! nspin=4: Bx, By, Bz are exactly zero; nonlocal SOC remains native QE.
  END SUBROUTINE

  SUBROUTINE read_field(dims,field)
    USE cell_base, ONLY : alat,at
    USE io_global, ONLY : stdout
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: dims(3)
    REAL(DP), ALLOCATABLE, INTENT(OUT) :: field(:,:)
    INTEGER :: iu,ios,header(3),source_spin
    INTEGER(int64) :: bytes,expected,n
    LOGICAL :: exists
    CHARACTER(80) :: magic,grid
    REAL(DP) :: source_alat,source_at(3,3)
    n=PRODUCT(INT(dims,int64))
    expected=28_int64+32_int64*n
    INQUIRE(FILE='V_HYC_dp',EXIST=exists,SIZE=bytes,IOSTAT=ios)
    IF (ios/=0 .OR. .NOT.exists) CALL errore('hyc_no_mag','missing V_HYC_dp',1)
    ! Contract: native endian, 4-byte record markers/integers, four float64 channels.
    IF (bytes/=expected .OR. 32_int64*n>HUGE(0)) &
       CALL errore('hyc_no_mag','V_HYC_dp size/layout mismatch (expected N x 4 float64)',1)
    OPEN(NEWUNIT=iu,FILE='V_HYC_dp',FORM='unformatted',STATUS='old',ACTION='read',IOSTAT=ios)
    IF (ios/=0) CALL errore('hyc_no_mag','cannot open V_HYC_dp',ABS(ios))
    READ(iu,IOSTAT=ios) header
    IF (ios/=0) CALL errore('hyc_no_mag','invalid V_HYC_dp header',ABS(ios))
    IF (ANY(header/=dims)) CALL errore('hyc_no_mag','V_HYC_dp must match the smooth FFT grid',1)
    ALLOCATE(field(INT(n),4))
    READ(iu,IOSTAT=ios) field
    IF (ios/=0) CALL errore('hyc_no_mag','invalid V_HYC_dp data record',ABS(ios))
    CLOSE(iu,IOSTAT=ios)
    IF (ios/=0) CALL errore('hyc_no_mag','cannot close V_HYC_dp',ABS(ios))
    IF (.NOT.ALL(ieee_is_finite(field))) CALL errore('hyc_no_mag','non-finite V_HYC_dp',1)
    INQUIRE(FILE='V_HYC_dp.meta',EXIST=exists)
    IF (exists) THEN
       OPEN(NEWUNIT=iu,FILE='V_HYC_dp.meta',STATUS='old',ACTION='read',IOSTAT=ios)
       IF (ios/=0) CALL errore('hyc_no_mag','cannot open metadata',ABS(ios))
       READ(iu,*,IOSTAT=ios) magic,grid,source_spin
       IF (ios==0) READ(iu,*,IOSTAT=ios) header
       IF (ios==0) READ(iu,*,IOSTAT=ios) source_alat,source_at
       CLOSE(iu)
       IF (ios/=0) CALL errore('hyc_no_mag','invalid metadata',ABS(ios))
       IF (magic/='HYC_QE76_LOCAL_POTENTIAL_RY' .OR. &
           (grid/='smooth' .AND. grid/='dense') .OR. &
           (source_spin/=1 .AND. source_spin/=2 .AND. source_spin/=4)) &
          CALL errore('hyc_no_mag','unknown metadata convention',1)
       IF (ANY(header/=dims)) CALL errore('hyc_no_mag','metadata grid mismatch',1)
       IF (.NOT.ieee_is_finite(source_alat) .OR. .NOT.ALL(ieee_is_finite(source_at))) &
          CALL errore('hyc_no_mag','non-finite metadata lattice',1)
       IF (ABS(source_alat-alat)>1.d-8*MAX(1.d0,ABS(alat)) .OR. MAXVAL(ABS(source_at-at))>1.d-8) &
          CALL errore('hyc_no_mag','metadata lattice mismatch',1)
    ELSE
       WRITE(stdout,'(5x,a)') 'HYC legacy field without metadata: cell/units require external provenance.'
    ENDIF
  END SUBROUTINE

  SUBROUTINE hyc_clear_scalar()
    IF (ALLOCATED(scalar)) DEALLOCATE(scalar)
    saved_grid=0
  END SUBROUTINE
END MODULE
