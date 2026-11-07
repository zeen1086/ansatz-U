!
! Copyright (C) 2021 Quantum ESPRESSO Foundation
! Distributed under the GNU General Public License; see License.
!
PROGRAM post_processing_example
  USE kinds, ONLY : DP
  USE io_global, ONLY : ionode, ionode_id, stdout
  USE io_files, ONLY : tmp_dir, prefix
  USE mp_global, ONLY : mp_startup
  USE mp_images, ONLY : intra_image_comm, nimage
  USE mp_pools, ONLY : npool
  USE mp_bands, ONLY : nbgrp
  USE mp, ONLY : mp_bcast
  USE environment, ONLY : environment_start, environment_end
  USE fft_base, ONLY : dfftp, dffts
  USE fft_types, ONLY : fft_type_descriptor
  USE scatter_mod, ONLY : gather_grid
  USE scf, ONLY : vrs, vltot, v, kedtau
  USE gvecs, ONLY : doublegrid
  USE lsda_mod, ONLY : nspin
  USE cell_base, ONLY : at, alat
  USE klist, ONLY : xk, nkstot
  USE, INTRINSIC :: ieee_arithmetic, ONLY : ieee_is_finite
  IMPLICIT NONE
  LOGICAL :: needwf=.FALSE., write_v_qe=.FALSE., exists
  INTEGER :: ios, iu, mode, is, ix, iy, iz, ir, ip, ik, ngrid, npadded
  CHARACTER(LEN=256) :: outdir
  CHARACTER(LEN=16) :: export_grid='smooth'
  CHARACTER(LEN=7) :: label
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  TYPE(fft_type_descriptor) :: grid
  REAL(DP), ALLOCATABLE :: local(:), gathered(:), field(:,:), raw(:,:)
  NAMELIST /inputpp/ outdir, prefix, write_v_qe, export_grid

  CALL mp_startup()
  CALL environment_start('HYC_EXPORT')
  IF (npool /= 1 .OR. nbgrp /= 1 .OR. nimage /= 1) &
     CALL errore('hyc_export', 'use one image, pool and band group: -ni 1 -nk 1 -nb 1', 1)
  IF (ionode) CALL input_from_file()
  prefix='pwscf'
  CALL get_environment_variable('ESPRESSO_TMPDIR',outdir)
  IF (LEN_TRIM(outdir)==0) outdir='./'
  ios=0
  IF (ionode) THEN
     READ(5,inputpp,IOSTAT=ios)
     tmp_dir=trimcheck(outdir)
  ENDIF
  CALL mp_bcast(ios,ionode_id,intra_image_comm)
  IF (ios/=0) CALL errore('hyc_export','reading inputpp',ABS(ios))
  CALL mp_bcast(tmp_dir,ionode_id,intra_image_comm)
  CALL mp_bcast(prefix,ionode_id,intra_image_comm)
  CALL mp_bcast(write_v_qe,ionode_id,intra_image_comm)
  CALL mp_bcast(export_grid,ionode_id,intra_image_comm)
  IF (export_grid/='smooth' .AND. export_grid/='dense') &
     CALL errore('hyc_export','export_grid must be smooth or dense',1)
  ! Old export requests remain valid; the old read-only/no-output branch does not.
  IF (ionode) THEN
     INQUIRE(FILE='HYC.in',EXIST=exists)
     IF (exists) THEN
        OPEN(NEWUNIT=iu,FILE='HYC.in',STATUS='old',ACTION='read',IOSTAT=ios)
        IF (ios==0) THEN
           READ(iu,'(A7,I1)',IOSTAT=ios) label,mode
           CLOSE(iu)
           IF (ios==0 .AND. mode/=2) ios=1
        ENDIF
     ENDIF
  ENDIF
  CALL mp_bcast(ios,ionode_id,intra_image_comm)
  IF (ios/=0) CALL errore('hyc_export','HYC.in must request export mode 2, or be absent',ABS(ios))
  CALL read_file_new(needwf)
  IF (export_grid=='smooth') THEN
     grid=dffts
     CALL set_vrs(vrs,vltot,v%of_r,kedtau,v%kin_r,dfftp%nnr,nspin,doublegrid)
  ELSE
     grid=dfftp
  ENDIF
  ngrid=grid%nr1*grid%nr2*grid%nr3
  npadded=grid%nr1x*grid%nr2x*grid%nr3x
  ! QE gather_grid requires a full receive buffer even on non-root ranks.
  ALLOCATE(local(grid%nnr),gathered(npadded))
  IF (ionode) THEN
     ALLOCATE(field(ngrid,4))
     field=0.0_DP
  ENDIF
  DO is=1,nspin
     IF (export_grid=='smooth') THEN
        local=vrs(1:grid%nnr,is)
     ELSE
        local=v%of_r(:,is)
        IF (nspin/=4 .OR. is==1) local=local+vltot
     ENDIF
#if defined(__MPI)
     CALL gather_grid(grid,local,gathered)
#else
     gathered=local(1:npadded)
#endif
     IF (ionode) THEN
        ir=0
        DO iz=1,grid%nr3
           DO iy=1,grid%nr2
              DO ix=1,grid%nr1
                 ir=ir+1
                 ip=ix+(iy-1)*grid%nr1x+(iz-1)*grid%nr1x*grid%nr2x
                 field(ir,is)=gathered(ip)
              ENDDO
           ENDDO
        ENDDO
     ENDIF
  ENDDO
  DEALLOCATE(local,gathered)
  IF (ionode) THEN
     IF (.NOT.ALL(ieee_is_finite(field))) CALL errore('hyc_export','non-finite potential',1)
     IF (write_v_qe) THEN
        CALL write_field('V_QE',field(:,1:nspin))
     ENDIF
     IF (nspin==2) THEN
        ! Preserve the HYC collinear convention: [V0,Bz,0,0].
        ALLOCATE(raw(ngrid,2))
        raw=field(:,1:2)
        field(:,1)=0.5_DP*(raw(:,1)+raw(:,2))
        field(:,2)=0.5_DP*(raw(:,1)-raw(:,2))
        DEALLOCATE(raw)
     ENDIF
     CALL write_field('V_HYC_dp',field)
     DEALLOCATE(field)
     OPEN(NEWUNIT=iu,FILE='V_HYC_dp.meta',STATUS='replace',IOSTAT=ios)
     IF (ios/=0) CALL errore('hyc_export','opening metadata',ABS(ios))
     WRITE(iu,*,IOSTAT=ios) 'HYC_QE76_LOCAL_POTENTIAL_RY ',TRIM(export_grid),nspin
     IF (ios==0) WRITE(iu,*,IOSTAT=ios) grid%nr1,grid%nr2,grid%nr3
     IF (ios==0) WRITE(iu,*,IOSTAT=ios) alat,at
     IF (ios/=0) CALL errore('hyc_export','writing metadata',ABS(ios))
     CLOSE(iu,IOSTAT=ios)
     IF (ios/=0) CALL errore('hyc_export','closing metadata',ABS(ios))
     OPEN(NEWUNIT=iu,FILE='k_list',STATUS='replace',IOSTAT=ios)
     IF (ios/=0) CALL errore('hyc_export','opening k_list',ABS(ios))
     DO ik=1,nkstot
        WRITE(iu,'(3F12.8)',IOSTAT=ios) MATMUL(TRANSPOSE(at),xk(:,ik))
        IF (ios/=0) CALL errore('hyc_export','writing k_list',ABS(ios))
     ENDDO
     CLOSE(iu,IOSTAT=ios)
     IF (ios/=0) CALL errore('hyc_export','closing k_list',ABS(ios))
     WRITE(stdout,*) 'HYC_EXPORT_DONE ',TRIM(export_grid),ngrid,nspin
  ENDIF
  CALL environment_end()
  CALL stop_pp()
CONTAINS
  SUBROUTINE write_field(filename,values)
    CHARACTER(LEN=*),INTENT(IN) :: filename
    REAL(DP),INTENT(IN) :: values(:,:)
    INTEGER :: unit,status
    OPEN(NEWUNIT=unit,FILE=filename,FORM='unformatted',STATUS='replace',ACTION='write',IOSTAT=status)
    IF (status/=0) CALL errore('hyc_export','opening '//filename,ABS(status))
    WRITE(unit,IOSTAT=status) grid%nr1,grid%nr2,grid%nr3
    IF (status==0) WRITE(unit,IOSTAT=status) values
    IF (status/=0) CALL errore('hyc_export','writing '//filename,ABS(status))
    CLOSE(unit,IOSTAT=status)
    IF (status/=0) CALL errore('hyc_export','closing '//filename,ABS(status))
  END SUBROUTINE
END PROGRAM
