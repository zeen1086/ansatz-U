PROGRAM export_labels
  USE kinds, ONLY: DP
  USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
  USE io_global, ONLY: ionode, ionode_id, stdout
  USE io_files, ONLY: tmp_dir, prefix
  USE mp_global, ONLY: mp_startup
  USE mp_images, ONLY: intra_image_comm, nproc_image
  USE mp, ONLY: mp_bcast
  USE environment, ONLY: environment_start, environment_end
  USE fft_base, ONLY: dfftp
  USE fft_rho, ONLY: rho_r2g
  USE lsda_mod, ONLY: nspin
  USE scf, ONLY: rho, rho_core, rhog_core, tau_core, v, vltot
  USE noncollin_module, ONLY: noncolin, lspinorb, i_cons, pointlist, factlist, lsign, ux
  USE ldaU, ONLY: eth
  USE ions_base, ONLY: nat
  USE cell_base, ONLY: alat, at
  IMPLICIT NONE

  INTERFACE
    SUBROUTINE v_of_rho(rho,rho_core,rhog_core,tau_core,ehart,etxc,vtxc,eth,etotefield,charge,v)
      USE kinds, ONLY: DP
      USE scf, ONLY: scf_type
      USE fft_base, ONLY: dfftp
      USE gvect, ONLY: ngm
      TYPE(scf_type), INTENT(INOUT) :: rho,v
      REAL(DP), INTENT(IN) :: rho_core(dfftp%nnr),tau_core(dfftp%nnr)
      COMPLEX(DP), INTENT(IN) :: rhog_core(ngm)
      REAL(DP), INTENT(OUT) :: ehart,etxc,vtxc,eth,charge
      REAL(DP), INTENT(INOUT) :: etotefield
    END SUBROUTINE
  END INTERFACE

  CHARACTER(LEN=256) :: outdir
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  LOGICAL :: source_lsign, needwf
  REAL(DP) :: source_ux(3), ehart, etxc, vtxc, etotefield, charge
  REAL(DP), ALLOCATABLE :: core_r(:,:), moments(:,:), charges(:), fields(:,:)
  COMPLEX(DP), ALLOCATABLE :: core_g(:,:)
  INTEGER :: ios, unit, ia, ix, iy, iz, source, point, npoint
  NAMELIST /inputpp/ outdir, prefix, source_lsign, source_ux

  CALL mp_startup()
  CALL environment_start('V1_LABEL_EXPORT')
  IF (nproc_image /= 1) CALL errore('export_labels','requires one MPI rank',1)
  IF (ionode) CALL input_from_file()
  outdir = './'
  prefix = 'pwscf'
  source_lsign = .FALSE.
  source_ux = 0._DP
  ios = 0
  IF (ionode) READ(5,inputpp,IOSTAT=ios)
  CALL mp_bcast(ios,ionode_id,intra_image_comm)
  IF (ios /= 0) CALL errore('export_labels','reading inputpp',ABS(ios))
  CALL mp_bcast(outdir,ionode_id,intra_image_comm)
  CALL mp_bcast(prefix,ionode_id,intra_image_comm)
  CALL mp_bcast(source_lsign,ionode_id,intra_image_comm)
  CALL mp_bcast(source_ux,ionode_id,intra_image_comm)
  IF (.NOT.ALL(ieee_is_finite(source_ux))) CALL errore('export_labels','nonfinite GGA axis',1)
  IF (source_lsign) THEN
    IF (ABS(SUM(source_ux**2)-1._DP)>1.e-10_DP) CALL errore('export_labels','GGA axis must be unit',1)
  ELSE
    IF (ANY(source_ux/=0._DP)) CALL errore('export_labels','unsigned GGA requires zero axis',1)
  END IF
  tmp_dir = trimcheck(outdir)
  needwf = .FALSE.
  CALL read_file_new(needwf)
  IF (nspin/=4 .OR. .NOT.noncolin .OR. lspinorb) &
    CALL errore('export_labels','requires noncollinear no-SOC SCF',1)
  IF (dfftp%nnr /= dfftp%nr1x*dfftp%nr2x*dfftp%nr3x) &
    CALL errore('export_labels','requires undistributed density',1)
  IF (MIN(dfftp%nr1,dfftp%nr2,dfftp%nr3)<=0) CALL errore('export_labels','invalid grid',1)

  CALL rho_r2g(dfftp,rho%of_r,rho%of_g)
  ALLOCATE(core_r(dfftp%nnr,1),core_g(SIZE(rhog_core),1))
  core_r(:,1) = rho_core
  CALL rho_r2g(dfftp,core_r,core_g)
  rhog_core = core_g(:,1)
  DEALLOCATE(core_r,core_g)
  IF (.NOT.ALLOCATED(pointlist)) ALLOCATE(pointlist(dfftp%nnr))
  IF (.NOT.ALLOCATED(factlist)) ALLOCATE(factlist(dfftp%nnr))
  CALL make_pointlists()
  ALLOCATE(moments(3,nat),charges(nat))
  CALL get_locals(charges,moments,rho%of_r)

  ! read_file_new does not restore the SCF GGA axis; the caller supplies its recorded branch.
  lsign = source_lsign
  ux = source_ux
  i_cons = 0
  etotefield = 0._DP
  CALL v_of_rho(rho,rho_core,rhog_core,tau_core,ehart,etxc,vtxc,eth,etotefield,charge,v)
  v%of_r(:,1) = v%of_r(:,1)+vltot
  npoint = dfftp%nr1*dfftp%nr2*dfftp%nr3
  ALLOCATE(fields(npoint,4))
  point = 0
  DO iz=1,dfftp%nr3
    DO iy=1,dfftp%nr2
      DO ix=1,dfftp%nr1
        point = point+1
        source = ix+dfftp%nr1x*((iy-1)+dfftp%nr2x*(iz-1))
        fields(point,:) = v%of_r(source,:)
      END DO
    END DO
  END DO
  IF (.NOT.ALL(ieee_is_finite(fields))) CALL errore('export_labels','nonfinite fields',1)
  IF (.NOT.ALL(ieee_is_finite(moments)) .OR. .NOT.ALL(ieee_is_finite(charges))) &
    CALL errore('export_labels','nonfinite local moments',1)

  OPEN(NEWUNIT=unit,FILE='V_LABELS.bin',STATUS='new',FORM='unformatted')
  WRITE(unit) dfftp%nr1,dfftp%nr2,dfftp%nr3
  WRITE(unit) fields
  CLOSE(unit)
  OPEN(NEWUNIT=unit,FILE='moments.dat',STATUS='new')
  DO ia=1,nat
    WRITE(unit,'(I8,4ES26.17)') ia,moments(:,ia),charges(ia)
  END DO
  CLOSE(unit)
  OPEN(NEWUNIT=unit,FILE='LABELS.meta',STATUS='new')
  WRITE(unit,'(A)') 'V1_LABEL_FIELDS_RY'
  WRITE(unit,'(3I10)') dfftp%nr1,dfftp%nr2,dfftp%nr3
  WRITE(unit,'(I2,3ES26.17)') MERGE(1,0,source_lsign),source_ux
  WRITE(unit,'(9ES26.17)') at*alat
  CLOSE(unit)
  WRITE(stdout,'(A)') 'V1_LABEL_EXPORT_DONE'
  CALL environment_end('V1_LABEL_EXPORT')
  CALL stop_pp()
END PROGRAM export_labels
