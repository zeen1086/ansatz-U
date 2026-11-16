PROGRAM fixed_density_branch
  USE read_input, ONLY: read_input_file
  USE kinds, ONLY: DP
  USE io_global, ONLY: ionode,ionode_id,stdout
  USE io_files, ONLY: tmp_dir,prefix
  USE mp_global, ONLY: mp_startup
  USE mp_images, ONLY: intra_image_comm,nproc_image
  USE mp, ONLY: mp_bcast
  USE environment, ONLY: environment_start,environment_end
  USE fft_base, ONLY: dfftp
  USE fft_rho, ONLY: rho_r2g
  USE lsda_mod, ONLY: nspin
  USE scf, ONLY: rho,rho_core,rhog_core,v,vltot,scf_type
  USE noncollin_module, ONLY: noncolin,lspinorb,i_cons,domag,lsign,ux
  USE xc_lib, ONLY: xclib_dft_is
  USE cell_base, ONLY: omega
  USE ener, ONLY: etxcc
  IMPLICIT NONE
  INTERFACE
    SUBROUTINE read_file_new(needwf)
      LOGICAL,INTENT(INOUT)::needwf
    END SUBROUTINE
    SUBROUTINE v_xc(rho,rho_core,rhog_core,etxc,vtxc,vout)
      USE kinds, ONLY: DP
      USE scf, ONLY: scf_type
      USE fft_base, ONLY: dfftp
      USE gvect, ONLY: ngm
      USE lsda_mod, ONLY: nspin
      TYPE(scf_type),INTENT(INOUT)::rho
      REAL(DP),INTENT(IN)::rho_core(dfftp%nnr)
      COMPLEX(DP),INTENT(IN)::rhog_core(ngm)
      REAL(DP),INTENT(OUT)::etxc,vtxc
      REAL(DP),INTENT(INOUT)::vout(dfftp%nnr,nspin)
    END SUBROUTINE
  END INTERFACE
  CHARACTER(512)::outdir,path
  CHARACTER(256),EXTERNAL::trimcheck
  CHARACTER(32)::tag
  CHARACTER(8)::branch
  INTEGER::ios,np,i,j,u,inp,ncase,g(3),mode,c0,c1,rate,branch_mode
  LOGICAL::needwf
  REAL(DP)::etxc,vtxc,ehart,charge,zero,axis(3)
  REAL(DP),ALLOCATABLE::original(:,:),fields(:,:),hc(:,:),cr(:,:)
  COMPLEX(DP),ALLOCATABLE::cg(:,:)
  NAMELIST/inputpp/outdir,prefix,branch_mode
  CALL mp_startup()
  CALL environment_start('GGA_BRANCH_REEXPORT')
  CALL get_command_argument(1,path)
  CALL read_input_file('PW',TRIM(path))
  CALL iosys()
  CALL setup()
  CALL init_run()
  branch_mode=1
  IF(nproc_image/=1)CALL errore('native_large','one MPI task required',1)
  IF(nspin/=4.OR..NOT.noncolin.OR.lspinorb.OR..NOT.domag.OR.xclib_dft_is('meta')) &
    CALL errore('gga_reexport','original noncollinear noSOC nonmeta contract',1)
  np=dfftp%nr1x*dfftp%nr2x*dfftp%nr3x
  IF(np/=dfftp%nnr)CALL errore('gga_reexport','serial grid',1)
  ALLOCATE(original(np,4),fields(np,4),hc(np,4),cr(np,1),cg(SIZE(rhog_core),1))
  cr(:,1)=rho_core
  CALL rho_r2g(dfftp,cr,cg)
  rhog_core=cg(:,1)
  OPEN(NEWUNIT=u,FILE='CORE_INPUT.bin',FORM='unformatted',STATUS='new')
  WRITE(u)[dfftp%nr1x,dfftp%nr2x,dfftp%nr3x]
  WRITE(u)rho_core,rhog_core
  CLOSE(u)
  OPEN(NEWUNIT=inp,FILE='input_list.txt',STATUS='old')
  OPEN(NEWUNIT=j,FILE='NATIVE_METADATA.dat',STATUS='new')
  READ(inp,*)ncase
  i_cons=0
  DO i=1,ncase
    READ(inp,*)tag,path,axis
    IF(ABS(SUM(axis**2)-1._DP)>1.e-12_DP)CALL errore('gga_reexport','source axis norm',i)
    OPEN(NEWUNIT=u,FILE=TRIM(path),FORM='unformatted',STATUS='old')
    READ(u)g
    IF(ANY(g/=[dfftp%nr1x,dfftp%nr2x,dfftp%nr3x]))CALL errore('gga_reexport','density grid',i)
    READ(u)original
    CLOSE(u)
    rho%of_r=original
    CALL rho_r2g(dfftp,rho%of_r,rho%of_g)
    hc=0._DP
    CALL v_h(rho%of_g(:,1),ehart,charge,hc)
    DO mode=branch_mode,branch_mode
      CALL system_clock(c0,rate)
      lsign=mode==1
      ux=0._DP
      branch='FALSE'
      IF(lsign)THEN
        ux=axis
        branch='TRUE'
      ENDIF
      v%of_r=0._DP
      CALL v_xc(rho,rho_core,rhog_core,etxc,vtxc,v%of_r)
      fields(:,1)=hc(:,1)+v%of_r(:,1)+vltot
      zero=SUM(fields(:,1))/np
      fields(:,1)=fields(:,1)-zero
      fields(:,2:4)=v%of_r(:,2:4)
      OPEN(NEWUNIT=u,FILE=TRIM(tag)//'_'//TRIM(branch)//'.bin',FORM='unformatted',STATUS='new')
      WRITE(u)g
      WRITE(u)fields
      CLOSE(u)
      OPEN(NEWUNIT=u,FILE=TRIM(tag)//'_ENERGY.dat',STATUS='new')
      WRITE(u,'(A)') '# EH Exc Exc_core charge vtxc integral_HXC integral_vion scalar_zero; Ry'
      WRITE(u,'(8ES26.17)') ehart,etxc,etxcc,charge,vtxc, &
        SUM(original*(hc+v%of_r))*omega/np,SUM(original(:,1)*vltot)*omega/np,zero
      CLOSE(u)
      OPEN(NEWUNIT=u,FILE=TRIM(tag)//'_HXC.bin',STATUS='new',FORM='unformatted')
      WRITE(u)g
      WRITE(u)hc+v%of_r
      CLOSE(u)
      CALL system_clock(c1)
      WRITE(j,'(A,1X,A,1X,ES26.17,3(1X,ES26.17),1X,ES26.17)') &
        TRIM(tag),TRIM(branch),zero,ux,REAL(c1-c0,DP)/rate
      WRITE(stdout,*)'GGA_STATE_COMPLETE ',TRIM(tag),' ',TRIM(branch)
    ENDDO
  ENDDO
  CLOSE(inp)
  CLOSE(j)
  WRITE(stdout,*)'GGA_BRANCH_REEXPORT_COMPLETE'
  CALL environment_end('GGA_BRANCH_REEXPORT')
  CALL stop_pp()
END PROGRAM
