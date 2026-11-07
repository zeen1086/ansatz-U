program fitting_B

        !This code fits the ansatz coefficients from the extracted potential 
        !Detail should be changed for different material and ansatz 
        !We consider a linear ansatz between V and the four Spins
        !By changing  act='plot', it can output a predicted potential 
        !with a given spin orientation to be compared using plt_compare.f90 

        use fit_mask, only: update_allowed
        implicit none

        integer :: i,ix,iy,iz,nx,ny,nz,j
        integer :: z_box
        integer :: n_file, n_kind,nR,i_file,i_kind
        integer :: i_dir,j_dir,i2_dir

        real(kind(0d0)),allocatable  :: read_VA(:,:)
        real(kind(0d0)),allocatable  :: VAs(:,:,:,:)
        real(kind(0d0)),allocatable  :: spin(:,:,:)
        real(kind(0d0)),allocatable  :: r_vec(:,:,:,:)
        real(kind(0d0)),allocatable  :: read_angle(:)

        character(len=20)       :: file_char,kind_char,formatname,spinor_char
        character(len=200)      :: path,read_angle_form
        real(kind(0d0)),allocatable  :: VA_ave(:,:),VA_abs_sum(:,:)
        real(kind(0d0)),allocatable  :: V0(:,:),Vn(:,:)
        real(kind(0d0)),allocatable  :: B0(:),BH(:,:),JH(:,:)
        real(kind(0d0)),allocatable  :: B_DM(:,:,:)
        real(kind(0d0)),allocatable  :: B0_aniso(:,:),BH_aniso(:,:,:)
        real(kind(0d0)),allocatable  :: B_general(:,:,:,:)
        real(kind(0d0)),allocatable  :: blocked(:,:,:,:)
        real(kind(0d0)) ::pi
        real(kind(0d0)) ::lambda,lambda_H,lambda_DM,lambda_B
        real(kind(0d0)),allocatable  :: DV(:,:),cal_V(:,:)
        integer :: i_iter,n_iter,i_count
        integer :: n_spinor
        real(kind(0d0)),allocatable:: Err(:),r_hat(:,:),sub_skip(:)
        real(kind(0d0)) :: skip,l_rate,wt
        character(len=4) :: act
        integer,allocatable:: one_zero(:)
        integer,allocatable:: sym_BZ(:,:),pt_TF(:)
        integer :: ni
        logical :: output_exists

        real(kind(0d0)) ::rand


  character(len=1024) :: control_file,fields_dir,spin_file,restart_file,output_file
  integer :: diagonal_iters,full_iters,ios
  namelist /fit_control/ n_spinor,n_kind,diagonal_iters,full_iters,lambda_DM, &
       fields_dir,spin_file,restart_file,output_file
  n_spinor=1; n_kind=30; diagonal_iters=0; full_iters=0; lambda_DM=-1d0
  fields_dir=''; spin_file=''; restart_file=''; output_file=''
  call get_command_argument(1,control_file)
  if(len_trim(control_file)==0) error stop 'Supply fit_control namelist path'
  open(newunit=ios,file=trim(control_file),status='old',action='read')
  read(ios,nml=fit_control); close(ios)
  if(n_spinor<1.or.n_kind<1.or.diagonal_iters<0.or.full_iters<0) error stop 'Invalid fit dimensions/stages'
  if(diagonal_iters+full_iters<1.or.lambda_DM<=0d0) error stop 'Set explicit iterations and positive learning rate'
  if(len_trim(fields_dir)==0.or.len_trim(spin_file)==0.or.len_trim(output_file)==0) error stop 'Explicit fit paths required'
  inquire(file=trim(output_file),exist=output_exists)
  if(output_exists) error stop 'Fit output exists; choose a new output_file'
  allocate(Err(4))
  allocate(r_hat(n_spinor,3))
  allocate(sub_skip(n_spinor))
  allocate(read_angle(3*n_spinor))
  call random_seed()

  pi=3.1415926d0
  n_file=1
  ! n_kind comes from fit_control
  n_iter=diagonal_iters+full_iters
  lambda_H=0.0d0
  lambda_B=100000.0d0
  ! lambda_DM comes from fit_control
  act='fit'
  wt=1.0



   open(unit=101,file=trim(fields_dir)//'/spin1/V_HYC_dp',access='sequential',action='read'&
                &,form='unformatted')
                read(101) nx,ny,nz
                ni=nx*ny*nz
close(101)
!  open(unit=101,file='./sym_BZ',action='read'&
!                &,access='sequential',form='unformatted')
!        read(101) ni
!        print*, ni
        allocate(read_VA(ni,4))
        allocate(VAs(n_file,n_kind,ni,4))
        allocate(pt_TF(ni))
!  close(101)

  pt_TF(:)=1

  allocate(spin(n_file,n_kind,3*n_spinor))

  allocate(r_vec(n_file,n_kind,n_spinor,3))
  allocate(VA_ave(ni,4))

  allocate(V0(ni,3))
  allocate(Vn(ni,8))
  allocate(JH(ni,n_spinor))
  JH=0d0

  allocate(DV(ni,4))
  allocate(cal_V(ni,4))

  allocate(B_general(ni,3,n_spinor,3)) ! nR, i_dir, j_dir, i2_dir
  allocate(blocked(ni,3,n_spinor,3))
  blocked(:,:,:,:)=1.0d0
  call make_block()

  if (act .eq. 'plot' )then
        call do_plot()
        stop
  endif



  Do i_file=1,n_file
        print*,i_file
        call int2char(i_file,file_char)
        path=trim(spin_file)

        call int2char(n_spinor*3,spinor_char)
        formatname='('//trim(spinor_char)//'F12.6)'
        open(unit=102,file=path,access='sequential',action='read')
        do i_kind=1,n_kind
                read(102,*) read_angle(:)
                spin(i_file,i_kind,:)=read_angle(:)

                do i =1,n_spinor
                        r_vec(i_file,i_kind,i,1)=spin(i_file,i_kind,(i-1)*3+1)
                        r_vec(i_file,i_kind,i,2)=spin(i_file,i_kind,(i-1)*3+2)
                        r_vec(i_file,i_kind,i,3)=spin(i_file,i_kind,(i-1)*3+3)
                enddo
                !print*,i_file,r_vec(i_file,i_kind,:,:)
        enddo
        close(102)
        Do i_kind=1,n_kind
                call int2char(i_kind,kind_char)
                path=trim(fields_dir)//'/spin'//trim(kind_char)//'/V_HYC_dp'
                open(unit=101,file=trim(path)&
                        &,access='sequential',action='read'&
                        &,form='unformatted')
                        read(101) i
                        read(101) read_VA
     ! print*,i,ni
                close(101)
                VAs(i_file,i_kind,:,:)=read_VA(:,:)
        enddo

  enddo
!stop

          VA_ave(:,:)=0.0d0

  Do  i_file=1,n_file
  Do  i_kind=1,n_kind
        VAs(i_file,i_kind,:,:)=VAs(i_file,i_kind,:,:)-VA_ave(:,:)
        print*,sum(abs(VAs(i_file,i_kind,:,2))),sum(abs(VAs(i_file,i_kind,:,3))),sum(abs(VAs(i_file,i_kind,:,4)))
  enddo
  enddo
  V0(:,:)=0.0001d0
  Vn(:,:)=0.0001d0
  B_general(:,:,:,:)=0.0000000000d0
  B_general(:,1,:,1)=-0.000010000d0
  B_general(:,2,:,2)=-0.000010000d0
  B_general(:,3,:,3)=-0.000010000d0

  !goto 300! remove if you want to restart from a existing one
  if(len_trim(restart_file)>0) then
  open(unit=102,file=trim(restart_file),action='read',access='sequential',form='unformatted')
        read(102) VA_ave
        read(102) V0
        read(102) Vn
        read(102) B_general
        read(102) JH
  close(102)

  endif
  300 continue




  Do i_iter=1,n_iter
        Err(:)=0.0d0
        i_count=0

        DO i_file=1,n_file
        DO i_kind=1,n_kind
             !   if (i_kind .gt. 7 .and. i_kind .lt. 9 )cycle
                i_count=i_count+1
          !      if (i_file .eq. 2 .or. i_file .eq. 3 .or. i_file .eq. 7 .or. i_file .eq. 8 .or. i_file .eq. 17&
          !              & .or. i_file .eq. 23 ) cycle
               ! call random_number(rand)
               ! if (rand .lt. 0.2) cycle
                r_hat(:,:)=r_vec(i_file,i_kind,:,:)
                call V_func_value(r_hat(:,:),cal_V)
                DV(:,:)=cal_V(:,:)-VAs(i_file,i_kind,:,:)

                Err(1)=Err(1)+sum(abs(DV(:,1)*pt_TF(:))**2)
                Err(2)=Err(2)+sum(abs(DV(:,2)*pt_TF(:))**2)
                Err(3)=Err(3)+sum(abs(DV(:,3)*pt_TF(:))**2)
                Err(4)=Err(4)+sum(abs(DV(:,4)*pt_TF(:))**2)
            !    call random_number(rand)
            !    if (rand .lt. 0.2) cycle

                do i_dir=1,3
                Do j_dir=1,n_spinor
                Do i2_dir=1,3
                        if (.not.update_allowed(i_dir,i2_dir,i_iter,diagonal_iters)) cycle
                        B_general(:,i_dir,j_dir,i2_dir)=&
                                &B_general(:,i_dir,j_dir,i2_dir)-&
                                &lambda_DM*DV(:,i_dir+1)*r_hat(j_dir,i2_dir)
                       ! print*,sum(abs(lambda_DM*DV(:,i_dir+1)*r_hat(j_dir,i2_dir)))
                ENDDO
                ENDDO
                ENDDO

                B_general(:,:,:,:)=B_general(:,:,:,:)*blocked(:,:,:,:)

                700 continue


        enddo
        enddo

        Err(:)=Err(:)/i_count*n_file*n_kind
        print'(I8,3E24.10)',i_iter,Err(2)**0.5/(sum(abs(VAs(:,:,:,2))**2))**0.5,&
                &Err(3)**0.5/(sum(abs(VAs(:,:,:,3))**2))**0.5,&
                &Err(4)**0.5/(sum(abs(VAs(:,:,:,4))**2))**0.5
        if (mod(i_iter,50) .eq. 49.or.i_iter==n_iter) then
                print*,'write'
        open(unit=102,file=trim(output_file),action='write',access='sequential',position='append',status='replace',form='unformatted')
        write(102) VA_ave
        write(102) V0
        write(102) Vn
        write(102) B_general
        write(102) JH
        close(102)
        endif


ENDDO




contains

subroutine int2char(in_int,out_char)
        integer,intent(in)              :: in_int
        character(len=20),intent(out)    :: out_char

IF (in_int .LT. 10 )            WRITE(out_char,'(I1)')in_int
IF ((in_int .LT. 100)  .and. (in_int .GE. 10)   )   WRITE(out_char,'(I2)')in_int
IF ((in_int .LT. 1000) .and. (in_int .GE. 100)  )   WRITE(out_char,'(I3)')in_int
IF ((in_int .LT. 10000) .and. (in_int .GE. 1000) )  WRITE(out_char,'(I4)')in_int
IF ((in_int .LT. 100000) .and. (in_int .GE. 10000) )WRITE(out_char,'(I5)')in_int

end subroutine int2char


subroutine V_func_value(angle_in,V_out)
        implicit none
        real(kind(0d0)), intent(in) :: angle_in(n_spinor,3)
        real(kind(0d0)), intent(out):: V_out(ni,4)
        real(kind(0d0)) ::r(n_spinor,3)
        integer ::i_site,i_dir,j_dir,i2_dir

r(:,:)=angle_in(:,:)



V_out(:,:)=0d0

do i_dir =1,3
do j_dir =1,n_spinor
do i2_dir=1,3
        V_out(:,i_dir+1)=V_out(:,i_dir+1)+B_general(:,i_dir,j_dir,i2_dir)*r(j_dir,i2_dir)
enddo
enddo
enddo


end subroutine V_func_value



subroutine do_plot()
        implicit none
        real(kind(0d0)):: V_plot(ni,4)
        real(kind(0d0))::angle_in(n_spinor,3)

        open(unit=102,file='V_func_tmp',action='read',access='sequential',form='unformatted')
        read(102) VA_ave
        read(102) V0
        read(102) Vn
        read(102) B_general
        close(102)




        call int2char(n_spinor*3,spinor_char)
        formatname='('//trim(spinor_char)//'F12.6)'
        open(unit=102,file='./angle_1/spin_info',access='sequential',action='read')
        do i_kind=1,1
                read(102,*) read_angle(:)
        enddo
        close(102)



        angle_in(1,:)=read_angle(1:3)
      !  angle_in(2,:)=read_angle(4:6)
     !   angle_in(3,:)=read_angle(7:9)
     !   angle_in(4,:)=read_angle(10:12)


        call  V_func_value(angle_in,V_plot)
        V_plot(:,:)= V_plot(:,:)+VA_ave(:,:)
        open(unit=103,file='V_fit',action='write',access='sequential',position='append'&
                &,form='unformatted',status='replace')
        write(103) ni
        write(103) V_plot
        close(103)





end subroutine do_plot

subroutine make_block()



        do ix=0,nx-1
        do iy=0,ny-1
        do iz=0,nz/2-10
                i=1+ix+iy*nx+ny*nx*iz
 !               blocked(i,:,:,:)=0.0d0
             !   blocked(i,:,11:12,:)=0.0d0
        enddo
        enddo
        enddo

        do ix=0,nx-1
        do iy=0,ny-1
        do iz=nz/2+10,nz-1
                i=1+ix+iy*nx+ny*nx*iz
  !              blocked(i,:,:,:)=0.0d0
             !   blocked(i,:,11:12,:)=0.0d0
        enddo
        enddo
        enddo




end subroutine make_block


end program


