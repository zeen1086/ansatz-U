program cal_RVR

        implicit none
        character(len=1024) :: the_input

        integer         ::iwan,jwan,iR,jR,nwan,nwan_sp,nR
        integer         :: nx,ny,nz,ix,iy,iz,j
        integer         :: z_cut,ni
        complex(kind(0d0)),allocatable  :: read_wan(:,:,:)
        complex(kind(0d0)),allocatable  :: ini_wan(:,:,:),fin_wan(:,:,:),tmp_wan(:,:)
        integer,allocatable     :: ini_TF(:),fin_TF(:)
        integer :: convert_size_x,convert_size_y,convert_size_z,d_ext,i_dir,j_dir,k_dir,i,i2_dir
        integer,allocatable ::dir_vec(:,:),dir_inv(:)
        integer,allocatable :: hex_vec(:,:)
        character(len=200) :: wan_dir,func_dir
        character(len=20)  :: iwan_char,jwan_char,ndir_char,column_char
        integer :: arg_status
        complex(kind(0d0)) :: cI
        !!
        integer         :: n_func,i_func,wan_min,wan_max
        complex(kind(0d0)),allocatable  :: RVR(:,:,:),RVR_local(:,:,:)
        complex(kind(0d0)),allocatable  ::RVR_tmp(:,:),RVR_ave(:,:),RVR_ave_local(:,:)
        complex(kind(0d0)),allocatable  :: RVR_check(:,:,:)
        !! func !!
        real(kind(0d0)),allocatable :: read_V(:,:)
        real(kind(0d0)),allocatable  :: B_general(:,:,:,:)
        complex(kind(0d0)),allocatable  :: B_sig(:,:,:,:)
        real(kind(0d0)),allocatable :: ave_read(:,:),V_back(:,:)
        integer :: nspin
        !! func !!

         !!!mpi!!!
        real(kind(0d0)) :: start, finish,start2, finish2
        complex(kind(0d0)),allocatable ::wan_fun_local(:,:,:)
        integer,allocatable     :: pool(:,:)
        integer                 :: muti_pool,imuti_pool,ierr,my_rank,size
        integer                 :: npool ,ipool
        integer         :: n_dep, s_comp, B_comp
        !!!mpi!!! 
        character(len=200)      :: useless_string,nscf_dir,THE_dir


include 'mpif.h'


cI=dcmplx(0.0d0,1.0d0)
call get_command_argument(1,the_input)
if (len_trim(the_input)==0) error stop 'Supply explicit THE.in path as argument 1'
open(unit=101,file=trim(the_input),action='read',access='sequential',status='old')
read(101,'(A9,A200)')useless_string,func_dir
read(101,'(A9,A200)')useless_string,nscf_dir
read(101,'(A8,A200)')useless_string,THE_dir
read(101,'(A7,I3)')  useless_string,n_func
read(101,'(A15,I2)') useless_string,convert_size_x
read(101,'(A15,I2)') useless_string,convert_size_y
read(101,'(A15,I2)') useless_string,convert_size_z
read(101,'(A5,I3)') useless_string, nwan
read(101,'(A6,I2)') useless_string, nspin
read(101,*) !'(A3,I3)') useless_string, nx
read(101,*) !'(A3,I3)') useless_string, ny
read(101,*) !'(A3,I3)') useless_string, nz
read(101,*) ! A6 axial
read(101,*) ! A9 period_xf
read(101,*) ! A9 period_y
read(101,*) ! A9 period_z
read(101,*)
read(101,*)
read(101,*)
read(101,*)
read(101,*)
read(101,'(A8,I3)') useless_string, wan_min
read(101,'(A8,I3)') useless_string, wan_max
close(101)

!nwan=1
!nspin=1
!nwan_sp=nwan*2/nspin
!n_func=24 !B_general j_dir*i2_dir

call get_command_argument(2,column_char)
read(column_char,*,iostat=arg_status) iwan
if(arg_status/=0) error stop 'Argument 2 must be a spinor-Wannier column index'
call int2char(iwan,column_char)


n_dep=1
s_comp=3
B_comp=3

if(n_func/=n_dep*s_comp) error stop 'Single-Ni RVR requires n_func=3; do not truncate extra channels'
if(nspin/=1.and.nspin/=2) error stop 'nspin must be 1 or 2'


if (nspin .eq. 1) then
        if (iwan .lt. wan_min*2-1 .or. iwan .gt. wan_max*2) then
                print*,'iwan out of min,max'
                error stop 'Column outside requested Wannier range'
        endif
else
        if (iwan .lt. wan_min .or. iwan .gt. wan_max) then
                print*,'iwan out of min,max'
                error stop 'Column outside requested Wannier range'
        endif
endif

wan_dir=trim(nscf_dir)//'./make_wan/'
d_ext=(2*convert_size_x+1)*(2*convert_size_y+1)*(2*convert_size_z+1)
call  int2char(d_ext,ndir_char)



allocate(dir_vec(d_ext,3))
allocate(dir_inv(d_ext))
i=1
do i_dir=-convert_size_x,convert_size_x
do j_dir=-convert_size_y,convert_size_y
do k_dir=-convert_size_z,convert_size_z
        dir_vec(i,:)=(/i_dir,j_dir,k_dir/)
        i=i+1
enddo
enddo
enddo




do i_dir=1,d_ext
do j_dir=1,d_ext
        if (all(dir_vec(i_dir,:) .eq. -dir_vec(j_dir,:))) then
                dir_inv(i_dir)=j_dir
                exit
        endif
enddo
enddo



open(unit=101,file=trim(func_dir)//'/reference_V_HYC_dp',access='sequential',form='unformatted',action='read')
read(101) nx,ny,nz
print*,nx,ny,nz
close(101)

allocate(read_V(nx*ny*nz,4))


open(unit=101,file=trim(func_dir)//'/reference_V_HYC_dp',access='sequential',form='unformatted',action='read')
read(101) nx,ny,nz
if (nspin .eq. 2) then
        read(101) read_V
else
        read_V(:,:)=0.0d0
endif
close(101)



ni=nx*ny*nz
allocate(V_back(ni,4))
allocate(B_general(ni,B_comp,n_dep,S_comp))
allocate(B_sig(ni,B_comp,n_dep,S_comp))
allocate(fin_wan(d_ext,ni,2))
if (nspin .eq. 1) allocate(tmp_wan(d_ext,ni))
allocate(hex_vec(ni,3))


do ix=0 ,nx-1
do iy=0 ,ny-1
do iz=0 ,nz-1
        i=1+(ix)+(iy)*nx+(iz)*nx*ny
        hex_vec(i,:)=(/ix,iy,iz/)
enddo
enddo
enddo




!open(unit=101,file=trim(func_dir)//'./R_vec',action='read'&
!                &,access='sequential',form='unformatted')
!read(101) ni
!allocate(V_back(ni,4))
!allocate(B_general(ni,B_comp,n_dep,S_comp))
!allocate(B_sig(ni,B_comp,n_dep,S_comp))
!allocate(fin_wan(d_ext,ni,2))
!if (nspin .eq. 1) allocate(tmp_wan(d_ext,ni))
!allocate(hex_vec(ni,3))
!read(101) hex_vec
!close(101)



!open(unit=101,file=trim(func_dir)//'./sym_BZ',action='read'&
!                &,access='sequential',form='unformatted')
!read(101) ni
!allocate(ave_read(ni,4))
!allocate(V_back(ni,4))
!allocate(B_general(ni,B_comp,n_dep,S_comp))
!allocate(B_sig(ni,B_comp,n_dep,S_comp))
!allocate(fin_wan(d_ext,ni,2))
!if (nspin .eq. 1) allocate(tmp_wan(d_ext,ni))
!allocate(hex_vec(ni,3))
!close(101)
!
!open(unit=101,file=trim(func_dir)//'./R_vec',action='read'&
!                &,access='sequential',form='unformatted')
!read(101) hex_vec
!close(101)



do i=1,ni
        ix=hex_vec(i,1)
        iy=hex_vec(i,2)
        iz=hex_vec(i,3)
        j=1+ix+nx*iy+nx*ny*iz
        V_back(i,:)=read_V(j,:)
enddo





open(unit=101,file=trim(func_dir)//'./V_func_tmp',&
                &access='sequential',action='read'&
                &,form='unformatted')

        read(101) !ave_read
        read(101) 
        read(101) 
        read(101) B_general

close(101)



if (B_comp .eq. 3) then
        B_sig(:,1,:,:)=B_general(:,1,:,:)-cI*B_general(:,2,:,:)
        B_sig(:,2,:,:)=B_general(:,1,:,:)+cI*B_general(:,2,:,:)
        B_sig(:,3,:,:)=B_general(:,3,:,:)
elseif(B_comp .eq. 2 ) then
        B_sig(:,1,:,:)=B_general(:,1,:,:)-cI*B_general(:,2,:,:)
        B_sig(:,2,:,:)=B_general(:,1,:,:)+cI*B_general(:,2,:,:)
       ! B_sig(:,3,:,:)=B_general(:,3,:,:)
endif

deallocate(B_general)

!ave_read(:,:)=ave_read(:,:)-V_back(:,:)

!ave_read(:,:)=dcmplx(0.0d0,0.0d0)
!B_general(:,:,:,:)=dcmplx(0.0001d0,0.0d0)


!allocate(ini_wan(d_ext,nx*ny*nz,2))
!allocate(fin_wan(d_ext,nx*ny*nz,2))
allocate(ini_TF(d_ext))
allocate(fin_TF(d_ext))



if  (nspin .eq. 1) then

        nwan_sp=2*(wan_max-wan_min+1)
        npool=nwan_sp
        allocate(pool(npool,1))
        ipool=1
        !do iwan=2*wan_min-1,2*wan_max
        do jwan=2*wan_min-1,2*wan_max
                pool(ipool,1)=jwan
                ipool=ipool+1
        enddo
        !enddo
else
        nwan_sp=(wan_max-wan_min+1)
        npool=nwan_sp
        allocate(pool(npool,1))
        ipool=1
        !do iwan=wan_min,wan_max
        do jwan=wan_min,wan_max
                pool(ipool,:)=jwan
                ipool=ipool+1
        enddo
        !enddo

endif





allocate(RVR(n_func,d_ext,nwan_sp*d_ext))
allocate(RVR_local(n_func,d_ext,nwan_sp*d_ext))
allocate(RVR_tmp(d_ext,d_ext))
!allocate(RVR_ave(nwan_sp*d_ext,nwan_sp*d_ext))
!allocate(RVR_ave_local(nwan_sp*d_ext,nwan_sp*d_ext))


call mpi_init(ierr)
call mpi_comm_rank(MPI_COMM_WORLD,my_rank,ierr)
call mpi_comm_size(MPI_COMM_WORLD,size,ierr)
if(size<2) then
 call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
 error stop 'RVR requires at least two MPI ranks (one dispatcher)'
endif


if ( my_rank.eq.1)  then
open(unit=103,file='./report_'//trim(column_char),action='write',access='sequential',position='append',status='replace')
write(103,*) 'start calculation for explicit spinor-Wannier column'
write(103,*) 'nx,ny,nz =', nx,ny,nz
close(103)
endif


muti_pool=0
DO
        muti_pool=muti_pool+1
        IF(muti_pool*(size-1) .gt. npool) EXIT
ENDDO
muti_pool=muti_pool-1


if (.not. my_rank.eq.0)  then
                allocate(read_wan(d_ext,ni,nspin))
                allocate(ini_wan(d_ext,ni,2))
                if (nspin .eq. 1) then
                        if (mod(iwan,2) .eq. 0) then
                                call int2char(iwan/2,iwan_char)
                        else
                                call int2char((iwan+1)/2,iwan_char)
                        endif
                else
                        call int2char(iwan,iwan_char)
                endif

  !              call sleep(my_rank)
                open(unit=102+my_rank,file=trim(wan_dir)//&
                        &'./wan_data/wan_'//trim(iwan_char),access='stream',form='unformatted',action='read')
                read(102+my_rank) read_wan
                close(102+my_rank)


                fin_wan(:,:,:)=dcmplx(0.0d0,0.0d0)


                if (nspin .eq. 1) then
                        if (mod(iwan,2) .eq. 0) then
                                fin_wan(:,:,2)=read_wan(:,:,1)
                        else
                                fin_wan(:,:,1)=read_wan(:,:,1)
                        endif
                else
                        fin_wan(:,:,:)=read_wan(:,:,:)
                endif
                if ( my_rank.eq.1)  then
                open(unit=103,file='./report_'//trim(column_char),action='write',access='sequential',position='append')
                write(103,*) 'read fin wan done'
                close(103)
                endif



        RVR_local(:,:,:)=dcmplx(0.0d0,0.0d0)
        !RVR_ave_local(:,:)=dcmplx(0.0d0,0.0d0)
        DO imuti_pool=0,muti_pool
                ipool=imuti_pool*(size-1)+my_rank
                if (ipool .gt. npool) cycle
                jwan=pool(ipool,1)


                if ( my_rank.eq.1)  then
                open(unit=103,file='./report_'//trim(column_char),action='write',access='sequential',position='append')
                write(103,*) 'imuti_pool / muti_pool' , imuti_pool,' /',muti_pool
                write(103,*) 'current iwan, jwan=' , iwan,jwan
                close(103)
                endif

 !               if ( jwan .gt. 40 ) cycle
 !               if (.not. (jwan-1)/2 .eq. (iwan-1)/2) cycle





                if (nspin .eq. 1) then
                        if (mod(jwan,2) .eq. 0) then
                                call int2char(jwan/2,jwan_char)
                        else
                                call int2char((jwan+1)/2,jwan_char)
                        endif
                else
                        call int2char(jwan,jwan_char)
                endif


     !           if (.not. jwan .eq. iwan) cycle


 !               call sleep(my_rank)
                open(unit=104+my_rank,file=trim(wan_dir)//&
                        &'./wan_data/wan_'//trim(jwan_char),access='stream',form='unformatted',action='read')
                read(104+my_rank) read_wan
                close(104+my_rank)

        
                ini_wan(:,:,:)=dcmplx(0.0d0,0.0d0)

                 if (nspin .eq. 1) then
                        if (mod(jwan,2) .eq. 0) then
                                ini_wan(:,:,2)=read_wan(:,:,1)
                        else
                                ini_wan(:,:,1)=read_wan(:,:,1)
                        endif
                else
                        ini_wan(:,:,:)=read_wan(:,:,:)
                endif



                if ( my_rank.eq.1)  then
                open(unit=103,file='./report_'//trim(column_char),action='write',access='sequential',position='append')
                write(103,*) 'read ini wan done'
                close(103)
                endif


                open(unit=102,file=trim(wan_dir)//&
                        &'./cal_TF_'//trim(iwan_char),access='sequential',action='read')
                read(102,'('//trim(ndir_char)//'I8)') fin_TF
                close(102)

                open(unit=102,file=trim(wan_dir)//&
                        &'./cal_TF_'//trim(jwan_char),access='sequential',action='read')
                read(102,'('//trim(ndir_char)//'I8)') ini_TF
                close(102)

                if ( my_rank.eq.1)  then
                open(unit=103,file='./report_'//trim(column_char),action='write',access='sequential',position='append')
                write(103,*) 'read done'
                close(103)
                endif



                do j_dir=1,n_dep
                do i2_dir=1,S_comp
                        i_func=1+(i2_dir-1)+S_comp*(j_dir-1)
                        !print*,my_rank,iwan,jwan,i_func
                        call cpu_time(start)
                        call cal_each_func(j_dir,i2_dir,RVR_tmp)
                        if (nspin .eq. 1) then

                                RVR_local(i_func,:,&
                                        &(jwan-(2*wan_min-1))*d_ext+1:(jwan-(2*wan_min-1)+1)*d_ext)&
                                                &=RVR_tmp(:,:)
                        else
                                RVR_local(i_func,:,&
                                        &(jwan-wan_min)*d_ext+1:(jwan-wan_min+1)*d_ext)&
                                                &=RVR_tmp(:,:)
                        endif





                        call cpu_time(finish)
                        if ( my_rank.eq.1)  then
                                open(unit=103,file='./report_'//trim(column_char),action='write',access='sequential',position='append')
                                write(103,*) j_dir,i2_dir,i_func, '-th function took', finish-start, 's'
                                close(103)
                        endif
                enddo
                enddo

              !  call cal_ave_part(RVR_tmp)
                !RVR_ave_local((iwan-1)*d_ext+1:iwan*d_ext,(jwan-1)*d_ext+1:jwan*d_ext)=RVR_tmp(:,:)


        enddo
        call MPI_Barrier(  MPI_COMM_WORLD, ierr)
        call MPI_SEND(RVR_local(:,:,:),n_func*(nwan_sp*d_ext)*d_ext&
                                &,MPI_DOUBLE_COMPLEX,0,my_rank,MPI_COMM_WORLD,ierr)
       ! call MPI_SEND(RVR_ave_local(:,:),(nwan_sp*d_ext)**2&
       !                         &,MPI_DOUBLE_COMPLEX,0,my_rank+2*size,MPI_COMM_WORLD,ierr)

endif !my_rank


if (my_rank.eq.0)  then
        RVR(:,:,:)=dcmplx(0.0d0,0.0d0)
        call MPI_Barrier(  MPI_COMM_WORLD, ierr)
       ! RVR_ave(:,:)=dcmplx(0.0d0,0.0d0)
        do ipool=1,size-1
                call MPI_RECV(RVR_local(:,:,:),n_func*(nwan_sp*d_ext)*d_ext,MPI_DOUBLE_COMPLEX,&
                    &MPI_ANY_SOURCE,ipool,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
                RVR(:,:,:)=RVR(:,:,:)+RVR_local(:,:,:)
       !         call MPI_RECV(RVR_ave_local(:,:),(nwan_sp*d_ext)**2,MPI_DOUBLE_COMPLEX,&
       !             &MPI_ANY_SOURCE,ipool+2*size,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
 !               RVR_ave(:,:)=RVR_ave(:,:)+RVR_ave_local(:,:)
        enddo

!1~3  center....

        open(unit=103,file='./tmp/B_general_'//trim(column_char),action='write',access='sequential',position='append',status='replace',form='unformatted')
                write(103) RVR(:,:,:)
        close(103)

       ! open(unit=103,file='./tmp/B_ave_HYC_change',action='write',access='sequential',position='append',status='replace',form='unformatted')
       !         write(103) RVR_ave(:,:)
       ! close(103)



       ! do iwan=1,nwan
       ! do jwan=1,nwan
       ! do i_dir=1,d_ext
       ! do j_dir=1,d_ext
       !         i=1+(i_dir-1)+d_ext*(iwan-1)
       !         j=1+(j_dir-1)+d_ext*(jwan-1)
       !         if (sum(abs(RVR(:,i,j))**2) .lt. 10.0**(-3)) cycle
       !         print*, iwan,jwan,i_dir,j_dir,sum(abs(RVR(:,i,j))**2)

       ! enddo
       ! enddo
       ! enddo
       ! enddo


!allocate(RVR_check(n_func,2*d_ext,2*d_ext))
!do iwan=1,2
!do jwan=1,2
!do i_dir=1,d_ext
!do j_dir=1,d_ext
!        i=1+(i_dir-1)+d_ext*(iwan-1)
!        j=1+(j_dir-1)+d_ext*(jwan-1)
!        RVR_check(:,i,j)=RVR(:,i,j)

!enddo
!enddo
!enddo
!enddo


!open(unit=103,file='./cmap',action='write',access='sequential',position='append',status='replace')
!        do iwan=1,1
!        do jwan=2,2        
!        do i_dir=1,d_ext
!        do j_dir=1,d_ext        
!                i=1+(i_dir-1)+d_ext*(iwan-1)
!                j=1+(j_dir-1)+d_ext*(jwan-1)
!                if (.not. abs(RVR_check(1,i,j)) .gt. 10.0**(-10)) cycle
!                write(103,*) i,j, abs(RVR_check(1:5,i,j))
!
!        enddo
!        enddo
!        enddo
!        enddo
!close(103)





endif


call mpi_finalize(ierr)


contains

subroutine int2char(in_int,out_char)
        integer,intent(in)              :: in_int
        character(len=20),intent(out)    :: out_char

IF (in_int .LT. 10 )   WRITE(out_char,'(I1)')in_int
IF ((in_int .LT. 100)  .and. (in_int .GE. 10)   )       WRITE(out_char,'(I2)')in_int
IF ((in_int .LT. 1000) .and. (in_int .GE. 100)  )       WRITE(out_char,'(I3)')in_int
IF ((in_int .LT. 10000) .and. (in_int .GE. 1000) )      WRITE(out_char,'(I4)')in_int
IF ((in_int .LT. 100000) .and. (in_int .GE. 10000) )    WRITE(out_char,'(I5)')in_int

end subroutine int2char


subroutine cal_each_func(j_in,i_in,RVR_out)
        implicit none   
        integer,intent(in) :: j_in,i_in
        complex(kind(0d0)),intent(out) :: RVR_out(d_ext,d_ext)
        complex(kind(0d0)) :: pot(ni,B_comp)
        integer ::iv,jv,i_dir_cal,j_dir_cal
        complex(kind(0d0)) ::V_psi(ni,2)


        RVR_out=dcmplx(0.0d0,0.0d0)

        pot(:,:)=dcmplx(0.0d0,0.0d0)
        !pot(:,1)=B_general(:,1,j_in,i_in)-cI*B_general(:,2,j_in,i_in)
        !pot(:,2)=B_general(:,1,j_in,i_in)+cI*B_general(:,2,j_in,i_in)
        !pot(:,3)=B_general(:,3,j_in,i_in)

        pot(:,:) =B_sig(:,:,j_in,i_in)


        do j_dir_cal=1,d_ext
                jv=dir_inv(j_dir_cal)        
                if ( .not.  ini_TF(jv) .eq. 1) cycle
                if (B_comp .eq. 3 ) then
                        V_psi(:,1)=pot(:,1)*ini_wan(jv,:,2)+pot(:,3)*ini_wan(jv,:,1)
                        V_psi(:,2)=pot(:,2)*ini_wan(jv,:,1)-pot(:,3)*ini_wan(jv,:,2)
                elseif (B_comp .eq. 2) then
                        V_psi(:,1)=pot(:,1)*ini_wan(jv,:,2)
                        V_psi(:,2)=pot(:,2)*ini_wan(jv,:,1)
                endif
               ! V_psi(:,3)=pot(:,3)*ini_wan(jv,:,1)
               ! V_psi(:,4)=pot(:,3)*ini_wan(jv,:,2)


        
        do i_dir_cal=1,d_ext
                iv=dir_inv(i_dir_cal)

                if ( .not. (fin_TF(iv) .eq. 1 .and. ini_TF(jv) .eq. 1)) cycle

                RVR_out(i_dir_cal,j_dir_cal)=&
                        !&    sum(dconjg(fin_wan(iv,:,1))*pot(:,1)*ini_wan(jv,:,2))&
                        !&+   sum(dconjg(fin_wan(iv,:,2))*pot(:,2)*ini_wan(jv,:,1))&
                        !&+   sum(dconjg(fin_wan(iv,:,1))*pot(:,3)*ini_wan(jv,:,1))&
                        !&-   sum(dconjg(fin_wan(iv,:,2))*pot(:,3)*ini_wan(jv,:,2))
                        &     sum(dconjg(fin_wan(iv,:,1))*V_psi(:,1))&
                        &+    sum(dconjg(fin_wan(iv,:,2))*V_psi(:,2))!&
                        !&+    sum(dconjg(fin_wan(iv,:,1))*V_psi(:,3))&
                        !&-    sum(dconjg(fin_wan(iv,:,2))*V_psi(:,4))


        enddo
        enddo



end subroutine cal_each_func



subroutine cal_ave_part(RVR_out)
        implicit none
        complex(kind(0d0)),intent(out) :: RVR_out(d_ext,d_ext)
        complex(kind(0d0)) :: pot(ni,B_comp)
        complex(kind(0d0)) ::V_psi(ni,4)
        integer ::iv,jv,i_dir_cal,j_dir_cal

        RVR_out=dcmplx(0.0d0,0.0d0)
        


        !pot(:,:)=B_general(:,:,j_in,i_in)
        if (B_comp .eq. 3) then
        pot(:,1)=ave_read(:,2)-cI*ave_read(:,3)
        pot(:,2)=ave_read(:,2)+cI*ave_read(:,3)
        pot(:,3)=ave_read(:,4)
        elseif (B_comp .eq. 2) then
        pot(:,1)=ave_read(:,2)-cI*ave_read(:,3)
        pot(:,2)=ave_read(:,2)+cI*ave_read(:,3)
        endif



        do i_dir_cal=1,d_ext
        do j_dir_cal=1,d_ext
                iv=dir_inv(i_dir_cal)
                jv=dir_inv(j_dir_cal)

                if ( .not. (fin_TF(iv) .eq. 1 .and. ini_TF(jv) .eq. 1)) cycle

                if (B_comp .eq. 3) then
                        V_psi(:,1)=pot(:,1)*ini_wan(jv,:,2)
                        V_psi(:,2)=pot(:,2)*ini_wan(jv,:,1)
                        V_psi(:,3)=pot(:,3)*ini_wan(jv,:,1)
                        V_psi(:,4)=pot(:,3)*ini_wan(jv,:,2)
                elseif(B_comp .eq. 2) then
                        V_psi(:,1)=pot(:,1)*ini_wan(jv,:,2)
                        V_psi(:,2)=pot(:,2)*ini_wan(jv,:,1)
                        V_psi(:,3)=0.0d0
                        V_psi(:,4)=0.0d0
                endif



                RVR_out(i_dir_cal,j_dir_cal)=&
                       ! & sum(dconjg(fin_wan(iv,:,1))*pot(:,1)*ini_wan(jv,:,2))&
                       ! &+sum(dconjg(fin_wan(iv,:,2))*pot(:,2)*ini_wan(jv,:,1))&
                       ! &+sum(dconjg(fin_wan(iv,:,1))*pot(:,3)*ini_wan(jv,:,1))&
                       ! &-sum(dconjg(fin_wan(iv,:,2))*pot(:,3)*ini_wan(jv,:,2))
                        
                       ! & sum(dconjg(fin_wan(iv,:,1))*V_psi(:,1))&
                       ! &+sum(dconjg(fin_wan(iv,:,2))*V_psi(:,2))&
                       ! &+sum(dconjg(fin_wan(iv,:,1))*V_psi(:,3))&
                       ! &-sum(dconjg(fin_wan(iv,:,2))*V_psi(:,4))
                        & dot_product(fin_wan(iv,:,1),V_psi(:,1))&
                        &+dot_product(fin_wan(iv,:,2),V_psi(:,2))&
                        &+dot_product(fin_wan(iv,:,1),V_psi(:,3))&
                        &-dot_product(fin_wan(iv,:,2),V_psi(:,4))



        enddo
        enddo




end subroutine cal_ave_part






end program

