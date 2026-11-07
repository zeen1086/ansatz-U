program cal_RVR

        implicit none
        character(len=1024) :: the_input

        integer         ::iwan,jwan,iR,jR,nwan,nwan_sp,nR
        integer         :: nx,ny,nz,ix,iy,iz,j,i1,i2
        integer         :: z_cut,ni
        complex(kind(0d0)),allocatable  :: read_wan(:,:,:)
        complex(kind(0d0)),allocatable  :: ini_wan(:,:,:),fin_wan(:,:,:),tmp_wan(:,:)
        integer,allocatable     :: ini_TF(:),fin_TF(:)
        integer :: convert_size_x,convert_size_y,convert_size_z,d_ext,i_dir,j_dir,k_dir,i,i2_dir
        integer,allocatable ::dir_vec(:,:),dir_inv(:)
        integer,allocatable :: hex_vec(:,:)
        character(len=200) :: wan_dir,func_dir
        character(len=20)  :: iwan_char,jwan_char,ndir_char,file_char
        complex(kind(0d0)) :: cI
        !!
        integer         :: n_func,i_func
        complex(kind(0d0)),allocatable  :: RVR(:,:,:),RVR_local(:,:,:)
        complex(kind(0d0)),allocatable  ::RVR_tmp(:,:),RVR_ave(:,:),RVR_ave_local(:,:)
        complex(kind(0d0)),allocatable  :: RVR_check(:,:,:)
        !! func !!
        real(kind(0d0)),allocatable :: read_V(:,:)
        real(kind(0d0)),allocatable  :: B_general(:,:,:,:)
        real(kind(0d0)),allocatable :: ave_read(:,:),V_back(:,:)
        integer :: nspin,wan_min,wan_max
        !! func !!

         !!!mpi!!!
        real(kind(0d0)) :: start, finish,start2, finish2
        complex(kind(0d0)),allocatable ::wan_fun_local(:,:,:)
        integer,allocatable     :: pool(:,:)
        integer                 :: muti_pool,imuti_pool,ierr,my_rank,size
        integer                 :: npool ,ipool

        !!!mpi!!! 
        character(len=200)      :: useless_string,nscf_dir,THE_dir




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
read(101,'(A3,I3)') useless_string, nx
read(101,'(A3,I3)') useless_string, ny
read(101,'(A3,I3)') useless_string, nz
read(101,*) ! A6 axial
read(101,*) ! A9 period_x
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

wan_dir=trim(nscf_dir)//'./make_wan/'
!convert_size=1
d_ext=(2*convert_size_x+1)*(2*convert_size_y+1)*(2*convert_size_z+1)



allocate(dir_vec(d_ext,3))
i=1
do i_dir=-convert_size_x,convert_size_x
do j_dir=-convert_size_y,convert_size_y
do k_dir=-convert_size_z,convert_size_z
        dir_vec(i,:)=(/i_dir,j_dir,k_dir/)
        print*,i,dir_vec(i,:)
        i=i+1
enddo
enddo
enddo


if  (nspin .eq. 1) then
        nwan_sp=2*(wan_max-wan_min+1)
        i1=wan_min*2-1
        i2=wan_max*2
else
        nwan_sp=(wan_max-wan_min+1)
        i1=wan_min
        i2=wan_max
endif







allocate(RVR(n_func,nwan_sp*d_ext,nwan_sp*d_ext))
allocate(RVR_local(n_func,d_ext,nwan_sp*d_ext))
!allocate(RVR_ave(nwan_sp*d_ext,nwan_sp*d_ext))
!allocate(RVR_ave_local(nwan_sp*d_ext,nwan_sp*d_ext))




RVR(:,:,:)=dcmplx(0.0d0,0.0d0)
!RVR_ave(:,:)=dcmplx(0.0d0,0.0d0)

do i =i1,i2
        print*,i
        call int2char(i,file_char)
        open(unit=103,file='./tmp/B_general_'//trim(file_char),action='read'&
                                &,access='sequential',status='old',form='unformatted')
                read(103) RVR_local(:,:,:)
        close(103)

        if (nspin .eq.1) then
                RVR(:,(i-(2*wan_min-1))*d_ext+1:(i-(2*wan_min-1)+1)*d_ext,:)&
                                                &=RVR_local(:,:,:)
        else
                RVR(:,(i-wan_min)*d_ext+1:(i-wan_min+1)*d_ext,:)&
                                                &=RVR_local(:,:,:)
        endif
enddo

        open(unit=103,file='./B_general',action='write',access='sequential',position='append',status='replace',form='unformatted')
                write(103) RVR(:,:,:)
        close(103)

print*, sum(abs(RVR(:,:,:))**2)

        !open(unit=103,file='./B_ave',action='write',access='sequential',position='append',status='replace',form='unformatted')
        !        write(103) RVR_ave(:,:)
        !close(103)


        !open(unit=103,file='check_cmap',action='write',access='sequential',position='append',status='replace')
        !do i=1,20
        !do j=1,20
        !        write(103,*)i,j,log(abs(RVR(6,(i-1)*d_ext+14,(j-1)*d_ext+14)))
        !enddo
        !enddo
        !close(103)
 !               print*,RVR(3,0*d_ext+14,(i-1)*d_ext+14)*7.0,&
 !                       &RVR(6,0*d_ext+14,(i-1)*d_ext+14)*7.0,&
 !                       &RVR(9,0*d_ext+14,(i-1)*d_ext+14)*7.0,&
 !                        &RVR(12,0*d_ext+14,(i-1)*d_ext+14)*7.0
 !       enddo
 !print*, RVR(9,(1-1)*d_ext+14,(21-1)*d_ext+14)
 !print*, RVR(9,(11-1)*d_ext+14,(21-1)*d_ext+14)
 !print*, RVR(9,(31-1)*d_ext+14,(21-1)*d_ext+14)
 !print*, RVR(9,(11-1)*d_ext+2,(21-1)*d_ext+14)
 !print*, RVR(9,(31-1)*d_ext+5,(21-1)*d_ext+14)
 !print*, RVR(9,(1-1)*d_ext+11,(21-1)*d_ext+14)

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
        complex(kind(0d0)) :: pot(ni,3)
        integer ::iv,jv,i_dir_cal,j_dir_cal
        complex(kind(0d0)) ::V_psi(ni,2)


        RVR_out=dcmplx(0.0d0,0.0d0)

        pot(:,:)=dcmplx(0.0d0,0.0d0)
        !pot(:,:)=B_general(:,:,j_in,i_in)
        pot(:,1)=B_general(:,1,j_in,i_in)-cI*B_general(:,2,j_in,i_in)
        pot(:,2)=B_general(:,1,j_in,i_in)+cI*B_general(:,2,j_in,i_in)
        pot(:,3)=B_general(:,3,j_in,i_in)




        do j_dir_cal=1,d_ext
                jv=dir_inv(j_dir_cal)        
                if ( .not.  ini_TF(jv) .eq. 1) cycle
                V_psi(:,1)=pot(:,1)*ini_wan(jv,:,2)+pot(:,3)*ini_wan(jv,:,1)
                V_psi(:,2)=pot(:,2)*ini_wan(jv,:,1)-pot(:,3)*ini_wan(jv,:,2)
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
        complex(kind(0d0)) :: pot(ni,3)
        complex(kind(0d0)) ::V_psi(ni,4)
        integer ::iv,jv,i_dir_cal,j_dir_cal

        RVR_out=dcmplx(0.0d0,0.0d0)
        


        !pot(:,:)=B_general(:,:,j_in,i_in)
        pot(:,1)=ave_read(:,2)-cI*ave_read(:,3)
        pot(:,2)=ave_read(:,2)+cI*ave_read(:,3)
        pot(:,3)=ave_read(:,4)




        do i_dir_cal=1,d_ext
        do j_dir_cal=1,d_ext
                iv=dir_inv(i_dir_cal)
                jv=dir_inv(j_dir_cal)

                if ( .not. (fin_TF(iv) .eq. 1 .and. ini_TF(jv) .eq. 1)) cycle

                V_psi(:,1)=pot(:,1)*ini_wan(jv,:,2)
                V_psi(:,2)=pot(:,2)*ini_wan(jv,:,1)
                V_psi(:,3)=pot(:,3)*ini_wan(jv,:,1)
                V_psi(:,4)=pot(:,3)*ini_wan(jv,:,2)


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
