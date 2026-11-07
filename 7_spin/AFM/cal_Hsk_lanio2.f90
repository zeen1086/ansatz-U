program spin
        implicit none
        character(len=1024) :: the_input

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
        integer :: int3(3),super_nx,super_ny,super_nz,n_func,convert_size_x,convert_size_y,convert_size_z
        integer :: irept,nrept,super_nrept,i_locat,j_locat
        integer :: max_x,max_y,max_z,min_x,min_y,min_z,super_i,super_j,d_ext
        integer :: nspin, nwan_sp, nwan, wan_max,wan_min,i_dir,j_dir,k_dir
        integer,allocatable     :: dir_vec(:,:),inv_dir(:)
        character(len=200)      :: useless_string,func_dir,nscf_dir,THE_dir
        character(len=200)      :: RVR_dir
        character(len=20)       :: rept_char
        complex(kind(0d0)),allocatable  :: RVR(:,:,:),RVR_POT_ave(:,:),RVR_tmp(:,:)
        complex(kind(0d0)),allocatable  :: H_sk(:,:,:)
        complex(kind(0d0))      :: cI

CALL execute_command_line('mkdir -p data')


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
read(101,*)!'(A6,F8.4)')  useless_string,axial
read(101,'(A9,I3)') useless_string,period_x
read(101,'(A9,I3)') useless_string, period_y
read(101,'(A9,I3)') useless_string, period_z
read(101,*)
read(101,*)
read(101,*)
read(101,*)
read(101,*)
read(101,'(A8,I3)') useless_string, wan_min
read(101,'(A8,I3)') useless_string, wan_max
read(101,'(A3,F16.8)') useless_string, a0
read(101,'(A2,F16.8)') useless_string, celldm3
close(101)

if(n_func/=3) error stop 'This single-Ni tensor assembly requires n_func=3; do not truncate extra channels'
if(nspin/=1.and.nspin/=2) error stop 'nspin must be 1 or 2'
if(wan_min<1.or.wan_max>nwan.or.wan_max<wan_min) error stop 'Invalid Wannier range'
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
RVR_dir=trim(THE_dir)//'./6_RVR/'



call read_RVR()


allocate(coefficient(n_func))

H_sk(:,:,:)=dcmplx(0.0d0,0.0d0)

do ix=0,nx-1
do iy=0,ny-1
do iz=0,nz-1
        i=1+iz+nz*iy+ny*nz*ix
        print*,i,nx*ny*nz
        coefficient(:)=0.0d0
        do i_dir=1,d_ext
!                if (.not. all(dir_vec(i_dir,:) .eq. 0)) cycle
                ini_x=cell(i,1)+dir_vec(i_dir,1)
                ini_y=cell(i,2)+dir_vec(i_dir,2)
                ini_z=cell(i,3)+dir_vec(i_dir,3)
               ! if (.not. (( ini_y .eq. 1  .and.  ini_x .eq. 0) .or.&
               !         &( ini_y .eq. 0  .and.  ini_x .eq. 1))) cycle
               ! if (.not. ini_x .eq. 0) cycle
  !              if (.not. ini_z .eq. 0) cycle
  !              if (.not. (all(cell(i,:) .eq. (/0,0,0/)) .or.all(cell(i,:) .eq.(/0,3,0/)))) cycle
                if (.not. (ini_x .ge.0 .and. ini_x/period_x .eq. 0)) cycle
                if (.not. (ini_y .ge.0 .and. ini_y/period_y .eq. 0)) cycle
                if (.not. (ini_z .ge.0 .and. ini_z/period_z .eq. 0)) cycle
                i_locat=1+ini_z+period_z*ini_y+period_z*period_y*ini_x

                coefficient(1:3)=spin_vec(i,:,1)
                !coefficient(4:6)=spin_vec(i,:,2)
                !coefficient(7:9)=spin_vec(i,:,3)
                !coefficient(10:12)=spin_vec(i,:,4)




!print*,coefficient(1:27)

                
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


                
                !print*,fin_x,fin_y,fin_z,super_x,super_y,super_z,jrept

                jx=fin_x-super_x*period_x
                jy=fin_y-super_y*period_y
                jz=fin_z-super_z*period_z
                j_locat=1+jz+period_z*jy+period_z*period_y*jx

       ! print*,j_locat

        do iwan=wan_min,wan_max
        do jwan=wan_min,wan_max
        do is=1,2/nspin
        do js=1,2/nspin
                iRVR=1+(i_dir-1)+d_ext*((iwan-wan_min)*2/nspin+is-1)
                jRVR=1+(j_dir-1)+d_ext*((jwan-wan_min)*2/nspin+js-1)
                i_dim=(is+2/nspin*(iwan-1))+nwan*2/nspin*(i_locat-1)
                j_dim=(js+2/nspin*(jwan-1))+nwan*2/nspin*(j_locat-1)
               ! print*,RVR_tmp(iRVR,jRVR)
                
                ! Retain all tensor RVR entries (complex128).
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





do irept=1,super_nrept
        call int2char(irept,rept_char)
        open(unit=102,file='./data/Hsk_'//trim(rept_char),action='write',access='sequential',position='append'&
                &,status='replace',form='unformatted')
        write(102)      H_sk(:,:,irept)
        close(102)
        !if (.not. mod(irept,10) .eq. 8) cycle
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
 character(len=1024):: texture
 integer::u,stat,site
 call get_command_argument(2,texture)
 if(len_trim(texture)==0) error stop 'Supply texture path as argument 2: direct coefficients, no implicit amplitude'
max_x=-99
max_y=-99
max_z=-99
min_x=99
min_y=99
min_z=99

open(unit=102,file=trim(THE_dir)//'6_RVR/super_H/ir_vec_HYC',action='read',access='sequential',status='old')
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

print*,max_x,max_y,max_z
print*,min_x,min_y,min_z



nx=(max_x-min_x+1)*period_x
ny=(max_y-min_y+1)*period_y
nz=(max_z-min_z+1)*period_z


print*, nx,ny,nz

 allocate(cell(nx*ny*nz,3),spin_vec(nx*ny*nz,3,1))
 open(newunit=u,file=trim(texture),status='old',action='read',iostat=stat)
 if(stat/=0) error stop 'Cannot open explicit legacy tensor texture'
 do ix=0,nx-1
 do iy=0,ny-1
 do iz=0,nz-1
 site=1+iz+nz*iy+nz*ny*ix
 cell(site,:)=[ix+min_x*period_x,iy+min_y*period_y,iz+min_z*period_z]
 read(u,*,iostat=stat) spin_vec(site,:,1)
 if(stat/=0) error stop 'Texture requires nx*ny*nz rows in x,y,z order (z fastest)'
 enddo
 enddo
 enddo
 close(u)
 cI=dcmplx(0d0,1d0)
end subroutine




subroutine load_super_size()
        irept=0
        open(unit=102,file=trim(THE_dir)//'6_RVR/super_H/ir_vec_HYC',action='read',access='sequential',status='old')
        do
        read(102,'(3I8)',iostat=io) int3(:)
        if (io/=0) exit
        irept=irept+1
        enddo
        close(102)
        super_nrept=irept
        allocate(super_ir_vec(super_nrept,3))

        irept=1
        open(unit=102,file=trim(THE_dir)//'6_RVR/super_H/ir_vec_HYC',action='read',access='sequential',status='old')
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

        !open(unit=103,file=trim(RVR_dir)//'./B_ave',action='read',access='sequential',form='unformatted')
        !        read(103) RVR_POT_ave
        !close(103)
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




        !do i_dir=1,d_ext
        !do j_dir=1,d_ext
        !if (.not. all(dir_vec(i_dir,:) .eq. 0)) cycle
        !if (.not. (dir_vec(j_dir,3) .eq. 0)) cycle
        
        !do iwan=6,6!wan_min,wan_max
        !do jwan=6,6!wan_min,wan_max
        !do is=2,2!2/nspin
        !do js=2,2!2/nspin
        !        iRVR=1+(i_dir-1)+d_ext*((iwan-wan_min)*2+is-1)
        !        jRVR=1+(j_dir-1)+d_ext*((jwan-wan_min)*2+js-1)
        !print*, abs(RVR(3,jRVR,iRVR)),dir_vec(j_dir,:)
        !enddo
        !enddo
        !enddo
        !enddo
        !enddo
        !enddo

end subroutine read_RVR



end program

