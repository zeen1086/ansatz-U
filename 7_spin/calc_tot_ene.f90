program wan_interpolat
        use mpi

        implicit none
        character(len=1024) :: the_input

        complex(kind(0d0)),allocatable  :: H_wan(:,:,:),H_k(:,:),H_buff(:,:),eig_vec(:,:,:)
        complex(kind(0d0)),allocatable  :: eig_vec_local(:,:,:),AAR(:,:,:,:)
        complex(kind(0d0)),allocatable  ::unmk(:,:),mat_tmp(:,:),unmk_Hermit(:,:)
        complex(kind(0d0)),allocatable  :: delHH(:,:,:),AAq_tmp(:,:,:),AA(:,:,:),D_h(:,:,:)
        integer :: ik,jk,nk,nkx,nky,nkz,irept,nrept,io,read3_int(3)
        integer :: i,j,wan_x,wan_y,wan_z,nstep,ix,iy,iz,iwan,jwan
        integer :: jx,jy,jz,ikxy,H_dim,ibnd,jbnd,i_tmp,j_tmp
        integer :: nspin, period_x, period_y, period_z,super_size
        integer :: c_ix,c_iy,c_iz,c_jx,c_jy,c_jz,n1,n2,i_dim,j_dim,is,js
        integer :: super_i,super_j,c_i,c_j,npt
        integer :: nE,iE,i1,i2,i3,i4,f1,f2,f3,f4,i_dir
        real(kind(0d0)) ,allocatable     :: E_list(:),eig_list(:,:),spin_list(:,:),spin_tmp(:)
        complex(kind(0d0)),allocatable  :: sig_H(:,:,:),sig_AH(:,:,:)
        complex(kind(0d0)),allocatable  :: sig_H_local(:,:,:),sig_AH_local(:,:,:)
        real(kind(0d0)) ,allocatable     ::k_list(:,:),eigen(:,:),eigen_local(:,:)
        real(kind(0d0)) ,allocatable     :: cr_vec(:,:),eig_tmp(:),Sig_kxy(:),eig_tmp_tot(:,:)
        real(kind(0d0)) ,allocatable     :: nb(:,:),R_list(:,:),pt_k(:)
        real(kind(0d0)) ,allocatable     :: fold_tmp(:,:),fold(:,:,:)
        real(kind(0d0)) :: a0,a(3,3),Ef,Sig_k,dE,fac,pt_ini_k(3),pt_fin_k(3)
        real(kind(0d0)) :: length,vec_tmp(3),b(3,3)
        integer :: n_cell,nbnd,nwan,i_cell
        integer,allocatable     :: ir_vec(:,:)
        character(len=20)       :: irept_char
        complex(kind(0d0))      :: cI,sig_xy_local,phase,phase2,ctmp1,ctmp2,ctmp3
        complex(kind(0d0))      :: I2pi
        complex(kind(0d0)),allocatable:: phase0(:,:),phase1(:,:),phase_k(:)
        complex(kind(0d0)),allocatable:: H_wan_tmp(:,:),RH_tmp(:,:,:),AA_tmp(:,:,:)
        real(kind(0d0)) :: pi,E_max,E_min       
        logical, allocatable    :: rept_TF(:)
        integer,allocatable     :: deg(:),inv_k(:)
        character(len=200)      :: useless_string,nscf_dir,THE_dir,path,title
        character(len=20)       :: kz_char,kxy_char,i_char,j_char
        character(len=20)       :: super_size_char
        integer,allocatable     :: pool(:,:)
        integer                 :: muti_pool,imuti_pool,ierr,my_rank,size
        integer                 :: npool ,ipool,v_max,nline,iline,size_it
        real(kind(0d0)) :: t1,t2,t3,t4,eta,w,celldm3
        complex(kind(0d0))      ::A2_tmp(2,2)
        real(kind(0d0)) :: min_E
        real(kind(0d0)) :: weight, nk_total
        double precision :: filling, occ_per_band, total_E
        integer :: i_occ, n_electron, count_below_fermi
        real(kind(0d0)) :: Ef_low, Ef_high, Ef_mid
        real(kind(0d0)) :: filling_low, filling_high, filling_mid
        real(kind(0d0)) :: tol
        integer :: max_iter, iter
        character(len=20) :: theta_char
        ! 角度ファイル引数 → angle_id を抽出して出力名に使用
        character(len=256) :: angles_file, angles_arg
        character(len=64)  :: angle_id
call get_command_argument(1,the_input)
if (len_trim(the_input)==0) error stop 'Supply explicit THE.in path as argument 1'
open(unit=101,file=trim(the_input),action='read',access='sequential',status='old')
read(101,*)!useless_string,func_dir
read(101,'(A9,A200)')useless_string,nscf_dir
read(101,'(A8,A200)')useless_string,THE_dir
read(101,*)!'(A7,I3)')  useless_string,n_func
read(101,*)!'(A13,I2)') useless_string,convert_size
read(101,*)
read(101,*)
read(101,'(A5,I3)') useless_string, nwan
read(101,'(A6,I2)') useless_string, nspin
read(101,'(A3,I3)') useless_string, nkx
read(101,'(A3,I3)') useless_string, nky
read(101,'(A3,I3)') useless_string, nkz
read(101,*) ! A6 axial
read(101,'(A9,I3)') useless_string,period_x
read(101,'(A9,I3)') useless_string, period_y
read(101,'(A9,I3)') useless_string, period_z
read(101,'(A)') useless_string
read(useless_string(index(useless_string,'=')+1:),*) n_electron
read(101,'(A4,F8.4)') useless_string, eta
read(101,'(A4,F16.8)') useless_string, fac
read(101,*)!'(A6,I3)') useless_string,b_min
read(101,*)!'(A6,I3)') useless_string,b_max
read(101,*)!'(A6,I3)') useless_string,wan_min
read(101,*)!'(A6,I3)') useless_string,wan_max
read(101,'(A3,F16.8)') useless_string, a0
read(101,'(A2,F16.8)') useless_string, celldm3
read(101,'(A3,F16.8)') useless_string, Ef
close(101)



               a(1,:) = (/   1.000000,   0.000000,   0.000000 /)
               a(2,:) = (/  0.000000,   1.000000,   0.000000 /)
               a(3,:) = (/   0.000000,   0.000000,   0.85231865  /)

               b(1,:) = (/  1.000000,  0.000000,  0.000000 /)
               b(2,:) = (/  0.000000,  1.000000,  0.000000 /)
               b(3,:) = (/  0.000000,  0.000000,  1.0/0.85231865 /)

a(:,:)=a(:,:)*a0*0.529177 !Ang

a(1,:)=a(1,:)!*period_x
a(2,:)=a(2,:)!*period_y
a(3,:)=a(3,:)!*period_z


E_max=3.0
E_min=0.0
nE=300
allocate(E_list(nE))
allocate(sig_H(2,2,nE))
allocate(sig_AH(2,2,nE))
allocate(sig_H_local(2,2,nE))
allocate(sig_AH_local(2,2,nE))
Do iE=1,nE
        E_list(iE)=E_min+(iE-1)*(E_max-E_min)/nE
ENDDO




cI=dcmplx(0.000,1.000)
pi=3.141592653500
I2pi=2*pi*cI
H_dim=nwan*2/nspin*period_x*period_y*period_z



if (.not. mod(nkx,period_x) .eq. 0) then
print*, 'not conmenserate'
stop
endif
if (.not. mod(nky,period_y) .eq. 0) then
print*, 'not conmenserate'
stop
endif
if (.not. mod(nkz,period_z) .eq. 0) then
print*, 'not conmenserate'
stop
endif
nkx=nkx/period_x
nky=nky/period_y
nkz=nkz/period_z

allocate(H_k(H_dim,H_dim))
allocate(H_buff(H_dim,H_dim))
!allocate(eig_vec_local(nkx*nky,H_dim,H_dim))
!allocate(eigen_local(nkx*nky,H_dim))
allocate(unmk(H_dim,H_dim))
allocate(eig_tmp(H_dim))
allocate(delHH(H_dim,H_dim,3))
!allocate(AA(H_dim,H_dim,3))
!allocate(AAq_tmp(H_dim,H_dim,3))
allocate(D_h(H_dim,H_dim,3))
allocate(mat_tmp(H_dim,H_dim))
allocate(unmk_Hermit(H_dim,H_dim))
allocate(spin_tmp(H_dim))
irept=0
open(unit=101,file=trim(THE_dir)//'6_RVR/super_H/ir_vec_HYC',action='read',access='sequential')
do
        read(101,'(3I8)',iostat=io) read3_int(:)
                if(io/=0) exit
                irept=irept+1
        enddo
close(101)
nrept=irept
allocate(ir_vec(nrept,3))
allocate(cr_vec(nrept,3))
allocate(H_wan(H_dim,H_dim,nrept))
allocate(H_wan_tmp(H_dim,H_dim))
allocate(RH_tmp(H_dim,H_dim,3))
!allocate(AAR(H_dim,H_dim,nrept,3))
!allocate(AA_tmp(H_dim,H_dim,3))
        !allocate(deg(nrept))
irept=1
open(unit=101,file=trim(THE_dir)//'6_RVR/super_H/ir_vec_HYC',action='read',access='sequential')
do
        read(101,'(3I8)',iostat=io) read3_int(:)
        if(io/=0) exit
        ir_vec(irept,:)=read3_int(:)
        irept=irept+1
enddo
close(101)

ir_vec(:,1)=ir_vec(:,1)*period_x
ir_vec(:,2)=ir_vec(:,2)*period_y
ir_vec(:,3)=ir_vec(:,3)*period_z


do irept=1,nrept
        cr_vec(irept,1)=sum(a(:,1)*ir_vec(irept,:))
        cr_vec(irept,2)=sum(a(:,2)*ir_vec(irept,:))
        cr_vec(irept,3)=sum(a(:,3)*ir_vec(irept,:))
enddo



call mpi_init(ierr)
call mpi_comm_rank(MPI_COMM_WORLD,my_rank,ierr)
call mpi_comm_size(MPI_COMM_WORLD,size,ierr)


nkx = 40
nky = 40
nkz = 40
nk = nkx * nky * nkz
npool = nk

allocate(k_list(nk,3))
allocate(pt_k(nk))
allocate(pool(nk,1))
allocate(eig_list(nk,H_dim))
allocate(spin_list(nk,H_dim))

ik = 0
do iz = 0, nkz-1
    do iy = 0, nky-1
        do ix = 0, nkx-1
            ik = ik + 1
            k_list(ik,1) = dble(ix) / nkx
            k_list(ik,2) = dble(iy) / nky
            k_list(ik,3) = dble(iz) / nkz
            pt_k(ik) = dble(ik)
            pool(ik,1) = ik
        end do
    end do
end do

k_list(:,1)=k_list(:,1)/period_x
k_list(:,2)=k_list(:,2)/period_y
k_list(:,3)=k_list(:,3)/period_z




super_size=period_x*period_y*period_z
print *, 'super_size ',super_size

call int2char(super_size,super_size_char)
allocate(nb(super_size,3))
allocate(phase0(super_size,super_size))
allocate(phase1(super_size,super_size))
allocate(phase_k(nrept))
allocate(R_list(super_size,3))
allocate(fold_tmp(H_dim,super_size))
allocate(fold(nk,H_dim,super_size))
do ix=0,period_x-1
do iy=0,period_y-1
do iz=0,period_z-1
        i=1+iz+period_z*iy+period_z*period_y*ix
        nb(i,:)=(/ix*1.000/period_x,iy*1.000/period_y,iz*1.000/period_z/)
        
        R_list(i,1)=sum(a(:,1)*(/ix,iy,iz/))
        R_list(i,2)=sum(a(:,2)*(/ix,iy,iz/))
        R_list(i,3)=sum(a(:,3)*(/ix,iy,iz/))


enddo
enddo
enddo

muti_pool=0
DO
        muti_pool=muti_pool+1
        IF(muti_pool*(size-1) .gt. npool) EXIT
ENDDO
muti_pool=muti_pool-1
write(*,*) 'multipool = ', muti_pool

allocate(eig_tmp_tot(muti_pool,H_dim))

if (my_rank.eq.0)  then

        open(unit=101,file='report_E',action='write',access='sequential',position='append',status='replace')
        write(101,*) 'start reading'
        close(101)
        
        do irept=1,nrept
        !print*,irept
                call int2char(irept, irept_char)
                open(unit=101,file='./tot_HR/hr_'//trim(irept_char),action='read'&
                        &,access='stream',form='unformatted')
                        read(101) H_wan(:,:,irept)
                close(101)
                call MPI_Bcast( H_wan(:,:,irept),H_dim**2 , MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD, ierr)
                call MPI_Barrier(  MPI_COMM_WORLD, ierr)
                open(unit=101,file='report_E',action='write',access='sequential',position='append',status='old')
                write(101,*) 'read HR irept=,',irept,'out of ',nrept, 'done'
                close(101)



        enddo

        eig_list(:,:) =0.000

        Do iz=1,1
                call int2char(iz-1,kz_char)
                do size_it=1,(size-1)

                        call MPI_RECV(eig_list(((muti_pool)*(size_it-1)+1):((muti_pool)*size_it),:),H_dim*(muti_pool),MPI_DOUBLE_PRECISION,&
                        &MPI_ANY_SOURCE,size_it,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
                        !call MPI_RECV(spin_list(ipool,:),H_dim,MPI_DOUBLE_PRECISION,&
                        !&MPI_ANY_SOURCE,ipool+npool,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
                        !call MPI_RECV(fold(ipool,:,:),H_dim*super_size,MPI_DOUBLE_PRECISION,&
                        !&MPI_ANY_SOURCE,ipool+2*npool,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
                enddo
                call MPI_Barrier(  MPI_COMM_WORLD, ierr)

                print *, 'eig_list(:,:) = ',eig_list(:,:)

                !Save eig_list results to a file for rank 0
                open(unit=103, file='eig_list_rank_0.dat', action='write', status='replace')
                write(103,*) 'Process rank: 0'
                write(103,*) 'eig_list dimensions: ', nk, 'x', H_dim
                do ik = 1, nk
                    do i_dim = 1, H_dim
                        write(103,*) ik, i_dim, eig_list(ik, i_dim)
                    end do
                end do
                close(103)

        enddo
        
        ! Sort all eigenvalues in eig_list in ascending order
        !call sort_all_eigenvalues(eig_list, nk, H_dim)
        
        !open(unit=101,file='./energy_out_'//trim(title),action='write',access='sequential',position='append',status='replace')
        !do i_dim=1,H_dim
        !do ik=1,nk
        !        if (abs(eig_list(ik,i_dim)-Ef) .gt. 0.25) cycle
        !        write(101,'(I8,6F18.6)') ik, k_list(ik,:),pt_k(ik),eig_list(ik,i_dim),spin_list(ik,i_dim)
        !enddo
        !        write(101,*)
        !enddo
        !close(101)


        !open(unit=101,file='./energy_out_'//trim(title)//'_fold_rank',action='write',access='sequential',position='append',status='replace')
        !do i_dim=1,H_dim
        !do ik=1,nk
               ! if (abs(eig_list(ik,i_dim)-Ef) .gt. 0.25) cycle
        !        write(101,'(I8,5F12.6,'//trim(super_size_char)//'F12.6)') ik, k_list(ik,:),&
        !                                &pt_k(ik),eig_list(ik,i_dim),fold(ik,i_dim,:)
        !enddo
        !        write(101,*)
        !enddo
        !close(101)
        
        occ_per_band = 1.0 
        ! n_electron is the TOTAL active count supplied by THE.in (legacy sharp filling).

        ! 1. Define search range based on eigenvalue extrema
        Ef_low = minval(eig_list)
        Ef_high = maxval(eig_list)

        print *, 'Ef_low = ', Ef_low
        print *, 'Ef_high = ', Ef_high
        ! 2. Tolerance and limit
        tol = 1e-4
        max_iter = 100
        iter = 0

        ! 3. Initial fill counts
        filling_low = count_filling(eig_list, nk, H_dim, Ef_low, occ_per_band)
        filling_high = count_filling(eig_list, nk, H_dim, Ef_high, occ_per_band)
  
        print *, 'filling_low  = ', filling_low 
        print *, 'filling_high = ', filling_high  

        if (n_electron < filling_low .or. n_electron > filling_high) then
        print *, 'ERROR: n_electron is outside the eigenvalue range.'
        stop
        end if

        ! 4. Bisection loop
        do while (abs(Ef_high - Ef_low) > tol .and. iter < max_iter)
        iter = iter + 1
        Ef_mid = 0.5 * (Ef_low + Ef_high)
        filling_mid = count_filling(eig_list, nk, H_dim, Ef_mid, occ_per_band)
        print *, 'Filling = ', filling_mid
        if (filling_mid < dble(n_electron)) then
                Ef_low = Ef_mid
        else
                Ef_high = Ef_mid
        end if
        end do

        Ef = 0.5 * (Ef_low + Ef_high)

        print *, 'Estimated Fermi level (Ef) = ', Ef
        print *, 'Filling = ', count_filling(eig_list, nk, H_dim, Ef, occ_per_band)

        open(unit=101, file='auto_Ef.dat', action='write', status='replace')
        write(101,'(A,F18.8)') 'Fermi energy (auto-determined) = ', Ef
        write(101,'(A,F18.8)') 'Filling = ', count_filling(eig_list, nk, H_dim, Ef, occ_per_band)
        write(101,'(A,I10)')   'Target n_electron = ', n_electron
        write(101,'(A,I5)')    'Bisection iterations = ', iter
        close(101)

        total_E = 0.0
        nk_total = nkx * nky * nkz
        weight = 1.0 / nk_total
        count_below_fermi = 0 

        do ik = 1, nk
        do i_dim = 1, H_dim
        if (eig_list(ik, i_dim) < Ef) then
                count_below_fermi = count_below_fermi + 1
                total_E = total_E + eig_list(ik, i_dim)
        end if 
        end do 
        end do
        print *, 'total_E before nomalization= ', total_E
        print *, 'counter ', count_below_fermi
        total_E = total_E * weight

        ! 結果をファイルに出力
        open(unit=101, file='sum_energy.dat', action='write', status='replace')
        write(101,'(A,F18.8)') 'Sum of eigenvalues over all energy below fermi level = ', total_E
        close(101)

endif



if (.not. my_rank.eq.0)  then

        do irept=1,nrept
                call MPI_Bcast( H_wan(:,:,irept),H_dim**2 ,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD, ierr)
                call MPI_Barrier(  MPI_COMM_WORLD, ierr)
        enddo

        Do iz=1,1!HYC_z,HYC_z
        DO imuti_pool=0,(muti_pool-1)
                ipool=imuti_pool*(size-1)+my_rank
                if (ipool .gt. npool) cycle
                call cpu_time(t1)
                ikxy=pool(ipool,iz)
                ik=ikxy!1+(iz-1)+(ikxy-1)*nkz

                do c_ix=0,period_x-1
                do c_iy=0,period_y-1
                do c_iz=0,period_z-1
                        i=1+c_iz+period_z*c_iy+period_z*period_y*c_ix
                do c_jx=0,period_x-1
                do c_jy=0,period_y-1
                do c_jz=0,period_z-1
                        j=1+c_jz+period_z*c_jy+period_z*period_y*c_jx
                        phase0(i,j)=exp(I2pi*&
                                &dot_product(1.000*(/c_jx-c_ix,c_jy-c_iy,c_jz-c_iz/),k_list(ik,:)))
                enddo
                enddo
                enddo
                enddo
                enddo
                enddo


                do n1=1,super_size
                do c_ix=0,period_x-1
                do c_iy=0,period_y-1
                do c_iz=0,period_z-1
                        i=1+c_iz+period_z*c_iy+period_z*period_y*c_ix
                        phase1(i,n1)=exp(-I2pi*dot_product(1.000*(/c_ix,c_iy,c_iz/),nb(n1,:)))
                enddo
                enddo
                enddo
                enddo

                
                H_wan_tmp(:,:)=dcmplx(0.000,0.000)
                do irept=1,nrept
                phase_k(irept)=exp(I2pi*sum(ir_vec(irept,:)*k_list(ik,:)))
                H_wan_tmp(:,:)=H_wan_tmp(:,:)+H_wan(:,:,irept)*phase_k(irept)
                enddo



                H_k(:,:)=dcmplx(0.000,0.000)
                delHH(:,:,:)=dcmplx(0.000,0.000)
                
                
                do i=1,super_size
                        i3=2/nspin*nwan*(i-1)+1
                        f3=2/nspin*nwan*i
                do j=1,super_size
                        ctmp1=phase0(i,j)
                        i4=2/nspin*nwan*(j-1)+1
                        f4=2/nspin*nwan*j
                do n1=1,super_size
                        ctmp2=phase1(i,n1)*ctmp1
                        i1=2/nspin*nwan*(n1-1)+1
                        f1=2/nspin*nwan*n1
                        
                do n2=1,super_size
                        ctmp3=conjg(phase1(j,n2))*ctmp2
                        i2=2/nspin*nwan*(n2-1)+1
                        f2=2/nspin*nwan*n2


                        H_k(i1:f1,i2:f2)=H_k(i1:f1,i2:f2)+&
                                &ctmp3*H_wan_tmp(i3:f3,i4:f4)

                enddo
                enddo
                enddo
                enddo
               
                H_k=H_k/super_size
                call cpu_time(t2)
                call  diasym(H_k,eig_tmp(:),H_dim,unmk(:,:))
        

                !spin_tmp(:)=0.000
                !do i_dim=1,H_dim
                !do j_dim=1,H_dim
                !        spin_tmp(i_dim)=spin_tmp(i_dim)&
                !                &+(-1.0)**(j_dim+1)*abs(unmk(j_dim,i_dim))**2
                !ENDDO
                !ENDDO
        
                !Fold_tmp(:,:)=0.0

               ! do i_dim=1,H_dim
              !  do i=1,super_size
                     !   i3=2/nspin*nwan*(i-1)+1
                     !  f3=2/nspin*nwan*i
                     !   Fold_tmp(i_dim,i)=sum(abs(unmk(i3:f3,i_dim))**2)
                !enddo
                !enddo

 
        !if (ipool .eq. nk) then
        !do i_dim=1,h_dim
        !print*, i_dim,abs(unmk(i_dim,1:4))**2 
        !enddo
        !endif                       

                
                call cpu_time(t4)

        if ( my_rank.eq.1)  then
                open(unit=101,file='report_E',action='write',access='sequential',position='append')
                        write(101,*) 'build', t2-t1, 'diag', t4-t2, 'all',t4-t1
                close(101)
        endif

       ! print*,'build', t2-t1, 'diag', t3-t2, 'dipole', t4-t3, 'all',t4-t1,sig_xy_local
         eig_tmp_tot(imuti_pool+1,:)=eig_tmp(:)

               ! call MPI_SEND(spin_tmp(:),H_dim&
                !                &,MPI_DOUBLE_PRECISION,0,ikxy+npool,MPI_COMM_WORLD,ierr)
                !call MPI_SEND(Fold_tmp(:,:),H_dim*super_size&
                !                &,MPI_DOUBLE_PRECISION,0,ikxy+2*npool,MPI_COMM_WORLD,ierr)



        enddo
        !print *, 'eig_tmp_tot(:,:) = ',eig_tmp_tot(:,:)

        ! Save eig_tmp_tot results to a file specific to this process
        call int2char(my_rank, i_char)
        !open(unit=102, file='eig_tmp_tot_rank_'//trim(i_char)//'.dat', action='write', status='replace')
        !write(102,*) 'Process rank: ', my_rank
        !write(102,*) 'eig_tmp_tot dimensions: ', muti_pool, 'x', H_dim
        !do i = 1, muti_pool
        !    do j = 1, H_dim
        !        write(102,*) i, j, eig_tmp_tot(i, j)
        !    end do
        !end do
        !close(102)

        call MPI_SEND(eig_tmp_tot(:,:),H_dim*muti_pool&
                        &,MPI_DOUBLE_PRECISION,0,my_rank,MPI_COMM_WORLD,ierr)

        call MPI_Barrier(  MPI_COMM_WORLD, ierr)
        enddo


endif



call mpi_finalize(ierr)

contains

subroutine sort_all_eigenvalues(eigenvals, nk, nband)
    implicit none
    real(kind(0d0)) , intent(inout) :: eigenvals(nk, nband)
    integer, intent(in) :: nk, nband
    real(kind(0d0)) , allocatable :: temp_array(:)
    integer :: i, j, idx, total_elements
    
    total_elements = nk * nband
    allocate(temp_array(total_elements))
    
    ! Convert 2D array to 1D for sorting
    idx = 1
    do i = 1, nk
        do j = 1, nband
            temp_array(idx) = eigenvals(i, j)
            idx = idx + 1
        end do
    end do
    
    ! Sort the 1D array using quick sort
    call quick_sort(temp_array, 1, total_elements)
    
    ! Convert back to 2D array
    idx = 1
    do i = 1, nk
        do j = 1, nband
            eigenvals(i, j) = temp_array(idx)
            idx = idx + 1
        end do
    end do
    
    deallocate(temp_array)
    
    print *, 'Eigenvalues sorted in ascending order'
    print *, 'Min eigenvalue: ', eigenvals(1, 1)
    print *, 'Max eigenvalue: ', eigenvals(nk, nband)
end subroutine sort_all_eigenvalues

recursive subroutine quick_sort(arr, left, right)
    implicit none
    real(kind(0d0)) , intent(inout) :: arr(:)
    integer, intent(in) :: left, right
    integer :: pivot_idx
    
    if (left < right) then
        pivot_idx = partition(arr, left, right)
        call quick_sort(arr, left, pivot_idx - 1)
        call quick_sort(arr, pivot_idx + 1, right)
    end if
end subroutine quick_sort

function partition(arr, left, right) result(pivot_idx)
    implicit none
    real(kind(0d0)) , intent(inout) :: arr(:)
    integer, intent(in) :: left, right
    integer :: pivot_idx, i, j
    real(kind(0d0)) :: pivot, temp
    
    pivot = arr(right)
    i = left - 1
    
    do j = left, right - 1
        if (arr(j) <= pivot) then
            i = i + 1
            temp = arr(i)
            arr(i) = arr(j)
            arr(j) = temp
        end if
    end do
    
    temp = arr(i + 1)
    arr(i + 1) = arr(right)
    arr(right) = temp
    
    pivot_idx = i + 1
end function partition

function count_filling(eigenvals, nk, nband, Ef, occ_per_band) result(filling)
    implicit none
    real(kind(0d0)) , intent(in) :: eigenvals(nk, nband), Ef
    integer, intent(in) :: nk, nband
    double precision :: filling, cnt, weight, occ_per_band
    integer :: ik, ib

    filling = 0.0
    cnt     = 0.0
    weight  = 1.0d0 / dble(nk)   ! ← 実数除算にする

    print *, 'nk = ', nk
    print *, 'weight = ', weight

    do ik = 1, nk
        do ib = 1, nband
            if (eigenvals(ik, ib) < Ef) then
                cnt     = cnt + occ_per_band
            end if
        end do
    end do
    filling = cnt * weight
    print *, 'cnt = ', cnt
end function count_filling

subroutine int2char(in_int,out_char)
        integer,intent(in)              :: in_int
        character(len=20),intent(out)    :: out_char

IF (in_int .LT. 10 )    WRITE(out_char,'(I1)')in_int
IF ((in_int .LT. 100)  .and. (in_int .GE. 10)   )    WRITE(out_char,'(I2)')in_int
IF ((in_int .LT. 1000) .and. (in_int .GE. 100)  )    WRITE(out_char,'(I3)')in_int
IF ((in_int .LT. 10000) .and. (in_int .GE. 1000) )   WRITE(out_char,'(I4)')in_int
IF ((in_int .LT. 100000) .and. (in_int .GE. 10000) ) WRITE(out_char,'(I5)')in_int

end subroutine int2char


subroutine diasym(a,eig,n,D)
 implicit none

 integer,intent(in)  :: n
 complex(kind(0d0)),intent(in)  :: a(n,n)
 complex(kind(0d0))              :: a_tmp(n,n)
 complex(kind(0d0)),intent(out)  :: D(n,n)

 integer  :: l,inf
 complex(kind(0d0))  :: Lwork(2*n-1)
 real(kind(0d0)) ,intent(out)      :: eig(n)
 real(kind(0d0)) :: Rwork(3*n-2)

 a_tmp(:,:)=a(:,:)

 l=2*n-1
 call zheev('V','U',n,a_tmp,n,eig,Lwork,l,Rwork,inf)
 if (inf/=0) error stop 'ZHEEV failed'
 D(:,:)=a_tmp(:,:)



end subroutine diasym



end program



