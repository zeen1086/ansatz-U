program eval_joint
 use joint_model
 implicit none
 type(joint_potential)::m
 character(4096)::control,model,texture,prefix,path,summary
 character(12),parameter::output_names(4)=[character(12)::'V_fit','dV_fit','N_VU_fit','dN_dVU_fit']
 integer::u,ios,i,derivative_site=1
 real(dp)::ueff=1._dp,occupation_tolerance=1e-8_dp,e(3,4),axes(3,4),de(3,4),lo,hi,herm
 logical::enforce_occupation_bounds=.true.,exists
 real(dp),allocatable::field(:,:),dfield(:,:)
 complex(dp)::n(10,10,4),dn(10,10,4),vu(10,10,4),dvu(10,10,4)
 namelist /joint_eval_control/ ueff,derivative_site,occupation_tolerance,enforce_occupation_bounds
 call get_command_argument(1,control);call get_command_argument(2,model)
 call get_command_argument(3,texture);call get_command_argument(4,prefix)
 call require(command_argument_count()>=3.and.command_argument_count()<=4.and. &
  min(len_trim(control),len_trim(model),len_trim(texture))>0, &
  'Usage: eval_joint eval.nml V_func_tmp four_spins_axes.dat [output_prefix]')
 summary='eval_report.txt'
 if(len_trim(prefix)>0) then
  summary=trim(prefix)//'.txt';prefix=trim(prefix)//'.'
 endif
 open(newunit=u,file=trim(control),status='old',action='read');read(u,nml=joint_eval_control);close(u)
 call require(derivative_site>=1.and.derivative_site<=4,'derivative_site is1..4 in DFT Ni order')
 call require(ieee_is_finite(occupation_tolerance).and.occupation_tolerance>=0,'Invalid occupation tolerance')
 call read_joint(trim(model),m);allocate(field(product(m%grid),4),dfield(product(m%grid),4))
 open(newunit=u,file=trim(texture),status='old',action='read')
 do i=1,4
  read(u,*,iostat=ios) e(:,i),axes(:,i);call require(ios==0,'Texture needs four rows: ex ey ez axis_x axis_y axis_z')
 enddo
 close(u)
 call require(all(ieee_is_finite(axes)).and.maxval(abs(sum(axes**2,dim=1)-1))<1e-12_dp,'Rotation axes must be unit vectors')
 de=0;de(:,derivative_site)=cross3(axes(:,derivative_site),e(:,derivative_site))
 call joint_value(m,e,de,field,n,dfield,dn);call hubbard_value(n,dn,ueff,vu,dvu)
 call occupation_bounds(n,lo,hi,herm)
 if(enforce_occupation_bounds) call require(lo>=-occupation_tolerance.and.hi<=1+occupation_tolerance, &
  'N leaves occupation bounds at requested texture; no clipping/output; inspect model extrapolation')
 ! Check the entire output set before creating its first file.
 do i=1,5
  if(i<=4) then
   path=trim(prefix)//trim(output_names(i))
  else
   path=summary
  endif
  inquire(file=trim(path),exist=exists);call require(.not.exists,'Output exists: '//trim(path))
 enddo
 open(newunit=u,file=trim(prefix)//trim(output_names(1)),status='new',form='unformatted',convert='little_endian')
 write(u) int(m%grid,int32);write(u) field/ry_ev;close(u)
 open(newunit=u,file=trim(prefix)//trim(output_names(2)),status='new',form='unformatted',convert='little_endian')
 write(u) int(m%grid,int32);write(u) dfield/ry_ev;close(u)
 open(newunit=u,file=trim(prefix)//trim(output_names(3)),status='new',form='unformatted',convert='little_endian')
 write(u) int(m%atom_ids,int32),ueff;write(u) n;write(u) vu;close(u)
 open(newunit=u,file=trim(prefix)//trim(output_names(4)),status='new',form='unformatted',convert='little_endian')
 write(u) int(m%atom_ids,int32),int(derivative_site,int32);write(u) dn;write(u) dvu;close(u)
 open(newunit=u,file=trim(summary),status='new')
 write(u,'(a,es24.15)') 'Ueff_eV=',ueff;write(u,'(a,i0)') 'derivative_site_1based=',derivative_site
 write(u,'(a,2es24.15)') 'N_eigenvalue_range=',lo,hi;write(u,'(a,es24.15)') 'N_Hermiticity_error=',herm
 write(u,'(a,4es24.15)') 'field_component_RMS_eV=',sqrt(sum(field**2,dim=1)/size(field,1))
 write(u,'(a,4es24.15)') 'field_derivative_component_RMS_eV_rad=',sqrt(sum(dfield**2,dim=1)/size(field,1))
 write(u,'(a,4es24.15)') 'all_site_dN_Frobenius=',sqrt(sum(sum(abs(dn)**2,dim=1),dim=1))
 write(u,'(a)') 'N is dimensionless; all10x10 complex entries retained. VU=Ueff*(I/2-N), dVU=-Ueff*dN.'
 write(u,'(a)') 'fields/dfields record order: scalar,Bx,By,Bz; grid x fastest. Hubbard is separate and counted once.'
 write(u,'(a)') 'For full-model torque include scalar derivative and dN on every affected Ni, then project in the common gauge.'
 close(u)
 print '(a)', 'Joint fields, occupations, Hubbard potentials and all-site derivatives written.'
end program
