program rvr
 use local_rvr
 implicit none
 character(4096)::command,path,arg
 integer::argc,row,k,ios
 logical::keep
 argc=command_argument_count()
 call require(argc>=2,'Usage: rvr prepare|project|combine projection.nml [row [kz_index] | --keep-intermediates]')
 call get_command_argument(1,command);call get_command_argument(2,path)
 select case(trim(command))
 case('prepare')
  call require(argc==2,'Usage: rvr prepare projection.nml')
  call prepare_fields(trim(path))
 case('project')
  call require(argc<=4,'Usage: rvr project projection.nml [row [kz_index]]')
  row=0;k=0
  if(argc>=3) then
   call get_command_argument(3,arg);read(arg,*,iostat=ios) row
   call require(ios==0,'Invalid row');call require(row>=1,'Row must be positive')
  endif
  if(argc==4) then
   call get_command_argument(4,arg);read(arg,*,iostat=ios) k
   call require(ios==0,'Invalid kz index');call require(k>=1,'kz index must be positive')
  endif
  call project_row(trim(path),row,k,k)
 case('combine')
  call require(argc<=3,'Usage: rvr combine projection.nml [--keep-intermediates]')
  keep=argc==3
  if(keep) then
   call get_command_argument(3,arg)
   call require(arg=='--keep-intermediates','Unknown combine option')
  endif
  call combine_rows(trim(path),keep)
 case default
  call require(.false.,'Unknown rvr command: '//trim(command))
 end select
end program
