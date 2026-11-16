program cal_RVR_split
 use local_rvr
 implicit none
 character(4096)::path,arg
 integer::column,k,ios
 call require(command_argument_count()==2.or.command_argument_count()==3, &
 'Usage: cal_RVR_split projection.nml row [kz_index]; row=0 prepares shared V1 spatial fields')
 call get_command_argument(1,path);call get_command_argument(2,arg)
 read(arg,*,iostat=ios) column;call require(ios==0.and.column>=0,'Invalid row')
 k=0
 if(command_argument_count()==3) then
  call get_command_argument(3,arg);read(arg,*,iostat=ios) k;call require(ios==0.and.k>=1,'Invalid kz index')
 endif
 if(column==0) then
  call prepare_fields(trim(path))
 else
  call project_row(trim(path),column,k,k)
 endif
end program
