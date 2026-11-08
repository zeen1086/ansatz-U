program cal_RVR_split
 use joint_rvr
 implicit none
 character(4096)::control,column
 integer::i,ios
 call get_command_argument(1,control);call get_command_argument(2,column)
 call require(command_argument_count()>=1.and.command_argument_count()<=2, &
  'Usage: joint/cal_RVR_split projection.nml [spatial_orbital_1_to_120]; omitted selects all')
 i=0
 if(len_trim(column)>0) then
  read(column,*,iostat=ios) i;call require(ios==0,'Spatial orbital must be an integer')
  call require(i>=1.and.i<=120,'Explicit spatial orbital must be 1..120')
 endif
 call project_joint_row(trim(control),i)
end program
