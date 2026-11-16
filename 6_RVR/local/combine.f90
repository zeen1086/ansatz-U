program combine
 use local_rvr
 implicit none
 character(4096)::path,arg
 logical::keep_intermediates
 call require(command_argument_count()==1.or.command_argument_count()==2, &
 'Usage: combine projection.nml [--keep-intermediates]')
 keep_intermediates=command_argument_count()==2
 if(keep_intermediates) then
  call get_command_argument(2,arg)
  call require(arg=='--keep-intermediates','Unknown combine option')
 endif
 call get_command_argument(1,path);call combine_rows(trim(path),keep_intermediates)
end program
