program combine
 use local_rvr
 implicit none
 character(4096)::path
 call require(command_argument_count()==1,'Usage: combine projection.nml')
 call get_command_argument(1,path);call combine_rows(trim(path))
end program
