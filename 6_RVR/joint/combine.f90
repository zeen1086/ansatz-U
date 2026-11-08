program combine
 use joint_rvr
 implicit none
 character(4096)::control
 call get_command_argument(1,control)
 call require(command_argument_count()==1,'Usage: joint/combine projection.nml')
 call combine_joint(trim(control))
end program
