module fit_mask
 implicit none
contains
 logical function update_allowed(i,j,iteration,diagonal_iters)
 integer,intent(in)::i,j,iteration,diagonal_iters
 update_allowed=iteration>diagonal_iters.or.i==j
 end function
end module
