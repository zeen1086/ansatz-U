program sample_local_vfunc
 use local_vfunc
 implicit none
 type(local_potential)::m
 character(4096)::path,points,output
 integer::n,i,u,v,ios
 real(dp)::point(3)
 real(dp),allocatable::a(:),b(:)
 call require(command_argument_count()==3,'Usage: sample_local_vfunc model points output')
 call get_command_argument(1,path);call get_command_argument(2,points);call get_command_argument(3,output)
 call read_vfunc(trim(path),m,.true.);allocate(a(2*m%np+1),b(2*m%np+1))
 open(newunit=u,file=trim(points),status='old');read(u,*) n
 open(newunit=v,file=trim(output),status='new',access='stream',form='unformatted',convert='little_endian')
 do i=1,n
  read(u,*,iostat=ios) point;call require(ios==0,'Truncated sample points')
  call local_value(m,point,a);call plane_value(m,point,b);write(v) a,b
 enddo
 close(u);close(v)
end program
