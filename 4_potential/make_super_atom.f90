program make_super_atom
  use parameter_mod
        implicit none
        character(len=20) :: atom_name
        real(kind(0d0))    :: read3(3)
        integer :: ns1,ns2,ns3,i_atom
        integer :: is1,is2,is3,io

        ns1=p_atom_super_x
        ns2=p_atom_super_y
        ns3=p_atom_super_z
        open(unit=102,file='super_atom',action='write',access='append',status='replace')
        open(unit=101,file='prim_atom',action='read',access='sequential',status='old')
        do
                read(101,*,iostat=io) atom_name,read3(:)
                if (io/=0) exit
                do is3=0,ns3-1
                do is2=0,ns2-1
                do is1=0,ns1-1
                        write(102,'(A14,3F20.10)') atom_name,&
                                        &read3(1)/ns1+1.0/ns1*is1,&
                                        &read3(2)/ns2+1.0/ns2*is2,&
                                        &read3(3)/ns3+1.0/ns3*is3
                enddo
                enddo
                enddo
        enddo
        close(101)
        close(102)

end program

