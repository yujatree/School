program initial_making

        implicit none
        integer :: nseed, i, j, num_atoms, num_of_types
        double precision :: rand, density, overlap_distance, box_size, dx, dy, dz 
	integer, dimension(:), allocatable :: types
        double precision, dimension(:), allocatable :: mass, x, y, z
	open(1, file='initial.input', status='old')
	open(2, file='initial', status='unknown')
        nseed=2*int(secnds(0.0))+1235 
        call srand(nseed)
        
        read(1,*) num_atoms, density, overlap_distance

	box_size=(dble(num_atoms)/density)**(1.d0/3.d0)

        allocate(x(num_atoms))
	allocate(y(num_atoms))
	allocate(z(num_atoms))
	allocate(types(num_atoms))   	 		! 각 atom들의 type
	
	! temporary input for types of atoms-----------------------------------------------------
	num_of_types=1		     	 		! type 갯수
	
	allocate(mass(num_of_types))	 		! 각 type들의 mass 값
	mass=1.d0
	
	types=1

        ! writing formats for lammps-------------------------------------------------------------

	write(2,'(A)') '# LAMMPS DATA FILE'
	write(2,*)

	write(2,'(I0, A)') num_atoms, ' atoms'
	write(2,'(I0, A)') num_of_types, ' atom types' 
	write(2,*)

	write(2,10) 0.d0, box_size, 'xlo xhi'
	write(2,10) 0.d0, box_size, 'ylo yhi'
	write(2,10) 0.d0, box_size, 'zlo zhi'
	write(2,*)
10 	format(F12.6, F12.6,' ', A)
	
	write(2,'(A)') 'Masses'
	write(2,*)

	do i=1, num_of_types
      	   write(2,'(I4,F12.6)') i, mass(i)
	enddo
	write(2,*)

	write(2,'(A)') 'Atoms'
	write(2,'(A)') '# id type   x   	y	    z   ' 

	! create new atoms (considering overlap)------------------------------------------------	
        do i=1, num_atoms
           
100        x(i)=rand()*box_size   		
           y(i)=rand()*box_size 
           z(i)=rand()*box_size
           
	   if (i .ne. 1) then
	      do j=1, i-1

		 dx=x(j)-x(i)	      	
		 dy=y(j)-y(i)
		 dz=z(j)-z(i)
	         dx=nint(dx/box_size)*box_size-dx
		 dy=nint(dy/box_size)*box_size-dy
		 dz=nint(dz/box_size)*box_size-dz
		 
		 if (dx*dx+dy*dy+dz*dz < overlap_distance*overlap_distance) goto 100
	   
              enddo
	   endif

	   write(2,'(I4, I4, F12.6,F12.6,F12.6)') i, types(i), x(i), y(i), z(i) 
print*, i, 'th atoms is done!'
        enddo
        
end program


