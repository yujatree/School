program reader

	implicit none
	integer :: i, total_timestep, tau, num_of_steps, num_of_atoms, line, step, id, atom_type 
	integer, dimension(:), allocatable :: timestep
	double precision, dimension(:,:), allocatable :: box_bounds_x, box_bounds_y, box_bounds_z, &
							 xu, yu, zu, vx, vy, vz, fx, fy, fz
	character(128) :: filename

	!--------------------------------------------------------------------------------------------	

!	write(*,*) 'The file you want to read : '
!	read(*,*) filename
!	write(*,*) 'Total time and tau'
!	read(*,*) total_timestep, tau

	filename='lj_NpT.lammpstrj'
	total_timestep=100000
	tau=100

	open(1, file='~/md/NPT/'//trim(filename), status='old')

	num_of_steps = total_timestep/tau+1

	! Read Number of Atoms-----------------------------------------------------------------------

	do i=1,4
	   if (i<4) read(1,*)
	   if (i==4) read(1,*) num_of_atoms
	enddo

	rewind(1)

	! Allocate arrays----------------------------------------------------------------------------

	allocate(timestep(num_of_steps))
	allocate(box_bounds_x(num_of_steps,2))
	allocate(box_bounds_y(num_of_steps,2))
	allocate(box_bounds_z(num_of_steps,2))

	allocate(xu(num_of_steps, num_of_atoms))
	allocate(yu(num_of_steps, num_of_atoms))
	allocate(zu(num_of_steps, num_of_atoms))
	allocate(vx(num_of_steps, num_of_atoms))
	allocate(vy(num_of_steps, num_of_atoms))
	allocate(vz(num_of_steps, num_of_atoms))
	allocate(fx(num_of_steps, num_of_atoms))
	allocate(fy(num_of_steps, num_of_atoms))
	allocate(fz(num_of_steps, num_of_atoms))

	
	! Write on arrays----------------------------------------------------------------------------

	do i=1, num_of_steps*(9+num_of_atoms) ! to read all lines

	   line = mod(i,117)

	   if (line==0) line=117
	   if (line==1) step = int(i/117)+1
	   if (line==2) read(1, *) timestep(step)
	
	   if (line==1 .or. line==3 .or. line==4 .or. line==5 .or. line==9) then
	      read(1,*)
	      goto 100
	   endif

	   if (line==6) read(1,*) box_bounds_x(step,1), box_bounds_x(step,2)
	   if (line==7) read(1,*) box_bounds_y(step,1), box_bounds_y(step,2)
	   if (line==8) read(1,*) box_bounds_z(step,1), box_bounds_z(step,2)
	   
	   if (line>9) then
	      read(1,*) id, atom_type, &
			xu(step,id), yu(step,id), zu(step,id), &
	      	        vx(step,id), vy(step,id), vz(step,id), &
			fx(step,id), fy(step,id), fz(step,id) 
	   endif      
100	enddo
	
	print*, 'complete reading all lines!'
	
	!--------------------------------------------------------------------------------------------	


end program	   
	
	
	     
	
