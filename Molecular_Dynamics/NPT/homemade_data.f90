program reader

	implicit none
	integer :: i, j, l, total_timestep, tau, num_of_steps, num_of_atoms, line, step, id, atom_type 
	double precision :: dx, dy, dz, distance, sum_of_v, fit_u, u, sum_of_u, &
			    temp, poteng, kineng, toteng, press, volume, density, &
			    box_length_x, box_length_y, box_length_z
	double precision, parameter :: cutoff=2.5d0, epsil=1.d0, sig=1.d0
	integer, dimension(:), allocatable :: timestep
	double precision, dimension(:,:), allocatable :: box_bounds_x, box_bounds_y, box_bounds_z, &
							 xu, yu, zu, vx, vy, vz, fx, fy, fz, &
							 v
	character(128) :: filename

	open(2, file='./homemade_log.lammps', status='unknown')

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

	allocate(v(num_of_steps, num_of_atoms))

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

	! Homemade log data--------------------------------------------------------------------------	
	
	write(2,'(A)') '#     Step      Time           Temp         PotEng         KinEng         &
			     &TotEng          Press           Volume        Density'

	
	fit_u=4.d0*epsil*((sig/cutoff)**12.d0-(sig/cutoff)**6.d0)
	print*, fit_u

	do i=1, num_of_steps

	   sum_of_v=0.d0
	   sum_of_u=0.d0

	   box_length_x=box_bounds_x(i,2)-box_bounds_x(i,1)
	   box_length_y=box_bounds_y(i,2)-box_bounds_y(i,1)
	   box_length_z=box_bounds_z(i,2)-box_bounds_z(i,1)

	   do j=1, num_of_atoms

	      v(i,j)=vx(i,j)*vx(i,j)+vy(i,j)*vy(i,j)+vz(i,j)*vz(i,j)
	      sum_of_v=sum_of_v+v(i,j)
	         
	      if (j .ne. 1) then
	         do l=1, j-1
		    dx=xu(i,l)-xu(i,j)
		    dy=yu(i,l)-yu(i,j)
	   	    dz=zu(i,l)-zu(i,j)
		    dx=nint(dx/box_length_x)*box_length_x-dx
		    dy=nint(dy/box_length_y)*box_length_y-dy
		    dz=nint(dz/box_length_z)*box_length_z-dz
		    distance=dsqrt(dx*dx+dy*dy+dz*dz)

		    if (distance<=cutoff) then
		       u=4.d0*epsil*((sig/distance)**12.d0-(sig/distance)**6.d0)-fit_u
		    else
		       u=0.d0
		    endif

		    sum_of_u=sum_of_u+u

	         enddo
	      endif

	   enddo

	   temp=sum_of_v/dble(3*num_of_atoms-3)
	   poteng=sum_of_u/dble(num_of_atoms)
!	   poteng=sum_of_u/dble(num_of_atoms)
	   kineng=sum_of_v/dble(num_of_atoms)/2.d0
	   toteng=poteng+kineng
	   volume=box_length_x*box_length_y*box_length_z
	   density=num_of_atoms/volume
	   press=0.d0

	   write(2,10) timestep(i), timestep(i)*0.005, temp, poteng, kineng, toteng, &
		       press, volume, density
10 	   format(I10, F10.1, F15.8, F15.8, F15.8, F15.8, F15.8, F17.6, F15.8) 

	enddo

	!--------------------------------------------------------------------------------------------	

end program	   
	
	
	     
	
