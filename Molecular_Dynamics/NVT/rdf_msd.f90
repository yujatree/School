program reader

	implicit none
	integer :: i, j, k, l, m, n, total_timestep, tau, num_of_steps, num_of_atoms, line, step, id, atom_type, &
		   bin, equilstep
	double precision :: dx, dy, dz, distance, sum_of_v, avg_of_v, u, sum_of_u, &
			    temp, poteng, kineng, toteng, press, volume, density, &
			    box_length_x, box_length_y, box_length_z, &
			    bin_size, median, r, rho, rdf_avg, rdf_std, &
			    displ, msd, square_msd, msd_avg, msd_std
	double precision, parameter :: cutoff=2.5d0, epsil=1.d0, sig=1.d0, pi=datan(1.d0)*4.d0
	integer, dimension(:), allocatable :: timestep
	double precision, dimension(:), allocatable :: g, rdf, square_rdf
	double precision, dimension(:,:), allocatable :: box_bounds_x, box_bounds_y, box_bounds_z, &
							 xu, yu, zu, vx, vy, vz, fx, fy, fz
	character(128) :: filename

	open(2, file='./RDF.dat', status='unknown')
	open(3, file='./MSD.dat', status='unknown')
	!--------------------------------------------------------------------------------------------	

	filename='lj_NVT.lammpstrj'
	total_timestep=100000
	tau=100

	open(1, file='~/md/NVT/'//trim(filename), status='old')
	
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
	 
	! Let's get RDF!-----------------------------------------------------------------------------	

	equilstep=1
	bin_size=0.02d0 ; bin=500

	allocate(g(bin))
	allocate(rdf(bin))
	allocate(square_rdf(bin))

	rdf=0.d0
	square_rdf=0.d0

	do i=equilstep, num_of_steps

	   g=0.d0	

	   box_length_x=box_bounds_x(i,2)-box_bounds_x(i,1)
	   box_length_y=box_bounds_y(i,2)-box_bounds_y(i,1)
	   box_length_z=box_bounds_z(i,2)-box_bounds_z(i,1)
 	   rho=dble(num_of_atoms)/(box_length_x*box_length_y*box_length_z)
	
   	   do j=1, num_of_atoms
	      do l=1, j-1

	         if (j==1) goto 110
	         dx=xu(i,j)-xu(i,l)
	         dy=yu(i,j)-yu(i,l)
	         dz=zu(i,j)-zu(i,l)
	         dx=nint(dx/box_length_x)*box_length_x-dx
	         dy=nint(dy/box_length_y)*box_length_y-dy
	         dz=nint(dz/box_length_z)*box_length_z-dz
	         distance=dsqrt(dx*dx+dy*dy+dz*dz)

	         m=int(distance/bin_size)+1
                 g(m)=g(m)+2.d0

	      enddo
110	   enddo

	   g=g/dble(num_of_atoms)/rho

	   do j=1, bin
	      rdf(j)=rdf(j)+g(j)
	      square_rdf(j)=square_rdf(j)+g(j)*g(j)
	   enddo

	enddo

	do i=1, bin
	   median=dble(i-0.5d0)*bin_size
	   r=dble(i-1)*bin_size
	   rdf_avg=rdf(i)/dble(num_of_steps-equilstep+1)/(4.d0*pi*r*r*bin_size)
	   rdf_std=square_rdf(i)/dble(num_of_steps-equilstep+1)/(4.d0*pi*r*r*bin_size)**2.d0-rdf_avg*rdf_avg
	   write(2,*) median, rdf_avg, rdf_std
	enddo
	 
	! Let's get MSD!-----------------------------------------------------------------------------	

	do i=1, int(num_of_steps/3)

	   msd=0.d0
	   square_msd=0.d0

	   do j=1, num_of_steps-i

	      displ=0.d0

	      do k=1, num_of_atoms
	
		 dx=xu(i+j,k)-xu(j,k)
	   	 dy=yu(i+j,k)-yu(j,k)
		 dz=zu(i+j,k)-zu(j,k)
		 displ=displ+dx*dx+dy*dy+dz*dz

	      enddo
	    
	      displ=displ/dble(num_of_atoms)
	      msd=msd+displ
	      square_msd=square_msd+displ*displ

	   enddo

	   msd_avg=msd/dble(num_of_steps-i)
	   msd_std=square_msd/dble(num_of_steps-i)-msd_avg*msd_avg

	   write(3,*) timestep(i)*0.005, msd_avg, msd_std

	enddo
   
	!--------------------------------------------------------------------------------------------	

end program	   	
