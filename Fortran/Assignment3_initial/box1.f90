program point_particles_in_a_box

        implicit none
        integer :: num, i, j, nseed, bond, overlap
        double precision, dimension(:), allocatable :: x, y, z
        double precision :: dx, dy, dz, distance, rand, length, rho
        open(1, file='box1', status='unknown')
	open(2, file='inp', status='old')
        nseed=2*int(secnds(0.0))+1235 
        call srand(nseed)
        
        read(2,*) length, rho	        	 ! inp : box 길이와 밀도
        
        bond=0
        overlap=0
        num=int(rho*length*length*length)        ! particle 수

        allocate(x(num))
	allocate(y(num))
	allocate(z(num))

        do i=1, num
           
	   x(i)=rand()*length              ! array에 random x, y, z 저장
           y(i)=rand()*length 
           z(i)=rand()*length
           
	   if (i .ne. 1) then
	      do j=1, i-1
	      
	         dx=x(j)-x(i)
		 dy=y(j)-y(i)
		 dz=z(j)-z(i)
		 distance=dsqrt(dx*dx+dy*dy+dz*dz)
		 
		 if (0.5d0<=distance .and. distance<1.5d0) then 
		    bond=bond+1
		 elseif (distance<0.5d0) then
		    overlap=overlap+1
		 endif
              enddo
	   endif

        enddo
        
        write(*,*) '# of bonds: ', bond, '# of overlaps: ', overlap

end
