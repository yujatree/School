program point_particles_in_a_box

        implicit none
        integer :: num, i, j, nseed, bond
        double precision, dimension(:), allocatable :: x,y,z
        double precision :: x, y, z, dx, dy, dz, distance, rand, length, rho
        open(1, file='box2', status='unknown')
	open(2, file='inp', status='old')
        nseed=2*int(secnds(0.0))+1235 
        call srand(nseed)
        
        read(2,*) length, rho		        ! inp : box 길이와 밀도
        
	bond=0
        num=int(rho*length*length*length)       ! particle 수

        allocate(x(num))
	allocate(y(num))
	allocate(z(num))

        do i=1, num
           
100        x(i)=rand()*length   		  ! data array에 random x, y, z 저장
           y(i)=rand()*length 
           z(i)=rand()*length
           
	   if (i .ne. 1) then
	      do j=1, i-1
	      
	         dx=x(j)-x(i)
		 dy=y(j)-y(i)
		 dz=z(j)-z(i)
		 distance=dsqrt(dx*dx+dy*dy+dz*dz)
		 
		 if (distance < 0.5d0) goto 100
	      	 if (distance < 1.5d0) bond=bond+1
              enddo
	   endif

        enddo
        
	write(*,*) '# of bond: ', bond

end


