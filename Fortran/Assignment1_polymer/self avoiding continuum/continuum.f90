program self_avoiding_walk

        implicit none
        double precision :: pi, theta, phi, r, rand, dx, dy, dz, ddx, ddy, ddz, distance, bond_distance, overlap_distance
        integer :: i, j, k, nseed, monomer
	double precision, dimension(:), allocatable :: x, y, z
        open(1, file='one_self.inp', status='old')
	open(2, file='data1', status='unknown')
        nseed=2*int(secnds(0.0))+1235
        call srand(nseed)
        
	read(1,*) monomer, bond_distance, overlap_distance   ! monomer 수, bond_distance, overlap되는 거리
        pi=datan(1.d0)*4.d0

	allocate(x(monomer))
	allocate(y(monomer))
	allocate(z(monomer))
	
	x=0.d0 ; y=0.d0 ; z=0.d0
	write(2,*) x(1), y(1), z(1)
           
	do j=1, monomer-1

10         theta=rand()*pi
	   phi=rand()*2.d0*pi
           dx=bond_distance*dsin(theta)*dcos(phi)
           dy=bond_distance*dsin(theta)*dsin(phi)
	   dz=bond_distance*dcos(theta)
	      
	   if (j == 1) then
	      x(j+1)=dx ; y(j+1)=dy ; z(j+1)=dz
	   endif

	   if (j .ne. 1) then
	      
	      do k=1, j
		 ddx=x(k)-x(j)-dx
		 ddy=y(k)-y(j)-dy
		 ddz=z(k)-z(j)-dz
		 distance=dsqrt(ddx*ddx+ddy*ddy+ddz*ddz)
		 if (distance < overlap_distance) goto 10
	      enddo
	      
	      x(j+1)=x(j)+dx
	      y(j+1)=y(j)+dy
	      z(j+1)=z(j)+dz
	      
	   endif
	 	     
	   write(2,*) x(j+1), y(j+1), z(j+1)
        
        enddo

end

