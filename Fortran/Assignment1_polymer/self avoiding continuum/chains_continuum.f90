program self_avoiding_walk

        implicit none
        double precision :: pi, theta, phi, r, rand, displ, dx, dy, dz, ddx, ddy, ddz, distance, mov_distance, overlap_distance
        integer :: i, j, k, nseed, num, monomer, n_overlap
	double precision, dimension(:), allocatable :: x, y, z
        open(1, file='self.inp', status='old')
	open(2, file='./50mer_data', status='unknown', position='append')
        nseed=2*int(secnds(0.0))+1235
        call srand(nseed)
        
	read(1,*) num, monomer, mov_distance, overlap_distance   ! chain 수, monomer 수, overlap 거리
        pi=datan(1.d0)*4.d0
        print*, num, monomer, mov_distance, overlap_distance

	allocate(x(monomer))
	allocate(y(monomer))
	allocate(z(monomer))
	
        do i=1, num

11	   x=0.d0 ; y=0.d0 ; z=0.d0
           n_overlap = 0

           do j=1, monomer-1
10            theta=rand()*pi
	      phi=rand()*2.d0*pi
              dx=mov_distance*dsin(theta)*dcos(phi)
              dy=mov_distance*dsin(theta)*dsin(phi)
	      dz=mov_distance*dcos(theta)

	      if (j == 1) then
		 x(j+1)=dx ; y(j+1)=dy ; z(j+1)=dz
	      endif

	      if (j .ne. 1) then
		 do k=1, j

		    ddx=x(k)-x(j)-dx
		    ddy=y(k)-y(j)-dy
		    ddz=z(k)-z(j)-dz
		    distance=(ddx*ddx+ddy*ddy+ddz*ddz)
		    
		    if (distance < overlap_distance*overlap_distance) then
                       n_overlap = n_overlap + 1
                       if ( n_overlap > 1000 ) then
			  print*, 'it got stuck...! retry~'
			  goto 11
		       endif
                       goto 10
                    endif
		    
	         enddo
	      
	         x(j+1)=x(j)+dx
	         y(j+1)=y(j)+dy
	         z(j+1)=z(j)+dz
	      endif 	     
           enddo

	   if (mod(i,10000)==0) write(*,*) i, 'th chain is made!'
	   
	   displ=dsqrt(x(monomer)*x(monomer)+y(monomer)*y(monomer)+z(monomer)*z(monomer))
	   write(2,*) displ
	
        enddo

end


