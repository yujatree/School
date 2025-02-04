program Markov_chain_algorithm

	implicit none
	integer :: nseed, i, j, mov, hit, miss
	double precision :: rand, pi, displacement, x, y, r, theta, dx, dy, distance
	open(1, file='m.inp', status='old') 		! mov, displacement: input
	open(2, file='tragectory', status='unknown')
	open(3, file='hit', status='unknown')
	open(4, file='miss', status='unknown')
	nseed=2*int(secnds(0.0))+1235
	call srand(nseed)

	pi=datan(1.d0)*4.d0
	read(1,*) mov, displacement
	hit=0
	miss=0
	x=rand()
	y=rand()   ! initial position
	
	do i=1, mov
           r=rand()*displacement ! 이동거리
	   theta=rand()*2.d0*pi  ! 이동방향
	   dx=r*dcos(theta)
	   dy=r*dsin(theta)
	   
	   if (((x+dx)<0 .or. (x+dx)>1) .or. ((y+dy)<0 .or. (y+dy)>1)) then   ! box 밖 이동 > reject
	      dx=0.d0
	      dy=0.d0
	   endif
	  
	   x=x+dx
	   y=y+dy

	   distance=dsqrt((x-0.5d0)*(x-0.5d0)+(y-0.5d0)*(y-0.5d0))            ! center에서부터의 거리

	   if (distance<=0.5) then
	      write(3,*) x, y
	      hit=hit+1
	   else
	      write(4,*) x, y
	      miss=miss+1
	   endif

	enddo

	write(*,*) '# of hit: ', hit, '# of miss: ', miss
	write(*,*) 'a ratio of hit to total trials: ', dble(hit)/dble(mov)	  
	write(*,*) 'pi: ', dble(hit)/dble(mov)*4.d0

end
