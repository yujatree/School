program Markov_chain_algorithm

	implicit none
	integer :: nseed, i, j, mov, hit, miss
	double precision :: rand, pi, displacement, x, y, r, theta, dx, dy, distance
	open(1, file='m.inp', status='old') ! box_size, mov, displacement
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
	y=rand() 		! initial position
	
write(2,*) 'initial position: (',x,',',y,')'
	
	do i=1, mov
           r=rand()*displacement ! 이동거리
	   theta=rand()*2.d0*pi  ! 이동방향
	   dx=r*dcos(theta)
	   dy=r*dsin(theta)
	   
	    
	   if (((x+dx)<0 .or. (x+dx)>1) .or. ((y+dy)<0 .or. (y+dy)>1)) then        ! box 밖 이동 > reject
write(2,*) x+dx, y+dy, '> oh, no... move is rejected!'
	      endif
	   if ((x+dx)<0) dx=dx+1.d0  		! box 나갈경우 pbc로 돌아오기
	   if ((x+dx)>1) dx=dx-1.d0
	   if ((y+dy)<0) dy=dy+1.d0
	   if ((y+dy)>1) dy=dy-1.d0

	   x=x+dx
	   y=y+dy

write(2,*) i,'th position: (', x,',', y,')'

	   distance=dsqrt((x-0.5d0)*(x-0.5d0)+(y-0.5d0)*(y-0.5d0))         ! center에서부터의 거리^2

write(2,*) 'distance: ', distance

	   if (distance<=0.5) then
	      write(3,*) x, y
	      hit=hit+1
write(2,*) '> hit!'

	   else
	      write(4,*) x, y
	      miss=miss+1
write(2,*) '> miss!'

	   endif

write(2,*) '..............'

	enddo

write(2,*) '# of hit: ', hit, '# of miss: ', miss
write(2,*) 'a ratio of hit to total trials: ', dble(hit)/dble(mov)	  
write(2,*) 'pi: ', dble(hit)/dble(mov)*4.d0
end
