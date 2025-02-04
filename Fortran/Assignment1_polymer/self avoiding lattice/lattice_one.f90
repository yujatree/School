program self_avoiding_walk

        implicit none
        double precision :: rand, r
        integer :: i, nseed, mov, x, y, dx, dy
	integer, dimension(:,:), allocatable :: lattice	! 2D lattice space
        open(1, file='lattice.inp', status='old')
	open(2, file='real_polymer', status='unknown')	
        nseed=2*int(secnds(0.0))+1235
        call srand(nseed)
        
	read(1,*) mov     	! movement 수 (monomer 수) input

	allocate(lattice(mov*2, mov*2))  ! 최대 이동거리로 space 형성

	x=mov ; y=mov    	! 초기 위치
	lattice=0

	do i=1, mov

	   write(2,*) x-mov, y-mov

	   lattice(x,y)=1
	   
10	   r=rand()
	   if (0.d0<=r .and. r<0.25d0) then
    	      dx=1 ; dy=0
	   elseif (0.25d0<=r .and. r<0.5d0) then
              dx=-1 ; dy=0
	   elseif (0.5d0<=r .and. r<0.75d0) then
	      dx=0 ; dy=1
	   else
	      dx=0 ; dy=-1
	   endif
	
	   if (lattice(x+dx,y+dy)==1) then
	      if (lattice(x+1,y)==1 .and. lattice(x-1,y)==1 .and. lattice(x,y+1) .and. lattice(x,y-1)) then
		 write(*,*) 'got stuck!' ; goto 100
	      endif
	      goto 10 
	   endif

	   x=x+dx ; y=y+dy

        enddo
		     
100     write(*,*) 'movement stopped!'

end

