program self_avoiding_walk
	use modules
        implicit none
        double precision :: rand, r, displ, displsum, squaresum, avg, var, std
        integer :: i, j, nseed, mov, num, x, y, dx, dy
	integer, dimension(:), allocatable :: histogram
	integer, dimension(:,:), allocatable :: lattice(:,:) 	! 2D lattice space
	double precision, dimension(:), allocatable :: displacement
        open(1, file='lattice.inp', status='old')
	open(2, file='lattice_displacement', status='unknown')
        nseed=2*int(secnds(0.0))+1235
        call srand(nseed)
 
	read(1,*) mov, num    	! movement 수 (monomer 수), polymer 수,  input

	allocate(lattice(mov*2, mov*2))  ! 최대 이동거리로 space 형성
	allocate(displacement(num))

	displsum=0.d0
	squaresum=0.d0

	do i=1, num

100	   x=mov ; y=mov ! 초기위치
	   lattice=0

	   do j=1, mov

	      lattice(x,y)=1
	   
10	      r=rand()
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
	         if (lattice(x+1,y)==1 .and. lattice(x-1,y)==1 .and. lattice(x,y+1) .and. lattice(x,y-1)) goto 100
	         goto 10 
	      endif

	      x=x+dx ; y=y+dy

           enddo
		     
	   displ=dsqrt(dble((x-mov)*(x-mov)+(y-mov)*(y-mov)))
	   displacement(i) = displ
write(2,*) displ
	   displsum=displsum+displ
	   squaresum=squaresum+displ*displ

	enddo

	avg=displsum/dble(num)
	var=squaresum/dble(num)-avg*avg
	std=dsqrt(var)

	write(*,*) '평균: ', avg, '표준 편차: ', std

	call make_histogram(num, displacement, histogram)

end

