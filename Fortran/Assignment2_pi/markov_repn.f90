program Markov_chain_algorithm

	implicit none
	integer :: nseed, i, j, mov, hit, rep
	double precision :: rand, p, displacement, x, y, r, theta, dx, dy, distance, pi, pisum, squaresum, avg, var, std
	open(1, file='mn.inp', status='old')	 ! mov, displacement, rep : input
	open(2, file='change_of_delta_accuracy', status='unknown', position='append')
	nseed=2*int(secnds(0.0))+1235
	call srand(nseed)

	p=datan(1.d0)*4.d0
	read(1,*) mov, displacement, rep
	pisum=0.d0
	squaresum=0.d0
	x=rand()
	y=rand() 	! initial position

	do i=1, rep	
	   hit=0
	
	   do j=1, mov
              r=rand()*displacement ! 이동거리
	      theta=rand()*2.d0*p   ! 이동방향
	      dx=r*dcos(theta)
	      dy=r*dsin(theta)
	   
	      if ((x+dx)<0) dx=dx+1.d0  		! box 밖으로 이동 > pbc 적용
	      if ((x+dx)>1) dx=dx-1.d0
	      if ((y+dy)<0) dy=dy+1.d0
	      if ((y+dy)>1) dy=dy-1.d0

	      x=x+dx
	      y=y+dy

	      distance=dsqrt((x-0.5d0)*(x-0.5d0)+(y-0.5d0)*(y-0.5d0))         ! center에서부터의 거리

	      if (distance<=0.5) hit=hit+1

	   enddo

	   pi=dble(hit)/dble(mov)*4.d0
	   pisum=pisum+pi
	   squaresum=squaresum+pi*pi
if (mod(i,10000)==0) print*, i, 'th rep of',displacement, 'pi is', pi

	enddo

	avg=pisum/dble(rep)
	var=squaresum/dble(rep)-avg*avg
	std=dsqrt(var)
	write(*,*) 'pi 평균: ', avg, '표준편차: ', std

	write(2,*) displacement, avg, std 
end
