program A_point_in_a_space
        
       implicit none
       double precision :: pi, theta, r, rand, x, y
       integer :: i, nseed
       open(1, file='data1', status='unknown')
       nseed = 2*int(secnds(0.0))+1235
       call srand(nseed)

       pi=datan(1.d0)*4.d0
       x=0.d0
       y=0.d0                        ! 점 처음 위치: (0,0)

       do i = 1, 10**6
          theta = rand()*2.d0*pi     ! 이동방향(각도: 2pi x (0~1))
          r = rand()                 ! 이동거리

          x = x+(r*dcos(theta))
          y = y+(r*dsin(theta))

          write(1, *) x, y

       enddo
 
end 





