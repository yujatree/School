program n_points_moving

        implicit none
        double precision :: pi, theta, r, rand, displ, x, y
        integer :: i, j, nseed
        open(1, file='data2', status='unknown')
        nseed=2*int(secnds(0.0))+1235
        call srand(nseed)
        
        pi=datan(1.d0)*4.d0
        x=0.d0
        y=0.d0

        do i=1, 10**6                    ! 점 i 개

           do j=1, 10**6                 ! j번 움직임

              theta=rand()*2.d0*pi
              r=rand()

              x=x+(r*dcos(theta))
              y=y+(r*dsin(theta))

           enddo

           displ=dsqrt(x*x+y*y)          ! displacement
           write(1, *) i, displ          ! data: (점 번호, 변위)
           x=0.d0                        ! 자리 리셋
           y=0.d0

        enddo

end


