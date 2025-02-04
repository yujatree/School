program n_points_moving

        implicit none
        double precision :: r, rand, displ
        integer :: i, j, nseed, x, y
        open(1, file='data2', status='unknown')
        nseed=2*int(secnds(0.0))+1235
        call srand(nseed)
        
        x=0
        y=0

        do i=1, 10**6                    ! 점 i 개

           do j = 1, 10**6               ! j 번 움직임
              r = rand()
              
              if (r>=0.d0 .and. r<0.25d0) then
                 x=x+1
              elseif (0.25d0<=r .and. r<0.5d0) then
                 x=x-1
              elseif (0.5d0<=r .and. r<0.75d0) then
                 y=y+1
              else
                 y=y-1
              endif

           enddo

           displ=dsqrt(dble(x*x+y*y))    ! displacement
           write(1, *) i, displ          ! data: (점 번호, 변위)
           x=0                           ! 자리 리셋
           y=0

        enddo

end


