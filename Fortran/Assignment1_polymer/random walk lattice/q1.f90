program A_point_in_a_space
        
       implicit none
       double precision :: j, rand
       integer :: i, nseed, x, y
       open(1, file='data1', status='unknown')
       nseed = 2*int(secnds(0.0))+1235
       call srand(nseed)

       x=0
       y=0      ! 점 처음 위치: (0,0)

       do i = 1, 10**6
          j = rand()
          if (j>=0.d0 .and. j<0.25d0) then
             x=x+1
          elseif (0.25d0<=j .and. j<0.5d0) then
             x=x-1
          elseif (0.5d0<=j .and. j<0.75d0) then
             y=y+1
          else
             y=y-1
          endif
          ! 1/4 확률로 (+x, -x, +y, -y) 이동

          write(1, *) x, y
      
       enddo
       
end   
