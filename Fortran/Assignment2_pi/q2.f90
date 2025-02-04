program mensuration_by_divisioin

        implicit none
        double precision :: s, pi
        integer :: num, i, j
        open(1, file='pi', status='unknown')
        
        do i=1, 9

           s=0.d0
           num=10**i            ! 사각형의 수

           do j=1, num
              
              num=dble(num)
              s=s+(1.d0/num)*dsqrt(1.d0-((1.d0/num)*dble(j))**2.d0)

           enddo
        
           pi=s*4.d0

           write(1,*) num,  pi

        enddo
end



