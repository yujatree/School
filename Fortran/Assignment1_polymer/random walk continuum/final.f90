program n_points_moving

        implicit none
        double precision, allocatable, dimension(:) :: Displacement  ! displacement       
        integer, allocatable, dimension(:) :: Histogram              ! histogram
        double precision :: pi, theta, r, rand, displ, x, y, bin_size, displsum, squaresum, avg, var, std
        integer :: i, j, nseed, num, mov, bin, p
        open(1, file='inx.inp', status='unknown')       ! 점 갯수, 움직임 수, 간격 정보 input
        open(2, file='datax', status='unknown')         ! histogram 데이터 (for graph)
        open(3, file='averagex', status='unknown')
        nseed=2*int(secnds(0.0))+1235
        call srand(nseed)

        read(1,*) num, mov, bin_size
        pi=datan(1.d0)*4.d0
        x=0.d0
        y=0.d0
        displsum=0.d0
        squaresum=0.d0

        allocate(Displacement(num))

        do i=1, num       ! 점 num 개
           do j=1, mov    ! mov번 움직임
              
              theta=rand()*2.d0*pi
              r=rand()
              x=x+(r*dcos(theta))
              y=y+(r*dsin(theta))

           enddo
           displ=dsqrt(x**2.d0+y**2.d0)                ! displacement
           Displacement(i)=displ                       ! data array에 displacement 저장 
           x=0.d0                                      ! 자리 리셋
           y=0.d0
        enddo

!-------------------------------------------------------------------------------------------------------
        ! 히스토그램

        bin=int(maxval(Displacement)/bin_size)+1                  ! histogram 칸 수

        allocate(Histogram(bin))
        Histogram=0

        do i=1,num
           
           p=int(Displacement(i)/bin_size)+1
           Histogram(p)=Histogram(p)+1
            
        
           displsum=displsum+Displacement(i)
           squaresum=squaresum+Displacement(i)*Displacement(i)

        enddo  
        
        do i=1, bin
           write(2,*) (2*i-1)*bin_size/2, Histogram(i)/dble(num)/bin_size   ! 그래프를 위한 histogram (x,y) 저장
        enddo
        
!--------------------------------------------------------------------------------------------------------
        ! 평균과 표준편차

        avg=displsum/dble(num)                                     ! 평균
        var=squaresum/dble(num)-avg*avg                            ! 분산
        std=sqrt(var)                                              ! 표준편차

        write(*,*) '평균: ', avg, '표준 편차: ', std

        write(3,*) avg, maxval(Histogram)/dble(num)/bin_size

end


