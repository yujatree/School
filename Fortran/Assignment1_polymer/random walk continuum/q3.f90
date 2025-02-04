program Distribution

        implicit none
        integer, dimension(:), allocatable:: histogram
        double precision :: n, bin_size, displ, displsum, squaresum, avg, var, std
        integer :: bin, i, j, m, num
        open(1, file='data2', status='unknown')         ! displacement file
        open(2, file='data3', status='unknown')         ! histogram 저장 file
        open(3, file='average', status='unknown')
        open(4, file='in.inp', status='unknown')        ! 'data 수=num, 간격=bin_size' input

        read(4,*) num, bin_size
        displsum=0.d0
        squaresum=0.d0
        displ=0.d0
        
        ! 가장 큰 displacement 값 받아서 Histogram 틀 만들기

        do i=1, num
           read(1, *) m, n
           if (n>=displ) then
              displ=n
           endif
        enddo    

        bin=int(displ/bin_size)+1                       ! bin: histogram 칸수 
                                                        ! (e.g. bin_size=5, then 0~5, 5~10, ... max-5~max)
        
        allocate(histogram(bin))
        histogram=0

        rewind(1)

        do i=1, num                                     ! num 개의 displacement 데이터 histogram화
           
           read(1, *) m, n
           displ=n
           
           n=int(n/bin_size)+1
           histogram(n)=histogram(n)+1
           
           displsum=displsum+displ                      ! displacement 합
           squaresum=squaresum+displ*displ              ! displacement 제곱의 합 
        
        enddo                                           ! histogram 완성
        
        do i=1, bin
           write(2,*) (2*i-1)*bin_size/2, histogram(i)/dble(num)/bin_size      ! probability distribution 저장
        enddo

        avg=displsum/dble(num)                          ! 평균
        var=squaresum/dble(num)-avg*avg                 ! 분산
        std=dsqrt(var)                                  ! 표준편차

        write(*,*) '평균: ', avg, '표준편차: ', std
        
        write(3,*) avg, maxval(histogram)/dble(num)/bin_size
        
end

