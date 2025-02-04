program calculation_of_pi

        implicit none
        integer :: nseed, k, i, j, rep, hitin
        double precision :: rand, pi, p, q, pisum, squaresum, avg, var, std
        open(1, file='in.inp', status='unknown')        ! 10^n개의 점으로 pi 값 받는 것을 몇회 반복할지
        open(2, file='pi_geometrical', status='unknown')
        nseed=2*int(secnds(0.0))+1235
        call srand(nseed)

        read(1,*) rep      
        
        do k=2,8
           
           pisum=0.d0
           squaresum=0.d0

           do i=1, rep
              hitin=0

              do j=1, 10**k
                 p=rand()
                 q=rand()

                 if (dsqrt(p*p+q*q)<=1) then            ! 원 안에 들어올 경우 hit in!
                    hitin=hitin+1
                 endif   
              enddo
          
              pi=4*dble(hitin)/(10**k)
           
              pisum=pisum+pi
              squaresum=squaresum+pi*pi
           enddo
           
           avg=pisum/dble(rep)

           std=dsqrt(var)

           write(2,*) 10**k, avg, std
            
        enddo
        
end
