program calculation_of_pi

        implicit none
        integer :: nseed, k, i, j, rep, hitin
        real(8) :: rand, pi, p, q, s, ss, avg, var, std
        open(1, file='in.inp', status='unknown')        ! 10^n개의 점으로 pi 값 받는 것을 몇회 반복할지
        open(2, file='pi_geometrical', status='unknown')
        nseed=2*int(secnds(0.0))+1235
        call srand(nseed)

        read(1,*) rep

        write(2,*) '#', 'average', 'standard deviation'

s=0.d0
           ss=0.d0
           avg=0.d0
           var=0.d0
           std=0.d0

           do i=1, 100

              hitin=0

              do j=1, 10**k

                 p=rand()
                 q=rand()

                 if (sqrt(p*p+q*q)<=1) then
                    hitin=hitin+1
                 endif

              enddo

              pi=4*dble(hitin)/(10**k)

              s=s+pi
              ss=ss+pi*pi

           enddo
           avg=s/dble(rep)
           var=ss/dble(rep)-avg*avg
           std=sqrt(var)

           write(2,*) 10**k, avg, std
end
