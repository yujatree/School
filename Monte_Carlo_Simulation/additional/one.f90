program ising_model

	implicit none
	integer :: nseed, lattice_size, sum_of_Sj, i, j, k, l, m, n, o, sum_of_spins, no_change
	double precision :: rand, r, energy_difference, magnetization, sum_of_some_M, E, sum_of_some_E
	integer, parameter :: rep=90000000, bin_size=10000000
	double precision, parameter :: kt=1.d0, coupling_constant= 1.d0, h=0.d0
	integer, dimension(:,:), allocatable :: lattice
	double precision, dimension(:), allocatable :: mean_of_M, mean_of_E
	character (len=128) :: dir, trial
	open(1, file='1000.input', status='old')
!	open(2, file='process', status='unknown')
!	open(3, file='M_of_t', status='unknown')	        ! M(t)
!	open(4, file='M_of_t_mean', status='unknown') 		! M(t)의 fluctuation 완화
!	open(5, file='E_of_t', status='unknown')
!	open(60, file='E_of_t_mean', status='unknown')
!	open(7, file='stop_by_M', status='unknown')
!	open(8, file='stop_by_E', status='unknown')
	nseed=2*int(secnds(0.0))+1235
	call srand(nseed)
	
	dir='/home/jina/monte_carlo/additional/data/animation/lowT/'
	
	read(1,*) lattice_size
	allocate(lattice(lattice_size,lattice_size))

	allocate(mean_of_M(rep/bin_size))
	allocate(mean_of_E(rep/bin_size))
	sum_of_some_M=0.d0
	sum_of_some_E=0.d0
	no_change=0

!------------------------------------------------------------------------
! Generate an initial configuration
	
	do i=1, lattice_size
	   do j=1, lattice_size
	      r=rand()
	      if (0.d0<=r .and. r<0.5d0) lattice(i,j)=1
	      if (0.5d0<=r) lattice(i,j)=-1    
	   enddo
	enddo
	
	open(10, file=trim(dir)//'00000000', status='unknown')
	write(10,'(1000(I2,X))', advance='no') lattice
!-------------------------------------------------------------------------
! Pick a random point

	do i=1, rep

	   ! 임의의 point = +1 or -1
	   m=int(rand()*lattice_size)+1
	   n=int(rand()*lattice_size)+1
   
	   ! Sum of Sj for first term-------------------------------------
	   
	   sum_of_Sj=0

	   ! up 
	   if (m==1) then
	      sum_of_Sj=sum_of_Sj+lattice(lattice_size, n)
	   else
	      sum_of_Sj=sum_of_Sj+lattice(m-1,n)
	   endif
	   ! down
	   if (m==lattice_size) then
	      sum_of_Sj=sum_of_Sj+lattice(1,n)
	   else
	      sum_of_Sj=sum_of_Sj+lattice(m+1,n)
	   endif
	   ! left
	   if (n==1) then
	      sum_of_Sj=sum_of_Sj+lattice(m,lattice_size)
	   else
	      sum_of_Sj=sum_of_Sj+lattice(m,n-1)
	   endif
	   ! right
	   if (n==lattice_size) then
	      sum_of_Sj=sum_of_Sj+lattice(m,1)
	   else
	      sum_of_Sj=sum_of_Sj+lattice(m,n+1)
	   endif	   

	   ! difference of Hamiltonian
	   energy_difference=2.d0*2.d0*dble(lattice(m,n))*(coupling_constant*dble(sum_of_Sj)+2.d0*h)

	   ! Acceptance of flipping--------------------------------------

	   r=rand()
	   if (r<exp(-energy_difference/kt)) then
	      lattice(m,n)=-lattice(m,n)   ! flipped!
	   endif

	   
	   ! Check Magnetization-----------------------------------------

	   if (mod(i,100)==0) then

	      sum_of_spins=0

	      do k=1, lattice_size
	         do l=1, lattice_size
	            sum_of_spins=sum_of_spins+lattice(k,l)
	         enddo
	      enddo

	      magnetization=dble(sum_of_spins)/dble(lattice_size*lattice_size)

	      if (mod(i,10000)==0) print*, i, 'th magnetization: ', magnetization
!	      write(3,*) i, magnetization		! 시간에 따른 Magnetization 변화 그래프

	   endif
	
	   ! Mean of Magnetization for check the equilibrium--------------
	   
	   sum_of_some_M=sum_of_some_M+magnetization
	   
	   if (mod(i,bin_size)==0) then
	      mean_of_M(i/bin_size)=sum_of_some_M/dble(bin_size)
	      print*, i, '까지의 평균 M: ', mean_of_M(i/bin_size)

!	      write(4,*) i-bin_size/2, mean_of_M(i/bin_size)		! 그래프를 위해 저장
	      
	      if (i .ne. bin_size) then
		 if (abs(mean_of_M(i/bin_size)-mean_of_M(i/bin_size-1))<0.0001) then
		    write(*,*) 'equilibration is completed by M!'
!		    write(7,*) i, mean_of_M(i/bin_size)
		    goto 10
	      	 endif
	      endif

	      sum_of_some_M=0.d0
	   endif
	   
	   ! Check the Energy of state--------------------------------------
	
	   if (mod(i,100)==0) then 

	      E=0.d0

	      do m=1, lattice_size
		 do n=1, lattice_size
	     	    
	   	    sum_of_Sj=0.d0

	   	    ! up 
	   	    if (m==1) then
	      	       sum_of_Sj=sum_of_Sj+lattice(lattice_size, n)
	   	    else
	   	       sum_of_Sj=sum_of_Sj+lattice(m-1,n)
	   	    endif
	   	    ! down
	   	    if (m==lattice_size) then
	      	       sum_of_Sj=sum_of_Sj+lattice(1,n)
	   	    else
	   	       sum_of_Sj=sum_of_Sj+lattice(m+1,n)
	   	    endif
	   	    ! left
	   	    if (n==1) then
	   	       sum_of_Sj=sum_of_Sj+lattice(m,lattice_size)
	   	    else
	   	       sum_of_Sj=sum_of_Sj+lattice(m,n-1)
	   	    endif
	   	    ! right
	   	    if (n==lattice_size) then
	   	       sum_of_Sj=sum_of_Sj+lattice(m,1)
	   	    else
	   	       sum_of_Sj=sum_of_Sj+lattice(m,n+1)
	   	    endif	   
		    
		    E=E+lattice(m,n)*(-0.5d0*coupling_constant*sum_of_Sj-h)
		 enddo
	      enddo
!	   if (mod(i,100000)==0) print*, i, 'th Energy: ', E
!	   write(5,*) i, E
	   endif
		    
	   ! Mean of Energy for check the equilibrium--------------------

	   sum_of_some_E=sum_of_some_E+E
	   
	   if (mod(i,bin_size)==0) then
	      mean_of_E(i/bin_size)=sum_of_some_E/dble(bin_size)
!	      print*, i, '까지의 평균 E: ', mean_of_E(i/bin_size)
	
!	      write(60,*) i-bin_size/2, mean_of_E(i/bin_size)
	
	      if (i .ne. bin_size) then
	         if (abs(mean_of_E(i/bin_size)-mean_of_E(i/bin_size-1))<1.d0) then
!		    print*, 'equilibration is completed by E!'
!		    write(8,*) i, mean_of_E(i/bin_size)
!		    goto 10
		 endif
	      endif
	
	      sum_of_some_E=0.d0
	   endif
	
	   ! Save for animation------------------------------------------
	   
	   if (mod(i,100000)==0) then
	      write(trial,1000) i
1000	      format(i8.8)      
!	      write(*,*) i, 'th lattice is saved for animation!'
	      open(11, file=trim(dir)//trim(trial), status='unknown')
	      write(11,'(1000(I2,X))', advance='no') lattice 
	   endif

		   
	enddo

!---------------------------------------------------------------------------------
! End of trials

10 	print*, 'Magnetization: ', magnetization

end
