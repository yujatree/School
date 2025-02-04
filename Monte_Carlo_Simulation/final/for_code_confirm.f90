program ising_model

	implicit none
	integer :: nseed, lattice_size, sum_of_Sj, i, j, k, l, m, n, o, p, sum_of_spins, no_change, ensemble_max(1), ensemble_min(1)
	double precision :: coupling_constant, h, rand, kt, r, energy_difference, magnetization, sum_of_ensemble, ensemble_average, sum_of_some_E, E
	integer , parameter :: rep=10000000, num_of_ensemble=100, bin_size=10000
	integer, dimension(:,:), allocatable :: lattice
	double precision, dimension(:), allocatable :: mean_of_E, M_of_ensemble
	open(1, file='input', status='old')   ! input of lattice_size, J, h
	open(2, file='/home/jina/monte_carlo/final/data/1_1.dat', status='unknown', position='append')
	
	nseed=2*int(secnds(0.0))+1235
	call srand(nseed)
	
	read(1,*) lattice_size, coupling_constant, h
	allocate(lattice(lattice_size,lattice_size))

	allocate(mean_of_E(rep/bin_size))
	allocate(M_of_ensemble(num_of_ensemble))
		
! For kt-M data----------------------------------------------------------

	do o=10, 50
	   
	   kt=0.1d0*dble(o)
           sum_of_ensemble=0.d0
    
	   do p=1, num_of_ensemble

	      sum_of_some_E=0.d0
	      !------------------------------------------------------------------------
	      ! Generate an initial configuration
	
100	      do i=1, lattice_size
		 do j=1, lattice_size
		    r=rand()
	      
		    if (0.d0<=r .and. r<0.7d0) lattice(i,j)=1
		    if (0.7d0<=r) lattice(i,j)=-1    
		 enddo
	      enddo

	      !-------------------------------------------------------------------------
	      ! Pick a random point, i

	      do i=1, rep

		 m=int(rand()*lattice_size)+1
		 n=int(rand()*lattice_size)+1
	  	 ! lattice(m,n): 임의의 point = +1 or -1

	         ! Sum of Sj for first term (considering pbc)
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

	   	 energy_difference=2.d0*dble(lattice(m,n))*(coupling_constant*dble(sum_of_Sj)+2.d0*h)

	   	 ! Acceptance of flipping----------------------------------------------

	   	 r=rand()
	   
	   	 if (r<exp(-energy_difference/kt)) then
	   	    lattice(m,n)=-lattice(m,n)    ! flipped!
	   	 endif

	   	 ! Check Magnetization------------------------------------------------

	   	 if (mod(i,10000)==0) then

	   	    sum_of_spins=0

	   	    do k=1, lattice_size
	    	       do l=1, lattice_size
	                  sum_of_spins=sum_of_spins+lattice(k,l)
	               enddo
	   	    enddo

	   	    magnetization=dble(sum_of_spins)/dble(lattice_size*lattice_size) 
		    print*, i, 'th M: ', magnetization

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
	   
	         endif
		    
	         ! Mean of Energy for check the equilibrium--------------------

	         sum_of_some_E=sum_of_some_E+E
	   
	         if (mod(i,bin_size)==0) then
	      	    mean_of_E(i/bin_size)=sum_of_some_E/dble(bin_size)
	
	            if (i .ne. bin_size) then
	               if (abs(mean_of_E(i/bin_size)-mean_of_E(i/bin_size-1))<0.01d0) then
		          print*, 'equilibration is completed by E!'
		          goto 10
		       endif
	            endif
	
	            sum_of_some_E=0.d0
	         endif
	      enddo

10	      print*, p, 'th ensemble의 M: ', magnetization

	      M_of_ensemble(p)=magnetization
	
           enddo
	   
	   ensemble_average=sum(M_of_ensemble)/(dble(num_of_ensemble))

           !---------------------------------------------------------------------------------
           ! End of trials

           print*, 'kt: ', kt,  'M: ', ensemble_average

           write(2,*) kt, ensemble_average

50	enddo

end program
