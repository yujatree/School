program ising_model

	implicit none
	integer :: nseed, lattice_size, sum_of_Sj, i, j, k, l, m, n, o, p, sum_of_spins, no_change, ensemble_max(1), ensemble_min(1)
	double precision :: coupling_constant, h, rand, kt, r, energy_difference, magnetization, sum_of_ensemble, ensemble_average, sum_of_some_M
	integer , parameter :: rep=10000000, num_of_ensemble=10, bin_size=10000
	integer, dimension(:,:), allocatable :: lattice
	double precision, dimension(:), allocatable :: mean_of_M, M_of_ensemble
	open(1, file='input', status='old')
	open(2, file='/home/jina/monte_carlo/multi_ising/data/m.dat', status='unknown')
	
	nseed=2*int(secnds(0.0))+1235
	call srand(nseed)
	
	read(1,*) lattice_size, coupling_constant, h
	allocate(lattice(lattice_size,lattice_size))

	allocate(mean_of_M(rep/bin_size))
	allocate(M_of_ensemble(num_of_ensemble))
		
! For kt-M data----------------------------------------------------------

	do o=1000, 3000
           kt=0.001d0*dble(o)
           sum_of_ensemble=0.d0
    
	   do p=1, num_of_ensemble

	      sum_of_some_M=0.d0
	      !------------------------------------------------------------------------
	      ! Generate an initial configuration
	
100	      do i=1, lattice_size
		 do j=1, lattice_size
		    r=rand()
	      
		    if (0.d0<=r .and. r<0.5d0) lattice(i,j)=1
		    if (0.5d0<=r) lattice(i,j)=-1    
		 enddo
	      enddo

	      !-------------------------------------------------------------------------
	      ! Pick a random point

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

	   	 if (mod(i,100)==0) then

	   	    sum_of_spins=0

	   	    do k=1, lattice_size
	    	       do l=1, lattice_size
	                  sum_of_spins=sum_of_spins+lattice(k,l)
	               enddo
	   	    enddo

	   	    magnetization=dble(sum_of_spins)/dble(lattice_size*lattice_size) 

!		 if (mod(i,10000)==0) print*, i, 'th M for equilibration: ', magnetization
	         endif
		 
		 ! Mean of Magnetization of check the equilibrium----------------------
 
		 sum_of_some_M=sum_of_some_M + magnetization
		 
		 if (mod(i,bin_size)==0) then
		
		    mean_of_M(i/bin_size)=sum_of_some_M/dble(bin_size)
		    if (i .ne. bin_size) then
		       if (abs(mean_of_M(i/bin_size)-mean_of_M(i/bin_size-1))<0.0001) then
			  write(*,*) 'equilibration is completed!'
			  goto 10
		       endif
		    endif

		    sum_of_some_M=0.d0

		 endif
	      enddo
	  
10	      print*,  p, 'th ensemble의 M: ', magnetization

	      M_of_ensemble(p)=abs(magnetization)	
           enddo
	   
print*, M_of_ensemble
	   ensemble_max=maxloc(M_of_ensemble)
	   ensemble_min=minloc(M_of_ensemble)
	   M_of_ensemble(ensemble_max)=0.d0
	   M_of_ensemble(ensemble_min)=0.d0
	   ensemble_average=sum(M_of_ensemble)/(dble(num_of_ensemble)-2.d0)
print*, M_of_ensemble

           !---------------------------------------------------------------------------------
           ! End of trials

           print*, 'kt: ', kt,  'M: ', ensemble_average

           write(2,*) kt, ensemble_average

	enddo

end
