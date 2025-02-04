program ising_model

	implicit none
	integer :: nseed, lattice_size, sum_of_Sj, i, j, k, l, m, n, o, p, sum_of_spins, no_change, ensemble_max(1), ensemble_min(1)
	double precision :: coupling_constant, h, rand, kt, r, energy_difference, magnetization, sum_of_ensemble, ensemble_average, sum_of_some_E, E, avg_of_Sj
	integer , parameter :: rep=1000000, num_of_ensemble=100, bin_size=10000
	integer, dimension(:,:), allocatable :: lattice
	double precision, dimension(:), allocatable :: mean_of_E, M_of_ensemble
	open(1, file='input', status='old')
	open(2, file='/home/jina/monte_carlo/mean_field/multi_data/upper', status='unknown')
	
	nseed=2*int(secnds(0.0))+1235
	call srand(nseed)
	
	read(1,*) lattice_size, coupling_constant, h
	allocate(lattice(lattice_size,lattice_size))

	allocate(mean_of_E(rep/bin_size))
	allocate(M_of_ensemble(num_of_ensemble))
		
! For kt-M data----------------------------------------------------------

	do o=10, 50
!	   if (20<=o .and. o<24) then:
!	      kt=0.01d0*dble(10*o)
!	      print*, 'kt: ', kt
!	   else
	      kt=0.1d0*dble(o)
	      print*, 'kt: ', kt, '----------------------------------------'
!	   endif

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
	      ! Pick a random point

	      do i=1, rep

		 m=int(rand()*lattice_size)+1
		 n=int(rand()*lattice_size)+1
	  	 ! lattice(m,n): 임의의 point = +1 or -1

	         ! Sum of Sj for avg of Sj
	         sum_of_Sj=0

	         do k=1, lattice_size
	            do l=1, lattice_size
		       sum_of_Sj=sum_of_Sj+lattice(k,l)
	            enddo
	         enddo
	
	         avg_of_Sj=dble(sum_of_Sj-lattice(m,n))/dble(lattice_size*lattice_size-1)

	         ! difference of Hamiltonian
	         energy_difference=dble(lattice(m,n))*(8.d0*coupling_constant*avg_of_Sj+2.d0*h)

	   	 ! Acceptance of flipping----------------------------------------------

	   	 r=rand()

	   	 if (r<exp(-energy_difference/kt)) then
	   	    lattice(m,n)=-lattice(m,n)    ! flipped!

	   	 endif

	   	 ! Check Magnetization------------------------------------------------

	   	 if (mod(i,1000)==0) then

	   	    sum_of_spins=0

	   	    do k=1, lattice_size
	    	       do l=1, lattice_size
	                  sum_of_spins=sum_of_spins+lattice(k,l)
	               enddo
	   	    enddo

	   	    magnetization=dble(sum_of_spins)/dble(lattice_size*lattice_size) 

		 if (mod(i,100000)==0) print*, i, 'th M for equilibration: ', magnetization
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
	      
!  	            if (mod(i,100000)==0) print*, i, 'th Energy: ', E
	   
	         endif
		    
	         ! Mean of Energy for check the equilibrium--------------------

	         sum_of_some_E=sum_of_some_E+E
	   
	         if (mod(i,bin_size)==0) then
	      	    mean_of_E(i/bin_size)=sum_of_some_E/dble(bin_size)
!	            print*, i, '까지의 평균 E: ', mean_of_E(i/bin_size)
	
	            if (i .ne. bin_size) then
	               if (abs(mean_of_E(i/bin_size)-mean_of_E(i/bin_size-1))<10.d0) then
!		          print*, 'equilibration is completed by E!'
		          goto 10
		       endif
	            endif
	
	            sum_of_some_E=0.d0
	         endif
	      enddo

10	      print*, p, 'th ensemble의 M: ', magnetization

	      M_of_ensemble(p)=magnetization
	
           enddo
	   
!	   print*, M_of_ensemble
!	   ensemble_max=maxloc(M_of_ensemble)
!	   ensemble_min=minloc(M_of_ensemble)
!	   M_of_ensemble(ensemble_max)=0.d0
!	   M_of_ensemble(ensemble_min)=0.d0
	   ensemble_average=sum(M_of_ensemble)/(dble(num_of_ensemble))
!	   print*, M_of_ensemble
	  

           !---------------------------------------------------------------------------------
           ! End of trials

           print*, 'kt: ', kt,  'M: ', ensemble_average

           write(2,*) kt, ensemble_average

	enddo

end
