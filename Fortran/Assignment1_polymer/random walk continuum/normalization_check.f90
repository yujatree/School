program normalization_check

	implicit none
	integer :: i, bin
	double precision :: x, p, total, bin_size

	open(1, file='data3', status='old')

	write(*,*) '# of bins, bin size: '
	read(*,*) bin, bin_size
	
	total=0.d0
	p= 0.d0
	do i=1, bin
	   read(1,*) x, p
	   total=total+p
	enddo
	total=total*bin_size
print*, 'sum of probability: ', total

end
	
