program main
	implicit none

	integer, parameter :: maxbinexp = 32, maxiterexp = 8
	integer :: i, j, it, jazz, k, fit, skip, nbins, sumint!, checksum
	double precision :: r, y, D, p, res
	!Using 2 bit integers works as long as less than 32,767
	!values per bin (roughly 1e9 iterations)
	integer, dimension(:), allocatable :: fin

	if (kind(fin) .eq. 2 .and. 10**maxiterexp > 1e8) then
		stop "Use 4 bit integers with more than 1e8 iterations"
	endif

	open (unit=23, file="infodim.out")

	nbins = 2**maxbinexp
	allocate(fin(nbins)); 
	do k = 6,maxiterexp
		
		!Number of iterations
		fit = 10**k
		fin=0.d0

		!Loop to get transformed trajectory
		y = 0.5d0
		r = 0.5d0
		do it = 1,fit
			!Randomly transform
			call random_number(r)
			if(r.lt.1.d0/3.d0) then
				y = (1.d0+y+y)/3.d0
			else if (r.gt.1.d0/3.d0) then
				y = y/3.d0
			endif
			!write(123,*) it, r,y
			!Binning operat	ion
			j = floor(y*nbins) + 1
			!print*, "after j", r, y, j, size(fin)
			fin(j) = fin(j) + 1
		enddo

		print*, "maximum bin count = ", maxval(fin)
		do i =1,size(fin,1)
			write(100+k,*) i, dble(fin(i))/dble(fit)
		enddo

		!If we store all bins we can keep averaging in powers of two  
		!to get fractal dimensions from one run at about 64 bins (maxbinexp-6)
		do j=0,maxbinexp-6 !-1
			skip = 2**j
			D = 0.d0
			res = dble(nbins)/dble(skip)
			!Loop over blocks of skip averaging 
			do i = 1,nbins,skip
				if (skip .gt. 2) then
					sumint = sum(int(fin(i:i+skip-1)))
				else
					sumint = sum(fin(i:i+skip-1))
				endif
				p = dble(sumint)/dble(fit)
				!checksum = checksum + sum(fin(i:i+skip-1))
				if(p.ne.0) D = D + p*dlog(p)/dlog(1.d0/res)
			enddo
			!print*, skip, nbins, checksum, fit
			print*, nbins/skip,fit,1.d0/dlog(1.d0/res),D
			write(23,*) k, j, D
		enddo

	enddo
	deallocate(fin)
	close(23,status='keep')

end program
