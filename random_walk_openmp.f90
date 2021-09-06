! This csn be compiled in gfortran with integers promoted to 8byte
! which we can do with -fdefault-integer-8
! and for intel fortran is might be -integer-size 64
program main
	use omp_lib
    implicit none

    integer, parameter :: maxbinexp = 27, maxiterexp = 8, pexp=2, blocks=3
    integer :: i, j, k, n, it, fit, skip, nbins, sumint!, checksum
	!Specify integers of type 4 for MPI and to keep fin size down
	!the maximum number in a bin is then 2^31 = 2,147,483,648
    integer(kind=4) :: ierr, nproc, myid, maxproc
    integer(kind=4), dimension(:), allocatable :: fin
    double precision :: r, y, D, p, nbins_s, t1, t2, t3

	open (unit=23, file="infodim.out")

	nbins = pexp**maxbinexp
	allocate(fin(nbins)); 
	fin=0.d0

	do n=1,blocks
		t1 = omp_get_wtime()

		!Number of iterations
		fit = 10**maxiterexp

		! Initialize OpenMP
		!$omp parallel shared ( fin ) &
		!$omp private ( myid, r, y, j )
			myid = omp_get_thread_num()    
			nproc = omp_get_num_threads()
			maxproc = omp_get_max_threads()

			!Loop to get transformed trajectory
			y = dble(myid)/dble(nproc)
			print*, "ID = ", myid, " of ", nproc, & 
					" max possible", maxproc, "starting y=", y
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
				!Binning operation
				j = int(floor(y*nbins)) + 1
				!This ensures values are not missed
				!but slows the code by 2x
				!!!!$omp CRITICAL (STACKPROT) 
					  fin(j) = fin(j) + int(1, kind=4)
				!!!!$omp END CRITICAL (STACKPROT)

			enddo
		!$omp end parallel

		t2 = omp_get_wtime()
		print*, "maximum bin count = ", maxval(fin), &
				"time taken iterations = ", t2-t1 

		!If we store all bins we can keep averaging in powers of two  
		!to get fractal dimensions from one run at about 64 bins (maxbinexp-6)
		!or in powers of 3 down to 81 bins (maxbinexp-4)
		!which means we use 12/pexp = 6 or 4 for 2 or 3 respectively
		do j=0,maxbinexp-12/pexp !-1
			skip = pexp**j
			D = 0.d0
			nbins_s = dble(nbins)/dble(skip)
			!Loop over blocks of skip averaging 
			do i = 1,nbins,skip
				if (skip .gt. 2) then
					sumint = sum(int(fin(i:i+skip-1)))
				else
					sumint = sum(fin(i:i+skip-1))
				endif
				p = dble(sumint)/dble(fit*nproc*n)
				!checksum = checksum + sum(fin(i:i+skip-1))
				if(p.ne.0) D = D + p*log(p)/log(1.d0/nbins_s)

				!write(10000*pexp+(maxbinexp-j)*100+k,*) p

			enddo
			print*, nbins/skip,fit*nproc*n, sum(fin),1.d0/log(1.d0/nbins_s),D
			write(23,*) nbins/skip,fit*nproc*n, sum(fin),1.d0/log(1.d0/nbins_s),D
		enddo

		t3 = omp_get_wtime()
		print*, "Time taken summing ", t3-t2
		print*, "Time taken total ", t3-t1

	enddo

    deallocate(fin)
    close(23,status='keep')


end program
