
! This csn be compiled in gfortran with integers promoted to 8byte
! which we can do with -fdefault-integer-8
! and for intel fortran is might be -integer-size 64
program main
	use omp_lib
    implicit none

    integer, parameter :: maxbinexp = 17, maxiterexp = 9, pexp=3, blocks=10
    integer :: i, j, k, n, it, fit, skip, nbins, sumint, sumtotal
	!Specify integers of type 4 for MPI and to keep fin size down
	!the maximum number in a bin is then 2^31 = 2,147,483,648
    integer(kind=4) :: ierr, nproc, myid, maxproc
    integer(kind=4), dimension(:), allocatable :: fin
    double precision :: q, p, qnew, pnew, d, DI, y, delta, t1, t2, t3

	open (unit=23, file="qp.out")

	call cpu_time(t1)

	nbins = pexp**maxbinexp
	allocate(fin(nbins)); 
	fin=0.d0
	d = dsqrt(1.0d0/72.0d0)

   CALL OMP_SET_NUM_THREADS(6)

	do n=1,blocks
		t1 = omp_get_wtime()

		!Number of iterations
		fit = 10**maxiterexp

		! Initialize OpenMP
		!$omp parallel shared ( fin ) &
		!$omp private ( myid, q, p, qnew, pnew, y, j )
			myid = omp_get_thread_num()    
			nproc = omp_get_num_threads()
			maxproc = omp_get_max_threads()

			q = 0.d0 !dble(myid)/dble(nproc)
			p = 0.d0 !dble(nproc-myid)/dble(nproc)
			print*, "ID = ", myid, " of ", nproc, & 
					" max possible", maxproc, "starting q=", q, "starting p=", p 
            
	        do it=1,fit

		        if(q-p.lt.-4*d) qnew = + (11.0d0/ 6)*q - ( 7.0d0/ 6)*p + 14*d
		        if(q-p.lt.-4*d) pnew = - ( 7.0d0/ 6)*q + (11.0d0/ 6)*p - 10*d
		        if(q-p.ge.-4*d) qnew = + (11.0d0/12)*q - ( 7.0d0/12)*p -  7*d
		        if(q-p.ge.-4*d) pnew = - ( 7.0d0/12)*q + (11.0d0/12)*p -  1*d

		        q = qnew; p = pnew

		        y = 0.5d0*((q+p)/sqrt(2.d0)+1.d0)
		        j = int(floor(y*nbins)) + 1

				!$omp CRITICAL (STACKPROT) 
					  fin(j) = fin(j) + int(1, kind=4)
				!$omp END CRITICAL (STACKPROT)

	        enddo
		!$omp end parallel
	    t2 = omp_get_wtime()

        !Get number of iterations
        sumtotal = 0
	    do j=1,size(fin,1)
            sumtotal = sumtotal + fin(j)
        enddo

		print*, "maximum bin count = ", maxval(fin), & 
                "sum values=", sumtotal, &
				"time taken iterations = ", t2-t1 

		!If we store all bins we can keep averaging in powers of two  
		!to get fractal dimensions from one run at about 64 bins (maxbinexp-6)
		!or in powers of 3 down to 81 bins (maxbinexp-4)
		!which means we use 12/pexp = 6 or 4 for 2 or 3 respectively
		do j=0,maxbinexp-12/pexp 
			skip = pexp**j
			DI = 0.d0
			delta = 1.d0/(dble(nbins)/dble(skip))

			!Loop over blocks of skip averaging 
			do i = 1,nbins,skip
                sumint = 0
                do k=i,i+skip-1
        			sumint = sumint + fin(k)
                enddo
				p = dble(sumint)/dble(sumtotal)
				if(p.ne.0) DI = DI + p*log(p)/log(delta)

				!write(10000*pexp+(maxbinexp-j)*100+k,*) p

			enddo
			print*, maxbinexp-j, nbins/skip,fit*n*nproc, sumtotal,1.d0/log(delta),DI
			write(23,*) nbins/skip,fit*n*nproc, sumtotal,1.d0/log(delta),DI
		enddo

	    t3 = omp_get_wtime()
        print*, "time taken average = ", t3-t2 
        print*, "time taken block = ", t3-t1 
    enddo
	close(23)


end program main
