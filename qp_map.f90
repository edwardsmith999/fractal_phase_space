
! This csn be compiled in gfortran with integers promoted to 8byte
! which we can do with -fdefault-integer-8
! and for intel fortran is might be -integer-size 64
program main
    implicit none

    integer, parameter :: maxbinexp = 5, maxiterexp = 5, pexp=2
    integer :: i, j, k, n, it, fit, skip, nbins, sumint!, checksum
	!Specify integers of type 4 for MPI and to keep fin size down
	!the maximum number in a bin is then 2^31 = 2,147,483,648
    integer(kind=4) :: ierr, nproc, myid, maxproc
    integer(kind=4), dimension(:,:), allocatable :: fin
    double precision :: q, p, qnew, pnew, d, prob, x,y

	open (unit=23, file="qp.out")
	open (unit=24, file="qpbins.out")
	open (unit=25, file="qptraj.out")

	nbins = pexp**maxbinexp
	allocate(fin(nbins, nbins)); 
	fin=0.d0

	fit = 10**maxiterexp
	d = dsqrt(1.0d0/72.0d0)
	q = 0.d0; p = 0.d0
	do it=1,fit

		if(q-p.lt.-4*d) qnew = + (11.0d0/ 6)*q - ( 7.0d0/ 6)*p + 14*d
		if(q-p.lt.-4*d) pnew = - ( 7.0d0/ 6)*q + (11.0d0/ 6)*p - 10*d
		if(q-p.ge.-4*d) qnew = + (11.0d0/12)*q - ( 7.0d0/12)*p -  7*d
		if(q-p.ge.-4*d) pnew = - ( 7.0d0/12)*q + (11.0d0/12)*p -  1*d

		q = qnew; p = pnew
		write(23,*) q, p

		i = int(floor(0.5*(q+sqrt(2.d0))*nbins)) + 1
		j = int(floor(0.5*(p+sqrt(2.d0))*nbins)) + 1

		x = sqrt(1.d0/2.d0)*(q-p)
		y = sqrt(1.d0/2.d0)*(q+p)

		i = int(floor(0.5*(x+1.d0)*nbins)) + 1
		j = int(floor(0.5*(y+1.d0)*nbins)) + 1

		!write(25,*), q, p
		!print*, x, y, i, j, nbins

		fin(i,j) = fin(i,j) + int(1, kind=4)
	enddo


	D = 0.d0
	do i =1,nbins
	do j =1,nbins
		prob = dble(fin(i,j))/dble(fit)
		!PRINT*, I,J,PROB,d,log(1.d0/(nbins*nbins))
		if(prob.ne.0) D = D + prob*log(prob)/log(1.d0/(nbins*nbins))
		write(24,*) i,j,prob*log(prob)/log(1.d0/(nbins*nbins))

	enddo
	enddo
	print*, "D=", D

	close(23)
	close(24)


end program main
