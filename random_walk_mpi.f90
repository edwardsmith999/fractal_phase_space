! This csn be compiled in gfortran with integers promoted to 8byte
! which we can do with -fdefault-integer-8
! and for intel fortran is might be -integer-size 64
program main
    use mpi
    implicit none

    integer, parameter :: maxbinexp = 28, maxiterexp = 8, pexp=2
    integer :: i, j, k, it, fit, skip, nbins, sumint!, checksum
	!Specify integers of type 4 for MPI and to keep fin size down
	!the maximum number in a bin is then 2^31 = 2,147,483,648
    integer(kind=4) :: ierr, nproc, myid, rootid
    integer(kind=4), dimension(:), allocatable :: fin
    double precision :: r, y, D, p, nbins_s, t1, t2

    ! Initialize MPI
    call MPI_init(ierr)
    call MPI_comm_size (MPI_COMM_WORLD, nproc, ierr)
    call MPI_comm_rank (MPI_COMM_WORLD, myid, ierr)
    rootid = 0

    if (myid .eq. rootid) open (unit=23, file="infodim.out")

    nbins = pexp**maxbinexp
    allocate(fin(nbins)); 
    do k = maxiterexp,maxiterexp

        t1 = MPI_wtime()
    
        !Number of iterations
        fit = 10**k
        fin=0.d0

        !Loop to get transformed trajectory
        y = dble(myid)/dble(nproc)
        print*, "starting values = ", myid, nproc, y
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
            fin(j) = fin(j) + int(1, kind=4)
        enddo

        !Collect results from all processes
        if (myid .eq. rootid) then
            call MPI_Reduce(MPI_IN_PLACE, fin, nbins, MPI_INTEGER, &
                            MPI_SUM, rootid, MPI_COMM_WORLD, ierr)
        else
            call MPI_Reduce(fin, fin, nbins, MPI_INTEGER, &
                            MPI_SUM, rootid, MPI_COMM_WORLD, ierr)
        endif

        !Only average then on root
        if (myid .eq. rootid) then
            !print*, "maximum bin count = ", maxval(fin)

            !If we store all bins we can keep averaging in powers of two  
            !to get fractal dimensions from one run at about 64 bins (maxbinexp-6)
            do j=0,maxbinexp-2 !-1
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
                    p = dble(sumint)/dble(fit*nproc)
                    !checksum = checksum + sum(fin(i:i+skip-1))
                    if(p.ne.0) D = D + p*log(p)/log(1.d0/nbins_s)

                    !write(10000*pexp+(maxbinexp-j)*100+k,*) p

                enddo
                !print*, skip, nbins, checksum, fit
                print*, nbins/skip,fit*nproc, sum(fin),1.d0/log(1.d0/nbins_s),D
                write(23,*) k, j, D
            enddo

            t2 = MPI_wtime()
            print*, "Time take is ", t2-t1

        endif

    enddo
    deallocate(fin)
    if (myid .eq. rootid) close(23,status='keep')

    call MPI_finalize(ierr)

end program
