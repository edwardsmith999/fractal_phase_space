! This csn be compiled in gfortran with integers promoted to 8byte
! which we can do with -fdefault-integer-8
! and for intel fortran is might be -integer-size 64
program main
    use mpi
    implicit none

    integer, parameter :: maxbinexp = 18, maxiterexp = 10, pexp=3
    integer :: i, j, k, it, fit, skip, nbins, sumint, sumtotal
    INTEGER, DIMENSION (1) :: seed
	!Specify integers of type 4 for MPI and to keep fin size down
	!the maximum number in a bin is then 2^31 = 2,147,483,648
    integer(kind=4) :: ierr, nproc, myid, rootid, chunk, split
    integer(kind=4), dimension(:), allocatable :: fin
    !integer(kind=4), dimension(:), allocatable :: fin
    double precision :: r, y, D, nbins_s, t1, t2, t3, t4, t5
    double precision :: q, p, qnew, pnew

    ! Initialize MPI
    call MPI_init(ierr)
    call MPI_comm_size (MPI_COMM_WORLD, nproc, ierr)
    call MPI_comm_rank (MPI_COMM_WORLD, myid, ierr)
    rootid = 0

    t1 = MPI_wtime()

    nbins = pexp**maxbinexp
    allocate(fin(nbins));
    print*, "Array allocated with nbins=", nbins, " with ", &
        dble(4*nbins)/1024**3, "Gb per proc and total ", dble(4*nbins*nproc)/1024**3, "Gb"

    t2 = MPI_wtime()
    if (myid .eq. rootid) print*, "Array Allocation time=", t2-t1

    !Number of iterations
    fit = 10**maxiterexp
    fin=0

    !Loop to get transformed trajectory
    seed(1) = myid
    call random_seed(put=seed)
    y = dble(myid)/dble(nproc)
    call random_number(r)
    call random_number(r)
    print*, "starting values proc ", myid+1, " of ", nproc, &
            " y= ", y, "rand=", r

    do it = 1,fit
        !Randomly transform
        call random_number(r)
        if(r.lt.1.d0/3.d0) then
            y = (1.d0+y+y)/3.d0
        else if (r.gt.1.d0/3.d0) then
            y = y/3.d0
        endif

        !Binning operation
        j = int(floor(y*nbins)) + 1
        fin(j) = fin(j) + int(1, kind=4)
    enddo

    t3 = MPI_wtime()
    if (myid .eq. rootid) print*, "Run mapping time=", t3-t2

    print*, "value on proc=", myid, " for fin(1)=", &
            fin(1), "fin(10)=", fin(10)

    !Collect results from all processes, need to do in batches
	!to prevent send array from being greater than 32 bit integer
    split = ceiling(nbins/2**31)
    chunk = nbins/split
    if (myid .eq. rootid) then
        do i=1,split
            call MPI_Reduce(MPI_IN_PLACE, fin((i-1)*chunk+1:i*chunk),chunk, MPI_INTEGER, &
                            MPI_SUM, rootid, MPI_COMM_WORLD, ierr)
        enddo
        print*, "value after reduce for fin(1)=", &
                      fin(1), "fin(10)=", fin(10)
    else
        do i=1,split
            call MPI_Reduce(fin((i-1)*chunk+1:i*chunk), &
                            fin((i-1)*chunk+1:i*chunk), chunk,MPI_INTEGER, &
                            MPI_SUM, rootid, MPI_COMM_WORLD, ierr)
        enddo
    endif

    t4 = MPI_wtime()
    if (myid .eq. rootid) print*, "MPI_reduce time=", t4-t3

    !Only average then on root
    if (myid .eq. rootid) then


        !Get number of iterations
        sumtotal = 0
        do j=1,size(fin,1)
            sumtotal = sumtotal + fin(j)
        enddo

        !If we store all bins we can keep averaging in powers of two
        !to get fractal dimensions from one run at about 64 bins
        !(maxbinexp-6)
        !or in powers of 3 down to 81 bins (maxbinexp-4)
        !which means we use 12/pexp = 6 or 4 for 2 or 3 respectively
        do j=0,maxbinexp-12/pexp !-1
            skip = pexp**j
            D = 0.d0
            nbins_s = dble(nbins)/dble(skip)
            !Loop over blocks of skip averaging
            do i = 1,nbins,skip
                sumint = 0
                do k=i,i+skip-1
                    sumint = sumint + fin(k)
                enddo
                p = dble(sumint)/dble(sumtotal)
                if(p.ne.0) D = D + p*log(p)/log(1.d0/nbins_s)

            enddo
            print'(4i12, 2f20.15)', maxbinexp-j, nbins/skip,fit*nproc,sumtotal,1.d0/log(1.d0/nbins_s),D

        enddo

        t5 = MPI_wtime()
        print*, "Bins averaging time ", t5-t4
        print*, "Total Time take is ", t5-t1

    endif
    deallocate(fin)

    call MPI_finalize(ierr)

end program