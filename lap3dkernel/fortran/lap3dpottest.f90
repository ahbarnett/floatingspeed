! test speed of Laplace 3D kernel eval in Fortran 90.
! Barnett 9/26/18 - 10/1/18
! ntest>1 added 3/4/20

! Issues: why is run time so variable, compared to other language versions??

! real*8:  
! gfortran-7 lap3dpottest.f90 -O3 -fopenmp -funroll-loops -march=native; ./a.out
! 1.1 Gpair/s,    @ 8 thr, or often at 4 thr too.
! or 1.3 if use _sep ordering (slightly better)
! This is same as numba jit in python3 intel version.


! -----------All this was single (real*4, which is real): ---------------------
! gfortran-8 lap3dpottest.f90 -O3 -fopenmp -march=native; ./a.out
! around 1.9 @ 8 thr

! gfortran-8 lap3dpottest.f90 -O3 -fopenmp; ./a.out
! only 0.7   @ 8 threads

! gfortran-8 lap3dpottest.f90 -Ofast -fopenmp -funroll-loops; ./a.out
! 1.2     @ 8 thr

! gfortran-8 lap3dpottest.f90 -O3 -funroll-loops; ./a.out
! gfortran-8 lap3dpottest.f90 -O3; ./a.out
! 0.5 Gpair/s

! gfortran lap3dpottest.f90 -O3 -funroll-loops; ./a.out
! 0.17

! gfortran-7 lap3dpottest.f90 -O3 -funroll-loops; ./a.out
! 0.25

! gfortran-7 lap3dpottest.f90 -O3; ./a.out
! 0.1


program lap3dpottest

  use omp_lib
  
  implicit none
  integer :: ns=10000, nt=10000, ntest=20, i
  real*8, allocatable :: x(:,:), y(:,:), q(:), pot(:)
  real :: t0,t1,t
  real*8 :: tot

  print *,'ns=',ns,'   nt=',nt
  
  allocate(x(3,nt))
  allocate(pot(nt))
  allocate(y(3,ns))
  allocate(q(ns))

  call random_number(x)
  call random_number(y)
  call random_number(q)
  
  !call cpu_time(t0)             ! use if no openmp
  !print *,'timer tick : ',omp_get_wtick()    ! 1 ns resolution
  print *,'ntest = ',ntest,' ...'
  t0 =omp_get_wtime()
  do i=1,ntest
     call lap3dpot(pot,y,q,x,ns,nt)
  end do
  !call cpu_time(t1)
  t = omp_get_wtime()-t0
  tot = 0.D0
  do i=1,nt
     tot = tot + pot(i)
  end do
  print *,'tot=',tot
  print *,ns*nt,'src-targ pairs in',t,'s:',ntest*ns*nt/t/1.E9,'Gpair/s'

  ! Amazingly, gfortran -Wall doesn't complain if we fail to flip the arrays
  ! but call a routine expecting flipped arrays! ... shameful ...
  
  if (0.eq.1) then
     t0 =omp_get_wtime()
     do i=1,ntest
        call lap3dpot_tr(pot,y,q,x,ns,nt)
     end do
     t = omp_get_wtime()-t0
     tot = 0.D0
     do i=1,nt
        tot = tot + pot(i)
     end do
     print *,'tot=',tot
     print *,'tr:',ns*nt,'src-targ pairs in',t,'s:',ntest*ns*nt/t/1.E9,'Gpair/s'
  endif

  deallocate(x)
  deallocate(y)
  deallocate(q)
  deallocate(pot)
  
end program lap3dpottest

! ----------------------------------------------------------------------
subroutine lap3dpot(pot, y,q,x,ns,nt)
  ! writes to pot the potential at targets x due to sources at y
  ! with charges q
  
  implicit none
  integer ns, nt, i, j
  real*8 :: x(3,nt), y(3,ns), q(ns), pot(nt)
  real*8 :: prefac, r2ij
  real*8 :: pi = 4*atan(1.D0)
  
  prefac = 1.D0/(4*pi)
  !$omp parallel
  !$omp do private(i,j,r2ij)
  do i=1,nt
     pot(i) = 0.0
     do j=1,ns
        r2ij = (x(1,i)-y(1,j))**2+(x(2,i)-y(2,j))**2+(x(3,i)-y(3,j))**2
        pot(i) = pot(i) + prefac * q(j) / sqrt(r2ij)
     end do
  end do
  !$omp end do
  !$omp end parallel
  
end subroutine lap3dpot

! ----------------------------------------------------------------------
subroutine lap3dpot_tr(pot, y,q,x,ns,nt)
  ! writes to pot the potential at targets x due to sources at y
  ! with charges q. Transposed arrays with x,y,z separated (n*3 not 3*n)
  
  implicit none
  integer ns, nt, i, j
  real*8 :: x(nt,3), y(ns,3), q(ns), pot(nt)
  real*8 :: prefac, r2ij
  real*8 :: pi = 4*atan(1.D0)
  
  prefac = 1.D0/(4*pi)
  !$omp parallel
  !$omp do private(i,j,r2ij)
  do i=1,nt
     pot(i) = 0.0
     do j=1,ns
        r2ij = (x(i,1)-y(j,1))**2+(x(i,2)-y(j,2))**2+(x(i,3)-y(j,3))**2
        pot(i) = pot(i) + prefac * q(j) / sqrt(r2ij)
     end do
  end do
  !$omp end do
  !$omp end parallel
  
end subroutine lap3dpot_tr

