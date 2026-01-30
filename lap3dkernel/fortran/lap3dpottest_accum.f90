! test speed of Laplace 3D kernel eval in Fortran 90.
! Variant: use local accumulator to aid SIMD/optimizer.
!
! Based on lap3dpottest.f90

program lap3dpottest_accum

  implicit none
  integer :: ns=10000, nt=10000, ntest=20, i
  real*8, allocatable :: x(:,:), y(:,:), q(:), pot(:)
  real*8 :: t0,t1
  real :: t
  real*8 :: tot,omp_get_wtime

  print *,'ns=',ns,'   nt=',nt

  allocate(x(3,nt))
  allocate(pot(nt))
  allocate(y(3,ns))
  allocate(q(ns))

  call random_number(x)
  call random_number(y)
  call random_number(q)

  print *,'ntest = ',ntest,' ...'
  t0 = omp_get_wtime()
  do i=1,ntest
     call lap3dpot_accum(pot,y,q,x,ns,nt)
  end do
  t1 = omp_get_wtime()
  t = t1-t0
  tot = 0.D0
  do i=1,nt
     tot = tot + pot(i)
  end do
  print *,'tot=',tot
  print *,ns*nt,'src-targ pairs in',t,'s:',ntest*ns*nt/t/1.E9,'Gpair/s'

  deallocate(x)
  deallocate(y)
  deallocate(q)
  deallocate(pot)

end program lap3dpottest_accum

! ----------------------------------------------------------------------
subroutine lap3dpot_accum(pot, y,q,x,ns,nt)
  ! writes to pot the potential at targets x due to sources at y
  ! with charges q

  implicit none
  integer ns, nt, i, j
  real*8 :: x(3,nt), y(3,ns), q(ns), pot(nt)
  real*8 :: prefac, r2ij
  real*8 :: pi = 4*atan(1.D0)
  real*8 :: accum

  prefac = 1.D0/(4*pi)
  !$omp parallel
  !$omp do private(i,j,r2ij,accum)
  do i=1,nt
     accum = 0.0D0
     !$omp simd reduction(+:accum)
     do j=1,ns
        r2ij = (x(1,i)-y(1,j))**2+(x(2,i)-y(2,j))**2+(x(3,i)-y(3,j))**2
        accum = accum + prefac * q(j) / sqrt(r2ij)
     end do
     pot(i) = accum
  end do
  !$omp end do
  !$omp end parallel

end subroutine lap3dpot_accum
