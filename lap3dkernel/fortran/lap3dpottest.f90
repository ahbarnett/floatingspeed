! test speed of Laplace 3D kernel eval in Fortran 90.
! Barnett 9/26/18 - 10/1/18
! ntest>1 added 3/4/20
! ChatGPT5.2 via codex CLI, added variants using SIMD, 1/29/26.

! build:
! gcc -O3 -march=native -c rsqrtps.c -o rsqrtps.o
! gfortran lap3dpottest.f90 rsqrtps.o -O3 -fopenmp -funroll-loops -march=native -o lap3d_compare_O3_rsqrt
! ./lap3d_compare_O3_rsqrt

! Ryzen 5700 U laptop (ross), all 16 threads (full HT) 1/29/26, reported
! in codex:
!  - orig: 2.3120 Gpair/s
!  - accum: 5.9146 Gpair/s
!  - accum+simd: 6.7645 Gpair/s
!  - soa: 7.5631 Gpair/s             (SoA = struct of arrays)
!  - soa+simd: 8.1021 Gpair/s
!  - soa+rsqrt: 3.3898 Gpair/s          (due to func call overhead?)
!  - soa+rsqrb (block): 7.2522 Gpair/s   (still can't catch up with soa+simd)


program lap3dpottest

  !use omp_lib
  use iso_c_binding, only: c_float, c_double

  implicit none
  integer :: ns=10000, nt=10000, ntest=20, i
  real*8, allocatable :: x(:,:), y(:,:), q(:), pot(:)
  real*8 :: t0,t1
  real :: t
  real*8 :: tot, tot_ref, omp_get_wtime, gpair
  
  interface
     subroutine rsqrtps_nr4(input, output) bind(C, name="rsqrtps_nr4")
       use iso_c_binding, only: c_float, c_double
       real(c_float), dimension(4), intent(in) :: input
       real(c_double), dimension(4), intent(out) :: output
     end subroutine rsqrtps_nr4
     subroutine rsqrtps_block1024(input, output) bind(C, name="rsqrtps_block1024")
       use iso_c_binding, only: c_float, c_double
       real(c_float), dimension(1024), intent(in) :: input
       real(c_double), dimension(1024), intent(out) :: output
     end subroutine rsqrtps_block1024
  end interface

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
  t1 = omp_get_wtime()
  t = (t1-t0)/ntest
  tot = 0.D0
  do i=1,nt
     tot = tot + pot(i)
  end do
  print *,'orig tot=',tot
  gpair = ns*nt/t/1.E9
  write(*,'(A12,1X,I12,1X,A,1X,F7.4,1X,A,1X,F7.4,1X,A)') &
       'orig:', ns*nt, 'pairs in', t, 's:', gpair, 'Gpair/s'

  t0 =omp_get_wtime()
  do i=1,ntest
     call lap3dpot_accum(pot,y,q,x,ns,nt)
  end do
  t1 = omp_get_wtime()
  t = (t1-t0)/ntest
  tot = 0.D0
  do i=1,nt
     tot = tot + pot(i)
  end do
  print *,'accum tot=',tot
  gpair = ns*nt/t/1.E9
  write(*,'(A12,1X,I12,1X,A,1X,F7.4,1X,A,1X,F7.4,1X,A)') &
       'accum:', ns*nt, 'pairs in', t, 's:', gpair, 'Gpair/s'

  t0 =omp_get_wtime()
  do i=1,ntest
     call lap3dpot_accum_simd(pot,y,q,x,ns,nt)
  end do
  t1 = omp_get_wtime()
  t = (t1-t0)/ntest
  tot = 0.D0
  do i=1,nt
     tot = tot + pot(i)
  end do
  print *,'accum+simd tot=',tot
  gpair = ns*nt/t/1.E9
  write(*,'(A12,1X,I12,1X,A,1X,F7.4,1X,A,1X,F7.4,1X,A)') &
       'accum+simd:', ns*nt, 'pairs in', t, 's:', gpair, 'Gpair/s'

  t0 =omp_get_wtime()
  do i=1,ntest
     call lap3dpot_soa(pot,y,q,x,ns,nt)
  end do
  t1 = omp_get_wtime()
  t = (t1-t0)/ntest
  tot = 0.D0
  do i=1,nt
     tot = tot + pot(i)
  end do
  print *,'soa tot=',tot
  gpair = ns*nt/t/1.E9
  write(*,'(A12,1X,I12,1X,A,1X,F7.4,1X,A,1X,F7.4,1X,A)') &
       'soa:', ns*nt, 'pairs in', t, 's:', gpair, 'Gpair/s'

  t0 =omp_get_wtime()
  do i=1,ntest
     call lap3dpot_soa_simd(pot,y,q,x,ns,nt)
  end do
  t1 = omp_get_wtime()
  t = (t1-t0)/ntest
  tot = 0.D0
  do i=1,nt
     tot = tot + pot(i)
  end do
  print *,'soa+simd tot=',tot
  gpair = ns*nt/t/1.E9
  write(*,'(A12,1X,I12,1X,A,1X,F7.4,1X,A,1X,F7.4,1X,A)') &
       'soa+simd:', ns*nt, 'pairs in', t, 's:', gpair, 'Gpair/s'
  tot_ref = tot

  t0 =omp_get_wtime()
  do i=1,ntest
     call lap3dpot_soa_simd_rsqrtps(pot,y,q,x,ns,nt)
  end do
  t1 = omp_get_wtime()
  t = (t1-t0)/ntest
  tot = 0.D0
  do i=1,nt
     tot = tot + pot(i)
  end do
  print *,'soa+simd rsqrtps tot=',tot
  print *,'soa+rsqrt tot diff=',tot-tot_ref,'rel=',(tot-tot_ref)/tot_ref
  gpair = ns*nt/t/1.E9
  write(*,'(A12,1X,I12,1X,A,1X,F7.4,1X,A,1X,F7.4,1X,A)') &
       'soa+rsqrt:', ns*nt, 'pairs in', t, 's:', gpair, 'Gpair/s'

  t0 =omp_get_wtime()
  do i=1,ntest
     call lap3dpot_soa_simd_rsqrtps_block(pot,y,q,x,ns,nt)
  end do
  t1 = omp_get_wtime()
  t = (t1-t0)/ntest
  tot = 0.D0
  do i=1,nt
     tot = tot + pot(i)
  end do
  print *,'soa+simd rsqrtps blk tot=',tot
  print *,'soa+rsqrt blk diff=',tot-tot_ref,'rel=',(tot-tot_ref)/tot_ref
  gpair = ns*nt/t/1.E9
  write(*,'(A12,1X,I12,1X,A,1X,F7.4,1X,A,1X,F7.4,1X,A)') &
       'soa+rsqrb:', ns*nt, 'pairs in', t, 's:', gpair, 'Gpair/s'

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
     do j=1,ns
        r2ij = (x(1,i)-y(1,j))**2+(x(2,i)-y(2,j))**2+(x(3,i)-y(3,j))**2
        accum = accum + prefac * q(j) / sqrt(r2ij)
     end do
     pot(i) = accum
  end do
  !$omp end do
  !$omp end parallel
  
end subroutine lap3dpot_accum

! ----------------------------------------------------------------------
subroutine lap3dpot_accum_simd(pot, y,q,x,ns,nt)
  ! writes to pot the potential at targets x due to sources at y
  ! with charges q (SIMD reduction on inner loop)

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
  
end subroutine lap3dpot_accum_simd

! ----------------------------------------------------------------------
subroutine lap3dpot_soa(pot, y,q,x, ns,nt)
  ! SoA (struct-of-arrays) layout for better contiguous access in the inner loop

  implicit none
  integer ns, nt, i, j
  real*8 :: x(3,nt), y(3,ns), q(ns), pot(nt)
  real*8 :: xs(ns), ys(ns), zs(ns)
  real*8 :: xt(nt), yt(nt), zt(nt)
  real*8 :: prefac, r2ij
  real*8 :: pi = 4*atan(1.D0)
  real*8 :: accum, dx, dy, dz

  prefac = 1.D0/(4*pi)
  xs = y(1,:)
  ys = y(2,:)
  zs = y(3,:)
  xt = x(1,:)
  yt = x(2,:)
  zt = x(3,:)
  !$omp parallel
  !$omp do private(i,j,r2ij,accum,dx,dy,dz)
  do i=1,nt
     accum = 0.0D0
     do j=1,ns
        dx = xt(i) - xs(j)
        dy = yt(i) - ys(j)
        dz = zt(i) - zs(j)
        r2ij = dx*dx + dy*dy + dz*dz
        accum = accum + prefac * q(j) / sqrt(r2ij)
     end do
     pot(i) = accum
  end do
  !$omp end do
  !$omp end parallel

end subroutine lap3dpot_soa

! ----------------------------------------------------------------------
subroutine lap3dpot_soa_simd(pot, y,q,x, ns,nt)
  ! SoA (struct-of-arrays) layout + SIMD reduction on inner loop

  implicit none
  integer ns, nt, i, j
  real*8 :: x(3,nt), y(3,ns), q(ns), pot(nt)
  real*8 :: xs(ns), ys(ns), zs(ns)
  real*8 :: xt(nt), yt(nt), zt(nt)
  real*8 :: prefac, r2ij
  real*8 :: pi = 4*atan(1.D0)
  real*8 :: accum, dx, dy, dz

  prefac = 1.D0/(4*pi)
  xs = y(1,:)
  ys = y(2,:)
  zs = y(3,:)
  xt = x(1,:)
  yt = x(2,:)
  zt = x(3,:)
  !$omp parallel
  !$omp do private(i,j,r2ij,accum,dx,dy,dz)
  do i=1,nt
     accum = 0.0D0
     !$omp simd reduction(+:accum)
     do j=1,ns
        dx = xt(i) - xs(j)
        dy = yt(i) - ys(j)
        dz = zt(i) - zs(j)
        r2ij = dx*dx + dy*dy + dz*dz
        accum = accum + prefac * q(j) / sqrt(r2ij)
     end do
     pot(i) = accum
  end do
  !$omp end do
  !$omp end parallel

end subroutine lap3dpot_soa_simd

! ----------------------------------------------------------------------
subroutine lap3dpot_soa_simd_rsqrtps(pot, y,q,x, ns,nt)
  ! SoA (struct-of-arrays) layout + rsqrtps (float) with Newton refinements

  use iso_c_binding, only: c_float, c_double
  implicit none
  integer ns, nt, i, j
  real*8 :: x(3,nt), y(3,ns), q(ns), pot(nt)
  real*8 :: xs(ns), ys(ns), zs(ns)
  real*8 :: xt(nt), yt(nt), zt(nt)
  real*8 :: prefac, r2ij
  real*8 :: pi = 4*atan(1.D0)
  real*8 :: accum
  real(c_float) :: r2f(4)
  real(c_double) :: rinvf(4)
  integer :: jmax

  interface
     subroutine rsqrtps_nr4(input, output) bind(C, name="rsqrtps_nr4")
       use iso_c_binding, only: c_float, c_double
       real(c_float), dimension(4), intent(in) :: input
       real(c_double), dimension(4), intent(out) :: output
     end subroutine rsqrtps_nr4
  end interface

  prefac = 1.D0/(4*pi)
  xs = y(1,:)
  ys = y(2,:)
  zs = y(3,:)
  xt = x(1,:)
  yt = x(2,:)
  zt = x(3,:)
  !$omp parallel
  !$omp do private(i,j,accum,r2ij,r2f,rinvf,jmax)
  do i=1,nt
     accum = 0.0D0
     jmax = ns - mod(ns,4)
     do j=1,jmax,4
        r2f(1) = real((xt(i)-xs(j  ))**2 + (yt(i)-ys(j  ))**2 + (zt(i)-zs(j  ))**2, c_float)
        r2f(2) = real((xt(i)-xs(j+1))**2 + (yt(i)-ys(j+1))**2 + (zt(i)-zs(j+1))**2, c_float)
        r2f(3) = real((xt(i)-xs(j+2))**2 + (yt(i)-ys(j+2))**2 + (zt(i)-zs(j+2))**2, c_float)
        r2f(4) = real((xt(i)-xs(j+3))**2 + (yt(i)-ys(j+3))**2 + (zt(i)-zs(j+3))**2, c_float)
        call rsqrtps_nr4(r2f, rinvf)
        accum = accum + prefac * q(j  ) * rinvf(1)
        accum = accum + prefac * q(j+1) * rinvf(2)
        accum = accum + prefac * q(j+2) * rinvf(3)
        accum = accum + prefac * q(j+3) * rinvf(4)
     end do
     do j=jmax+1,ns
        r2ij = (xt(i)-xs(j))**2 + (yt(i)-ys(j))**2 + (zt(i)-zs(j))**2
        accum = accum + prefac * q(j) / sqrt(r2ij)
     end do
     pot(i) = accum
  end do
  !$omp end do
  !$omp end parallel

end subroutine lap3dpot_soa_simd_rsqrtps

! ----------------------------------------------------------------------
subroutine lap3dpot_soa_simd_rsqrtps_block(pot, y,q,x, ns,nt)
  ! SoA (struct-of-arrays) layout + rsqrtps in C over 1024-element blocks

  use iso_c_binding, only: c_float, c_double
  implicit none
  integer, parameter :: block = 1024
  integer ns, nt, i, j, k, jmax
  real*8 :: x(3,nt), y(3,ns), q(ns), pot(nt)
  real*8 :: xs(ns), ys(ns), zs(ns)
  real*8 :: xt(nt), yt(nt), zt(nt)
  real*8 :: prefac, r2ij
  real*8 :: pi = 4*atan(1.D0)
  real*8 :: accum
  real(c_float) :: r2f(block)
  real(c_double) :: rinvf(block)

  interface
     subroutine rsqrtps_block1024(input, output) bind(C, name="rsqrtps_block1024")
       use iso_c_binding, only: c_float, c_double
       real(c_float), dimension(1024), intent(in) :: input
       real(c_double), dimension(1024), intent(out) :: output
     end subroutine rsqrtps_block1024
  end interface

  prefac = 1.D0/(4*pi)
  xs = y(1,:)
  ys = y(2,:)
  zs = y(3,:)
  xt = x(1,:)
  yt = x(2,:)
  zt = x(3,:)
  !$omp parallel
  !$omp do private(i,j,k,jmax,accum,r2ij,r2f,rinvf)
  do i=1,nt
     accum = 0.0D0
     jmax = ns - mod(ns,block)
     do j=1,jmax,block
        !$omp simd
        do k=1,block
           r2f(k) = real((xt(i)-xs(j+k-1))**2 + (yt(i)-ys(j+k-1))**2 + &
                         (zt(i)-zs(j+k-1))**2, c_float)
        end do
        call rsqrtps_block1024(r2f, rinvf)
        !$omp simd reduction(+:accum)
        do k=1,block
           accum = accum + prefac * q(j+k-1) * rinvf(k)
        end do
     end do
     do j=jmax+1,ns
        r2ij = (xt(i)-xs(j))**2 + (yt(i)-ys(j))**2 + (zt(i)-zs(j))**2
        accum = accum + prefac * q(j) / sqrt(r2ij)
     end do
     pot(i) = accum
  end do
  !$omp end do
  !$omp end parallel

end subroutine lap3dpot_soa_simd_rsqrtps_block
