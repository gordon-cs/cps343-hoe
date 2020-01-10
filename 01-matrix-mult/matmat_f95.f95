! **********************************************************************
! * Jonathan Senning <jonathan.senning@gordon.edu>                     *
! * Department of Mathematics and Computer Science                     *
! * Gordon College, 255 Grapevine Road, Wenham MA 01984-1899           *
! *                                                                    *
! * FORTRAN 95 FREE-FORMAT VERSION                                     *
! *                                                                    *
! * Benchmark ijk, jki, and ikj matrix-matrix products.                *
! **********************************************************************
!
! $Smake: gfortran -std=f95 -O3 -o %F %f
!

program MatMat
  implicit none
  integer n
  parameter (n=500)
  double precision a(n,n),b(n,n)
  double precision c1(n,n),c2(n,n),c3(n,n)
  integer icmppr
  integer i,j
  real t1,t2
  real timijk,timjki,timikj
  real gflpct

  ! initialize matrices

  do j=1,n
     do i=1,n
        call random_number(a(i,j))
        call random_number(b(i,j))
     end do
  end do

  ! ijk product

  call cpu_time(t1)
  call matijk(c1,a,b,n)
  call cpu_time(t2)
  timijk=t2-t1

  ! jki product

  call cpu_time(t1)
  call matjki(c2,a,b,n)
  call cpu_time(t2)
  timjki=t2-t1

  ! ikj product

  call cpu_time(t1)
  call matikj(c3,a,b,n)
  call cpu_time(t2)
  timikj=t2-t1

  ! output results

  gflpct = 2.0*n**3/1.0e9

  write (*,"(' Fortran 95     (',i3,') ijk: ',f6.3,' gflops, jki: ',&
       &f6.3,' gflops, ikj: ',f6.3,' gflops')") &
       n,gflpct/timijk, gflpct/timjki, gflpct/timikj
  
  ! compare products

  if (icmppr(c1,c2,n) > 0) then
     write(*,"('c',i1,' <> c',i1,': Validation error')") 1, 2
  endif
  if (icmppr(c1,c3,n) > 0) then
     write(*,"('c',i1,' <> c',i1,': Validation error')") 1, 3
  endif
  if (icmppr(c2,c3,n) > 0) then
     write(*,"('c',i1,' <> c',i1,': Validation error')") 2, 3
  endif

  ! all done

end program MatMat

!-----------------------------------------------------------------------------
! ijk matrix-matrix product
!-----------------------------------------------------------------------------

subroutine matijk(c,a,b,n)
  integer n
  double precision c(n,n),a(n,n),b(n,n)
  do i=1,n
     do j=1,n
        c(i,j)=0.0
        do k=1,n
           c(i,j)=c(i,j)+a(i,k)*b(k,j)
        end do
     end do
  end do
  return
end subroutine matijk

!-----------------------------------------------------------------------------
! jki matrix-matrix product
!-----------------------------------------------------------------------------

subroutine matjki(c,a,b,n)
  integer n
  double precision c(n,n),a(n,n),b(n,n)
  do j=1,n
     do i=1,n
        c(i,j)=0.0
     end do
     do k=1,n
        do i=1,n   
           c(i,j)=c(i,j)+a(i,k)*b(k,j)
        end do
     end do
  end do
  return
end subroutine matjki

!-----------------------------------------------------------------------------
! ikj matrix-matrix product
!-----------------------------------------------------------------------------

subroutine matikj(c,a,b,n)
  integer n
  double precision c(n,n),a(n,n),b(n,n)
  do i=1,n
     do j=1,n
        c(i,j)=0.0
     end do
     do k=1,n
        do j=1,n   
           c(i,j)=c(i,j)+a(i,k)*b(k,j)
        end do
     end do
  end do
  return
end subroutine matikj

!-----------------------------------------------------------------------------
! compare products
!-----------------------------------------------------------------------------

integer function icmppr(c,d,n)
  integer n
  double precision c(n,n),d(n,n)
  double precision eps
  integer icount
  icount=0
  eps=1.0e-12
  do i=1,n
     do j=1,n
        if (abs(c(i,j)-d(i,j))>eps) then
           icount=icount+1
           write (*,*) '****warning****'
           write (*,"('(',i3,',',i3,'): ',f15.5,' <> ',f15.5)") &
                i,j,c(i,j),d(i,j)
        endif
     end do
  end do
  icmppr=icount
end function icmppr
