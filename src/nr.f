* This file contains routines adapted from
* Numerical Recipes in Fortran 90, Vol.2
* by Press, Teukoslky et al.

************************************************************************
      subroutine select(k,arr_in,n,outval)
************************************************************************
      !USE nrtype; USE nrutil, ONLY : assert,swap
      IMPLICIT NONE
      !INTEGER, INTENT(IN) :: k
      !REAL, DIMENSION(:), INTENT(INOUT) :: arr
      !REAL :: select
      INTEGER, INTENT(IN) :: k
      integer, intent(in) :: n
      REAL, DIMENSION(n), INTENT(IN) :: arr_in
      real, intent(out) :: outval
*     Returns the kth smallest value in the array arr. The input array
*     will be rearranged to have this value in location arr(k), with
*     all smaller elements moved to arr(1:k-1) (in arbitrary order)
*     and all larger elements in arr(k+1:) (also in arbitrary order).
      INTEGER :: i,r,j,l
      REAL :: a
      REAL, DIMENSION(size(arr_in)) :: arr
      arr=arr_in
      !call assert(k >= 1, k <= n, ’select args’)
      l=1
      r=n
      do
       if (r-l <= 1) then
        if (r-l == 1) call masked_swap_rs(arr(l),arr(r),arr(l)>arr(r))
        outval=arr(k)
        RETURN
       else
        i=(l+r)/2
        call swap_r(arr(i),arr(l+1))
        call masked_swap_rs(arr(l),arr(r),arr(l)>arr(r))
        call masked_swap_rs(arr(l+1),arr(r),arr(l+1)>arr(r))
        call masked_swap_rs(arr(l),arr(l+1),arr(l)>arr(l+1))
        i=l+1
        j=r
        a=arr(l+1)
        do
         do
          i=i+1
          if (arr(i) >= a) exit
         end do
         do
          j=j-1
          if (arr(j) <= a) exit
         end do
         if (j < i) exit
         call swap_r(arr(i),arr(j))
        end do
        arr(l+1)=arr(j)
        arr(j)=a
        if (j >= k) r=j-1
        if (j <= k) l=i
       end if
      end do
      END subroutine

*** swap
      SUBROUTINE swap_r(a,b)
      REAL, INTENT(INOUT) :: a,b
      REAL :: dum
      dum=a
      a=b
      b=dum
      END SUBROUTINE swap_r

      SUBROUTINE masked_swap_rs(a,b,mask)
      REAL, INTENT(INOUT) :: a,b
      LOGICAL, INTENT(IN) :: mask
      REAL :: swp
      if (mask) then
       swp=a
       a=b
       b=swp
      end if
      END SUBROUTINE masked_swap_rs
