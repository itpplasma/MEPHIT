module inside_of_wall
  implicit none
! this is a code taken from the web-site
! http://alienryderflex.com/polygon/
! which is in c but has been translated to f90
!! and slightly modified
  integer, private :: npc                        !! number of corners of the polygon (no repeats)
  real(kind=8), dimension(:), allocatable, private :: xpoly, ypoly  !! polygon data
  real(kind=8), dimension(:), allocatable, private  :: c, m          !! precalculated data 
  logical, private :: call1st=.true.
contains

  function PointIsInside(x, y) result(inside)
    implicit none
    logical :: inside                                 !! answer as logical
    real(kind=8), intent(in) :: x, y                  !! point to be tested
! local
    integer :: i, j

    j=npc
    inside=.false.

    do i=1,npc
       if(( ypoly(i) < y .and. ypoly(j) >= y ) .or. &
            ( ypoly(j) < y .and. ypoly(i) >= y ) ) inside = (inside .neqv. y*m(i)+c(i) < x)
       j=i
    enddo

  end function PointIsInside

  subroutine precalc_wall
    implicit none
! local
    integer :: i, j
    integer, parameter :: npc_max=10000
    real(kind=8), allocatable, dimension(:) :: xpoly_max, ypoly_max
    real(kind=8) :: dummy_r8

    if(call1st) then
       call1st = .false.
       allocate(xpoly_max(npc_max), ypoly_max(npc_max))
       npc = 0
!       open(1,file='FROMQQ/Dave/wall_true.dat')
       open(1,file='17151.wall_full')
       do i=1, npc_max
!          read(1,*,END=10) xpoly_max(i), ypoly_max(i)
          read(1,*,END=10) xpoly_max(i), dummy_r8, ypoly_max(i)
          npc = npc + 1
       enddo
10     continue
       close(1)
!       npc = npc-1 ! the last point should not coincide with the 1st
       allocate(xpoly(npc), ypoly(npc), c(npc), m(npc))
       xpoly(1:npc) = xpoly_max(1:npc)
       ypoly(1:npc) = ypoly_max(1:npc)
       deallocate(xpoly_max, ypoly_max)
    endif

    j=npc
    do i=1,npc
       if(ypoly(j) == ypoly(i))then
          c(i)=xpoly(i)
          m(i)=0
       else
          c(i)=xpoly(i)-(ypoly(i)*xpoly(j))/(ypoly(j)-ypoly(i))+(ypoly(i)*xpoly(i))/(ypoly(j)-ypoly(i))
          m(i)=(xpoly(j)-xpoly(i))/(ypoly(j)-ypoly(i))
       endif
       j=i
    enddo
  end subroutine precalc_wall
end module inside_of_wall
