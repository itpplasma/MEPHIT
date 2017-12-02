module points_in_polygon
implicit none
! this is a code taken from the web-site
! http://alienryderflex.com/polygon/
! which is in c but has been translated to f90
!! and slightly modified 
contains

function PointIsInside(x, y, npc, xpoly, ypoly, c, m) result(inside)
implicit none
logical :: inside                                 !! answer as logical
real, intent(in) :: x, y                          !! point to be tested
integer, intent(in) :: npc                        !! number of corners of the polygon (no repeats)
real, dimension(npc), intent(in) :: xpoly, ypoly  !! polygon data
real, dimension(npc), intent(in) :: c, m          !! precalculated data
! local
integer :: i, j

j=npc
inside=.false.

do i=1,npc
  if(( ypoly(i) < y .and. ypoly(j) >= y ) .or. &
     ( ypoly(j) < y .and. ypoly(i) >= y ) ) inside = (inside .xor. y*m(i)+c(i) < x)
  j=i
enddo

end function PointIsInside

subroutine precalc_values(npc, xpoly, ypoly, c, m)
implicit none
integer, intent(in) :: npc                       !! number of corners of the polygon (no repeats)
real, dimension(npc), intent(in) :: xpoly, ypoly !! polygon points
real, dimension(npc), intent(out) :: c, m        !! pre-calculated arrays for speed-up
! local
integer :: i, j

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

end subroutine precalc_values

end module points_in_polygon
