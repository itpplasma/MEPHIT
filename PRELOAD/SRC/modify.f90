program modify
  use from_nrtype
  implicit none
  integer(i4b), parameter :: nmax=100, num_el=27
  integer(i4b) :: i, j, k, n(num_el), num_wall, j_next, i_next 
  integer(i4b) :: npoi_add_1, nele_add_1, nbetween_1, npoi_add_2, nele_add_2, nbetween_2 
  real(dp), dimension(nmax,num_el) :: r, z
  real(dp) :: r_2, z_2, hr, hz
  character*2 dummy
  type :: wall_element
!     sequence
     real(dp) :: r_start 
     real(dp) :: r_fin
     real(dp) :: z_start
     real(dp) :: z_fin
     integer(i4b) :: material
  end type wall_element
  type(wall_element), dimension(:), allocatable :: wall

! 15 artificial points should be added to element # 15 between points 15 & 16
  npoi_add_1 = 15
  nele_add_1 = 15
  nbetween_1 = 15
! and 31 artificial points should be added to element # 20 between points 1 & 2
  npoi_add_2 = 31
  nele_add_2 = 20
  nbetween_2 = 1

  num_wall = 0
  open(1,file='FROMQQ/Dave/17151.wall')
  do i=1, num_el
     read(1,100)dummy, n(i)
     print *,i,n(i)

     if(i .eq. nele_add_1) then
        do j=1, nbetween_1
           read(1,*)r(j,i), z(j,i)
        enddo
        read(1,*) r_2, z_2
        hr = (r_2 - r(nbetween_1,i))/dfloat(npoi_add_1 + 1)
        hz = (z_2 - z(nbetween_1,i))/dfloat(npoi_add_1 + 1)
        do k=1, npoi_add_1
           r(nbetween_1+k,i) = r(nbetween_1,i) + hr*k 
           z(nbetween_1+k,i) = z(nbetween_1,i) + hz*k
        enddo
        r(nbetween_1+npoi_add_1+1,i) = r_2
        z(nbetween_1+npoi_add_1+1,i) = z_2

        do j=nbetween_1+2,n(i)
           read(1,*)r(nbetween_1+npoi_add_1+j-nbetween_1,i), z(nbetween_1+npoi_add_1+j-nbetween_1,i)
        enddo
        n(i) = n(i) + npoi_add_1
     elseif(i .eq. nele_add_2) then
        do j=1, nbetween_2
           read(1,*)r(j,i), z(j,i)
        enddo
        read(1,*) r_2, z_2
        hr = (r_2 - r(nbetween_2,i))/dfloat(npoi_add_2 + 1)
        hz = (z_2 - z(nbetween_2,i))/dfloat(npoi_add_2 + 1)
        do k=1, npoi_add_2
           r(nbetween_2+k,i) = r(nbetween_2,i) + hr*k 
           z(nbetween_2+k,i) = z(nbetween_2,i) + hz*k
        enddo
        r(nbetween_2+npoi_add_2+1,i) = r_2
        z(nbetween_2+npoi_add_2+1,i) = z_2

        do j=nbetween_2+2,n(i)
           read(1,*)r(nbetween_2+npoi_add_2+j-nbetween_2,i), z(nbetween_2+npoi_add_2+j-nbetween_2,i)
        enddo
        n(i) = n(i) + npoi_add_2
     else
        do j=1,n(i)
           read(1,*)r(j,i), z(j,i)
        enddo
     endif
     num_wall = num_wall + n(i)

  enddo
  close(1)

  allocate(wall(1:num_wall))
  k=0
  do i=1, num_el
     do j=1,n(i)
        k = k+1
        if(j .lt.n(i)) then
           j_next = j+1
           i_next = i
           wall(k)%material = 74 ! Tungsten
        else
           wall(k)%material = 0  ! Absent boundary
           if(i .lt. num_el) then
              i_next = i+1
              j_next = 1
           else
              i_next = 1
              j_next = 1
           endif
        endif
        wall(k)%r_start = r(j,i)
        wall(k)%r_fin   = r(j_next,i_next)
        wall(k)%z_start = z(j,i)
        wall(k)%z_fin   = z(j_next,i_next)
     enddo
  enddo
  open(2,file='FROMQQ/Dave/17151.wall_full')
  do k=1,num_wall
     write(2,101) wall(k)
  enddo
  close(2)
100 format(a2,i4)
101 format(4(ES15.8E2,1X),I4)

end program modify
