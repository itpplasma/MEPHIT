program merge_qq_tri
  implicit none
  integer, parameter :: legs=3
  integer :: i, j, k, m, ntri, ntri_qq, idummy, npoints
  real*8, dimension(:), allocatable :: R_coord, Z_coord
  integer, dimension(:,:), allocatable :: list_tri, list_tri_qq

  open(1, file='17151.aug_qq.ele')
  read(1,*) ntri_qq
  allocate(list_tri_qq(legs,ntri_qq))
  do m=1, ntri_qq
     read(1,*) idummy, (list_tri_qq(j,m), j=1,legs)
  enddo
  close(1)

  open(2, file='17151.aug.1.ele')
  read(2,*) ntri
  ntri = ntri + ntri_qq
  allocate(list_tri(legs,ntri))
  list_tri(:,1:ntri_qq) = list_tri_qq(:,:)
  do m=ntri_qq+1, ntri
     read(2,*) idummy, (list_tri(j,m), j=1,legs)
  enddo
  close(2)

  open(3, file='17151.aug.1.node')
  read(3,*) npoints
  allocate(R_coord(npoints), Z_coord(npoints))
  do m=1, npoints
     read(3,*) idummy, R_coord(m), Z_coord(m)
  enddo
  close(3)

  do m=1, ntri
     do k=1,legs+1
        j = modulo(k,legs) + 1
        write(210,*) R_coord(list_tri(j,m)), Z_coord(list_tri(j,m))
     enddo
     write(210,*) 
  enddo
  close(210)

end program merge_qq_tri
