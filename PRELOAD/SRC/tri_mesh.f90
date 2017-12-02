program tri_mesh
  implicit none
  integer, parameter :: legs=3, nwall=1000
  integer :: i, j, k, m, ntri, nfree, nbou_qq, ji, j1i, ji1, j1i1, nwall_true, iyqq_prev
  integer :: number_of_points, nr_core, nt_core, nr_pfr, nt_pfr, nr_sol, nt_sol, n_therm
  real*8, dimension(:), allocatable :: R_coord, Z_coord, psif, bmod, B_p
  real*8, dimension(nwall) :: R_wall, Z_wall
  real*8, dimension(:,:), allocatable :: xleg, yleg
  real*8, parameter ::  eps_dist=1.d-6
  real*8 :: dist, dummy_r8, btf, rtf
  integer, dimension(:,:), allocatable :: ipoi_core, ipoi_pfr, ipoi_sol ! number of point on the list if inside
  integer, dimension(:,:), allocatable :: list_tri, l_bou_qq, neighbour, list_tri_m
! list_tri- pointer from triangles to list of points (1st -- index of vertex (1-2-3), 2nd -- index of triangle)
! l_bou_qq- pointer from boundary points of the quadrangle mesh to list of points 
! (1st index -- start vertex, end vertex; 2nd -- index on the list of quadrangle boundary)
! neighbours -- pointer from triangles to neighbours (1st -- index of vertex (1-2-3), 2nd -- index of triangle
! sharing the leg 1-2-3)
  integer, dimension(:), allocatable :: list_free, list_bou_qq, ix_qq, iy_qq, i_therm 
! list of points do not belong to qadrangle mesh, to be used in outer triangulation routine
! the same as l_bou_qq, but without double points 
! (start point of the next boundary interval coinsides with the end point of the previous one) 
  integer, dimension(:), allocatable :: ind_free_s, list_free_s

  open(211,file='index_core.pnt')
  read(211,*) number_of_points, nr_core, nt_core 
  allocate(ipoi_core(nt_core,nr_core))
  do i=1, nr_core
     read(211,*) (ipoi_core(j,i), j=1,nt_core)
  enddo
  close(211)
  n_therm = (nr_core-1)*(nt_core-1)
  open(211,file='index_pfr.pnt')
  read(211,*) number_of_points, nr_pfr, nt_pfr 
  allocate(ipoi_pfr(-nt_pfr:nt_pfr, nr_pfr))
  do i=1, nr_pfr
     read(211,*) (ipoi_pfr(j,i), j=-nt_pfr,nt_pfr)
  enddo
  close(211)
  n_therm = n_therm + (nr_pfr-1)*nt_pfr*2
  open(211,file='index_sol.pnt')
  read(211,*) number_of_points, nr_sol, nt_sol 
  allocate(ipoi_sol(nt_sol,nr_sol))
  do i=1, nr_sol
     read(211,*) ipoi_sol(1:nt_sol,i)
  enddo
  close(211)
  n_therm = n_therm + (nr_sol-1)*(nt_sol-1)
  allocate(R_coord(number_of_points), Z_coord(number_of_points))
  allocate(psif(number_of_points), bmod(number_of_points), B_p(number_of_points))
  open(212,file='points.fmt')
  do i=1, number_of_points
     read(212,*)R_coord(i), Z_coord(i), psif(i), bmod(i), B_p(i)
  enddo
  read(212,*) btf, rtf
  close(212)

  allocate(list_tri(legs, nt_core+2*( (nr_core-1)*nt_core + (nr_pfr-1)*2*nt_pfr + (nr_sol-1)*(nt_sol-1))))
  allocate(list_free(2*((nr_pfr-1)*2*nt_pfr + (nr_sol-1)*(nt_sol-1))))
  allocate(l_bou_qq(2,2*((nr_pfr-1)*2*nt_pfr + (nr_sol-1)*(nt_sol-1))))
  allocate(ix_qq(nt_core+2*( (nr_core-1)*nt_core + (nr_pfr-1)*2*nt_pfr + (nr_sol-1)*(nt_sol-1))),      &
           iy_qq(nt_core+2*( (nr_core-1)*nt_core + (nr_pfr-1)*2*nt_pfr + (nr_sol-1)*(nt_sol-1))),      &
         i_therm(nt_core+2*( (nr_core-1)*nt_core + (nr_pfr-1)*2*nt_pfr + (nr_sol-1)*(nt_sol-1))))
  open(210,file='triangles.fmt')
  ntri = 0
  nfree = 0
  nbou_qq = 0
  do j=1, nt_core-1 ! around O-point
     ntri = ntri + 1
     list_tri(1,ntri) = ipoi_core(j,2)
     list_tri(2,ntri) = ipoi_core(j+1,2)
     list_tri(3,ntri) = ipoi_core(j,1)
     ix_qq(ntri) = j
     iy_qq(ntri) = 1
     i_therm(ntri) = 1
  enddo
  do i=2,nr_core-1
     do j=1, nt_core-1
        ntri = ntri + 1
        list_tri(1,ntri) = ipoi_core(j,i+1)
        list_tri(2,ntri) = ipoi_core(j+1,i+1)
        list_tri(3,ntri) = ipoi_core(j,i) 
        ix_qq(ntri) = j
        iy_qq(ntri) = i
        i_therm(ntri) = i_therm(ntri-1) + 1 
        ntri = ntri + 1
        list_tri(1,ntri) = ipoi_core(j+1,i+1)
        list_tri(2,ntri) = ipoi_core(j+1,i)
        list_tri(3,ntri) = ipoi_core(j,i)
        ix_qq(ntri) = j
        iy_qq(ntri) = i
        i_therm(ntri) = i_therm(ntri-1)
     enddo
  enddo
  iyqq_prev = iy_qq(ntri)
  do i=2,nr_pfr-1
     do j=-nt_pfr+1, nt_pfr-1
        ji = ipoi_pfr(j,i)
        j1i = ipoi_pfr(j+1,i)
        ji1 = ipoi_pfr(j,i+1)
        j1i1 = ipoi_pfr(j+1,i+1)
        if(ji .gt. 0) then
           if( j1i.gt.0 .and. ji1.gt.0 .and. j1i1.gt.0 ) then
              ntri = ntri + 1
              list_tri(1,ntri) = ji1
              list_tri(2,ntri) = j1i1
              list_tri(3,ntri) = ji
              ix_qq(ntri) = ix_qq(ntri-1) + 1 
              iy_qq(ntri) = iyqq_prev + i
              i_therm(ntri) = i_therm(ntri-1) + 1  
              ntri = ntri + 1
              list_tri(1,ntri) = j1i1
              list_tri(2,ntri) = j1i
              list_tri(3,ntri) = ji
              ix_qq(ntri) = ix_qq(ntri-1)
              iy_qq(ntri) = iyqq_prev + i
              i_therm(ntri) = i_therm(ntri-1)
           endif
        endif
     enddo
  enddo
  iyqq_prev = iy_qq(ntri)
  do i=1,nr_sol-1
     do j=1, nt_sol-1
        ji = ipoi_sol(j,i)
        j1i = ipoi_sol(j+1,i)
        ji1 = ipoi_sol(j,i+1)
        j1i1 = ipoi_sol(j+1,i+1)
        if(ji .gt. 0) then
           if( j1i.gt.0 .and. ji1.gt.0 .and. j1i1.gt.0 ) then
              ntri = ntri + 1
              list_tri(1,ntri) = ji1
              list_tri(2,ntri) = j1i1
              list_tri(3,ntri) = ji
              ix_qq(ntri) = ix_qq(ntri-1) + 1 
              iy_qq(ntri) = iyqq_prev + i
              i_therm(ntri) = i_therm(ntri-1) + 1  
              ntri = ntri + 1
              list_tri(1,ntri) = j1i1
              list_tri(2,ntri) = j1i
              list_tri(3,ntri) = ji
              ix_qq(ntri) = ix_qq(ntri-1)
              iy_qq(ntri) = iyqq_prev + i
              i_therm(ntri) = i_therm(ntri-1)
           endif
        endif
     enddo
  enddo
! compare lists of points which belong to quadrangle mesh and the whole list 
  do i=1,nr_pfr
     th_pfr: do j=-nt_pfr, nt_pfr
        ji = ipoi_pfr(j,i)
        if(ji .gt. 0) then
           do m=1, ntri
              do k=1,legs
                 if(ji .eq. list_tri(k,m)) cycle th_pfr 
              enddo
           enddo
           nfree = nfree + 1
           list_free(nfree) = ji
        endif
     enddo th_pfr
  enddo
  do i=1,nr_sol
     th_sol: do j=1, nt_sol
        ji = ipoi_sol(j,i)
        if(ji .gt. 0) then
           do m=1, ntri
              do k=1,legs
                 if(ji .eq. list_tri(k,m)) cycle th_sol 
              enddo
           enddo
           nfree = nfree + 1
           list_free(nfree) = ji
        endif
     enddo th_sol
  enddo

! -----------------------------17.07.16 begin
  allocate(ind_free_s(nfree), list_free_s(nfree))
  call sortin(list_free, ind_free_s, nfree)
  list_free_s(1:nfree) = list_free(ind_free_s(1:nfree))
  open(212,file='points_m.fmt')
  write(212,*) number_of_points - nfree 
  points: do i=1, number_of_points
     do j=1,nfree
        if(list_free_s(j) .eq. i) then
           cycle points
        endif
     enddo
     write(212,*)R_coord(i), Z_coord(i), psif(i), bmod(i), B_p(i)
  enddo points
  write(212,*) btf, rtf
  close(212)
  allocate(list_tri_m(legs, ntri))

  print *,list_free_s(1:nfree)
  do m=1, ntri
     do k=1,legs+1
        j = modulo(k,legs) + 1
        list_tri_m(j,m) = list_tri(j,m)
        free: do i= nfree,1,-1
           if(list_tri(j,m) .gt. list_free_s(i)) then
              list_tri_m(j,m) = list_tri(j,m) - i
              exit free
           endif
        enddo free
     enddo
  enddo
  open(212,file='points_m.fmt')
  read(212,*) number_of_points
  do i=1, number_of_points
     read(212,*)R_coord(i), Z_coord(i)
  enddo
  close(212)
! -----------------------------17.07.16 end



  do m=1, ntri
     do k=1,legs+1
        j = modulo(k,legs) + 1
        write(210,*) R_coord(list_tri_m(j,m)), Z_coord(list_tri_m(j,m))
     enddo
     write(210,*) 
  enddo
  close(210)

  open(2, file='17151.aug_qq.ele')
  write(2,*) ntri, legs, 0
  do m=1, ntri
     write(2,*) m, (list_tri_m(j,m), j=1,legs)
  enddo
  close(2)

  do i=1, nfree
     write(211,*) R_coord(list_free_s(i)), Z_coord(list_free_s(i))
  enddo
  close(211)
! search for the neighboring tr.'s sharing the same leg:
  allocate(xleg(legs,ntri),yleg(legs,ntri)) ! coords of midpoints of all legs
  do i=1,ntri
     do m=1,legs
        if(m .lt. legs) then
!          factor 0.5 omitted everywhere:
           xleg(m,i) = R_coord(list_tri_m(m,i)) + R_coord(list_tri_m(m+1,i))
           yleg(m,i) = Z_coord(list_tri_m(m,i)) + Z_coord(list_tri_m(m+1,i))
        else
           xleg(m,i) = R_coord(list_tri_m(m,i)) + R_coord(list_tri_m(1,i))
           yleg(m,i) = Z_coord(list_tri_m(m,i)) + Z_coord(list_tri_m(1,i))
        end if
     enddo
  enddo
  allocate(neighbour(legs, ntri))
  do k=1,ntri
     check_legs: do m=1,legs
        do i=1,ntri
           if( k .ne. i) then 
              do j=1,legs
                 if(abs(xleg(m,k)-xleg(j,i)).le.eps_dist .and. abs(yleg(m,k)-yleg(j,i)).le.eps_dist) then
                    neighbour(m,k) = i
                    cycle check_legs
                 endif
              enddo
           endif
        enddo
        nbou_qq = nbou_qq + 1
        neighbour(m,k) = -nbou_qq
        l_bou_qq(1,nbou_qq) = list_tri_m(m,k)
        if(m.lt.legs) then
           l_bou_qq(2,nbou_qq) = list_tri_m(m+1,k)
        else
           l_bou_qq(2,nbou_qq) = list_tri_m(1,k)
        endif
     enddo check_legs
  enddo
  allocate(list_bou_qq(nbou_qq+1))
  print *,'nbou_qq= ', nbou_qq
  list_bou_qq(1) = l_bou_qq(1,1)
  list_bou_qq(2) = l_bou_qq(2,1)
  l_bou_qq(1,1) = 0
  l_bou_qq(2,1) = 0
  write(212,*) R_coord(list_bou_qq(1)), Z_coord(list_bou_qq(1)),1
  write(212,*) R_coord(list_bou_qq(2)), Z_coord(list_bou_qq(2)),2
  do i=2, nbou_qq
     intervals: do j=1,nbou_qq
        do k=1,2
           if( l_bou_qq(k,j) .eq. list_bou_qq(i)) then
              if(k .eq. 1) then
                 list_bou_qq(i+1) = l_bou_qq(2,j)
              else
                 list_bou_qq(i+1) = l_bou_qq(1,j)
              endif
              l_bou_qq(1,j) = 0
              l_bou_qq(2,j) = 0
              exit intervals
           endif
        enddo
     enddo intervals
  enddo
  list_bou_qq(nbou_qq+1) = list_bou_qq(1)

  open(2, file='boundary_qq.fmt')
  do j=nbou_qq+1,1,-1
     write(2,*) list_bou_qq(j), R_coord(list_bou_qq(j)), Z_coord(list_bou_qq(j)), nbou_qq+1-j+1 
  enddo
  close(2)

  open(2, file='neighbour.fmt')
  write(2,*) nbou_qq
  do i=1, ntri
     write(2,999) neighbour(1:legs,i), ix_qq(i), iy_qq(i), i_therm(i)
  enddo
  close(2)
999 format(1000(i6,2x))
 
!  do i=1, nbou_qq+1
!     write(212,*) R_coord(list_bou_qq(1,i)), Z_coord(list_bou_qq(1,i))
!     write(212,*) R_coord(list_bou_qq(2,i)), Z_coord(list_bou_qq(2,i))
!     write(212,*) 
!     write(212,*) R_coord(list_bou_qq(i)), Z_coord(list_bou_qq(i))
!  enddo
!  close(212)

!!$  nwall_true = 0
!!$  open(1, file='FIELD/17151.wall_full')
!!$  do i=1,nwall
!!$     read(1,*,END=10) R_wall(i), dummy_r8, Z_wall(i)
!!$     nwall_true = nwall_true + 1
!!$  enddo
!!$10 continue
!!$  close(1)
!!$!  nwall_true = nwall_true - 1 ! last point = 1st, exclude it
!!$  open(2, file='17151.aug.node')
!!$  write(2,*) nwall_true + number_of_points, 2, 0, 1
!!$  write(2,*)'# volume, ', number_of_points, ' points'
!!$  do i=1, number_of_points
!!$     k = 0
!!$     do j=1,nbou_qq
!!$        if(i .eq. list_bou_qq(j)) then
!!$           k = 1
!!$           exit
!!$        endif
!!$     enddo
!!$     write(2,*) i, R_coord(i), Z_coord(i), k 
!!$  enddo  
!!$  write(2,*)'# the wall, ', nwall_true, ' points'
!!$  do i= number_of_points+1, number_of_points + nwall_true
!!$     write(2,*) i, R_wall(i-number_of_points), Z_wall(i-number_of_points), 2 
!!$  enddo
!!$  close(2)
!!$
!!$print *,nbou_qq, nwall_true
!!$
!!$
!!$  open(2, file='17151.aug.poly')
!!$  write(2,*)'# all vertices listed in .poly file: '
!!$  write(2,*) 0, 2, 0, 1
!!$  write(2,*)'#', nwall_true+nbou_qq,'  segments, each with a boundary marker'
!!$  write(2,*) nwall_true+nbou_qq,  1
!!$  write(2,*)'# boundary of the mesh:'
!!$  do j=1,nbou_qq
!!$     write(2,*) j, list_bou_qq(j), list_bou_qq(j+1), 1
!!$  enddo
!!$!  write(2,*) nbou_qq+1, list_bou_qq(nbou_qq), list_bou_qq(1), 1
!!$  write(2,*)'# the wall:'
!!$  do i=1, nwall_true-1
!!$     write(2,*) i+nbou_qq, i+number_of_points, i+number_of_points+1, 2 
!!$  enddo
!!$  write(2,*) nwall_true+nbou_qq, nwall_true + number_of_points, number_of_points+1, 2
!!$  write(2,*)'# the hole:'
!!$  write(2,*) 1
!!$  write(2,*) 1, R_coord(1), Z_coord(1)  ! O-point
end program tri_mesh
! --------------------------------
subroutine sortin(iarray, ipoi, n)
  implicit none
  integer :: n, iarray(n), ipoi(n), i, j, isave

  do  i = 1, n
     ipoi(i) = i
  enddo

  do i = 1, n-1
     do j = i+1, n
        if(iarray(ipoi(j)) .lt. iarray(ipoi(i))) then
           isave = ipoi(i)
           ipoi(i) = ipoi(j)
           ipoi(j) = isave
        endif
     enddo
  enddo
  return
end subroutine sortin
