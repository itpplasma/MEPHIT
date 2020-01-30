module check_tri
  integer, dimension(:,:), allocatable :: in_bc
end module check_tri
!--------------------------------------------------------------------------------
function check_in_tri(ind_tri, R, Z)
  use constants,   only : one3rd
  use mesh_mod,    only : mesh_point, mesh_element, ntri
  use check_tri,   only : in_bc
  implicit none
  double precision :: R, Z, r_bc, z_bc
  double precision, dimension(3) :: rv, zv
  integer :: ind_tri, i, j, in_point
  logical :: firststep=.true., check_in_tri
  save firststep

  if(firststep) then
     firststep = .false.
     allocate(in_bc(3,ntri))
     do i=1,ntri
        do j=1,3
           rv(j) = mesh_point(mesh_element(i)%i_knot(j))%rcoord
           zv(j) = mesh_point(mesh_element(i)%i_knot(j))%zcoord
        enddo
        r_bc = sum(rv)*one3rd
        z_bc = sum(zv)*one3rd
        in_bc(1,i) = nint( sign( 1.d0, ((rv(2)-rv(1))*(z_bc - zv(1)) - (r_bc - rv(1))*(zv(2)-zv(1))) ) )
        in_bc(2,i) = nint( sign( 1.d0, ((rv(3)-rv(2))*(z_bc - zv(2)) - (r_bc - rv(2))*(zv(3)-zv(2))) ) )
        in_bc(3,i) = nint( sign( 1.d0, ((rv(1)-rv(3))*(z_bc - zv(3)) - (r_bc - rv(3))*(zv(1)-zv(3))) ) )
     enddo
  endif

  do j=1,3
     rv(j) = mesh_point(mesh_element(ind_tri)%i_knot(j))%rcoord
     zv(j) = mesh_point(mesh_element(ind_tri)%i_knot(j))%zcoord
  enddo
  in_point = nint( sign( 1.d0, ((rv(2)-rv(1))*(Z - zv(1)) - (R - rv(1))*(zv(2)-zv(1))) ) )
  if(in_point .ne. in_bc(1,ind_tri)) then
     check_in_tri = .false.
     return
  endif
  in_point = nint( sign( 1.d0, ((rv(3)-rv(2))*(Z - zv(2)) - (R - rv(2))*(zv(3)-zv(2))) ) )
  if(in_point .ne. in_bc(2,ind_tri)) then
     check_in_tri = .false.
     return
  endif
  in_point = nint( sign( 1.d0, ((rv(1)-rv(3))*(Z - zv(3)) - (R - rv(3))*(zv(1)-zv(3))) ) )
  if(in_point .ne. in_bc(3,ind_tri)) then
     check_in_tri = .false.
     return
  endif
  check_in_tri = .true.
  return

end function check_in_tri
!-----------------------------------------------------------------------------
function look_4_tri(itri_fin, R_fin, Z_fin, ierr)
  use mesh_mod,    only : mesh_point, mesh_element, ntri_inbou, inbou_list
  implicit none
  integer :: look_4_tri, itri_fin, i, ip, ivert, j, k, itri, ierr
  integer, dimension(1) :: indb
  double precision :: R_fin, Z_fin, dist(3)
  logical :: check_in_tri, in_tri
  external check_in_tri
  ierr = 0
  look_4_tri = -8000000

  do i = 1, 3
     dist(i) =  abs(R_fin - mesh_point(mesh_element(itri_fin)%i_knot(i))%rcoord)   &
              + abs(Z_fin - mesh_point(mesh_element(itri_fin)%i_knot(i))%zcoord)
  enddo
  if(maxval(dist) .gt. 3.d0*mesh_element(itri_fin)%sizmaxtri) then
     print *, 'wrong triangle in look_4_tri'
     print *, mesh_element(itri_fin)%sizmaxtri, itri_fin
     print *,R_fin,Z_fin
     print *,mesh_point(mesh_element(itri_fin)%i_knot(:))%rcoord
     print *,mesh_point(mesh_element(itri_fin)%i_knot(:))%zcoord
     stop
  endif
  indb = minloc(dist)
  ip = mesh_element(itri_fin)%i_knot(indb(1))
  in_tri=.false.
  if(ip .ne. 1) then ! not the O-point
     do i=1, mesh_point(ip)%n_owners
        in_tri = check_in_tri(mesh_point(ip)%i_owner_tri(i), R_fin, Z_fin)
        if(in_tri) then
           look_4_tri = mesh_point(ip)%i_owner_tri(i)
           return
        endif
     enddo
  else
     do i=1, ntri_inbou
        in_tri = check_in_tri(inbou_list(i)%i_tri_all, R_fin, Z_fin)
        if(in_tri) then
           look_4_tri = inbou_list(i)%i_tri_all
           return
        endif
     enddo
  endif

!  print *, 'triangle not found 1', indb(1), ip
!!$  do k=1, mesh_point(ip)%n_owners
!!$     itri = mesh_point(ip)%i_owner_tri(k)
!!$     do i = 1, 4
!!$        j = modulo(i,3) + 1
!!$        write(2222,*) mesh_point(mesh_element(itri)%i_knot(j))%rcoord,    &
!!$             mesh_point(mesh_element(itri)%i_knot(j))%zcoord
!!$     enddo
!!$     write(2222,*)
!!$  enddo
!!$  write(2223,*)R_fin, Z_fin
!!$  close(2222)
!!$  close(2223)
!!$  pause

  ivert = indb(1) + 1  
  if(ivert.gt.3) ivert = 1
  ip = mesh_element(itri_fin)%i_knot(ivert) 
  if(ip .ne. 1) then ! not the O-point
     do i=1, mesh_point(ip)%n_owners
        in_tri = check_in_tri(mesh_point(ip)%i_owner_tri(i), R_fin, Z_fin)
        if(in_tri) then
           look_4_tri = mesh_point(ip)%i_owner_tri(i)
           return
        endif
     enddo
  else
     do i=1, ntri_inbou
        in_tri = check_in_tri(inbou_list(i)%i_tri_all, R_fin, Z_fin)
        if(in_tri) then
           look_4_tri = inbou_list(i)%i_tri_all
           return
        endif
     enddo
  endif
!  print *, 'triangle not found 2', ivert, ip
!!$  do k=1, mesh_point(ip)%n_owners
!!$     itri = mesh_point(ip)%i_owner_tri(k)
!!$     do i = 1, 4
!!$        j = modulo(i,3) + 1
!!$        write(2222,*) mesh_point(mesh_element(itri)%i_knot(j))%rcoord,    &
!!$             mesh_point(mesh_element(itri)%i_knot(j))%zcoord
!!$     enddo
!!$     write(2222,*)
!!$  enddo
!!$  write(2223,*)R_fin, Z_fin
!!$  close(2222)
!!$  close(2223)
!!$  pause

  ivert = ivert + 1  
  if(ivert.gt.3) ivert = 1
  ip = mesh_element(itri_fin)%i_knot(ivert) 
  if(ip .ne. 1) then ! not the O-point
     do i=1, mesh_point(ip)%n_owners
        in_tri = check_in_tri(mesh_point(ip)%i_owner_tri(i), R_fin, Z_fin)
        if(in_tri) then
           look_4_tri = mesh_point(ip)%i_owner_tri(i)
           return
        endif
     enddo
  else
     do i=1, ntri_inbou
        in_tri = check_in_tri(inbou_list(i)%i_tri_all, R_fin, Z_fin)
        if(in_tri) then
           look_4_tri = inbou_list(i)%i_tri_all
           return
        endif
     enddo
  endif
!!$  print *, 'triangle not found 3', ivert, ip
!!$  do k=1, mesh_point(ip)%n_owners
!!$     itri = mesh_point(ip)%i_owner_tri(k)
!!$     do i = 1, 4
!!$        j = modulo(i,3) + 1
!!$        write(2222,*) mesh_point(mesh_element(itri)%i_knot(j))%rcoord,    &
!!$             mesh_point(mesh_element(itri)%i_knot(j))%zcoord
!!$     enddo
!!$     write(2222,*)
!!$  enddo
!!$  write(2223,*)R_fin, Z_fin
!!$  close(2222)
!!$  close(2223)
  ierr = -1
  return
end function look_4_tri
!-----------------------------------------------------------------------------
!!$subroutine reflect_inbou(R, Z, Rnew, Znew, ind_tri, leg_exit)
!!$  use mesh_mod,    only : mesh_point, mesh_element, ntri
!!$  use check_tri,   only : in_bc
!!$  implicit none
!!$  double precision :: R, Z, Rnew, Znew, rv1, rv2, zv1, zv2
!!$  integer :: ind_tri, leg_exit, lineside, i1, i2
!!$  
!!$  i1 = leg_exit
!!$  i2 = modulo(leg_exit,3) + 1
!!$  rv1 = mesh_point(mesh_element(ind_tri)%i_knot(i1))%rcoord
!!$  zv1 = mesh_point(mesh_element(ind_tri)%i_knot(i1))%zcoord
!!$  rv2 = mesh_point(mesh_element(ind_tri)%i_knot(i2))%rcoord
!!$  zv2 = mesh_point(mesh_element(ind_tri)%i_knot(i2))%zcoord
!!$
!!$  lineside = nint( sign( 1.d0, ((rv2 - rv1)*(Znew - zv1) - (Rnew - rv1)*(zv2 - zv1)) ) )
!!$  if(lineside .eq. in_bc(leg_exit,ind_tri) ) then
!!$     R = Rnew
!!$     Z = Znew
!!$  else if(lineside .eq. 0) then
!!$     print *,'lineside = 0'
!!$     print *,rv2, rv1, Rnew, zv2, zv1, Znew
!!$     stop
!!$  else
!!$     
!!$  endif
!!$     
!!$end subroutine reflect_inbou
