program modify
!  use constants,     only : pi, nsorts, charge, one3rd 
  use mesh_mod,      only :  npoint, mesh_point, mesh_point_new
  implicit none
!
  integer          :: i
!
! read 
  open(1,file='RESULTS/points.dat',form='unformatted')
  read(1) npoint
  allocate(mesh_point(npoint))
  read(1) mesh_point
  close(1)


  allocate(mesh_point_new(npoint))

  do i=1,npoint
    mesh_point_new(i)%rcoord = mesh_point(i)%rcoord 
     mesh_point_new(i)%zcoord = mesh_point(i)%zcoord
     mesh_point_new(i)%psi_pol = mesh_point(i)%psi_pol
     mesh_point_new(i)%b_mod = mesh_point(i)%b_mod
     mesh_point_new(i)%b_phi = mesh_point(i)%b_phi
     mesh_point_new(i)%PhiovB2 = mesh_point(i)%PhiovB2
     mesh_point_new(i)%pl_parms_knot(1,:) = mesh_point(i)%pl_parms_knot(1,:)
     mesh_point_new(i)%pl_parms_knot(2,:) = mesh_point(i)%pl_parms_knot(2,:)
     mesh_point_new(i)%pl_parms_knot(3,:) = mesh_point(i)%pl_parms_knot(3,:)
     mesh_point_new(i)%pl_parms_knot(4,:) = 0.d0
     mesh_point_new(i)%n_owners = mesh_point(i)%n_owners
     mesh_point_new(i)%i_owner_tri(:) = mesh_point(i)%i_owner_tri(:)
     mesh_point_new(i)%weight_intp(:) = mesh_point(i)%weight_intp(:)
  enddo
  open(1,file='RESULTS/points_new.dat',form='unformatted')
  write(1) npoint
  write(1) mesh_point_new
  close(1)
end program modify
