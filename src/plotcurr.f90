!
  use from_nrtype
  use mesh_mod
!
  implicit none
!
  integer :: i,j,k,nboupoi
  double precision, dimension(3) :: R_vert,Z_vert
  double precision, dimension(6) :: reimcur
  double complex :: cur1,cur2,delta,cur_R,cur_Z
!
  open(1,file='START_PRFS/points.dat',form='unformatted')
  read (1) npoint
  allocate( mesh_point(npoint))
  read (1) mesh_point
  close(1)
!
  open(1,file='START_PRFS/triangles.dat',form='unformatted')
  read(1) ntri
  allocate(mesh_element(ntri))
  read(1) mesh_element
  close(1)
!
  open(1,file='RMP_EXCHANGE/currents.dat')
  open(2,file='PLOTTING/curr_toplot.dat')
  do i=1,ntri
    read(1,*) reimcur
    cur1=cmplx(reimcur(1),reimcur(4))
    cur2=cmplx(reimcur(2),reimcur(5))
    R_vert(:) = mesh_point( mesh_element(i)%i_knot(:) )%rcoord
    Z_vert(:) = mesh_point( mesh_element(i)%i_knot(:) )%zcoord
    delta=(Z_vert(2)-Z_vert(1))*(R_vert(2)-R_vert(3)) &
         -(Z_vert(2)-Z_vert(3))*(R_vert(2)-R_vert(1))
    cur_R=(cur1*(R_vert(2)-R_vert(3))+cur2*(R_vert(2)-R_vert(1)))/delta
    cur_Z=(cur1*(Z_vert(2)-Z_vert(3))+cur2*(Z_vert(2)-Z_vert(1)))/delta
    write (2,*) real(cur_R),real(cur_Z),dimag(cur_R),dimag(cur_Z)
  enddo
  close(1)
  close(2)
!
  open(1,file='PLOTTING/dens_temp_ePhi.dat')
  do i=1,ntri
    write (1,*) mesh_element(i)%D_part,mesh_element(i)%T_part,mesh_element(i)%ePhi_tri, &
                mesh_element(i)%thermforces
  enddo
  close(1)
!
  end
