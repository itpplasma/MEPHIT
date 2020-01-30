!
  use from_nrtype
  use mesh_mod
!
  implicit none
!
  integer :: i,j,k,nboupoi
  integer, dimension(:), allocatable :: notused,ipoinew,indboupoi
  double precision, dimension(:,:), allocatable :: rzpoi
!
  open(1,file='points.dat',form='unformatted')
  read (1) npoint
  allocate( mesh_point(npoint))
  read (1) mesh_point
  close(1)
!
  print *,'points:',npoint
!
  allocate(notused(npoint),ipoinew(npoint))
  notused=1
!
  open(1,file='triangles.dat',form='unformatted')
  read(1) ntri
  allocate(mesh_element(ntri))
  read(1) mesh_element
  close(1)
!
  print *,'triangles:',ntri
!
! write out all mesh points
!
  open(1,file='points_asc.dat')
  do i=1,npoint
    write (1,*) mesh_point(i)%rcoord,mesh_point(i)%zcoord
  enddo
  close(1)
!
!
! write out all mesh points ordered in triangles, find unused points
!
  open(1,file='mesh_all.dat')
  do i=1,ntri
    do j=1,3
      k=mesh_element(i)%i_knot(j)
      notused(k)=0
      write (1,*) mesh_point(k)%rcoord,mesh_point(k)%zcoord
    enddo
    k=mesh_element(i)%i_knot(1)
    write (1,*) mesh_point(k)%rcoord,mesh_point(k)%zcoord
    write (1,*) ' '
  enddo
  close(1)
!
!
! print indices of unused points
!
  do i=1,npoint
    if(notused(i).eq.1) print *,i
  enddo
!
!
! mapping to squezzed enumeration of points
!
  j=0
  do i=1,npoint
    if(notused(i).eq.1) j=j+1
    ipoinew(i)=i-j
  enddo
!
  print *,j,'unused points'
!
!
! point coordinate array with squezzed enumeration
!
  allocate(rzpoi(2,npoint-j))
  do i=1,npoint
    if(notused(i).eq.0) then
      rzpoi(1,ipoinew(i))=mesh_point(i)%rcoord
      rzpoi(2,ipoinew(i))=mesh_point(i)%zcoord
    endif
  enddo
!
  npoint=npoint-j
!
!
! remap pointers in triangles
!
  do i=1,ntri
    do j=1,3
      k=mesh_element(i)%i_knot(j)
      mesh_element(i)%i_knot(j)=ipoinew(k)
    enddo
  enddo
!
!
! test: new squezzed point array
!
  open(1,file='points_asc_sqz.dat')
  do i=1,npoint
    write(1,*) rzpoi(:,i)
  enddo
  close(1)
!
!
! test: write out all mesh points ordered in triangles again (for "diff" with previous file)
!
  open(1,file='mesh_all_sqz.dat')
  do i=1,ntri
    do j=1,3
      k=mesh_element(i)%i_knot(j)
      write (1,*) rzpoi(:,k)
    enddo
    k=mesh_element(i)%i_knot(1)
    write (1,*) rzpoi(:,k)
    write (1,*) ' '
  enddo
  close(1)
!
!
! read and remap pointers for boundary array
!
  open(1,file='boundary_qq.fmt')
  nboupoi=0
  do
    read(1,*,end=1) i
    nboupoi=nboupoi+1
  enddo
1 continue
  close(1)
!
  allocate(indboupoi(nboupoi))
  open(1,file='boundary_qq.fmt')
  do i=1,nboupoi
    read(1,*) indboupoi(i)
    indboupoi(i)=ipoinew(indboupoi(i))
  enddo
  close(1)
!
!
! test: write out boundary via squezzed enumeration
!
  open(1,file='boundary_sqz.dat')
  do i=1,nboupoi
    write (1,*) rzpoi(:,indboupoi(i))
  enddo
  close(1)
!
!
! mark boundary points (array "notused" is reused)
!
  deallocate(notused)
  allocate(notused(npoint))
  notused=0
  do i=1,nboupoi
    notused(indboupoi(i))=1
  enddo
!
!
! form the input file for Maxwell solver
!
  open(1,file='inputformaxwell.msh')
  write (1,*) npoint,ntri,nboupoi-1
!
  do i=1,npoint
    write (1,*) rzpoi(:,i),notused(i)
  enddo
!
  do i=1,ntri
    write (1,*) mesh_element(i)%i_knot(:),0
  enddo
!
  do i=1,nboupoi-1
    write (1,*) indboupoi(i),indboupoi(i+1),1
  enddo
!
  close(1)
  !

!!
! write out bmod and psipol at each triangle node
!
  open(1,file='bmod_psipol.dat')
  do i=1,ntri
     write (1,*) mesh_point(mesh_element(i)%i_knot(1))%b_mod, &
     mesh_point(mesh_element(i)%i_knot(2))%b_mod, & 
     mesh_point(mesh_element(i)%i_knot(3))%b_mod, &
     mesh_point(mesh_element(i)%i_knot(1))%psi_pol, &
     mesh_point(mesh_element(i)%i_knot(2))%psi_pol, &
     mesh_point(mesh_element(i)%i_knot(3))%psi_pol 
  enddo
  close(1)
  
end
