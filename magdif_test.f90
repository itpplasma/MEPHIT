program magdif_test
  
  use from_nrtype
  use constants,     only : pi, clight, amass, echarge, ev2erg, erg2ev, nsorts, charge,       &
                            one3rd, isort1, isort2
  use mesh_mod,      only : n_owners_max, legs, ntri, npoint, ntri_inbou, bphicovar,          &
                            mesh_point, mesh_element, inbou_list, grad_PhiovB2,               &
                            i_pfz, i_dpr, i_sol, i_dpl, i_inb, R_bou, Z_bou, S_bou,           &
                            mesh_element_rmp
  use sparse_mod,    only : sparse_example
  use magdif,        only : assemble_system

  implicit none

  integer :: k, npmin, npmax, n
  complex(8), allocatable, dimension(:) :: a, b, c, q, d, du, Mq
  complex(8) :: alpha
  
  print *, "MAGDIF test start"
  
  open(1,file='points.dat',form='unformatted')
  read(1) npoint
  print *, npoint
  allocate(mesh_point(npoint))
  read(1) mesh_point
  close(1)
  open(1,file='triangles.dat',form='unformatted')
  read(1) ntri
  print *, ntri
  allocate(mesh_element(ntri))
  read(1) mesh_element
  close(1)

  npmin = 52
  npmax = 101
  open(1,file='points.out')
  do k=npmin,npmax
     write(1,*) mesh_point(k)%rcoord, mesh_point(k)%zcoord
  end do
  close(1)

  n = 2
  
  allocate(a(npmax-npmin+1), b(npmax-npmin+1), c(npmax-npmin+1), q(npmax-npmin+1))
  allocate(d(npmax-npmin+1), du(npmax-npmin), Mq(npmax-npmin+1))
  do k=npmin, npmax
     a(k-npmin+1) = 1d0/sqrt(mesh_point(k)%rcoord**2 + mesh_point(k)%zcoord**2)
  end do

  b = (0d0,1d0)*n
  c = 1d0
  q = 1d0

  call assemble_system(npmax-npmin+1, a, b, c, q, d, du, alpha, Mq)

  !print *, d
  !print *, du
  !print *, alpha
  !print *, Mq

  call sparse_example(1)
  
  deallocate(a,b,c,q,d,du,Mq)
  print *, "MAGDIF test end"
end program magdif_test
