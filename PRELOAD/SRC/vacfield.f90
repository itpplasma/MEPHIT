! vacfield.f90
! writes out quantities from vacuum perturbation field
! into ascii files
! 

  use from_nrtype
  use mesh_mod
!
  implicit none
!
  integer :: i,j,k
  double precision :: R, Z, dum
  double precision :: BR_re,BR_im,BZ_re,BZ_im,bmod0
  double precision, dimension(3) :: p0
  complex(dp) :: hnorm
  type(triangle) :: elem
  type(knot), dimension(3) :: knots
  real(dp), dimension(3) :: rr, zz, lr, lz
  complex(dp), dimension(3) :: Bnflux
  complex(dp) :: Bnphiflux
!
  double complex :: B_Rn, B_Zn
!
  integer :: n_tor_out=2
!
  open(1,file='points.dat',form='unformatted')
  read (1) npoint
  allocate( mesh_point(npoint))
  read (1) mesh_point
  close(1)
!
  open(1,file='triangles.dat',form='unformatted')
  read(1) ntri
  allocate(mesh_element(ntri))
  read(1) mesh_element
  close(1)
!
  open(10,file='hpsi_vac.dat')
  open(11,file='BR_re.dat')
  open(12,file='BR_im.dat')
  open(13,file='BZ_re.dat')
  open(14,file='BZ_im.dat')
  open(15,file='bmod_psipol_vac.dat')
  open(16,file='p0.dat')
  open(17,file='Bn_flux.dat')
!
  R = mesh_point(1)%rcoord
  Z = mesh_point(1)%zcoord
  call field(R, 0d0, Z, dum, dum, dum, dum, dum, dum, dum, dum, dum, dum, dum, dum)
!
  do i=1,ntri
    R=0.d0
    Z=0.d0
    bmod0=0.d0
    do j=1,3
      k=mesh_element(i)%i_knot(j)
      R=R+mesh_point(k)%rcoord
      Z=Z+mesh_point(k)%zcoord
      bmod0=bmod0+mesh_point(k)%b_mod
    enddo
    R=R/3.d0
    Z=Z/3.d0
    bmod0=bmod0/3.d0
!
    call spline_bpol_n(n_tor_out, R, Z, B_Rn, B_Zn)
!
    hnorm=(mesh_element(i)%dPsi_dR*B_Rn+mesh_element(i)%dPsi_dZ*B_Zn)/bmod0
!
    write (10,*) real(hnorm),dimag(hnorm)
!
    BR_re=real(B_Rn)
    BR_im=dimag(B_Rn)
    BZ_re=real(B_Zn)
    BZ_im=dimag(B_Zn)
!
    write (11,*) i,BR_re
    write (12,*) i,BR_im
    write (13,*) i,BZ_re
    write (14,*) i,BZ_im

    write (15,*) mesh_point(mesh_element(i)%i_knot(1))%b_mod, &
     mesh_point(mesh_element(i)%i_knot(2))%b_mod, & 
     mesh_point(mesh_element(i)%i_knot(3))%b_mod, &
     mesh_point(mesh_element(i)%i_knot(1))%psi_pol, &
     mesh_point(mesh_element(i)%i_knot(2))%psi_pol, &
     mesh_point(mesh_element(i)%i_knot(3))%psi_pol

    write (16,*) mesh_point(mesh_element(i)%i_knot)%pl_parms_knot(1,1)&
         *mesh_point(mesh_element(i)%i_knot)%pl_parms_knot(3,1) &
         / mesh_point(mesh_element(i)%i_knot)%pl_parms_knot(2,1), &
         mesh_point(mesh_element(i)%i_knot)%pl_parms_knot(1,2)&
         *mesh_point(mesh_element(i)%i_knot)%pl_parms_knot(3,2) &
         / mesh_point(mesh_element(i)%i_knot)%pl_parms_knot(2,2)
    !
    
    elem = mesh_element(i)
    knots = mesh_point(elem%i_knot)
       
    rr(1) = (knots(1)%rcoord + knots(2)%rcoord)/2d0 ! edge midpoints
    rr(2) = (knots(2)%rcoord + knots(3)%rcoord)/2d0
    rr(3) = (knots(3)%rcoord + knots(1)%rcoord)/2d0
    
    zz(1) = (knots(1)%zcoord + knots(2)%zcoord)/2d0
    zz(2) = (knots(2)%zcoord + knots(3)%zcoord)/2d0
    zz(3) = (knots(3)%zcoord + knots(1)%zcoord)/2d0
    
    lr(1) = knots(2)%rcoord - knots(1)%rcoord ! edge vector components
    lr(2) = knots(3)%rcoord - knots(2)%rcoord ! (counterclockwise)
    lr(3) = knots(1)%rcoord - knots(3)%rcoord
    
    lz(1) = knots(2)%zcoord - knots(1)%zcoord
    lz(2) = knots(3)%zcoord - knots(2)%zcoord
    lz(3) = knots(1)%zcoord - knots(3)%zcoord
    
    do k=1,3
       call spline_bpol_n(n_tor_out, rr(k), zz(k), B_Rn, B_Zn)
       Bnflux(k) = lz(k)*rr(k)*B_Rn - lr(k)*rr(k)*B_Zn ! project vector onto edge normal
    end do
    Bnphiflux = -1d0/((0d0,1d0)*n_tor_out)*sum(Bnflux) ! make divergence free
    write(17,*) real(Bnflux(1)), aimag(Bnflux(1)),&
         real(Bnflux(2)), aimag(Bnflux(2)),&
         real(Bnflux(3)), aimag(Bnflux(3)),&
         real(Bnphiflux), aimag(Bnphiflux)
  enddo
!
  close(10)
  close(11)
  close(12)
  close(13)
  close(14)
  close(15)
  close(16)
  close(17)
!
  end
