! -------------------------------------------------------------------------------
subroutine born_particle(isort, R, phi, Z, vpar, vperp, ind_tri, cg_rng)
  use constants,   only : sq2m1, pi, ev2erg, amass
!  use mesh_mod,    only : mesh_point, mesh_element                   !<= SRR
  use mesh_mod,    only : mesh_point, mesh_element, mesh_element_rmp  !<= SRR
  use parmesh_mod, only : D_therm, V_therm, T_therm
  implicit none

  integer, intent(in) :: isort
  double precision, intent(out) :: R, phi, Z, vperp, vpar
  integer, intent(out) :: ind_tri
  double precision, intent(inout) :: cg_rng(0:5)
  real*8, dimension(:,:), allocatable :: vol_core, ind_core_tri 
  double precision :: vt, V0 
  double precision :: vx, vy, xii, vmod, Rvert(3), Zvert(3), r_bc, z_bc
  integer :: irho, itht, nmn, nmx, i, j, itri, nr_core, nt_core
!  double precision, dimension(:), allocatable :: prb_rho, prb_tht
!  double precision, dimension(:,:), allocatable :: prb_rho_tht
  double precision, dimension(:), allocatable :: prb_rho
  double precision, dimension(:,:), allocatable :: prb_rho_tht, prb_tht
  double precision ::  ran_u_01, gauss_openmp, dummy
  logical :: firstcall=.true.
  save firstcall, vol_core, ind_core_tri, prb_rho, prb_tht, prb_rho_tht, nr_core, nt_core  
  external  ran_u_01, gauss_openmp

  if(firstcall) then
     firstcall=.false.
     open(1,file='START_PRFS/source.dat',form='unformatted')
     read(1) nr_core, nt_core
     allocate(vol_core(nr_core, 2*nt_core), ind_core_tri(nr_core, 2*nt_core)) 
     read(1) vol_core
     read(1) ind_core_tri
     close(1)
!     allocate(prb_rho(0:nr_core), prb_tht(0:2*nt_core))
!     allocate(prb_rho_tht(nr_core, 2*nt_core))
     allocate(prb_rho(0:nr_core))
     allocate(prb_rho_tht(nr_core, 2*nt_core), prb_tht(0:2*nt_core,nr_core))
     prb_rho(:) = 0.d0  
!     prb_tht(:) = 0.d0  
     prb_tht(:,:) = 0.d0  
     prb_rho_tht(:,:) = 0.d0

    do i=1, nr_core
       prb_rho(i) = prb_rho(i-1)
       dummy=0.d0
       do j=1, 2*nt_core
          itri = ind_core_tri(i,j)
          if( itri .ne. 0) then
!             prb_rho_tht(i,j) = vol_core(i,j)*D_therm(mesh_element(itri)%i_therm,isort)      !<= SRR
             prb_rho_tht(i,j) = vol_core(i,j)*D_therm(mesh_element(itri)%i_therm,isort)    &  !<= SRR
                              / mesh_element_rmp(itri)%wsrr_rec_fix                           !<= SRR
             prb_rho(i) = prb_rho(i) + prb_rho_tht(i,j)
             dummy = dummy + prb_rho_tht(i,j)
          endif
          prb_tht(j,i) = dummy
       enddo
    enddo
    do i=1,nr_core
      if(prb_tht(2*nt_core,i).gt.0.d0) prb_tht(:,i)=prb_tht(:,i)/prb_tht(2*nt_core,i)
    enddo
    nmn = 0
    nmx = nr_core
    prb_rho(nmn:nmx) = prb_rho(nmn:nmx)/prb_rho(nmx)
  endif
!
  nmn = 0
  nmx = nr_core
  xii = ran_u_01(cg_rng)

  call binsrc(prb_rho, nmn, nmx, xii, irho)

!  do j=1, 2*nt_core
!     prb_tht(j) = prb_tht(j-1) + prb_rho_tht(irho,j)
!  enddo
  nmn = 0
  nmx = 2*nt_core
!  prb_tht(nmn:nmx) = prb_tht(nmn:nmx)/prb_tht(nmx)
  xii = ran_u_01(cg_rng)

!  call binsrc(prb_tht, nmn, nmx, xii, itht)
  call binsrc(prb_tht(:,irho), nmn, nmx, xii, itht)

  ind_tri = ind_core_tri(irho,itht)
  Rvert(:) = mesh_point( mesh_element(ind_tri)%i_knot(:) )%rcoord
  Zvert(:) = mesh_point( mesh_element(ind_tri)%i_knot(:) )%zcoord


           r_bc = ( mesh_point(mesh_element(ind_tri)%i_knot(1))%rcoord +                &
                    mesh_point(mesh_element(ind_tri)%i_knot(2))%rcoord +                &
                    mesh_point(mesh_element(ind_tri)%i_knot(3))%rcoord  )/3.d0
           z_bc = ( mesh_point(mesh_element(ind_tri)%i_knot(1))%zcoord +                &
                    mesh_point(mesh_element(ind_tri)%i_knot(2))%zcoord +                &
                    mesh_point(mesh_element(ind_tri)%i_knot(3))%zcoord  )/3.d0
!write(331,*) r_bc, z_bc

  vt = sqrt(2.d0*T_therm(mesh_element(ind_tri)%i_therm,isort)*ev2erg/amass(isort))
  V0 = V_therm(mesh_element(ind_tri)%i_therm,isort)

  call tri_ran(Rvert, Zvert, R, Z, cg_rng)
  phi = ran_u_01(cg_rng)*2.d0*pi
  vpar = vt*gauss_openmp(cg_rng)*sq2m1 ! /sq2 because of exp(-x**2/2) in GASDEV
  vpar = vpar + V0
  vx = vt*gauss_openmp(cg_rng)*sq2m1
  vy = vt*gauss_openmp(cg_rng)*sq2m1
  vperp = sqrt(vx**2 + vy**2)
!  vmod = sqrt(vperp**2 + vpar**2)
  return
end subroutine born_particle
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine binsrc(p, nmin, nmax, xi, i)
!
! Finds the index  i  of the array of increasing numbers   p  with dimension  n 
! which satisfies   p(i-1) <  xi  <  p(i) . Uses binary search algorithm.
!
  implicit none
!
  integer, intent(in)                                :: nmin, nmax
  double precision, dimension(nmin:nmax), intent(in) :: p
  integer, intent(out)                               :: i
  integer              :: imin, imax, k, n
  double precision     :: xi
!
  imin=nmin
  imax=nmax
  n=nmax-nmin
!
  do k=1,n
    i=(imax-imin)/2+imin
    if(p(i).gt.xi) then
      imax=i
    else
      imin=i
    endif
    if(imax.eq.imin+1) exit
  enddo
!
  i=imax
!
  return
end subroutine binsrc
!===================================================================================
subroutine tri_ran(rv, zv, r0, z0, cg_rng)
! generates random point (R0,Z0) uniformly distributed within a given 2D triangle
  use from_nrtype
  implicit none
  real(dp), dimension(3), intent(in)  :: rv, zv ! vertices  
  real(dp), intent(out) :: r0, z0  ! random point
  real(dp), intent(inout) :: cg_rng(0:5)
  real(dp) :: ran_u_01, xi1, xi3, r12, z12, r32, z32
  integer(i4b) :: iside0, iside2
  external ran_u_01

1 xi1 = ran_u_01(cg_rng)
  xi3 = ran_u_01(cg_rng)
  r12 = (rv(1) - rv(2))*xi1 
  z12 = (zv(1) - zv(2))*xi1 
  r32 = (rv(3) - rv(2))*xi3 
  z32 = (zv(3) - zv(2))*xi3 

  r0 = rv(2) + r12 + r32
  z0 = zv(2) + z12 + z32
  iside0 = nint( sign( 1.d0, ((rv(3)-rv(1))*(z0   -zv(1)) - (r0   -rv(1))*(zv(3)-zv(1))) ) )
  iside2 = nint( sign( 1.d0, ((rv(3)-rv(1))*(zv(2)-zv(1)) - (rv(2)-rv(1))*(zv(3)-zv(1))) ) )
  if (iside0 .eq. iside2) then
     return
  else if (iside0 .eq. 0) then
     go to 1
  else
     r12 = (rv(1) - rv(2))*(1.0_dp - xi1) 
     z12 = (zv(1) - zv(2))*(1.0_dp - xi1)
     r32 = (rv(3) - rv(2))*(1.0_dp - xi3)
     z32 = (zv(3) - zv(2))*(1.0_dp - xi3) 
     r0 = rv(2) + r12 + r32
     z0 = zv(2) + z12 + z32
     return
  end if

end subroutine tri_ran
! ----------------------------------------------------------------------------