module mesh_mod
  use from_nrtype
  use constants, only : nsorts

  integer(i4b), parameter :: n_owners_max=16, legs=3, npqq=4 ! 5 for version with barycenter
  integer(i4b) :: ntri, npoint ! number of trs. & vertices
  integer(i4b) :: ntri_inbou ! number of trs. building inner boundary 
  real(dp)     :: bphicovar  ! Bphi = btf*rtf/rrr; bphicovar=btf*rtf=const.
  integer(i4b), dimension(2) :: i_pfz, i_dpr, i_sol, i_dpl, i_inb 
  real(dp), dimension(:), allocatable :: R_bou, Z_bou, S_bou
!------------------------------------------------------------------------------------
  type :: knot
     sequence
     real(dp) :: rcoord 
     real(dp) :: zcoord
     real(dp) :: psi_pol
     real(dp) :: b_mod
     real(dp) :: b_phi
     real(dp) :: PhiovB2
     real(dp), dimension(4,nsorts) :: pl_parms_knot ! 1-n, 2-V, 3-T, 4-Tperp (interpolated to mesh points)
     integer(i4b) :: n_owners    !true number of triangles own the point    
     integer(i4b), dimension(n_owners_max) :: i_owner_tri  ! indexes of triangles own the point
     real(dp), dimension(n_owners_max) :: weight_intp ! weights of tr.'s for interpolation
!     logical :: qq_bc  ! barycenter of qq-cell, [y|n]
  end type knot
  type(knot), dimension(:), allocatable :: mesh_point
!------------------------------------------------------------------------------------
  type :: triangle
     sequence
     integer(i4b) :: ix_qq
     integer(i4b) :: iy_qq
!     integer(i4b) :: iqq_gc ! guard cell, obsolete (0-reg.cell, -1-left DP, -2-right DP, 1-core, 2-out.bound.GC), 3-PFZ  
     integer(i4b), dimension(legs) :: i_knot ! index of the vertice in the "knot" list
     integer(i4b), dimension(legs) :: neighbour ! index of the next tr. (sharing the borderline # LEG)
     integer(i4b), dimension(legs) :: neighbour_edge ! index of the next tr. edge sharing the borderline # LEG
     real(dp) :: sizmaxtri ! max. size of the triangle
     real(dp) :: det_3  ! (r2-r1)(z3-z1) - (r3-r1)(z2-z1) 
     real(dp) :: oneoverD3  ! 1/det_3
     real(dp) :: V_tri  ! triangular cell volume  
! derivatives (constant within each triangle)
     real(dp) :: dBinv_dR
     real(dp) :: dBinv_dZ
     real(dp) :: dPsi_dR
     real(dp) :: dPsi_dZ
     real(dp) :: dPhiovB2_dR
     real(dp) :: dPhiovB2_dZ
! test particles parameters:
     real(dp), dimension(nsorts) :: Dm_part  ! marker density
     real(dp), dimension(nsorts) :: D_part   ! particle density
     real(dp), dimension(nsorts) :: V_part
     real(dp), dimension(nsorts) :: T_part
     real(dp), dimension(nsorts) :: Vt2_part ! for C-G-L pressure
     real(dp) :: ePhi_tri ! potential related to triangle
! index of the triangle in the list of thermostat cells:
     integer(i4b) :: i_therm
!
     double precision, dimension(2,nsorts)    :: thermforces
  end type triangle
  type(triangle), dimension(:), allocatable :: mesh_element
!------------------------------------------------------------------------------------
  type :: bou_list
! used for interpolation in O-point now
     sequence
     integer(i4b) :: i_tri_all ! index of boundary tr. in the full list 
     real(dp) :: vol_norm ! volume normalized to 1 
  end type bou_list
  type(bou_list), dimension(:), allocatable :: inbou_list
!------------------------------------------------------------------------------------
  type :: triangle_rmp
     sequence
     double complex                           :: bnorm_vac
     double complex,   dimension(2,nsorts)    :: bnorm_times_thermforces
     double complex,   dimension(legs,nsorts) :: currents
  end type triangle_rmp
  type(triangle_rmp), dimension(:), allocatable :: mesh_element_rmp
!------------------------------------------------------------------------------------
contains
!
! A*(x-x1) + B*(y-y1) + C*(z-z1) = 0
!
!      |y2-y1  z2-z1|         |x2-x1  z2-z1|       |x2-x1  y2-y1| 
!  A = |            |,  B = - |            |,  C = |            |
!      |y3-y1  z3-z1|         |x3-x1  z3-z1|       |x3-x1  y3-y1|
!
!  z = z1 - A/C*(x-x1) - B/C*(y-y1)
!  dz/dx = -A/C,  dz/dy = -B/C
!
  function cell_linint(icell, r, z, f)
    use from_nrtype
    implicit none
    real(dp) :: cell_linint
    real(dp), dimension(legs), intent(in) :: f
    real(dp), intent(in) :: r, z  ! coords of the point in tr.
! index of triangle, index of plasma parameter to interpolate: 
    integer(i4b), intent(in) :: icell
    real(dp) :: r2r1, r3r1, z2z1, z3z1, f2f1, f3f1, f1

    r2r1 = mesh_point(mesh_element(icell)%i_knot(2))%rcoord -                           & 
           mesh_point(mesh_element(icell)%i_knot(1))%rcoord

    r3r1 = mesh_point(mesh_element(icell)%i_knot(3))%rcoord -                           & 
           mesh_point(mesh_element(icell)%i_knot(1))%rcoord

    z2z1 = mesh_point(mesh_element(icell)%i_knot(2))%zcoord -                           & 
           mesh_point(mesh_element(icell)%i_knot(1))%zcoord

    z3z1 = mesh_point(mesh_element(icell)%i_knot(3))%zcoord -                           & 
           mesh_point(mesh_element(icell)%i_knot(1))%zcoord

    f2f1 = f(2) - f(1)
    f3f1 = f(3) - f(1)

    cell_linint = f(1) + (z3z1*f2f1 - z2z1*f3f1)*mesh_element(icell)%oneoverD3*           &
                           (r - mesh_point(mesh_element(icell)%i_knot(1))%rcoord)         &
                       + (r2r1*f3f1 - r3r1*f2f1)*mesh_element(icell)%oneoverD3*           &
                           (z - mesh_point(mesh_element(icell)%i_knot(1))%zcoord)
    return
  end function cell_linint
!------------------------------------------------------------------------------------
  subroutine grad_PhiovB2(icell)
    use from_nrtype
    implicit none
    integer(i4b), intent(in) :: icell
    real(dp) :: r2r1, r3r1, z2z1, z3z1, f2f1, f3f1, PhiovB2_1

    r2r1 = mesh_point(mesh_element(icell)%i_knot(2))%rcoord -                      & 
           mesh_point(mesh_element(icell)%i_knot(1))%rcoord

    r3r1 = mesh_point(mesh_element(icell)%i_knot(3))%rcoord -                      & 
           mesh_point(mesh_element(icell)%i_knot(1))%rcoord

    z2z1 = mesh_point(mesh_element(icell)%i_knot(2))%zcoord -                      & 
           mesh_point(mesh_element(icell)%i_knot(1))%zcoord

    z3z1 = mesh_point(mesh_element(icell)%i_knot(3))%zcoord -                      & 
           mesh_point(mesh_element(icell)%i_knot(1))%zcoord

    PhiovB2_1 = mesh_point(mesh_element(icell)%i_knot(1))%PhiovB2
    f2f1 =  mesh_point(mesh_element(icell)%i_knot(2))%PhiovB2 - PhiovB2_1
    f3f1 =  mesh_point(mesh_element(icell)%i_knot(3))%PhiovB2 - PhiovB2_1
    mesh_element(icell)%dPhiovB2_dR = (z3z1*f2f1 - z2z1*f3f1)*mesh_element(icell)%oneoverD3
    mesh_element(icell)%dPhiovB2_dZ = (r2r1*f3f1 - r3r1*f2f1)*mesh_element(icell)%oneoverD3

    return
  end subroutine grad_PhiovB2
end module mesh_mod
