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
     real(dp) :: rcoord = 0d0
     real(dp) :: zcoord = 0d0
     real(dp) :: psi_pol = 0d0
     real(dp) :: b_mod = 0d0
     real(dp) :: b_phi = 0d0
     real(dp) :: PhiovB2 = 0d0
     real(dp), dimension(4,nsorts) :: pl_parms_knot = 0d0 ! 1-n, 2-V, 3-T, 4-Tperp (interpolated to mesh points)
     integer(i4b) :: n_owners = 0 !true number of triangles own the point
     integer(i4b), dimension(n_owners_max) :: i_owner_tri = 0 ! indexes of triangles own the point
     real(dp), dimension(n_owners_max) :: weight_intp = 0d0 ! weights of tr.'s for interpolation
!     logical :: qq_bc  ! barycenter of qq-cell, [y|n]
  end type knot
  type(knot), dimension(:), allocatable :: mesh_point
!------------------------------------------------------------------------------------
  type :: triangle
     sequence
     integer(i4b) :: ix_qq = 0
     integer(i4b) :: iy_qq = 0
!     integer(i4b) :: iqq_gc ! guard cell, obsolete (0-reg.cell, -1-left DP, -2-right DP, 1-core, 2-out.bound.GC), 3-PFZ  
     integer(i4b), dimension(legs) :: i_knot = 0 ! index of the vertice in the "knot" list
     integer(i4b), dimension(legs) :: neighbour = 0 ! index of the next tr. (sharing the borderline # LEG)
     integer(i4b), dimension(legs) :: neighbour_edge = 0! index of the next tr. edge sharing the borderline # LEG
     real(dp) :: sizmaxtri = 0d0 ! max. size of the triangle
     real(dp) :: det_3 = 0d0  ! (r2-r1)(z3-z1) - (r3-r1)(z2-z1)
     real(dp) :: oneoverD3 = 0d0  ! 1/det_3
     real(dp) :: V_tri = 0d0  ! triangular cell volume
! derivatives (constant within each triangle)
     real(dp) :: dBinv_dR = 0d0
     real(dp) :: dBinv_dZ = 0d0
     real(dp) :: dPsi_dR = 0d0
     real(dp) :: dPsi_dZ = 0d0
     real(dp) :: dPhiovB2_dR = 0d0
     real(dp) :: dPhiovB2_dZ = 0d0
     ! for linear interpolation of bphicovar
     integer(i4b) :: knot_h = 0
     real(dp) :: dbphicovdpsi = 0d0
! test particles parameters:
     real(dp), dimension(nsorts) :: Dm_part = 0d0 ! marker density
     real(dp), dimension(nsorts) :: D_part = 0d0  ! particle density
     real(dp), dimension(nsorts) :: V_part = 0d0
     real(dp), dimension(nsorts) :: T_part = 0d0
     real(dp), dimension(nsorts) :: Vt2_part = 0d0 ! for C-G-L pressure
     real(dp) :: ePhi_tri = 0d0 ! potential related to triangle
! index of the triangle in the list of thermostat cells:
     integer(i4b) :: i_therm = 0
!
     double precision, dimension(2,nsorts)    :: thermforces = 0d0
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
     real(dp)                                 :: wsrr_rec_fix            !<= SRR
! bnorm - scalar product of perturbed unit vector along the manetic field with gradient
! of unperturbed poloidal flux (bnorm_vac - for vacuum field perturbation):
     complex(dp)                              :: bnorm_vac
! products of bnorm with thermodynamic forces comuted using poloidal flux as radial variable:
     complex(dp),      dimension(2,nsorts)    :: bnorm_times_thermforces
! perturbed currents through prism faces ("triangle edges"):
     complex(dp),      dimension(legs,nsorts) :: currents
! perturbations of density, perpendicular and parallel stress tensor components
! (for slow rotations these are pressure tensor components) and parallel current density:
     complex(dp),      dimension(nsorts)      :: denspert,prespert_perp,prespert_par,parcurrpert
     !> area of the triangle
     real(dp) :: area = 0d0
     !> indices of nodes in #mesh_point for edge i, going counter-clockwise
     integer :: li(2) = [0, 0]
     !> indices of nodes in #mesh_point for edge o, going counter-clockwise
     integer :: lo(2) = [0, 0]
     !> indices of nodes in #mesh_point for edge f, going counter-clockwise
     integer :: lf(2) = [0, 0]
     !> index of edge i in #mesh_element
     integer :: ei = 0
     !> index of edge o in #mesh_element
     integer :: eo = 0
     !> index of edge f in #mesh_element
     integer :: ef = 0
     !> orientation of the triangle: if true, edge f lies on the outer flux surface
     logical :: orient = .true.
     !> \f$ R \f$ coordinate of triangle 'center'
     real(dp) :: R_Omega = 0d0
     !> \f$ Z \f$ coordinate of triangle 'center'
     real(dp) :: Z_Omega = 0d0
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
