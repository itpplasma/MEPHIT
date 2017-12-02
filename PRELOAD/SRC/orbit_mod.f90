module constants
  use from_nrtype
  double precision, parameter :: amp=1.6726d-24, ame=9.1094d-28, clight=2.9979d10,  &
       ev2erg=1.6022d-12,  erg2ev=1.d0/ev2erg, echarge=4.8032d-10,              &
       pi=3.141592653589793238462643383d0, Cs0=9.79d5*9.12870929d-1
  integer(i4b), parameter :: nsorts=2
  double precision, parameter, dimension(nsorts) :: amass=(/ ame, amp*2.d0 /), charge =(/-1.d0, 1.d0/)
  double precision, dimension(nsorts, nsorts) :: rmass, c_c_2 ! mass relation, charge**2
  double precision, parameter :: one3rd=1.d0/3.d0, sq2=1.414213562373095049d0,     &
                                 sq2m1=1.d0/sq2, sq3=1.73205080756887729353d0, sq3m1=1.d0/sq3
  double precision :: sigma_random=1.d0
  integer(i4b), parameter :: isort1 = 1, isort2 = 2
end module constants
! ------------------------------------------------------------------------
module accuracy_mod
  use from_nrtype
  integer(i4b), parameter :: nbisec=10
  double precision, parameter :: relerr=1.d-6, eps_newt=1.d-8, eps_cell = 1.d-3
  double precision, parameter :: slowstep=0.1d0
  double precision, parameter :: eps_dp = epsilon(1.d0), tinydp = tiny(1.d0)
end module accuracy_mod
!
module parmot_mod
  use from_nrtype
  integer(i4b) :: ichsvr, ichsvz
  double precision :: ro0, eeff, vgr, vgz, vlam, dtaumin, rmn, rmx, zmn, zmx
  double precision :: rmu=1.d35 ! inverse relativistic temperature
end module parmot_mod
!
!  module backpar_mod
!    double precision :: R0, Z0, Phi_R0, Phi_Z0, Phi_RZ0
!  end module backpar_mod
!
module renorm_mod
  use constants, only : nsorts
  double precision, dimension(nsorts) :: vref, rhoref, elpref
end module renorm_mod
!
module parmesh_mod
  use from_nrtype
  integer(i4b) :: n_therm  ! number of thermostat cells \le ntri
  double precision, dimension(:,:), allocatable :: w_cell_V, w_cell_T, D_therm, V_therm, T_therm
  double precision, dimension(:),   allocatable :: vol_therm
  double precision, dimension(:,:),   allocatable :: time_coll, sum_tr, v_tr, en_mean_tr, wsrr_rec, t4Dm
  double precision, dimension(:,:),   allocatable :: vt2_tr 
  double precision, dimension(:,:),   allocatable :: j_tr, q_tr, tm_tr
end module parmesh_mod
!
module for_macrostep
  use from_nrtype
  use constants, only : nsorts
  double precision :: dt0, wdamp_rate
  double precision, dimension(nsorts) :: dt0_factor
  integer(i4b)          :: nptrace
  double precision, dimension(nsorts) :: time_l, vt_s, wcv, wct
  double precision :: source0, relax_tau, D_perp, fulltime, wpartsrr
  double precision :: t_min=2.d1, d_min=1.d11 
  double precision :: perp_jump = 1.0d-1, par_jump = 4.0d-1
  double complex :: wpart_rmp
  integer :: n_tormode
  logical :: do_averages
end module for_macrostep
!
!!$  module contrbound_mod
!!$    logical :: checkb=.false.
!!$    integer :: icoord
!!$    double precision :: rzbound,godir
!!$  end module contrbound_mod
!
module inbou_mod
  use from_nrtype
  integer(i4b) :: npoi
  double precision, dimension(:),   allocatable :: bstart, volstart 
  double precision, dimension(:,:), allocatable :: xstart
end module inbou_mod

module ranstate
  use from_nrtype
  integer(i4b), parameter :: ntab=32
  integer(i4b) :: idum, idum2=123456789, iv(1:ntab)=0, iy=0
end module ranstate

module hujnja
  use from_nrtype
  use constants, only : nsorts
  integer(i4b) :: numstep, numstepmax, num_see_saw
  logical :: firstpart=.true.
!  logical :: prev_perp
!  double precision :: ubila_rr, dobavila_rr, istochnik, stenka, proebal
  integer :: laststeps=1
!  double precision, dimension(:), allocatable :: dens_avr_2d, temp_avr_2d
  double precision, dimension(:), allocatable :: rel_therm_max
end module hujnja
!
!!$module parmesh_fast_mod
!!$  use from_nrtype
!!$  integer(i4b) :: npoi_mag
!!$  integer(i4b),          dimension(:,:),   allocatable :: ipoi_mag
!!$  double precision :: bphicovar
!!$  double precision, dimension(:),     allocatable :: psi_arr,binv_arr
!!$  double precision, dimension(:),     allocatable :: phiovb2_arr
!!$end module parmesh_fast_mod
!
module reg_step_tria_mod
  use from_nrtype
  double precision :: abr,abz,aphir,aphiz,apsir,apsiz
  double precision :: binv_vert1,phiovb2_vert1,Rcell,sizmax
  double precision :: deltatri
  double precision, dimension(3) :: R_vert,Z_vert
end module reg_step_tria_mod

module for_mpi
  use from_nrtype
  use constants, only : nsorts
  integer(i4b) :: mpi_p_root=0, mype, npes
  double precision, dimension(nsorts) :: time_local
  double precision, dimension(:,:), allocatable :: Dm_local, Dp_local, Vp_local, Tp_local, dummy8, dummy8_2, Vt2_local 
  double precision, dimension(:), allocatable :: ephi_sheath, ephi_local, dummy8_1
  double precision, dimension(:,:), allocatable :: flux_q, flux_j, tor_mom, dens_phi

  character*2 mype_sym
  character*15 evolname
end module for_mpi
