program result_spectrum

  use from_nrtype
  use constants,     only : pi, clight, amass, echarge, ev2erg, erg2ev, nsorts, charge,       &
                            one3rd, isort1, isort2
  use renorm_mod,    only : vref, rhoref, elpref
!  use backpar_mod,   only : R0, Z0, Phi_R0, Phi_Z0, Phi_RZ0
  use mesh_mod,      only : n_owners_max, legs, ntri, npoint, ntri_inbou, bphicovar,           &
                            mesh_point, mesh_element, inbou_list, grad_PhiovB2,                &
                            i_pfz, i_dpr, i_sol, i_dpl, i_inb, R_bou, Z_bou, S_bou,            &
                            mesh_element_rmp
  use parmesh_mod,   only : n_therm, w_cell_V, w_cell_T, vol_therm, D_therm, V_therm, T_therm, &
                            time_coll, sum_tr, v_tr, en_mean_tr, t4Dm, j_tr, q_tr, tm_tr,      &
                            vt2_tr
  use for_macrostep, only : dt0, nptrace, time_l, vt_s, source0,                               &
                            relax_tau, wcv, wct, D_perp, fulltime, d_min, t_min,               &
                            n_tormode
  use ranstate,      only : idum, idum2, iv, iy
  use hujnja, only  : numstepmax, rel_therm_max !, v3_sum_in, v3_sum_out, v3_in, v3_out, part_num
!, ubila_rr, dobavila_rr, istochnik, stenka, proebal
  use for_mpi, only : mype, npes, mpi_p_root, time_local, Dm_local, Dp_local, Vp_local,        & 
                      Tp_local, ephi_local, mype_sym, evolname,              &
                      flux_j, flux_q, tor_mom, Vt2_local, ephi_sheath, dens_phi 
!#ifdef PARALLEL
!  use mpi
!!  include 'mpif.h'
!#endif
!
  implicit none
!
  integer(i4b)          :: ierr, mpi_triangle, mpi_bound, mpi_mesh_point, nboucell
  logical               :: restart
  integer(i4b)          :: i, j, iptrace, k, ind_therm, n_outbou, is
  integer(i4b)          :: isort, icount, npgroup, macrosteps, n_prfs, n_evol
  double precision :: bmod_ref, bmod00, oneovernp
  double precision :: R, phi, Z, T_source, ePhi_knot, r_bc, z_bc
  double precision, dimension(nsorts) :: T_knot, V_knot, D_knot,                         &
       dp_j_l, dp_j_r, dp_q_l, dp_q_r, dp_tm_l, dp_tm_r,                                 &
       j_outer_wall, j_pfz_wall, q_outer_wall, q_pfz_wall, tm_outer_wall, tm_pfz_wall
  double precision :: dens_tp_avr, temp_tp_avr, v2_tp_avr, v2_avr, temp_avr
! random numbers:
  double precision :: ran_u_01, gauss_openmp, ran2
  integer(i4b), parameter :: max_openmp_threads=64
  double precision, dimension(0:5), save :: seed, cg_rng, saverandom
  double precision, dimension(0:5,0:max_openmp_threads-1) :: cg_rng_sav
!
  double precision, dimension(:,:), allocatable :: D_Dm, D_Dp, D_Vp, D_Tp, D_Dt, D_Vt, D_Tt ! dispersions
  double precision, dimension(:), allocatable :: D_ephi
  double precision, dimension(:,:), allocatable :: D_j, D_q, D_tm
  double precision, dimension(nsorts) :: D_time_l, part_sink, temp_out
  double precision :: dummy_re,dummy_im
  double complex   :: bnorm_vac,bnorm_plas

  ! TODO: get this from axis.f90
  integer :: nr_core=1000
  double precision :: R_o=172.74461503599866, Z_o=10.313710665604722
  double precision :: R_x=142.27640211007648, Z_x=-95.294670798084724
  double precision :: delta_xo

  double precision :: dt, vpar, vperp
  integer :: ind_tri, ind_tri_orig, iedge_ent
  logical :: intri, check_in_tri
  
  double precision :: q
  double complex, allocatable :: bpsi(:)
  integer, parameter :: mmax = 24
  double complex :: bpsim(2*mmax+1)

  character(len=1024) :: temp

  external  ran_u_01, check_in_tri
!#ifdef PARALLEL
!  call MPI_INIT(ierr)
!  if(ierr .ne. MPI_SUCCESS) then 
!     write(*,*)'MPI init. err.'
!     stop 1
!  endif
!  call MPI_COMM_SIZE(MPI_COMM_WORLD,npes,ierr)
!  call MPI_COMM_RANK(MPI_COMM_WORLD,mype,ierr)
!#else
  npes = 1
  mype = 0
!#endif

  oneovernp = 1.d0/dble(npes)
  if(max_openmp_threads .lt. npes) stop 'too few threads'

  if (mype .eq. mpi_p_root) then
! read everything
     open(1,file='kin2d.inp')
     read (1,*) bmod_ref                   ! reference field, G, for Boozer $B_{00}$
     read (1,*) nptrace, source0, T_source ! number of particles to trace, source intencity  
     read (1,*) relax_tau                  ! underrelaxation factor 
     read (1,*) wcv(1), wct(1)             ! weights of cells, electrons
     read (1,*) wcv(2), wct(2)             ! weights of cells, ions
     read (1,*) D_perp                     ! anomalous diff.coeff., cm^2/s
     read (1,*) n_evol, n_prfs             ! number of steps between outputs (average vals. and profiles)
     read (1,*) npgroup                    ! number of particles per macrostep
     read (1,*) restart                    ! restart from last output
     close(1)

     open(1,file='points.dat',form='unformatted')
     read(1) npoint
     allocate(mesh_point(npoint))
     read(1) mesh_point
     close(1)
     open(1,file='triangles.dat',form='unformatted')
     read(1) ntri
     allocate(mesh_element(ntri))
     read(1) mesh_element
     read(1) bphicovar, fulltime, time_l
     close(1)
     allocate(mesh_element_rmp(ntri))
!
     n_tormode=2
     dt0=2d-4
!
     do i=1,ntri
       mesh_element_rmp(i)%currents(:,:) = (0.d0,0.d0)
       mesh_element_rmp(i)%bnorm_times_thermforces=mesh_element(i)%thermforces &
                                                  *(0.d0,0.d0)  
     enddo

!print *,sum(flux_j),sum(flux_q),sum(tor_mom)
!stop
     open(1,file='thermostat.dat',form='unformatted')
     read(1) n_therm
     allocate(D_therm(n_therm,nsorts),V_therm(n_therm,nsorts),T_therm(n_therm,nsorts))
     allocate(vol_therm(n_therm))
     read(1) D_therm
     read(1) V_therm
     read(1) T_therm
     read(1) vol_therm
     close(1)

!!$allocate(v3_sum_in(n_therm), v3_sum_out(n_therm), v3_in(n_therm), v3_out(n_therm), part_num(n_therm))
!!$v3_sum_in(:) = 0.d0
!!$v3_sum_out(:) = 0.d0
!!$part_num = 0.d0

     open(1,file='boundary.dat',form='unformatted')
     read(1) ntri_inbou
     allocate(inbou_list(ntri_inbou))
     read(1) inbou_list
!!$     read(1) i_pfz, i_dpr, i_sol, i_dpl, i_inb
!!$     nboucell = i_inb(2) - i_pfz(1) + 1
!!$     allocate(R_bou(nboucell), Z_bou(nboucell), S_bou(nboucell))
!!$     read(1) R_bou, Z_bou, S_bou
     close(1)

     open(1,file='wall_loads.dat',form='unformatted')
     read(1) nboucell
     allocate(flux_j(nboucell,nsorts), flux_q(nboucell,nsorts), tor_mom(nboucell,nsorts))
     read(1) flux_j
     read(1) flux_q
     read(1) tor_mom
     close(1)

  endif

!#ifdef PARALLEL
!  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!
!  call broadcast_all
!
!  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!#endif

  allocate(ephi_sheath(nboucell))

  do isort=1,nsorts
     mesh_element(1:ntri)%Dm_part(isort) = 0.d0 
  enddo
!!$  write(mype_sym,'(i2)') mype
!!$  if(mype .le. 9) then
!!$     evolname = 'RESULTS/evol.0' // trim(adjustl(mype_sym))
!!$  else
!!$     evolname = 'RESULTS/evol.'  // trim(adjustl(mype_sym))
!!$  endif
!!$  open(1001,file=evolname,status='replace')
!!$  close(1001)

  allocate(Dm_local(ntri,nsorts), Dp_local(ntri,nsorts), Vp_local(ntri,nsorts), Tp_local(ntri,nsorts), ephi_local(ntri), &
           Vt2_local(ntri,nsorts) )
  allocate(D_Dt(n_therm,nsorts), D_Vt(n_therm,nsorts), D_Tt(n_therm,nsorts))
  allocate(D_Dm(ntri,nsorts), D_Dp(ntri,nsorts), D_Vp(ntri,nsorts), D_Tp(ntri,nsorts), D_ephi(ntri))
  allocate(D_j(nboucell,nsorts), D_tm(nboucell,nsorts), D_q(nboucell,nsorts))
  allocate(dens_phi(ntri,nsorts))

  Dm_local = 0.d0
  Dp_local = 0.d0
  Vp_local = 0.d0
  Tp_local = 0.d0
  Vt2_local = 0.d0
  ephi_local = 0.d0
!  dummy8 = 0.d0
!  dummy8_2 = 0.d0
  D_Dt = 0.d0
  D_Vt = 0.d0
  D_Tt = 0.d0
  D_Dm = 0.d0
  D_Dp = 0.d0
  D_Vp = 0.d0
  D_Tp = 0.d0
  D_j = 0.d0
  D_q = 0.d0
  D_tm = 0.d0
  D_ephi = 0.d0
  dens_phi = 0.d0

  macrosteps = nptrace/npgroup
! source particle thermal velocity, cm/s
  vt_s(:) = sqrt(2.d0*T_source*ev2erg/amass(:))
!
! Normamization factors:
! Larmor raidus corresponding to the field strength egual to $B_{00}$ harmonic
! in Boozer coordinates:
  bmod00 = bmod_ref
! Larmor radius:
!  rlarm=v0*4.d0*p_mass*c/(echarge*bmod_ref)
  do is=1,nsorts
     vref(is) = vt_s(is)
     rhoref(is) = vt_s(is)*amass(is)*clight/(echarge*charge(is)*bmod_ref)*bmod00
     elpref(is) = 2.d0*charge(is)/((amass(is))*vref(is)**2)*ev2erg
  enddo

  allocate( w_cell_V(n_therm,nsorts), w_cell_T(n_therm,nsorts) )
  w_cell_V = 0.d0
  w_cell_T = 0.d0

  allocate(sum_tr(ntri,nsorts), time_coll(ntri,nsorts), v_tr(ntri,nsorts), en_mean_tr(ntri,nsorts), t4Dm(ntri,nsorts),   &
           j_tr(nboucell,nsorts), q_tr(nboucell,nsorts), tm_tr(nboucell,nsorts), vt2_tr(ntri,nsorts))


  ! Read data file
  allocate(bpsi(ntri))
  call get_command_argument(1,temp)
  open(1,file=trim(temp))
  do i=1,ntri
     read (1,*) dummy_re,dummy_im
     bpsi(i)=cmplx(dummy_re,dummy_im)
  enddo
  close(1)
  
  ! Q profile generation
  delta_xo = 1.d0/dfloat(nr_core-1)
  
  do i=2,nr_core-1
     R = R_x + sqrt(1.d0-delta_xo*dfloat(i-1))*(R_o - R_x)
     Z = Z_x + sqrt(1.d0-delta_xo*dfloat(i-1))*(Z_o - Z_x)

     !R = R_o + delta_xo*dfloat(i-1)*(R_x-R_o)
     !Z = Z_o + delta_xo*dfloat(i-1)*(Z_x-Z_o)
     
     do ind_tri=1,ntri
        intri = check_in_tri(ind_tri, R, Z)
        if (intri) then
           !print *, 'found triangle: ', ind_tri
           exit
        end if
     end do

     isort = 1
     vpar = 1.d9
     vperp = 0.d0
     phi = 0.d0
     ind_tri_orig = ind_tri
     iedge_ent = 0

     ! first shot until edge
     dt = 1.d0
     call reg_step_fast(isort, ind_tri, iedge_ent, R, phi, Z, vperp, vpar, dt, ierr)
     do
        dt = 1.d0
        call reg_step_fast(isort, ind_tri, iedge_ent, R, phi, Z, vperp, vpar, dt, ierr)
        if(ind_tri == ind_tri_orig) then ! last shot until edge
           dt = 1.d0
           call reg_step_fast(isort, ind_tri, iedge_ent, R, phi, Z, vperp, vpar, dt, ierr)
           exit
        end if
     end do

     q = phi/(2d0*pi)

     isort = 1
     vpar = 1.d9
     vperp = 0.d0
     phi = 0.d0
     iedge_ent = 0
     
     R = R_x + sqrt(1.d0-delta_xo*dfloat(i-1))*(R_o - R_x)
     Z = Z_x + sqrt(1.d0-delta_xo*dfloat(i-1))*(Z_o - Z_x)
     
     call spectrum
     
     !print *, delta_xo*dfloat(i-1), q, real(bpsim), aimag(bpsim)
     ! calculate average poloidal flux
     dummy_re = 1.d0/3.d0*(mesh_point(mesh_element(ind_tri)%i_knot(1))%psi_pol + &
&      mesh_point(mesh_element(ind_tri)%i_knot(2))%psi_pol + &
&      mesh_point(mesh_element(ind_tri)%i_knot(3))%psi_pol)
     print *, dummy_re, q, real(bpsim), aimag(bpsim)
     
     phi = 0d0
  end do

contains

  subroutine spectrum_old()

    integer :: m, ind_tri_old
    double precision :: phi_old
    double precision :: theta, dtheta
    
    bpsim = (0d0,0d0)

    do ind_tri=1,ntri
        intri = check_in_tri(ind_tri, R, Z)
        if (intri) then
           exit
        end if
     end do

     theta = 0d0
     phi = 0d0
     phi_old = 0d0

     ! first shot until edge
     dt = 1.d0
     call reg_step_fast(isort, ind_tri, iedge_ent, R, phi, Z, vperp, vpar, dt, ierr)
     do
        phi_old = phi
        ind_tri_old = ind_tri
        dt = 1.d0
        call reg_step_fast(isort, ind_tri, iedge_ent, R, phi, Z, vperp, vpar, dt, ierr)
        theta = (phi+phi_old)/(2d0*q)
        dtheta = (phi-phi_old)/q

        do m=-mmax,0,mmax
           bpsim(m+mmax+1) = bpsim(m+mmax+1) + dtheta/(2d0*pi)*bpsi(ind_tri_old)*exp((0d0,-1d0)*m*theta)
           !bpsim(m+1) = bpsim(m+1) + dtheta/(2d0*pi)*1*exp((0d0,-1d0)*m*theta)
        end do

        if(ind_tri_old == ind_tri_orig) exit ! last shot
     end do

   end subroutine spectrum_old

   subroutine spectrum()
    integer :: m, ntheta, kt, ind_tri_old
    double precision :: theta, dtheta
    
    bpsim = (0d0,0d0)
    ntheta = 1024*mmax
    dtheta = 2d0*pi/ntheta

    do ind_tri=1,ntri
        intri = check_in_tri(ind_tri, R, Z)
        if (intri) then
           exit
        end if
     end do

     theta = 0d0
     phi = 0d0

     do kt=0,ntheta-1
        theta = (kt+5d-1)*dtheta
        do
           if(phi/q > theta) exit
           dt = 1.d0
           ind_tri_old = ind_tri
           call reg_step_fast(isort, ind_tri, iedge_ent, R, phi, Z, vperp, vpar, dt, ierr)
        end do
        do m=-mmax,mmax
           bpsim(mmax+m+1) = bpsim(mmax+m+1) + dtheta/(2d0*pi)*bpsi(ind_tri_old)*exp((0d0,-1d0)*m*theta)
           !bpsim(m+1) = bpsim(m+1) + dtheta/(2d0*pi)*1d0*exp((0d0,-1d0)*m*theta)
        end do
     end do
   end subroutine spectrum
   
  
end program result_spectrum
