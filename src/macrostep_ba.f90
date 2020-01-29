subroutine macrostep_ba(cg_rng, npgroup, iptrace)
  use from_nrtype
  use accuracy_mod,  only : tinydp
  use constants,     only : nsorts, one3rd, erg2ev, amass, isort1, isort2, &
                            echarge, charge, pi, ev2erg
  use parmesh_mod,   only : w_cell_V, w_cell_T, D_therm, V_therm, T_therm, vol_therm,         &
                            time_coll, sum_tr, v_tr, en_mean_tr, n_therm, t4Dm,               &
                            j_tr, q_tr, tm_tr, vt2_tr
  use for_macrostep, only : dt0, nptrace, vt_s, source0, relax_tau, wpartsrr,                 & 
                            wcv, wct, D_perp, t_min, d_min, fulltime, perp_jump, par_jump,    &
                            wpart_rmp,dt0_factor,do_averages,wdamp_rate,n_tormode
  use mesh_mod,      only : ntri, ntri_inbou, mesh_point, inbou_list,                         &
                            cell_linint, npoint, grad_PhiovB2, mesh_element,                  &
                            mesh_element_rmp                                                        !<=SRR
  use for_mpi,       only : time_local, Dm_local, Dp_local, Vp_local, Tp_local, ephi_local,   &
                            mype, npes, mype_sym, evolname, flux_j, flux_q, tor_mom,          &
                            Vt2_local, dens_phi
  use hujnja,         only: numstep, numstepmax, firstpart, rel_therm_max, laststeps

  implicit none
!
! begin uncomment
  double precision, parameter :: fsplit_min=0.5d0, fsplit_max=2.d0
! end uncomment

!  double precision, parameter :: fsplit_min=0.25d0, fsplit_max=4.d0
!  double precision, parameter :: fsplit_max=3.d0, fsplit_min=1.d0/fsplit_max
!  double precision, parameter :: fsplit_max=29.d0, fsplit_min=1.d0/fsplit_max

! begin comment
!  double precision, parameter :: fsplit_max=3.d10, fsplit_min=1.d0/fsplit_max
! end comment

!  double precision, parameter :: fsplit_max=3.d0, fsplit_min=0.d0
!
  integer :: ierr, isort, ipart, i, j, npgroup, ind_tri, itri_sav, m, ind_therm, is
  integer :: ipart_max, nsplit, iptrace !, ipart_locmax
  double precision, dimension(nsorts) :: dens_cl, temp_cl, vel_cl, w_c_V, w_c_T, deltav, deltat
  double precision :: R, phi, Z, vpar, vperp, wpart, xii, ran_u_01
  double precision :: R_sav, phi_sav, Z_sav, vperp_sav, vpar_sav, dens_tp_avr, temp_tp_avr
  double precision :: dtstep, dt_rest, dt_coll, dtstep_sav, col_freq, col_freq_e
  double precision :: dtstep_perp,  dtstep_par, tau_l_1st
  double precision :: wsrr_rec_min, wsrr, fsplit, wsrr_min=1.d-5, dummy_vts = 0.d0
  double precision :: dt0_loc

  integer :: look_4_tri,iedge_ent
  double precision, dimension(:), allocatable :: wsrr_rec, dummy4Dm, dummy4ww
  double precision :: dummy1, dummy2, ePhi_knot , vt
  double precision, dimension(nsorts) :: T_knot, V_knot, D_knot, time_step, time_l_new
  double precision, dimension(3) :: dummy_verts
  double precision :: cg_rng(0:5)
  logical :: firststep=.true.
  type :: particle
     real(dp) :: R_p 
     real(dp) :: phi_p
     real(dp) :: Z_p
     real(dp) :: v_ll
     real(dp) :: v_t
     real(dp) :: w_p
     real(dp) :: w_srr
     integer(i4b) :: i_tri
     double complex :: w_rmp
  end type particle
  type(particle), dimension(:), allocatable :: part_data
  logical :: check_in_tri
  save wsrr_rec, firststep,            R, phi, Z, vpar, vperp, ind_tri, wpart, wsrr
  external check_in_tri, look_4_tri
double precision ::  gauss_openmp
external  gauss_openmp
!
  integer :: nstartp,npitch,k,nx,ix
  double precision :: hpitch,vmod,smt_jacobian,hz,xmax,x,f0norm
  integer,          dimension(:), allocatable :: itribeg
  double precision, dimension(:), allocatable :: Rbeg,Zbeg
!
  nstartp=0
  open(314,file='startpoints.dat')
  do
    read (314,*,end=1) R
    nstartp=nstartp+1
  enddo
1 close(314)
  allocate(Rbeg(nstartp),Zbeg(nstartp),itribeg(nstartp))
  open(314,file='startpoints.dat')
  do i=1,nstartp
    read (314,*) Rbeg(i),Zbeg(i),itribeg(i)
  enddo
  close(314)
!
  nx=10
  xmax=3.d0
  npitch=200 !100
  hpitch=1.d0/dfloat(npitch)
  npgroup=2*npitch*nstartp*nx


  if(firststep) then
     firststep = .false.
     allocate(wsrr_rec(ntri), rel_therm_max(ntri))
     rel_therm_max(:) = 0.d0
  endif
  wpart = sum(mesh_element(:)%V_tri*mesh_element(:)%D_part(1)) 
!  print *,'total number of particles = ',wpart
!  wpart = wpart*echarge/(2.d0*pi*dt0*npgroup)
!<=SMT  wpart = wpart/(2.d0*pi*npgroup)
   wpart = 2.d0*(Rbeg(nstartp)-Rbeg(1))*xmax/(sqrt(pi)*npgroup)  !<=SMT

  do is=1,nsorts
     w_cell_V(:,is) = wpart*wcv(is)   
     w_cell_T(:,is) = wpart*wct(is)
  enddo
!

  time_step(:) = 0.d0
  time_coll(:,:) = 0.d0
  sum_tr(:,:) = 0.d0
  v_tr(:,:) = 0.d0
  en_mean_tr(:,:) = 0.d0
  vt2_tr(:,:) = 0.d0
  t4Dm(:,:) = 0.d0
  j_tr(:,:) = 0.d0
  q_tr(:,:) = 0.d0
  tm_tr(:,:) = 0.d0
!
!
!  isort = isort1 ! attention 1 - electrons, 2- ions
!
laststeps=0
isort = 1
!!!ENABLE AGAIN>>  species: do isort=1, nsorts 
!
    dt0_loc=dt0*dt0_factor(isort)
    wdamp_rate=1.d0/dt0_loc         !<= Krook model
!    wdamp_rate=0.1d0/dt0_loc         !<= Krook model
!    wdamp_rate=10.d0/dt0_loc         !<= Krook model
!
!  wsrr_rec(:) = mesh_element(:)%D_part(isort)/maxval(mesh_element(:)%D_part(isort)) !<=SRR
  wsrr_rec(:) = mesh_element_rmp(:)%wsrr_rec_fix                                     !<=SRR
  
  allocate( dummy4ww(n_therm) )
  do i=1,ntri
     ind_therm = mesh_element(i)%i_therm
     dummy4ww(ind_therm) = wsrr_rec(i)
  enddo
  do i=1,n_therm
     w_cell_V(i,:) =  w_cell_V(i,:)*dummy4ww(i)
     w_cell_T(i,:) =  w_cell_T(i,:)*dummy4ww(i)
  enddo
  deallocate( dummy4ww )

  wsrr_rec_min = 1.d0
  wsrr_rec_min=max(wsrr_rec_min,wsrr_min)
  ipart_max = npgroup*int(2.d0/wsrr_rec_min)

  allocate(part_data(ipart_max))
!
  ipart=0
  do i=1,nstartp
    R = Rbeg(i)
    Z = Zbeg(i)
    ind_tri = itribeg(i)
    vt = sqrt(2.d0*T_therm(mesh_element(ind_tri)%i_therm,isort)*ev2erg/amass(isort))
    hz=mesh_element(ind_tri)%dPsi_dR*3.d0                     &  !<=SMT
      /sum(mesh_point(mesh_element(ind_tri)%i_knot(:))%b_mod)    !<=SMT
    do ix=1,nx
      x=xmax*(dfloat(ix)-0.5)/dfloat(nx)
      f0norm=2.d0*x**3*exp(-x**2)
      vmod=vt*x
      do k=-1,1,2
        do j=1,npitch
          if(k.eq.1) then
            vpar=(dfloat(j)-0.5d0)*hpitch*vmod
          else
            vpar=(dfloat(j-npitch)-0.5d0)*hpitch*vmod
          endif
          vperp=sqrt(vmod**2-vpar**2)
          phi=0.d0                    !<= start from $\varphi=0$
          wpart_rmp=(0.d0,0.d0)
          wsrr = mesh_element_rmp(ind_tri)%wsrr_rec_fix

          ipart=ipart+1
          part_data(ipart)%R_p = R
          part_data(ipart)%phi_p = phi
          part_data(ipart)%Z_p = Z
          part_data(ipart)%v_ll = vpar
          part_data(ipart)%v_t = vperp
!         part_data(ipart)%w_p = wpart/dt0_factor(isort)
!<=SMT        part_data(ipart)%w_p = wpart
          smt_jacobian=abs(vpar*hz)/vmod                               !<=SMT
          part_data(ipart)%w_p = wpart*smt_jacobian*vt*f0norm       &  !<=SMT
                               * mesh_element(ind_tri)%D_part(1)       !<=SMT
!<=SMT        part_data(ipart)%w_srr = wsrr
          part_data(ipart)%w_srr = 1.d0                                !<=SMT
          part_data(ipart)%i_tri = ind_tri
          part_data(ipart)%w_rmp = wpart_rmp
        enddo
      enddo
    enddo
  enddo
!
  print *,ipart,npgroup
!  ipart = npgroup
!
! -----------------------------------------------------------------------------------------------
  particles: do while(ipart.gt.0)
     if(ipart/1000*1000.eq.ipart) print *,'isort = ',isort,'  ipart= ',ipart
     R       = part_data(ipart)%R_p
     phi     = part_data(ipart)%phi_p
     Z       = part_data(ipart)%Z_p
     vpar    = part_data(ipart)%v_ll
     vperp   = part_data(ipart)%v_t
     wpart   = part_data(ipart)%w_p
     wsrr    = part_data(ipart)%w_srr
     ind_tri = part_data(ipart)%i_tri
     wpart_rmp = part_data(ipart)%w_rmp
!
     itri_sav = ind_tri
     ierr = 0
!
     iedge_ent=0
     do_averages=.false.
!      
! leave the starting triangle:      
     numstep=0
     do while (ind_tri.eq.itri_sav)
        dtstep =  R/sqrt(vpar**2 + vperp**2)*par_jump 
        call reg_step_fast(isort, ind_tri, iedge_ent, R, phi, Z, vperp, vpar, dtstep, ierr)
        if(ierr.eq.-1) then
! exit to the wall, skip this orbit
          ipart = ipart-1
          cycle particles
        elseif(ierr.ne.0) then
          print *,'loop 1/1 ierr = ',ierr
        endif
        numstep=numstep+1
        if(numstep.gt.100) then
! elliptic orbit within the triangle, skip it
          print *,'elliptic'
          ipart = ipart-1
          cycle particles
        endif
     enddo
!
! starting velocity at the triangle edge
     vpar_sav = vpar
!
! First integration: search for bounce time, toroidal shift and initial weight:
!
     dt_rest=0.d0  !<=bounce time
     wpart_rmp = (0.d0,0.d0)
     phi=0.d0
     numstep=0
!
     do
!
! re-enter the starting triangle
       do while (ind_tri.ne.itri_sav)
          dtstep =  R/sqrt(vpar**2 + vperp**2)*par_jump    
          call reg_step_fast(isort, ind_tri, iedge_ent, R, phi, Z, vperp, vpar, dtstep, ierr)
          dt_rest=dt_rest+dtstep
          if(ierr.eq.-1) then
! exit to the wall, skip this orbit
            ipart = ipart-1
            cycle particles
          elseif(ierr.ne.0) then
            print *,'loop 1/2 ierr = ',ierr
          endif
       enddo
!
! leave again the starting triangle:      
       do while (ind_tri.eq.itri_sav)
          dtstep =  R/sqrt(vpar**2 + vperp**2)*par_jump 
          call reg_step_fast(isort, ind_tri, iedge_ent, R, phi, Z, vperp, vpar, dtstep, ierr)
          dt_rest=dt_rest+dtstep
          if(ierr.eq.-1) then
! exit to the wall, skip this orbit
            ipart = ipart-1
            cycle particles
          elseif(ierr.ne.0) then
            print *,'loop 1/3 ierr = ',ierr
          endif
       enddo
       if(abs(vpar-vpar_sav).lt.abs(vpar)*1.d-5) exit
       numstep=numstep+1
       if(numstep.gt.1) then
! Inaccurate orbit, skip it
         print *,'low accuracy'
         ipart = ipart-1
         cycle particles
       endif
!
     enddo
!if(numstep.eq.0) then
!
! Second integration: contribution to averages:
!
     do_averages=.true.
write(2222,*) R,sqrt(vperp**2+vpar**2),vpar/sqrt(vperp**2+vpar**2),n_tormode*phi, dt_rest, &
              real(wpart_rmp),dimag(wpart_rmp),wdamp_rate*dt_rest
     wpart_rmp = wpart_rmp*exp(wdamp_rate*dt_rest) &
               /(exp(wdamp_rate*dt_rest+cmplx(0.d0,n_tormode*phi))-(1.d0,0.d0))
!write(2223,*) vpar/sqrt(vperp**2+vpar**2),real(wpart_rmp),dimag(wpart_rmp)
!<=SMT     wpartsrr = wpart*wsrr/dt_rest   !<- SRR
     wpartsrr = wpart*wsrr/dfloat(1+numstep)   !<- SRR <=SMT: numstep=0 for passing, 1 for trapped
     phi=0.d0
!
     dt_rest = dt_rest*(1.d0+1.d-6)   !<=fix - cross the last boundary (do not stop at its vicinity)
!
     do while(dt_rest .gt. tinydp)
       dtstep =  R/sqrt(vpar**2 + vperp**2)*par_jump 
       dtstep = min(dt_rest,dtstep)
       call reg_step_fast(isort, ind_tri, iedge_ent, R, phi, Z, vperp, vpar, dtstep, ierr) 
       dt_rest = dt_rest - dtstep
       if(ierr.ne.0) exit
     enddo
!endif
!
     ipart = ipart-1
  enddo particles
!print *, mype, R, phi, Z, vpar, vperp, wpart_rmp
! -----------------------------------------------------------------------------------------------
!
  deallocate(part_data)
!
205 FORMAT(1000(ES15.8E2,1X))
206 FORMAT(7(ES21.14E2,2X),i4,2x,i4)

  if(.not. allocated(dummy4Dm) ) allocate (dummy4Dm(n_therm))
!  enddo species
  return
end subroutine macrostep_ba
