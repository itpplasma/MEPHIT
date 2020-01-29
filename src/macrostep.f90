subroutine macrostep(cg_rng, npgroup, iptrace)
  use from_nrtype
  use accuracy_mod,  only : tinydp
  use constants,     only : nsorts, one3rd, erg2ev, amass, isort1, isort2, &
                            echarge, charge, pi
  use parmesh_mod,   only : w_cell_V, w_cell_T, D_therm, V_therm, T_therm, vol_therm,         &
                            time_coll, sum_tr, v_tr, en_mean_tr, n_therm, t4Dm,               &
                            j_tr, q_tr, tm_tr, vt2_tr
  use for_macrostep, only : dt0, nptrace, vt_s, source0, relax_tau, wpartsrr,                 & 
                            wcv, wct, D_perp, t_min, d_min, fulltime, perp_jump, par_jump,    &
                            wpart_rmp,dt0_factor,do_averages,wdamp_rate
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
  double precision :: dummy1, dummy2, ePhi_knot 
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

do_averages=.true.


!dtstep_par =  1.6d2/vt_s(isort)*par_jump  ! attention

  if(firststep) then
     firststep = .false.
     allocate(wsrr_rec(ntri), rel_therm_max(ntri))
     rel_therm_max(:) = 0.d0
  endif
!  wpart = sum(mesh_element(:)%V_tri)*source0*time_local(isort2)
  wpart = sum(mesh_element(:)%V_tri*mesh_element(:)%D_part(1)) 
!  print *,'total number of particles = ',wpart
!  wpart = wpart*echarge/(2.d0*pi*dt0*npgroup)
  wpart = wpart/(2.d0*pi*dt0*npgroup)

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
  species: do isort=1, nsorts 
!print *,'sort =',isort
!
  dt0_loc=dt0*dt0_factor(isort)
wdamp_rate=1.d0/dt0_loc
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
  do ipart=1,npgroup

!BEGIN ZOOM THE BAND
do     
     call born_particle(isort, R, phi, Z, vpar, vperp, ind_tri, cg_rng)
if( maxval(mesh_point( mesh_element(ind_tri)%i_knot(:) )%psi_pol) .gt. 8.35d6 .and. &
minval(mesh_point( mesh_element(ind_tri)%i_knot(:) )%psi_pol) .lt. 1.07d7) exit
enddo
!END ZOOM THE BAND
     wpart_rmp=(0.d0,0.d0)
     wsrr = mesh_element_rmp(ind_tri)%wsrr_rec_fix

     part_data(ipart)%R_p = R
     part_data(ipart)%phi_p = phi
     part_data(ipart)%Z_p = Z
     part_data(ipart)%v_ll = vpar
     part_data(ipart)%v_t = vperp
     part_data(ipart)%w_p = wpart/dt0_factor(isort)
     part_data(ipart)%w_srr = wsrr
     part_data(ipart)%i_tri = ind_tri
     part_data(ipart)%w_rmp = wpart_rmp

!istochnik = istochnik + wsrr*wpart

  enddo
!
  ipart = npgroup
!  ipart_locmax = npgroup
  dtstep_perp = (mesh_element(100)%sizmaxtri*perp_jump)**2/(2.d0*D_perp) 

! -----------------------------------------------------------------------------------------------
  do while(ipart.gt.0)
     numstep = 0
 !if(ipart/100*100.eq.ipart) print *,'isort = ',isort,'  ipart= ',ipart
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
     time_step(isort)=0.d0
     iedge_ent=0

     micro_i: do while (time_step(isort) .le. dt0_loc) 
!print *,time_step(isort),ipart
!      
        numstep = numstep + 1

!if((numstep/1000000)*1000000 .eq. numstep) then
!print *,'time= ',time_step(isort)
!write(111,*) R, Z
!endif

        ind_therm = mesh_element(ind_tri)%i_therm

!!$v3_sum_in(ind_therm) = v3_sum_in(ind_therm) + vpar
!!$part_num(ind_therm) = part_num(ind_therm) + 1.d0

         dens_cl=mesh_element(ind_tri)%D_part
         vel_cl=mesh_element(ind_tri)%V_part
         temp_cl=mesh_element(ind_tri)%T_part
!        do i=1,nsorts
! enough for the estimation of step length: 
!begin comment
!           dens_cl(i) = D_therm(ind_therm, i)  ! no interpolation of collis. freq.
!           vel_cl(i)  = V_therm(ind_therm, i)
!           temp_cl(i) = max(T_therm(ind_therm, i), t_min)
!end comment
!begin comment
!!$           dummy_verts(:) = mesh_point( mesh_element(ind_tri)%i_knot(:) )%pl_parms_knot(1,i)
!!$           dens_cl(i) =  max(cell_linint(ind_tri, R, Z, dummy_verts), d_min)
!!$           dummy_verts(:) = mesh_point( mesh_element(ind_tri)%i_knot(:) )%pl_parms_knot(2,i)
!!$           vel_cl(i)  = cell_linint(ind_tri, R, Z, dummy_verts)
!!$           dummy_verts(:) = mesh_point( mesh_element(ind_tri)%i_knot(:) )%pl_parms_knot(3,i)
!!$           temp_cl(i) = max(cell_linint(ind_tri, R, Z, dummy_verts), t_min)
!end uncomment
!        enddo
if(dens_cl(2).eq.0.d0) then
print *,R,Z
stop
endif
        call collis_time(dens_cl, temp_cl, vel_cl, vpar, vperp, isort, dt_coll, col_freq, col_freq_e)

dtstep_par =  R/sqrt(vpar**2 + vperp**2)*par_jump 

        dt_rest = min( dt_coll, dtstep_perp, dtstep_par )

!print *,dt_coll, dtstep_perp, dtstep_par, dt_rest
!pause

        dtstep_sav = dt_rest

        do i=1,nsorts
           w_c_V(i) = w_cell_V(ind_therm,i)
           w_c_T(i) = w_cell_T(ind_therm,i)
        enddo
        wpartsrr = wpart*wsrr   !<- SRR
        deltat = 0.d0
        deltav = 0.d0
        
!              if((wpartsrr/w_c_T(isort2)) .gt. rel_therm_max(ind_tri)) rel_therm_max(ind_tri) = wpartsrr/w_c_T(isort2)

        call collis_step(vel_cl, vpar, vperp, wpartsrr, w_c_V, w_c_T, isort, dtstep_sav, cg_rng, deltat, deltav)   !<- SRR

        do i=1,nsorts
           V_therm(ind_therm,i) = V_therm(ind_therm,i) + deltav(i)
           T_therm(ind_therm,i) = T_therm(ind_therm,i) + deltat(i)
        enddo

!v3_sum_out(ind_therm) = v3_sum_out(ind_therm) + vpar

        parallel: do while(dt_rest .gt. tinydp) 

           R_sav = R
           phi_sav = phi
           Z_sav = Z
           vperp_sav = vperp 
           vpar_sav = vpar
           itri_sav = ind_tri
           ierr = 0
           dtstep = dt_rest

laststeps=laststeps+1
           call reg_step_fast(isort, ind_tri, iedge_ent, R, phi, Z, vperp, vpar, dtstep, ierr)

           select case (ierr)

           case (0)

                    call account_trajectory(itri_sav, isort)
                 case (1)
                    print *,'0 or negative dt', mype
                    print *, itri_sav, R, phi, Z, vperp, vpar, dtstep
                    print *, mesh_point( mesh_element(itri_sav)%i_knot(1:3) )%rcoord
                    print *, mesh_point( mesh_element(itri_sav)%i_knot(:) )%zcoord
                    stop
                 case (2)
                    print *,'enormous parallel jump', mype
                    print *, vperp, vpar, dtstep, time_step
                    print *, itri_sav, R_sav, Z_sav
                    print *, ind_tri, R, Z
                    do i = 1,4
                       j = modulo(i,3) + 1
                       write(222,*) mesh_point( mesh_element(itri_sav)%i_knot(j) )%rcoord, &
                            mesh_point( mesh_element(itri_sav)%i_knot(j) )%zcoord
                    enddo
                    write(223,*) R_sav, Z_sav
                    write(223,*) R, Z
                    close(222)
                    close(223)
                    stop
                 case (-1)
                    call account_trajectory(itri_sav, isort)
!stenka = stenka + wsrr*wpart
              exit micro_i ! to load a new particle
!!$           case (-2)
!!$              call account_trajectory(itri_sav, isort)
!!$              reissue = .true.
!!$              call born_particle(vt_s(isort), R, phi, Z, vpar, vperp, ind_tri, cg_rng, reissue)
           case (-3)
              print *,'particle got lost, || ', mype
              call account_trajectory(itri_sav, isort)
!proebal = proebal + wsrr*wpart
              exit micro_i ! to load a new particle
           end select

           dt_rest = dt_rest - dtstep

        enddo parallel
!write (1001,*) time_step(isort),R,phi,Z,vperp,vpar,real(wpart_rmp),dimag(wpart_rmp)

        ierr = 0
!        call perp_step(dtstep_sav, R, phi, Z, vperp, vpar, ind_tri, ierr, cg_rng)

        if(ierr .eq. -1) then
!stenka  = stenka + wsrr*wpart   
           exit micro_i ! new particle        
        elseif(ierr .eq. -3) then
!proebal  = proebal + wsrr*wpart   
           exit micro_i ! new particle  
        endif

! Begin SRR
        fsplit = wsrr/wsrr_rec(ind_tri)
!
        if(fsplit .gt. fsplit_max) then
! splitting
           nsplit = nint(fsplit)
           nsplit = min(nsplit,ipart_max-ipart)  !< - safety
           if(nsplit .le. 1) cycle
           wsrr = wsrr/dfloat(nsplit)
           nsplit = nsplit - 1
           do while(nsplit .gt. 0)
              part_data(ipart)%R_p   = R
              part_data(ipart)%phi_p = phi
              part_data(ipart)%Z_p   = Z
              part_data(ipart)%v_ll  = vpar
              part_data(ipart)%v_t   = vperp
              part_data(ipart)%w_p   = wpart
              part_data(ipart)%w_srr = wsrr
              part_data(ipart)%i_tri = ind_tri
              part_data(ipart)%w_rmp = wpart_rmp
              nsplit = nsplit - 1
              ipart = ipart + 1
           enddo
        elseif(fsplit .lt. fsplit_min) then
! roulette
           xii = ran_u_01(cg_rng)
           if(xii .lt. fsplit) then
              wsrr = wsrr_rec(ind_tri)
           else
              exit micro_i
           endif
        endif
! End SRR
     enddo micro_i
     ipart = ipart-1
  enddo
!print *, mype, R, phi, Z, vpar, vperp, wpart_rmp
! -----------------------------------------------------------------------------------------------
  time_step(isort) = time_step(isort)/dfloat(npgroup)
  time_coll(:,isort) = time_coll(:,isort)/dfloat(npgroup)
  sum_tr(:,isort) = sum_tr(:,isort)/dfloat(npgroup)
  v_tr(:,isort) = v_tr(:,isort)/dfloat(npgroup)
  en_mean_tr(:,isort) = en_mean_tr(:,isort)/dfloat(npgroup)
  vt2_tr(:,isort) = vt2_tr(:,isort)/dfloat(npgroup)
  j_tr(:,isort) = j_tr(:,isort)/dfloat(npgroup)
  q_tr(:,isort) = q_tr(:,isort)/dfloat(npgroup)
  tm_tr(:,isort) = tm_tr(:,isort)/dfloat(npgroup)

!!$v3_sum_in = v3_sum_in/dfloat(npgroup)  
!!$v3_sum_out = v3_sum_out/dfloat(npgroup)
!!$part_num = part_num/dfloat(npgroup)


!print *,'j_tr=', sum(j_tr)
!
if(.false.) then
!  if(iptrace .eq. 1) then
     time_local(isort) = time_step(isort)
!  else
     time_l_new(isort) = relax_tau*time_local(isort) + (1.d0 - relax_tau)*time_step(isort)

     Dm_local(:,isort) = relax_tau*time_local(isort)/time_l_new(isort)*Dm_local(:,isort) +                       &
          (1.d0 - relax_tau)*time_step(isort)*t4Dm(:,isort)/time_l_new(isort)/(sum(time_coll(:,isort))*mesh_element(:)%V_tri)

     Dp_local(:,isort) = relax_tau*time_local(isort)/time_l_new(isort)*Dp_local(:,isort) +                       &
          (1.d0 - relax_tau)*time_step(isort)*sum_tr(:,isort)/time_l_new(isort)/(sum(time_coll(:,isort))*mesh_element(:)%V_tri)
!     where(Dp_local(:,isort) .ge. d_min)
     where(Dp_local(:,isort) .gt. 0.d0)
        Vp_local(:,isort) =  relax_tau*time_local(isort)/time_l_new(isort)*Vp_local(:,isort) +                   &
             (1.d0 - relax_tau)*v_tr(:,isort)/time_l_new(isort)/mesh_element(:)%V_tri/Dp_local(:,isort)

        Tp_local(:,isort) =  relax_tau*time_local(isort)/time_l_new(isort)*Tp_local(:,isort) +                   & 
             (1.d0 - relax_tau)*en_mean_tr(:,isort)*erg2ev*amass(isort)/time_l_new(isort)                        & 
             / mesh_element(:)%V_tri / Dp_local(:,isort)*one3rd 
        Vt2_local(:,isort) =  relax_tau*time_local(isort)/time_l_new(isort)*Vt2_local(:,isort) +                 & 
             (1.d0 - relax_tau)*vt2_tr(:,isort)*erg2ev*amass(isort)/time_l_new(isort)                            & 
             / mesh_element(:)%V_tri / Dp_local(:,isort)*5.d-1 

     elsewhere
        Dp_local(:,isort) = d_min
        Vp_local(:,isort) = 0.d0
        Tp_local(:,isort) = t_min
        Vt2_local(:,isort) = 0.d0
     end where
     flux_j(:,isort) =  relax_tau*time_local(isort)/time_l_new(isort)*flux_j(:,isort) +    &
          (1.d0 - relax_tau)*j_tr(:,isort)/time_l_new(isort)

     flux_q(:,isort) =  relax_tau*time_local(isort)/time_l_new(isort)*flux_q(:,isort) +    &
          (1.d0 - relax_tau)*q_tr(:,isort)/time_l_new(isort)

     tor_mom(:,isort) =  relax_tau*time_local(isort)/time_l_new(isort)*tor_mom(:,isort) +  &
          (1.d0 - relax_tau)*tm_tr(:,isort)/time_l_new(isort)

     time_local(isort) = time_l_new(isort)
endif
  fulltime = fulltime + time_step(isort)*dfloat(npgroup)
!
  deallocate(part_data)
!
205 FORMAT(1000(ES15.8E2,1X))
206 FORMAT(7(ES21.14E2,2X),i4,2x,i4)

  if(.not. allocated(dummy4Dm) ) allocate (dummy4Dm(n_therm))
if(.false.) then   !========>this brach is turned off<============
  dummy4Dm = 0.0d0
  D_therm = 0.0d0
  do i=1,ntri
     ind_therm = mesh_element(i)%i_therm
     D_therm(ind_therm,isort) = D_therm(ind_therm,isort) +       &
          Dp_local(i,isort)*mesh_element(i)%V_tri/vol_therm(ind_therm)
     dummy4Dm(ind_therm) = dummy4Dm(ind_therm) +                 &
          Dm_local(i,isort)*mesh_element(i)%V_tri/vol_therm(ind_therm)
  enddo

  Dp_local(:,isort) = D_therm(mesh_element(:)%i_therm,isort)
  Dm_local(:,isort) = dummy4Dm(mesh_element(:)%i_therm)
endif
  enddo species
  return
!----------------------------------------------------------------------------------------
contains

  subroutine account_trajectory(i_tri, isort)
    integer :: i_tri, isort
    double precision :: u, w
    time_step(isort) = time_step(isort) + dtstep*wsrr   !<- SRR
    t4Dm(i_tri,isort) = t4Dm(i_tri,isort) + dtstep
    time_coll(i_tri,isort) = time_coll(i_tri,isort) + dtstep*wsrr   !<- SRR
    sum_tr(i_tri,isort) = sum_tr(i_tri,isort) + wpart*wsrr*dtstep   !<- SRR
    u = (vpar + vpar_sav)*0.5d0
    w = (vperp + vperp_sav)*0.5d0
    v_tr(i_tri,isort) = v_tr(i_tri,isort) + dtstep*u*wpart*wsrr   !<- SRR
    en_mean_tr(i_tri,isort) = en_mean_tr(i_tri,isort) + dtstep*(u**2 + w**2)*wpart*wsrr    !<- SRR
    vt2_tr(i_tri,isort) = vt2_tr(i_tri,isort) + dtstep*w**2*wpart*wsrr    !<- SRR

!    dens_phi(i_tri,isort) = dens_phi(i_tri,isort) + dtstep*wsrr*wtr

  end subroutine account_trajectory

end subroutine macrostep
