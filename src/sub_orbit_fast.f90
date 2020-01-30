module esche_hujnja
  double precision, dimension(4) :: varphi_out
  integer :: line_type
end module esche_hujnja
!
!
!
subroutine reg_step_fast(isort, ind_tri, iedge_ent, R, phi, Z, vperp, vpar, dt, ierr)
  use constants,         only : one3rd, charge, amass, echarge
  use accuracy_mod,      only : eps_newt, eps_cell
  use mesh_mod,          only : legs, mesh_point, mesh_element, cell_linint,        &      
                                mesh_element_rmp, bphicovar
  use reg_step_tria_mod, only : abr,abz,aphir,aphiz,apsir,apsiz,                    &
                                binv_vert1,phiovb2_vert1,Rcell,sizmax,              &
                                deltatri,R_vert,Z_vert
  use for_macrostep, only :  wpartsrr,wpart_rmp,n_tormode,do_averages,wdamp_rate
use esche_hujnja
  implicit none
  integer :: isort, ierr, ierr_tri, ind_tri, leg_exit, ind_tri_new, iedge_ent
  double precision :: dt, R, phi, Z, vperp, vpar, R_sav, Z_sav, R_new, Z_new        &
       , vperp_sav, vpar_sav

!  integer :: itri_div
  double precision :: tor_momentum, b_abs, b_tor, a_phi, t_m_part, rad
  double precision :: damp_prob
  double precision, dimension(legs) :: fun_tmp
  double precision, dimension(2) :: t_m_leg, r_leg, z_leg
  integer :: m, i, i_st, i_fn, np_st, np_fn, np_3rd
  integer, parameter :: lim_see_saw = 10
  double precision, dimension(lim_see_saw) :: R_write, Z_write, dist_write

  logical :: check_in_tri !, found_div
  integer :: look_4_tri, itri_sav, n_bou, leg_entry
  external check_in_tri, look_4_tri
!
  double complex :: cdummy,expon,countden,countpres_perp,countpres_par,countparcurr
!
  double precision :: bphicovar_interp, psival(3), psi, psi_h
  
  ierr = 0
  ierr_tri = 0
!
  abr = mesh_element(ind_tri)%dBinv_dR
  abz = mesh_element(ind_tri)%dBinv_dZ
  aphir = mesh_element(ind_tri)%dPhiovB2_dR
  aphiz = mesh_element(ind_tri)%dPhiovB2_dZ
  apsir = mesh_element(ind_tri)%dPsi_dR
  apsiz = mesh_element(ind_tri)%dPsi_dZ
  binv_vert1 = 1.0d0/mesh_point(mesh_element(ind_tri)%i_knot(1))%b_mod
  phiovb2_vert1 = mesh_point(mesh_element(ind_tri)%i_knot(1))%PhiovB2
  R_vert(:) = mesh_point( mesh_element(ind_tri)%i_knot(:) )%rcoord
  Z_vert(:) = mesh_point( mesh_element(ind_tri)%i_knot(:) )%zcoord
  deltatri = mesh_element(ind_tri)%det_3
  sizmax = mesh_element(ind_tri)%sizmaxtri
  Rcell = sum(R_vert)*one3rd
  R_sav = R
  Z_sav = Z
  itri_sav = ind_tri
  vpar_sav = vpar
  vperp_sav = vperp
!
  expon=exp(cmplx(0.d0,n_tormode*phi))
  cdummy=-vpar*expon*(mesh_element_rmp(ind_tri)%bnorm_times_thermforces(1,isort)        &
        +(vperp**2+vpar**2)*mesh_element_rmp(ind_tri)%bnorm_times_thermforces(2,isort))
  if(do_averages) then
    countden=wpart_rmp*conjg(expon)/R
    countpres_perp=wpart_rmp*conjg(expon)*vperp**2/R
    countpres_par=wpart_rmp*conjg(expon)*vpar**2/R
    countparcurr=wpart_rmp*conjg(expon)*vpar/R
  endif
!
  psival(:) = mesh_point(mesh_element(ind_tri)%i_knot(:))%psi_pol
  psi = cell_linint(ind_tri, R, Z, psival)
  psi_h = psival(mesh_element(ind_tri)%knot_h)
  bphicovar_interp = bphicovar + mesh_element(ind_tri)%dbphicovdpsi*(psi-psi_h)
       
  call reg_step_tria(isort, leg_exit, iedge_ent, R, phi, Z, vperp, vpar, &
       bphicovar_interp, dt, ierr)
!
  damp_prob=0.5d0*wdamp_rate*dt            !damping of the weight (wdamp_rate=$\nu_0$)
!  wpart_rmp=wpart_rmp+0.5d0*dt*cdummy
  wpart_rmp=wpart_rmp*(1.d0-damp_prob)+0.5d0*dt*cdummy
  expon=exp(cmplx(0.d0,n_tormode*phi))
  cdummy=-vpar*expon*(mesh_element_rmp(ind_tri)%bnorm_times_thermforces(1,isort)        &
        +(vperp**2+vpar**2)*mesh_element_rmp(ind_tri)%bnorm_times_thermforces(2,isort))
!  wpart_rmp=wpart_rmp+0.5d0*dt*cdummy
  wpart_rmp=(wpart_rmp+0.5d0*dt*cdummy)/(1.d0+damp_prob)
  if(do_averages) then
    countden=countden+wpart_rmp*conjg(expon)/R
    countpres_perp=countpres_perp+wpart_rmp*conjg(expon)*vperp**2/R
    countpres_par=countpres_par+wpart_rmp*conjg(expon)*vpar**2/R
    countparcurr=countparcurr+wpart_rmp*conjg(expon)*vpar/R
!
    countden=countden*0.5d0*dt*wpartsrr*echarge*charge(isort)
    countpres_perp=countpres_perp*amass(isort)*0.25d0*dt*wpartsrr
    countpres_par=countpres_par*amass(isort)*0.5d0*dt*wpartsrr
    countparcurr=countparcurr*0.5d0*dt*wpartsrr*echarge*charge(isort)
    mesh_element_rmp(ind_tri)%denspert(isort)                             &
        = mesh_element_rmp(ind_tri)%denspert(isort) + countden
    mesh_element_rmp(ind_tri)%prespert_perp(isort)                        &
        = mesh_element_rmp(ind_tri)%prespert_perp(isort) + countpres_perp
    mesh_element_rmp(ind_tri)%prespert_par(isort)                         &
        = mesh_element_rmp(ind_tri)%prespert_par(isort) + countpres_par
    mesh_element_rmp(ind_tri)%parcurrpert(isort)                          &
        = mesh_element_rmp(ind_tri)%parcurrpert(isort) + countparcurr
  endif
!
  if(leg_exit .eq. 0) then
    iedge_ent=0
  else
    iedge_ent=mesh_element(ind_tri)%neighbour_edge(leg_exit)
  endif
!
  if(leg_exit .eq. 0) then ! step is over?
!-----------------------------------------------------------------
     if(.not. check_in_tri(ind_tri, R, Z) ) then ! not yet
ierr=-3
return
print *,'lost inside'
write(1002,*) ' '
write(1002,*) R,Z
write(1002,*) ' '
do i=1,3
write(1002,*) mesh_point(mesh_element(ind_tri)%i_knot(i) )%rcoord, &
              mesh_point(mesh_element(ind_tri)%i_knot(i) )%zcoord
enddo
i=1
write(1002,*) mesh_point(mesh_element(ind_tri)%i_knot(i) )%rcoord, &
              mesh_point(mesh_element(ind_tri)%i_knot(i) )%zcoord
stop
        ind_tri_new = look_4_tri(ind_tri, R, Z, ierr_tri) 
        if(ierr_tri .eq. -1) then
           print *, 'parallel, dt = ',dt,' particle got lost, leg exit=0!'
           print *, check_in_tri(itri_sav, R_sav, Z_sav)
           ierr = -3
           call write_tri
print *,line_type, varphi_out(:)
print *,R_sav, Z_sav, itri_sav, vpar_sav, vperp_sav
!           pause
        else
!!$           n_bou = mesh_element(ind_tri_new)%iqq_gc
!!$           if ( n_bou .ne. 0 ) &
!!$                call check_boundary(isort, ind_tri_new, n_bou, vpar, vperp, R, Z, ierr)
! if tr. does not belong to GC, just continue tracing
           ind_tri = ind_tri_new
           print *, 'parallel, particle got lost but found'
        endif
     endif
!-----------------------------------------------------------------
  else
!-----------------------------------------------------------------
     ind_tri_new = mesh_element(ind_tri)%neighbour(leg_exit)
     if(ind_tri_new .lt. 0) then
        ierr = -1
        return
     endif
!
     leg_entry = mesh_element(ind_tri)%neighbour_edge(leg_exit)
     if(do_averages) then
!       cdummy=charge(isort)*wpartsrr*wpart_rmp*conjg(expon)/R
       cdummy=echarge*charge(isort)*wpartsrr*wpart_rmp*conjg(expon)/R
       mesh_element_rmp(ind_tri)%currents(leg_exit,isort)                    &
          = mesh_element_rmp(ind_tri)%currents(leg_exit,isort) + cdummy
       mesh_element_rmp(ind_tri_new)%currents(leg_entry,isort)               &
          = mesh_element_rmp(ind_tri_new)%currents(leg_entry,isort) - cdummy
     endif
!
     ind_tri = ind_tri_new
!-----------------------------------------------------------------
  endif

205 FORMAT(1000(ES15.8E2,1X))
  return
!-----------------------------------------------------------------------------
contains

  subroutine write_tri
  integer :: i, j
    do i = 1,4
       j = modulo(i,3) + 1
       write(222,*) mesh_point( mesh_element(itri_sav)%i_knot(j) )%rcoord, &
            mesh_point( mesh_element(itri_sav)%i_knot(j) )%zcoord
    enddo
    close(222)
    write(223,*) R_sav, Z_sav, itri_sav
    write(223,*) R, Z, ind_tri
    close(223)
!    do i = 1,4
!       j = modulo(i,3) + 1
!       write(224,*) mesh_point( mesh_element(ind_tri)%i_knot(j) )%rcoord, &
!            mesh_point( mesh_element(ind_tri)%i_knot(j) )%zcoord
!    enddo
!    close(224)
!    pause
  end subroutine write_tri

end subroutine reg_step_fast
!=============================================================================
!!$subroutine check_boundary(isort, ind_tri, n_bou, vpar, vperp, R, Z, ierr)
!!$
!!$  use mesh_mod,      only : legs, mesh_point, mesh_element, cell_linint,  &
!!$                            bphicovar, i_pfz, i_dpr, i_sol, i_dpl, i_inb
!!$  use parmesh_mod,   only :  j_tr, q_tr, tm_tr
!!$  implicit none
!!$  integer :: isort, ind_tri, n_bou, ierr
!!$  double precision ::  dt, vpar, vperp, b_abs, R, Z
!!$  double precision, dimension(legs) :: fun_tmp
!!$
!!$  if(n_bou.ge.i_dpl(1) .and. n_bou.le.i_dpl(2)) then  ! left plate, absorbtion
!!$     j_tr(n_bou,isort) = j_tr(n_bou,isort) + wpartsrr
!!$     q_tr(n_bou,isort) = q_tr(n_bou,isort) + (vpar**2 + vperp**2)*wpartsrr
!!$     fun_tmp(:) = mesh_point(mesh_element(ind_tri)%i_knot(:))%b_mod
!!$     b_abs = cell_linint(ind_tri, R, Z, fun_tmp)
!!$     tm_tr(n_bou,isort) = tm_tr(n_bou,isort) + vpar*wpartsrr*bphicovar/b_abs
!!$     ierr = -1 
!!$  else if (n_bou.ge.i_dpr(1) .and. n_bou.le.i_dpr(2)) then  ! right plate, absorbtion
!!$     j_tr(n_bou,isort) = j_tr(n_bou,isort) + wpartsrr
!!$     q_tr(n_bou,isort) = q_tr(n_bou,isort) + (vpar**2 + vperp**2)*wpartsrr
!!$     fun_tmp(:) = mesh_point(mesh_element(ind_tri)%i_knot(:))%b_mod
!!$     b_abs = cell_linint(ind_tri, R, Z, fun_tmp)
!!$     tm_tr(n_bou,isort) = tm_tr(n_bou,isort) + vpar*wpartsrr*bphicovar/b_abs
!!$     ierr = -1
!!$  else if (n_bou.ge.i_sol(1) .and. n_bou.le.i_sol(2)) then ! outer boundary (SOL) crossed, absorbtion
!!$     j_tr(n_bou,isort) = j_tr(n_bou,isort) + wpartsrr
!!$     q_tr(n_bou,isort) = q_tr(n_bou,isort) + (vpar**2 + vperp**2)*wpartsrr
!!$     fun_tmp(:) = mesh_point(mesh_element(ind_tri)%i_knot(:))%b_mod
!!$     b_abs = cell_linint(ind_tri, R, Z, fun_tmp)
!!$     tm_tr(n_bou,isort) = tm_tr(n_bou,isort) + vpar*wpartsrr*bphicovar/b_abs
!!$     ierr = -1 
!!$  else if (n_bou.ge.i_pfz(1) .and. n_bou.le.i_pfz(2)) then  ! outer boundary PFZ crossed, absorbtion
!!$     j_tr(n_bou,isort) = j_tr(n_bou,isort) + wpartsrr
!!$     q_tr(n_bou,isort) = q_tr(n_bou,isort) + (vpar**2 + vperp**2)*wpartsrr
!!$     fun_tmp(:) = mesh_point(mesh_element(ind_tri)%i_knot(:))%b_mod
!!$     b_abs = cell_linint(ind_tri, R, Z, fun_tmp)
!!$     tm_tr(n_bou,isort) = tm_tr(n_bou,isort) + vpar*wpartsrr*bphicovar/b_abs
!!$     ierr = -1 
!!$  else if (n_bou.ge.i_inb(1) .and. n_bou.le.i_inb(2)) then  ! inner boundary crossed, reissue
!!$     ierr = -2 
!!$  end if
!!$  return
!!$end subroutine check_boundary
!=============================================================================
!
subroutine reg_step_tria(isort, leg_exit, iedge_ent, R, phi, Z, vperp, vpar, &
     bphicovar_interp, dt, ierr)
!
  use constants,         only : pi
  use mesh_mod,          only : bphicovar
  use accuracy_mod,      only : eps_newt
  use renorm_mod,        only : vref, rhoref, elpref
  use reg_step_tria_mod, only : abr,abz,aphir,aphiz,apsir,apsiz,        &
                                binv_vert1,phiovb2_vert1,Rcell,sizmax,  &
                                deltatri,R_vert,Z_vert
  implicit none
!
  logical :: lowtri,ellips
!
  integer :: isort,ierr,leg_exit,ind_vert,iedge_ent
  integer, dimension(1) :: indb
  double precision :: dt,R,phi,Z,vperp,vpar
  double precision :: binv,Rnew,Znew
  double precision :: alpha,coefphi,aphir_pt,aphiz_pt,apsir_pt,apsiz_pt
  double precision :: binv0,phiovb20,abar_rr,abar_rz,abar_zz,wbar,ubar
  double precision :: bbar_r,bbar_z,Delta,R_a,Z_a,R_c,Z_c,R_s,Z_s
  double precision :: rate,R_bou,Z_bou,c_c,c_s,c,discr,x,den,phase,co,si
  double precision :: sipl,simn,copl,comn,shpl,shmn,chpl,chmn,phipl,phimn
  double precision :: acoef_r,acoef_z,bcoef_r,bcoef_z,acoef,bcoef,ccoef
  double precision :: expon,R_side,Z_side
  double precision :: coint,siint,phinew
  double precision :: expon1,expon2,sc
  double precision, dimension(4) :: varphi
  double precision :: bphicovar_interp
integer :: jjj, j
!
  ierr=0
!
  alpha=vref(isort)/rhoref(isort)/bphicovar
  coefphi=elpref(isort)*vref(isort)**2
!
  aphir_pt=aphir*coefphi
  aphiz_pt=aphiz*coefphi
  apsir_pt=apsir*alpha
  apsiz_pt=apsiz*alpha
!
  binv0=binv_vert1+abr*(R-R_vert(1))+abz*(Z-Z_vert(1))
  phiovb20=phiovb2_vert1*coefphi &
          +aphir_pt*(R-R_vert(1))+aphiz_pt*(Z-Z_vert(1))
!
  wbar=vperp**2+vpar**2+phiovb20/binv0**2
  ubar=(vperp**2+2.d0*vpar**2)*binv0+2.d0*phiovb20/binv0
!
  abar_rr=(wbar*abr**2-apsir_pt**2)/(alpha*Rcell)
  abar_rz=(wbar*abr*abz-apsir_pt*apsiz_pt)/(alpha*Rcell)
  abar_zz=(wbar*abz**2-apsiz_pt**2)/(alpha*Rcell)
!
  bbar_r=(0.5d0*(ubar*abr-aphir_pt)+vpar*binv0*apsir_pt)/(alpha*Rcell)
  bbar_z=(0.5d0*(ubar*abz-aphiz_pt)+vpar*binv0*apsiz_pt)/(alpha*Rcell)
!
  Delta=abar_rr*abar_zz-abar_rz**2
  acoef_r=bbar_r*abar_zz-bbar_z*abar_rz
  acoef_z=bbar_z*abar_rr-bbar_r*abar_rz
!
  if(max(abs(acoef_r),abs(acoef_z))*eps_newt.gt.abs(delta)*sizmax  &
     .or.                                                         &
     (bbar_r**2+bbar_z**2)*eps_newt.gt.sizmax**2*abs(delta)) then
! parabola
    bcoef_r=-bbar_z
    bcoef_z=bbar_r
    varphi(4)=dt
    varphi(1:3)=2.d0*varphi(4)
! 
    do leg_exit=1,3
      ind_vert=modulo(leg_exit,3)+1
      R_side=(R_vert(ind_vert)-R_vert(leg_exit))
      Z_side=(Z_vert(ind_vert)-Z_vert(leg_exit))
      acoef=acoef_r*Z_side-acoef_z*R_side
      bcoef=bcoef_r*Z_side-bcoef_z*R_side
      ccoef=((R_vert(leg_exit)-R)*Z_side       &
           - (Z_vert(leg_exit)-Z)*R_side)*2.d0
      discr=bcoef**2-acoef*ccoef
      if(discr.gt.0.d0) then
        discr=sqrt(discr)
        if(leg_exit.eq.iedge_ent) then
          if(abs(bcoef+discr).gt.abs(bcoef-discr)) then
            dt=(bcoef+discr)/acoef
            if(dt.gt.0.d0) varphi(leg_exit)=min(varphi(leg_exit),dt)
          else
            dt=(bcoef-discr)/acoef
            if(dt.gt.0.d0) varphi(leg_exit)=min(varphi(leg_exit),dt)
          endif
        else
          dt=(bcoef+discr)/acoef
          if(dt.gt.0.d0) varphi(leg_exit)=min(varphi(leg_exit),dt)
          dt=(bcoef-discr)/acoef
          if(dt.gt.0.d0) varphi(leg_exit)=min(varphi(leg_exit),dt)
        endif
      endif
    enddo
! 
    indb=minloc(varphi)
    dt=varphi(indb(1))
    leg_exit=modulo(indb(1),4)
!
    Rnew=R-0.5d0*acoef_r*dt**2+bcoef_r*dt
    Znew=Z-0.5d0*acoef_z*dt**2+bcoef_z*dt
!   
    phinew=(vpar*binv0+(apsir_pt*acoef_r+apsiz_pt*acoef_z)*dt**2/6.d0 &
          -0.5d0*(apsir_pt*bcoef_r+apsiz_pt*bcoef_z)*dt)*dt
    phinew=phi+phinew*bphicovar_interp*(2d0/(R+Rnew))**2
!
    binv=binv_vert1+abr*(Rnew-R_vert(1))+abz*(Znew-Z_vert(1))
!
    vperp=vperp*sqrt(binv0/binv)
    vpar=(vpar*binv0+apsir_pt*(R-Rnew)+apsiz_pt*(Z-Znew))/binv
    R=Rnew
    Z=Znew
    phi=phinew
    return
  endif
!
  R_c=acoef_r/Delta
  Z_c=acoef_z/Delta
  R_a=R-R_c
  Z_a=Z-Z_c
  rate=sqrt(abs(Delta))
  R_s=-bbar_z/rate
  Z_s=bbar_r/rate
!
  if(Delta.gt.0.d0) then
! eliptic orbit
    ellips=.true.
  else
! hyperbolic orbit
    ellips=.false.
  endif
!  
  varphi(4)=dt*rate
  varphi(1:3)=2.d0*varphi(4)
!
  do leg_exit=1,3
    ind_vert=modulo(leg_exit,3)+1
    R_side=(R_vert(ind_vert)-R_vert(leg_exit))
    Z_side=(Z_vert(ind_vert)-Z_vert(leg_exit))
    c_c=R_c*Z_side-Z_c*R_side
    c_s=R_s*Z_side-Z_s*R_side
    c=(R_a-R_vert(leg_exit))*Z_side-(Z_a-Z_vert(leg_exit))*R_side
    if(ellips) then
      discr=c_s**2+c_c**2-c**2
      if(discr.ge.0.d0) then
        discr=sqrt(discr)
        den=c_s**2+c_c**2
        sipl=(-c*c_s+c_c*discr)/den
        simn=(-c*c_s-c_c*discr)/den
        copl=(-c*c_c-c_s*discr)/den
        comn=(-c*c_c+c_s*discr)/den
        if(leg_exit.eq.iedge_ent) then
          phipl=atan2(sipl,copl)
          phimn=atan2(simn,comn)
          if(abs(phipl).gt.abs(phimn)) then
            varphi(leg_exit)=modulo(phipl,2.d0*pi)
          else
            varphi(leg_exit)=modulo(phimn,2.d0*pi)
          endif
        else
          phipl=modulo(atan2(sipl,copl),2.d0*pi)
          phimn=modulo(atan2(simn,comn),2.d0*pi)
          varphi(leg_exit)=min(phipl,phimn)
        endif
      endif
    else
      discr=c**2-c_c**2+c_s**2
      if(discr.ge.0.d0) then
        discr=sqrt(discr)
        if(leg_exit.eq.iedge_ent) then
          sc=-sign(1.d0,c)
          den=abs(c)+discr
          expon1=sc*den/(c_c+c_s)
          expon2=sc*(c_c-c_s)/den
          if(abs(expon1-1.d0).gt.abs(expon2-1.d0)) then
            if(expon2.gt.0.d0.and.expon1.gt.1.d0) varphi(leg_exit)=log(expon1)
          else
            if(expon1.gt.0.d0.and.expon2.gt.1.d0) varphi(leg_exit)=log(expon2)
          endif
        else
          den=-1.d0
          expon=-sign(1.d0,c)*(abs(c)+discr)/(c_c+c_s)
          if(expon.gt.1.d0) den=expon
          expon=-sign(1.d0,c)*(c_c-c_s)/(abs(c)+discr)
          if(expon.gt.1.d0) then
            if(den.gt.0.d0) then
              den=min(expon,den)
            else
              den=expon
            endif
          endif
          if(den.gt.0.d0) varphi(leg_exit)=log(den)
        endif
      endif
    endif
  enddo
!
  indb=minloc(varphi)
  phase=varphi(indb(1))
  leg_exit=modulo(indb(1),4)
  dt=phase/rate
!
  if(ellips) then
    co=cos(phase)
    si=sin(phase)
    coint=si/rate
    siint=(1.d0-co)/rate
  else
    co=cosh(phase)
    si=sinh(phase)
    coint=si/rate
    siint=(co-1.d0)/rate
  endif
  Rnew=R_a+R_c*co+R_s*si
  Znew=Z_a+Z_c*co+Z_s*si
!
  phinew=(vpar*binv0+apsir_pt*(R-R_a)+apsiz_pt*(Z-Z_a))*dt   &
        -(apsir_pt*R_c+apsiz_pt*Z_c)*coint-(apsir_pt*R_s+apsiz_pt*Z_s)*siint
  !phinew=(vpar*binv0*dt)
  !phinew=phi+2d0*phase

  !do k=1,3
  !   phival(k) = 
  !end do
  
  phinew=phi+phinew*bphicovar_interp*(2d0/(R+Rnew))**2
!
  binv=binv_vert1+abr*(Rnew-R_vert(1))+abz*(Znew-Z_vert(1))
  vperp=vperp*sqrt(binv0/binv)
  vpar=(vpar*binv0+apsir_pt*(R-Rnew)+apsiz_pt*(Z-Znew))/binv
!
  if(dt.le.0.0d0) then
     print *,'reg_step_fast, dt<0',dt 
     ierr=1
     return
  elseif(abs(Rnew-R)+abs(Znew-Z).gt.2.00001d0*sizmax) then
     print *,'reg_step_fast, long jump', sizmax, abs(Rnew-R)+abs(Znew-Z)
     R=Rnew
     Z=Znew
     ierr=2
     return
  endif
!
  R=Rnew
  Z=Znew
  phi=phinew
!
end subroutine reg_step_tria
!----------------------------------------------------
function asinh_4lf(x)
  double precision :: x, asinh_4lf
  asinh_4lf = log(x + sqrt(x**2+1))
  return
end function asinh_4lf
!===============================================================================
!!$subroutine ibinsrc(ip, nmin, nmax, ixi, i)
!!$!
!!$! Finds the index  i  of the integer array of increasing numbers ip(nmin:nmax) 
!!$! which satisfies   ip(i-1)=<ixi. Uses binary search algorithm.
!!$!
!!$  implicit none
!!$!
!!$  integer                       :: n, nmin, nmax, i, imin, imax, k
!!$  integer                       :: ixi
!!$  integer, dimension(nmin:nmax) :: ip
!!$!
!!$  imin = nmin
!!$  imax = nmax
!!$  n = nmax - nmin
!!$  do k = 1, n
!!$     i = (imax - imin)/2 + imin
!!$     if(ip(i) .gt. ixi) then
!!$        imax = i
!!$     else
!!$        imin = i
!!$     endif
!!$     if(imax .eq. imin+1) exit
!!$  enddo
!!$  i = imin
!!$  return
!!$end subroutine ibinsrc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!$!
!!$  subroutine reg_step_plot(isort, leg_exit, R, phi, Z, vperp, vpar, dt,ierr)
!!$!
!!$  use constants,         only : pi
!!$  use mesh_mod,  only : bphicovar
!!$  use accuracy_mod,      only : eps_newt
!!$  use renorm_mod,        only : vref, rhoref, elpref
!!$  use reg_step_tria_mod, only : abr,abz,aphir,aphiz,apsir,apsiz,        &
!!$                                binv_vert1,phiovb2_vert1,Rcell,sizmax,  &
!!$                                deltatri,R_vert,Z_vert
!!$!
!!$  implicit none
!!$!
!!$  logical :: lowtri,ellips
!!$!
!!$  integer :: isort,ierr,leg_exit,ind_vert
!!$  integer, dimension(1) :: indb
!!$  double precision :: dt,R,phi,Z,vperp,vpar
!!$  double precision :: binv,Rnew,Znew
!!$  double precision :: alpha,coefphi,aphir_pt,aphiz_pt,apsir_pt,apsiz_pt
!!$  double precision :: binv0,phiovb20,abar_rr,abar_rz,abar_zz,wbar,ubar
!!$  double precision :: bbar_r,bbar_z,Delta,R_a,Z_a,R_c,Z_c,R_s,Z_s
!!$  double precision :: rate,R_bou,Z_bou,c_c,c_s,c,discr,x,den,phase,co,si
!!$  double precision :: sipl,simn,copl,comn,shpl,shmn,chpl,chmn,phipl,phimn
!!$  double precision :: acoef_r,acoef_z,bcoef_r,bcoef_z,acoef,bcoef,ccoef
!!$  double precision :: expon,R_side,Z_side
!!$  double precision, dimension(4) :: varphi
!!$!
!!$  ierr=0
!!$!
!!$  alpha=pi/3.d0
!!$  dt=0.1d0*sizmax
!!$!
!!$  bcoef_r=cos(alpha)
!!$  bcoef_z=sin(alpha)
!!$  varphi(4)=dt
!!$  varphi(1:3)=2.d0*varphi(4)
!!$! 
!!$  do leg_exit=1,3
!!$    ind_vert=modulo(leg_exit,3)+1
!!$    R_side=(R_vert(ind_vert)-R_vert(leg_exit))
!!$    Z_side=(Z_vert(ind_vert)-Z_vert(leg_exit))
!!$    bcoef=bcoef_r*Z_side-bcoef_z*R_side
!!$    ccoef=((R_vert(leg_exit)-R)*Z_side       &
!!$         - (Z_vert(leg_exit)-Z)*R_side)      &
!!$         + eps_newt*deltatri
!!$    dt=ccoef/(bcoef+sign(abs(ccoef)*epsilon(1.d0),bcoef))
!!$    if(dt.gt.0.d0) varphi(leg_exit)=min(varphi(leg_exit),dt)
!!$  enddo
!!$! 
!!$  indb=minloc(varphi)
!!$  dt=varphi(indb(1))
!!$  leg_exit=modulo(indb(1),4)
!!$!
!!$  Rnew=R+bcoef_r*dt
!!$  Znew=Z+bcoef_z*dt
!!$!
!!$  R=Rnew
!!$  Z=Znew
!!$!
!!$  end subroutine reg_step_plot
