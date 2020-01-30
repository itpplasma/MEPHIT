!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine perp_step(dtstep, R, phi, Z, vperp, vpar, ind_tri, ierr, cg_rng)
!
  use constants,     only : isort2 ! for ions only
  use accuracy_mod,  only : eps_newt, tinydp
  use for_macrostep, only : D_perp
  use mesh_mod,      only : mesh_point, mesh_element
  use hujnja, only: numstep
!
  implicit none
!
  integer :: ierr, ierr_tri, ind_tri, ind_vert, i, ind_tri_new, leg_exit, itri_sav, j
!  integer ::  numrefl, numreflmax=10
  integer, dimension(1) :: indb
  double precision :: dtstep, R, phi, Z, cg_rng(0:5), dltot, dl, sqxi
  double precision :: Rnew, Znew, ran_u_sq3, deltatri            , dltotsav
  double precision :: bcoef_r, bcoef_z, bcoef, ccoef
  double precision :: R_side, Z_side, vperp, vpar
  double precision, dimension(3) :: R_vert, Z_vert
  double precision, dimension(4) :: varphi
  logical :: check_in_tri, reissue
  integer :: look_4_tri

double precision :: polpot_do, polpot_p, polpot_m, xi1, xi2, Rin, Zin
integer :: istep, ind_tri_in 
  external ran_u_sq3, look_4_tri
!
  Rin = R
  Zin = Z
  ind_tri_in = ind_tri
  istep=1
  xi1 = ran_u_sq3(cg_rng)
  xi2 = ran_u_sq3(cg_rng)

111 continue
  ierr = 0
  ierr_tri = 0
!if(istep.eq.1) then
  bcoef_r = -sqrt(2.d0*D_perp*dtstep)*xi1 + D_perp/Rin*dtstep
  bcoef_z = -sqrt(2.d0*D_perp*dtstep)*xi2
!else
!  bcoef_r = sqrt(2.d0*D_perp*dtstep)*xi1 + D_perp/Rin*dtstep
!  bcoef_z = sqrt(2.d0*D_perp*dtstep)*xi2
!  R = Rin
!  Z = Zin
!  ind_tri = ind_tri_in
!endif
  dltot = sqrt(bcoef_r**2 + bcoef_z**2)
  bcoef_r = bcoef_r/dltot
  bcoef_z = bcoef_z/dltot
  countdown: do while(dltot .gt. tinydp) 
     deltatri = mesh_element(ind_tri)%det_3
     R_vert(:) = mesh_point( mesh_element(ind_tri)%i_knot(:) )%rcoord
     Z_vert(:) = mesh_point( mesh_element(ind_tri)%i_knot(:) )%zcoord
     varphi(4) = dltot
     varphi(1:3) = 2.d0*varphi(4)
! 
     do i=1,3
        ind_vert = modulo(i,3) + 1
        R_side = (R_vert(ind_vert) - R_vert(i))
        Z_side = (Z_vert(ind_vert) - Z_vert(i))
        bcoef = bcoef_r*Z_side - bcoef_z*R_side
        ccoef = ( (R_vert(i)-R)*Z_side       &
                - (Z_vert(i)-Z)*R_side )     
        ccoef =  ccoef + sign(eps_newt, ccoef)
        dl = ccoef/( bcoef + sign(abs(ccoef)*eps_newt,bcoef) )
        if(dl .gt. 0.d0) then
           varphi(i) = min(varphi(i),dl)
        endif
     enddo
!  
     indb = minloc(varphi)
     dl = varphi(indb(1))
     dltot = dltot - dl
     leg_exit = modulo(indb(1),4)
     if(abs(bcoef_r*dl)+abs(bcoef_z*dl).gt.2.00001d0*mesh_element(ind_tri)%sizmaxtri) then
        ierr=2
        print *,'long jump',bcoef_r*dl,bcoef_z*dl,mesh_element(ind_tri)%sizmaxtri
        stop
     endif

     if(leg_exit .eq. 0) then ! step is over?
        Rnew = R + bcoef_r*dl
        Znew = Z + bcoef_z*dl
        ind_tri_new = ind_tri
        if(.not. check_in_tri(ind_tri, Rnew, Znew) ) then
           ind_tri_new = look_4_tri(ind_tri, Rnew, Znew, ierr_tri)
           if(ierr_tri .eq. -1 ) then
              ierr = -100 ! particle got lost
              return
!!$           elseif(mesh_element(ind_tri_new)%iqq_gc .ne. 0) then ! boundary
!!$              call check_boundary(isort2, ind_tri_new, mesh_element(ind_tri_new)%iqq_gc, vpar, vperp, R, Z, ierr)
!!$              return               
           endif
! otherwise it was regular step to next triangle
        endif
     else 
        ind_tri_new = mesh_element(ind_tri)%neighbour(leg_exit)
        if(ind_tri_new .lt. 0) then ! wall
           ierr = -1
           return
        endif
! regular step to next triangle
        Rnew = R + bcoef_r*dl
        Znew = Z + bcoef_z*dl
        if(.not. check_in_tri(ind_tri_new, Rnew, Znew) ) then
           ind_tri_new = look_4_tri(ind_tri_new, Rnew, Znew, ierr_tri)
           if(ierr_tri .eq. -1) then
              ierr = -3 ! particle got lost
!print *,varphi, leg_exit, mesh_element(ind_tri)%neighbour(leg_exit)
              call write_tri
              return
!!$           elseif(mesh_element(ind_tri_new)%iqq_gc .ne. 0) then ! boundary
!!$              call check_boundary(isort2, ind_tri_new, mesh_element(ind_tri_new)%iqq_gc, vpar, vperp, R, Z, ierr)
!!$              return 
           endif
        endif
     endif

     R = Rnew
     Z = Znew
     ind_tri = ind_tri_new
  end do countdown

  return
!-----------------------------------------------------------------------------
contains

  subroutine write_tri
  integer :: i, j, itn
    do i = 1,4
       j = modulo(i,3) + 1
       write(222,*) mesh_point( mesh_element(ind_tri)%i_knot(j) )%rcoord, &
            mesh_point( mesh_element(ind_tri)%i_knot(j) )%zcoord
    enddo
    close(222)
    write(223,*) R, Z, ind_tri
    write(223,*) Rnew, Znew, ind_tri
    close(223)
    itn = mesh_element(ind_tri)%neighbour(leg_exit)
    do i = 1,4
       j = modulo(i,3) + 1
       write(224,*) mesh_point( mesh_element(itn)%i_knot(j) )%rcoord, &
            mesh_point( mesh_element(itn)%i_knot(j) )%zcoord
    enddo
    close(224)
!    pause
  end subroutine write_tri
end subroutine perp_step
