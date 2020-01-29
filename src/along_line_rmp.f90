program along_line
  use constants,     only : pi, nsorts, charge, one3rd 
  use mesh_mod,      only : ntri, npoint, ntri_inbou, mesh_point, mesh_element, cell_linint &
                            , bphicovar
  use for_macrostep, only : time_l, fulltime
  use parmesh_mod,   only : n_therm, w_cell_V, w_cell_T, D_therm, V_therm, T_therm 
  use reg_step_tria_mod, only : sizmax, deltatri, R_vert, Z_vert
  implicit none
!
  integer          :: i, j, k, isort, ind_tri, ind_tri_new, leg_exit, itri_sav, l_psi_start, np_layer
  double precision :: R, Z, R_start, Z_start, alpha, ddt, bmod, psipol, dummy, potential
  double precision :: Dm, R_sav, Z_sav, ePhi_knot, r_bc, z_bc, rho,  R_l, Z_l
  double precision, dimension(3) :: f_vert, B_vert, P_vert    
  double precision, dimension(:,:), allocatable :: Re_jR,Re_jZ,Aim_jR,Aim_jZ    
  double precision, dimension(:,:), allocatable :: Re_dn,Aim_dn,Re_dpp,Aim_dpp,Re_dpc,Aim_dpc
  logical          :: check_in_tri, disp_yes, inside_tri
  external check_in_tri
!
  open(1,file='along_line.inp')
  read (1,*) R_start, Z_start ! start point 
  read (1,*) alpha, ddt       ! angle/(2pi), subdivision factor 
  read (1,*) isort
  close(1)
  alpha = alpha*2.d0*pi

! read everything
  open(1,file='START_PRFS/points.dat',form='unformatted')
  read(1) npoint
  allocate(mesh_point(npoint))
  read(1) mesh_point
  close(1)
  open(1,file='START_PRFS/triangles.dat',form='unformatted')
  read(1) ntri
  allocate(mesh_element(ntri))
  read(1) mesh_element
  read(1) bphicovar, fulltime, time_l
  close(1)
  allocate(Re_jR(2,ntri),Re_jZ(2,ntri),Aim_jR(2,ntri),Aim_jZ(2,ntri))
  allocate(Re_dn(2,ntri),Re_dpp(2,ntri),Aim_dn(2,ntri),Aim_dpp(2,ntri))
  allocate(Re_dpc(2,ntri),Aim_dpc(2,ntri))
!
  open(1,file='PLOTTING/curr_toplot.dat')
  do i=1, ntri
    read(1,*) Re_jR(:,i),Re_jZ(:,i),Aim_jR(:,i),Aim_jZ(:,i)
  enddo
  close(1)
!
  open(1,file='PLOTTING/dens_toplot.dat')
  do i=1, ntri
    read(1,*) Re_dn(:,i),Aim_dn(:,i)
  enddo
  close(1)
!
  open(1,file='PLOTTING/ppres_toplot.dat')
  do i=1, ntri
    read(1,*) Re_dpp(:,i),Aim_dpp(:,i)
  enddo
  close(1)
!
  open(1,file='PLOTTING/parcurr_toplot.dat')
  do i=1, ntri
    read(1,*) Re_dpc(:,i),Aim_dpc(:,i)
  enddo
  close(1)
!
  open(1010,file='along_line_equ.dat',status='replace')
  open(1011,file='along_line_rmp.dat',status='replace')

  R = R_start
  Z = Z_start
  ind_tri = 0
! look for the start triangle:
  do i=1, ntri
     if( check_in_tri(i, R, Z) ) then
        ind_tri = i
        exit
     endif
  enddo
  if(ind_tri .eq. 0) stop 'point out of domain'
! beam tracing:
  do i = 1,20000
     R_sav = R
     Z_sav = Z
     itri_sav = ind_tri

     R_vert(:) = mesh_point( mesh_element(ind_tri)%i_knot(:) )%rcoord
     Z_vert(:) = mesh_point( mesh_element(ind_tri)%i_knot(:) )%zcoord
     Dm = mesh_element(ind_tri)%Dm_part(isort)
     B_vert(:) = mesh_point( mesh_element(ind_tri)%i_knot(:) )%b_mod
     bmod = cell_linint( ind_tri, R, Z, B_vert)
     P_vert(:) = mesh_point( mesh_element(ind_tri)%i_knot(:) )%psi_pol
     psipol = cell_linint(ind_tri, R, Z, P_vert)
     P_vert(:) = mesh_point( mesh_element(ind_tri)%i_knot(:) )%PhiovB2
     potential = cell_linint(ind_tri, R, Z, P_vert)*bmod**2

     rho = sqrt((R-R_start)**2 + (Z-Z_start)**2) 
     write(1010,204) rho, R, Z, Dm, bmod, psipol, potential
     j=ind_tri
     write(1011,204) rho, R, Z, Re_jR(:,j),Re_jZ(:,j),Aim_jR(:,j),Aim_jZ(:,j)   &
                              , Re_dn(:,j),Aim_dn(:,j),Re_dpp(:,j),Aim_dpp(:,j) &
                              , Re_dpc(:,j),Aim_dpc(:,j)

     deltatri = mesh_element(ind_tri)%det_3
     sizmax = mesh_element(ind_tri)%sizmaxtri

     call reg_step_plot(alpha, ddt, leg_exit, R, Z)

     if(leg_exit .eq. 0) then  ! step is over
        ind_tri_new = ind_tri
     else
        ind_tri_new = mesh_element(ind_tri)%neighbour(leg_exit)
        if (ind_tri_new .lt. 0) stop      ! outer boundary, death
     endif
     ind_tri = ind_tri_new
     if( .not. check_in_tri(ind_tri, R, Z) ) then
        print *, ind_tri, R, Z
        do k = 1,4
           j = modulo(k,3) + 1
           write(222,*) mesh_point( mesh_element(itri_sav)%i_knot(j) )%rcoord, &
                        mesh_point( mesh_element(itri_sav)%i_knot(j) )%zcoord
        enddo
        do k = 1,4
           j = modulo(k,3) + 1
           write(224,*) mesh_point( mesh_element(ind_tri)%i_knot(j) )%rcoord, &
                        mesh_point( mesh_element(ind_tri)%i_knot(j) )%zcoord
        enddo
        write(223,*) R_sav, Z_sav
        write(223,*) R, Z
        close(222)
        close(223)
        close(224)   
        close(1010)    
        stop 'out of triangle'
     endif
  enddo
204 format(1000(ES12.5E2,1X))
205 format(1000(ES23.16E2,1X))
  close(1010)
  close(1011)
  stop 'too many steps'
end program along_line
!---------------------------------------------------------------------------
subroutine reg_step_plot(alpha, ddt, leg_exit, R, Z)
!
  use accuracy_mod,      only : eps_newt
  use reg_step_tria_mod, only : sizmax, deltatri, R_vert, Z_vert
  implicit none
  double precision,    intent(in) :: alpha, ddt
  integer,            intent(out) :: leg_exit
  double precision, intent(inout) :: R, Z
  integer          :: ind_vert
  double precision :: dt, binv
  double precision :: bcoef_r, bcoef_z, bcoef, ccoef
  double precision :: R_side, Z_side
  integer, dimension(1) :: indb
  double precision, dimension(4) :: varphi
!
  dt = sizmax*ddt
!
  bcoef_r = cos(alpha)
  bcoef_z = sin(alpha)
  varphi(4) = dt
  varphi(1:3) = 2.d0*varphi(4)
! 
  do leg_exit = 1,3
     ind_vert = modulo(leg_exit,3) + 1
     R_side = (R_vert(ind_vert) - R_vert(leg_exit))
     Z_side = (Z_vert(ind_vert) - Z_vert(leg_exit))
     bcoef = bcoef_r*Z_side - bcoef_z*R_side
     ccoef = ((R_vert(leg_exit) - R)*Z_side       &
            - (Z_vert(leg_exit) - Z)*R_side)      &
             + eps_newt*deltatri
     dt = ccoef/(bcoef + sign(abs(ccoef)*epsilon(1.d0),bcoef))
     if(dt .gt. 0.d0) varphi(leg_exit) = min(varphi(leg_exit), dt)
  enddo
! 
  indb = minloc(varphi)
  dt = varphi(indb(1))
  leg_exit = modulo(indb(1),4)
!
  R = R + bcoef_r*dt
  Z = Z + bcoef_z*dt
end subroutine reg_step_plot

