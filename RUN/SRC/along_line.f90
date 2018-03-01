program along_line
  use constants,     only : pi, nsorts, charge, one3rd 
  use mesh_mod,      only : ntri, npoint, ntri_inbou, mesh_point, mesh_element, inbou_list, cell_linint &
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
  double precision, dimension(nsorts) :: T_knot, V_knot, D_knot, Tp_knot 
  double precision, dimension(:,:), allocatable :: D_Dm, D_Dp, D_Vp, D_Tp, D_Dt, D_Vt, D_Tt ! dispersions
  double precision, dimension(:), allocatable :: R_layer, Z_layer, N_layer, V_layer, T_layer, Dm_layer, Tp_layer
  double precision, dimension(:), allocatable ::  rel_therm_max
  logical          :: check_in_tri, disp_yes, inside_tri
  external check_in_tri
!
  open(1,file='along_line.inp')
  read (1,*) R_start, Z_start ! start point 
  read (1,*) alpha, ddt       ! angle/(2pi), subdivision factor 
  read (1,*) isort
  read (1,*) R_l, Z_l ! start point for magnetic surface
  close(1)
  alpha = alpha*2.d0*pi

! read everything
  open(1,file='RESULTS/points.dat',form='unformatted')
  read(1) npoint
  allocate(mesh_point(npoint))
  read(1) mesh_point
  close(1)
  open(1,file='RESULTS/triangles.dat',form='unformatted')
  read(1) ntri
  allocate(mesh_element(ntri))
  read(1) mesh_element
  read(1) bphicovar, fulltime, time_l
  close(1)
  open(1,file='RESULTS/thermostat.dat',form='unformatted')
  read(1) n_therm
  allocate(D_therm(n_therm,nsorts),V_therm(n_therm,nsorts),T_therm(n_therm,nsorts))
  read(1) D_therm
  read(1) V_therm
  read(1) T_therm
  close(1)
ntri_inbou=0
  open(1,file='RESULTS/inbou.dat',form='unformatted')
print *,ntri_inbou
  read(1) ntri_inbou
print *,ntri_inbou
  allocate(inbou_list(ntri_inbou))
  read(1) inbou_list
  close(1)

  open(1,file='RESULTS/wpartoverwcell.dat')
  allocate(rel_therm_max(ntri))
  do i=1,ntri
     read(1,*) dummy, dummy, rel_therm_max(i)
  enddo
  close(1)

  allocate( R_layer(ntri_inbou*2), Z_layer(ntri_inbou*2), N_layer(ntri_inbou*2),  & 
            V_layer(ntri_inbou*2), T_layer(ntri_inbou*2), Dm_layer(ntri_inbou*2), Tp_layer(ntri_inbou*2))
  inside_tri = .false.
  do i=1,ntri
     if(check_in_tri(i, R_l, Z_l) ) then
        l_psi_start = mesh_element(i)%layer_psi
        inside_tri = .true.
        exit
     endif
  enddo
  if(.not. inside_tri) stop 'no triangle found, psi_layer'
  np_layer = 0
  do i=1,ntri
     if( mesh_element(i)%layer_psi .eq. l_psi_start) then
        np_layer = np_layer + 1
        R_layer(np_layer) = sum(mesh_point(mesh_element(i)%i_knot(:))%rcoord)*one3rd
        Z_layer(np_layer) = sum(mesh_point(mesh_element(i)%i_knot(:))%zcoord)*one3rd
        N_layer(np_layer) = mesh_element(i)%D_part(isort)
        V_layer(np_layer) = mesh_element(i)%V_part(isort)
        T_layer(np_layer) = mesh_element(i)%T_part(isort)
        Dm_layer(np_layer) = mesh_element(i)%Dm_part(isort)
        Tp_layer(np_layer) = mesh_element(i)%Vt2_part(isort)
     endif
  enddo
  do i=1,np_layer
     write(1111,204) R_layer(i), Z_layer(i), N_layer(i), V_layer(i), T_layer(i), Dm_layer(i), Tp_layer(i) 
  enddo
  close(1111)

!do i=1,n_therm,10
!print *, sngl(T_therm(i:(i+9),isort))
!enddo
!print *, sngl(T_therm(25:72,isort))

!  inquire (file='RESULTS/dispersion.dat', exist=disp_yes)

  disp_yes = .false.

  if(disp_yes) then
     allocate(D_Dt(ntri,nsorts), D_Vt(ntri,nsorts), D_Tt(ntri,nsorts))
     allocate(D_Dm(ntri,nsorts), D_Dp(ntri,nsorts), D_Vp(ntri,nsorts), D_Tp(ntri,nsorts))
     open(1,file='RESULTS/dispersion.dat')
     do i=1,ntri
!!$     r_bc = (mesh_point(mesh_element(i)%i_knot(1))%rcoord + mesh_point(mesh_element(i)%i_knot(2))%rcoord  &
!!$          + mesh_point(mesh_element(i)%i_knot(3))%rcoord)*one3rd
!!$     z_bc = (mesh_point(mesh_element(i)%i_knot(1))%zcoord + mesh_point(mesh_element(i)%i_knot(2))%zcoord  &
!!$          + mesh_point(mesh_element(i)%i_knot(3))%zcoord)*one3rd
!     j = mesh_element(i)%i_therm
        read(1,204) r_bc, z_bc, D_Dm(i,isort), D_Dp(i,isort), D_Vp(i,isort), D_Tp(i,isort), &
                                               D_Dt(i,isort), D_Vt(i,isort), D_Tt(i,isort)
     enddo
     close(1)
  endif

  open(1010,file='along_line.dat',status='replace')
!  close(1010)

  do i=1,npoint
     ePhi_knot = 0.d0
     D_knot = 0.d0
     V_knot = 0.d0
     T_knot = 0.d0
     Tp_knot = 0.d0
     do j=1,mesh_point(i)%n_owners
        ePhi_knot = ePhi_knot + mesh_point(i)%weight_intp(j)                             &
             * mesh_element(mesh_point(i)%i_owner_tri(j))%ePhi_tri
!!$        do k=1,nsorts
!!$           D_knot(k) = D_knot(k) + mesh_point(i)%weight_intp(j)                          &
!!$                * D_therm(mesh_element(mesh_point(i)%i_owner_tri(j))%i_therm,k)
!!$           V_knot(k) = V_knot(k) + mesh_point(i)%weight_intp(j)                          &
!!$                * V_therm(mesh_element(mesh_point(i)%i_owner_tri(j))%i_therm,k)
!!$           T_knot(k) = T_knot(k) + mesh_point(i)%weight_intp(j)                          &
!!$                * T_therm(mesh_element(mesh_point(i)%i_owner_tri(j))%i_therm,k)
        do k=1,nsorts
           D_knot(k) = D_knot(k) + mesh_point(i)%weight_intp(j)                          &
                * mesh_element(mesh_point(i)%i_owner_tri(j))%D_part(k)
           V_knot(k) = V_knot(k) + mesh_point(i)%weight_intp(j)                          &
                * mesh_element(mesh_point(i)%i_owner_tri(j))%V_part(k)
           T_knot(k) = T_knot(k) + mesh_point(i)%weight_intp(j)                          &
                * mesh_element(mesh_point(i)%i_owner_tri(j))%T_part(k)
           Tp_knot(k) = T_knot(k) + mesh_point(i)%weight_intp(j)                          &
                * mesh_element(mesh_point(i)%i_owner_tri(j))%Vt2_part(k)
        enddo
     enddo
     do k=1,nsorts
        mesh_point(i)%pl_parms_knot(1,k) = D_knot(k)
        mesh_point(i)%pl_parms_knot(2,k) = V_knot(k)
        mesh_point(i)%pl_parms_knot(3,k) = T_knot(k)
        mesh_point(i)%pl_parms_knot(4,k) = Tp_knot(k)
     enddo
  enddo

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


!print *,mesh_element(ind_tri)%i_therm,T_therm(mesh_element(ind_tri)%i_therm,isort),ind_tri
!if(i/35*35.eq.i) pause
     rho = sqrt((R-R_start)**2 + (Z-Z_start)**2) 
!     write(1011,204) R,Z, mesh_element(ind_tri)%T_part(isort), D_Tp(ind_tri,isort)
     do k=1,nsorts
        f_vert(:) = mesh_point( mesh_element(ind_tri)%i_knot(:))%pl_parms_knot(1,k) 
        D_knot(k) = cell_linint(ind_tri, R, Z, f_vert)
        f_vert(:) = mesh_point( mesh_element(ind_tri)%i_knot(:))%pl_parms_knot(2,k)
        V_knot(k) = cell_linint(ind_tri, R, Z, f_vert)
        f_vert(:) = mesh_point( mesh_element(ind_tri)%i_knot(:))%pl_parms_knot(3,k)
        T_knot(k) = cell_linint(ind_tri, R, Z, f_vert)
        f_vert(:) = mesh_point( mesh_element(ind_tri)%i_knot(:))%pl_parms_knot(4,k)
        Tp_knot(k) = cell_linint(ind_tri, R, Z, f_vert)
     enddo
     if(disp_yes) then
        write(1011,204) R, Z, D_knot(isort), D_Dt(ind_tri,isort),   &
                              V_knot(isort), D_Vt(ind_tri,isort),   &
                              T_knot(isort), D_Tt(ind_tri,isort),   &
                              mesh_element(ind_tri)%D_part(isort),  D_Dp(ind_tri,isort),   &
                              mesh_element(ind_tri)%V_part(isort),  D_Vp(ind_tri,isort),   &
                              mesh_element(ind_tri)%T_part(isort),  D_Tp(ind_tri,isort)
     endif
!     write(1010,204) R, Z, cell_linint(ind_tri, 1, isort, R, Z), cell_linint(ind_tri, 2, isort, R, Z),  &
!                           cell_linint(ind_tri, 3, isort, R, Z), mesh_element(ind_tri)%D_part(isort),   &
!                           mesh_element(ind_tri)%V_part(isort), mesh_element(ind_tri)%T_part(isort)
     write(1010,204) R, Z,   &
          D_therm(mesh_element(ind_tri)%i_therm,isort),V_therm(mesh_element(ind_tri)%i_therm,isort) ,  &
          T_therm(mesh_element(ind_tri)%i_therm,isort),  &
          D_knot(isort),V_knot(isort),T_knot(isort),Tp_knot(isort),   &
          mesh_element(ind_tri)%Dm_part(isort), rel_therm_max(ind_tri),bmod,psipol,bphicovar, potential

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
  if(disp_yes) close(1011)
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

