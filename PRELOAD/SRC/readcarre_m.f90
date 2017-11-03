program readcarre_m
  use from_nrtype
  use constants,     only : pi, one3rd, nsorts, charge, &
                            ev2erg,amass
!  use field_eq_mod,  only : btf, rtf
  use mesh_mod
  use for_macrostep, only : fulltime, time_l, t_min, d_min
  use for_mpi,       only : flux_q, flux_j, tor_mom
  use parmesh_mod,   only :  n_therm, D_therm, V_therm, T_therm, vol_therm
  implicit none

  real(dp) :: rbtor, xv(npqq), yv(npqq), psimin, psimax, r_bc, z_bc, sumD, field_ratio
  real(dp), parameter ::  eps_dist=1.d-8, ti0=3.d3, di0 = 5.d13, ePhi0=2d3
  real(dp) :: btf, rtf, psi_loc, bmod_loc
  character(len=7) :: dummy7
  character(len=6) :: dummy6
  character(len=5) :: dummy5
  character(len=4) :: dummy4
  character(len=2) :: dummy2
  character(len=1) :: dummy1
  character(len=28) :: dummy28
  character(len=17) :: dummy17
  complex(dpc) :: vert1, vert2, vert3, vert4, vert_bc
  integer(i4b) :: nx, ny, ncut, niso, i, j, k, m, nnn, nlayer, ntriqq=2 !4 with bc
  integer(i4b) :: j1, j2, j3, idummy, nbou_qq, nr_core, nt_core, nbou, i_tht
  integer(i4b) :: inbou_count, i_c_test, numcellqq, numcell_gc, numtri_gc
  integer(i4b), dimension(:), allocatable :: nxcut, nycut, icell_qq !, ix_qq, iy_qq 
  real(dp), dimension(:,:), allocatable :: xleg, yleg, xvertice, yvertice
  real(dp), dimension(n_owners_max) :: area_one_tri, dist_bc
  real(dp) :: r2r1, r3r1, z2z1, z3z1, r3r2, z3z2, det_tmp, vol_sum=0.0d0, R1, Z1, R2, Z2
! guess of potential (look Rozhansky paper)
!  real(dp), parameter :: psi2r=9.0d0, d_right=4.5d0, d_left=10.0d0, alpha_r=300.d0, alpha_l=-600.d0, phi_r0=-1167.3d0  
!  real(dp) :: phi_guess, r_guess
  integer, dimension(:), allocatable :: ix_qq_s, iy_qq_s, i_therm_s 
  integer, dimension(:,:), allocatable :: neighbour_s

  real(dp) :: ePhi_knot
  real(dp), dimension(nsorts) :: D_knot, V_knot, T_knot
  real(dp), dimension(:,:), allocatable :: vol_core, ind_core_tri 
!
  real(dp) :: ePhi_prof,dePhi_dpsi
  real(dp), dimension(nsorts) :: dens_prof,temp_prof,ddens_dpsi,dtemp_dpsi,vt2inv
!
  integer :: ntri_core                                     !<= SRR
  integer,  dimension(:), allocatable :: indrho            !<= SRR
  real(dp) :: dummy                                        !<= SRR
  real(dp), dimension(:), allocatable :: psiofrho,rho      !<= SRR
!
  open(211,file='index_core.pnt')
  read(211,*) idummy, nr_core, nt_core 
  close(211)
!
  ntri_core=nt_core-1+2*(nt_core-1)*(nr_core-2)            !<= SRR
  allocate(indrho(ntri_core))                              !<= SRR
  ntri=0                                                   !<= SRR
  do j=1, nt_core-1 ! around O-point                       !<= SRR
    ntri = ntri + 1                                        !<= SRR
    indrho(ntri) = 1                                       !<= SRR
  enddo                                                    !<= SRR
  do i=2,nr_core-1                                         !<= SRR
     do j=1, nt_core-1                                     !<= SRR
        ntri = ntri + 1                                    !<= SRR
        indrho(ntri) = i                                   !<= SRR
        ntri = ntri + 1                                    !<= SRR
        indrho(ntri) = i                                   !<= SRR
     enddo                                                 !<= SRR
  enddo                                                    !<= SRR
  allocate(psiofrho(0:nr_core),rho(0:nr_core))             !<= SRR
  open(211,file='RZpsirho.dat')                            !<= SRR
  do i=0,nr_core                                           !<= SRR
    read(211,*) dummy,dummy,psiofrho(i),rho(i)             !<= SRR
  enddo                                                    !<= SRR
  close(211)                                               !<= SRR
!
  nt_core = nt_core - 1
  open(212,file='points_m.fmt')
  read(212,*) npoint
  allocate( mesh_point(npoint) )
  do i=1, npoint
     read(212,*) mesh_point(i)%rcoord, mesh_point(i)%zcoord, mesh_point(i)%psi_pol,    &
          mesh_point(i)%b_mod, mesh_point(i)%b_phi
     mesh_point(i)%PhiovB2 = 0.d0
     mesh_point(i)%pl_parms_knot = 0.d0
     mesh_point(i)%n_owners = 0
     mesh_point(i)%i_owner_tri(1:n_owners_max) = -1000000
     mesh_point(i)%weight_intp(1:n_owners_max) = 0.d0   
  enddo
  read(212,*) btf, rtf
  close(212)

  open(2, file='17151.aug_qq.ele')
  read(2,*) ntri
  allocate(mesh_element(ntri))
  do i=1,ntri
     mesh_element(i)%ix_qq = 0
     mesh_element(i)%iy_qq = 0
     mesh_element(i)%i_knot(1:legs) = 0
     mesh_element(i)%neighbour(1:legs) = -1000000
     mesh_element(i)%sizmaxtri = 0.d0
     mesh_element(i)%det_3 = 0.d0
     mesh_element(i)%oneoverD3 = 0.d0
     mesh_element(i)%V_tri = 0.d0
     mesh_element(i)%dBinv_dR = 0.d0
     mesh_element(i)%dBinv_dZ = 0.d0
     mesh_element(i)%dPsi_dR = 0.d0
     mesh_element(i)%dPsi_dZ = 0.d0
     mesh_element(i)%dPhiovB2_dR = 0.d0
     mesh_element(i)%dPhiovB2_dZ = 0.d0
     mesh_element(i)%Dm_part(1:nsorts) = 0.d0
     mesh_element(i)%D_part(1:nsorts) = 0.d0
     mesh_element(i)%V_part(1:nsorts) = 0.d0
     mesh_element(i)%T_part(1:nsorts) = 0.d0
     mesh_element(i)%Vt2_part(1:nsorts) = 0.d0
     mesh_element(i)%ePhi_tri = 0.d0
     mesh_element(i)%i_therm = 0
!
     mesh_element(i)%thermforces=0.d0
  end do
  do m=1, ntri
     read(2,*) idummy, (mesh_element(m)%i_knot(j), j=1,legs)
  enddo
  close(2)
!
  open(2,file='wsrr_rec_fix.dat')                          !<= SRR
  do i=1,ntri_core                                         !<= SRR
    m=indrho(i)                                            !<= SRR
    if(m.lt.nr_core-2) then                                !<= SRR
      write(2,*) (rho(m)-rho(m-1))/(rho(1)-rho(0))         !<= SRR
    else                                                   !<= SRR
      write(2,*) 1.d0                                      !<= SRR
    endif                                                  !<= SRR
  enddo                                                    !<= SRR
  do i=ntri_core+1,ntri                                    !<= SRR
    write(2,*) 1.d0                                        !<= SRR
  enddo                                                    !<= SRR
  close(2)                                                 !<= SRR
  deallocate(indrho,psiofrho,rho)                          !<= SRR
!
  allocate(ix_qq_s(ntri), iy_qq_s(ntri), i_therm_s(ntri))
  allocate(neighbour_s(legs,ntri))
  open(2, file='neighbour.fmt')
  read(2,*) nbou_qq
  do i=1, ntri
     read(2,999) neighbour_s(1:legs,i), ix_qq_s(i), iy_qq_s(i), i_therm_s(i)
  enddo
  close(2)
999 format(1000(i6,2x))
  do i=1,ntri
     mesh_element(i)%ix_qq = ix_qq_s(i)
     mesh_element(i)%iy_qq = iy_qq_s(i)
     mesh_element(i)%neighbour(1:legs) = neighbour_s(1:legs,i)
     mesh_element(i)%i_therm = i_therm_s(i)
  end do

  do i=1,ntri
     r2r1 = mesh_point(mesh_element(i)%i_knot(2))%rcoord -  &
            mesh_point(mesh_element(i)%i_knot(1))%rcoord
     r3r1 = mesh_point(mesh_element(i)%i_knot(3))%rcoord -  &
            mesh_point(mesh_element(i)%i_knot(1))%rcoord
     z2z1 = mesh_point(mesh_element(i)%i_knot(2))%zcoord -  &
            mesh_point(mesh_element(i)%i_knot(1))%zcoord
     z3z1 = mesh_point(mesh_element(i)%i_knot(3))%zcoord -  &
            mesh_point(mesh_element(i)%i_knot(1))%zcoord
     det_tmp = r2r1*z3z1 - r3r1*z2z1
     mesh_element(i)%oneoverD3 = 1.0_dp/det_tmp
     if(det_tmp .le. 0.0_dp) then
        det_tmp = -det_tmp
     endif
     mesh_element(i)%det_3 = det_tmp
  enddo


! search for triangles owning given knot:
! .or. mesh_element(mesh_point(nnn)%i_owner_tri(j))%iqq_gc.ne.0

  do i=2,npoint ! O-point excluded
     owner_tri: do k=1,ntri
        do m=1,legs
           if(i .eq. mesh_element(k)%i_knot(m)) then
              mesh_point(i)%n_owners = mesh_point(i)%n_owners + 1
              if(mesh_point(i)%n_owners .gt. n_owners_max) then
                 write(*,*) 'too many triangles own the point', i, mesh_point(i)%rcoord, mesh_point(i)%zcoord
                 stop
              else
                 mesh_point(i)%i_owner_tri(mesh_point(i)%n_owners) = k
              endif
              cycle owner_tri
           endif
        enddo
     enddo owner_tri
  enddo
! weights for averaging over tr.'s:
  do i=2,npoint ! O-point excluded
     do j=1,n_owners_max
        area_one_tri(j) = 0.d0
        dist_bc(j) = 1.d0
        if(mesh_point(i)%i_owner_tri(j) .gt. 0) then
           area_one_tri(j) = mesh_element(mesh_point(i)%i_owner_tri(j))%det_3
           r_bc = ( mesh_point(mesh_element(mesh_point(i)%i_owner_tri(j))%i_knot(1))%rcoord +                &
                    mesh_point(mesh_element(mesh_point(i)%i_owner_tri(j))%i_knot(2))%rcoord +                &
                    mesh_point(mesh_element(mesh_point(i)%i_owner_tri(j))%i_knot(3))%rcoord  )*one3rd
           z_bc = ( mesh_point(mesh_element(mesh_point(i)%i_owner_tri(j))%i_knot(1))%zcoord +                &
                    mesh_point(mesh_element(mesh_point(i)%i_owner_tri(j))%i_knot(2))%zcoord +                &
                    mesh_point(mesh_element(mesh_point(i)%i_owner_tri(j))%i_knot(3))%zcoord  )*one3rd
           dist_bc(j) = sqrt( (mesh_point(i)%rcoord - r_bc)**2 + (mesh_point(i)%zcoord - z_bc)**2 )
        endif
     enddo
     do j=1,mesh_point(i)%n_owners  
        area_one_tri(j) = 0.d0
     enddo
     if( sum(area_one_tri(:)/dist_bc(:)) .le. 1.d-8 ) then
        mesh_point(i)%weight_intp(:) = 0.d0
     else
        mesh_point(i)%weight_intp(:) = (area_one_tri(:)/dist_bc(:)) / sum(area_one_tri(:)/dist_bc(:))
     endif
  enddo    

  psimin = minval(mesh_point%psi_pol)
  psimax = maxval(mesh_point%psi_pol)
  print *,'psimin = ',psimin,'  psimax = ',psimax

! length measure and volume of cells:
  do i=1,ntri
     r2r1 = mesh_point(mesh_element(i)%i_knot(2))%rcoord - mesh_point(mesh_element(i)%i_knot(1))%rcoord
     r3r1 = mesh_point(mesh_element(i)%i_knot(3))%rcoord - mesh_point(mesh_element(i)%i_knot(1))%rcoord
     r3r2 = mesh_point(mesh_element(i)%i_knot(3))%rcoord - mesh_point(mesh_element(i)%i_knot(2))%rcoord
     z2z1 = mesh_point(mesh_element(i)%i_knot(2))%zcoord - mesh_point(mesh_element(i)%i_knot(1))%zcoord
     z3z1 = mesh_point(mesh_element(i)%i_knot(3))%zcoord - mesh_point(mesh_element(i)%i_knot(1))%zcoord
     z3z2 = mesh_point(mesh_element(i)%i_knot(3))%zcoord - mesh_point(mesh_element(i)%i_knot(2))%zcoord
     mesh_element(i)%sizmaxtri = max( abs(r2r1)+ abs(z2z1), abs(r3r1)+ abs(z3z1),  abs(r3r2)+ abs(z3z2) )
     r_bc = ( mesh_point(mesh_element(i)%i_knot(1))%rcoord +                &
              mesh_point(mesh_element(i)%i_knot(2))%rcoord +                &
              mesh_point(mesh_element(i)%i_knot(3))%rcoord  )*one3rd
     mesh_element(i)%V_tri = mesh_element(i)%det_3 * r_bc * pi ! *2/2


     do m=1,4 
        j=modulo(m,3)+1
        write(110,204) mesh_point(mesh_element(i)%i_knot(j))%rcoord, mesh_point(mesh_element(i)%i_knot(j))%zcoord, &
             mesh_point(mesh_element(i)%i_knot(j))%psi_pol
     enddo
     write(110,*)
  enddo

!write(99,204) mesh_point(2134)%rcoord, mesh_point(2134)%zcoord
!write(99,204) mesh_point(2927)%rcoord, mesh_point(2927)%zcoord

! neighbourhood of the O-point
  nbou = 0
  vol_sum = 0.0d0
  allocate(inbou_list(nt_core))
  do i=1,ntri
     if( mesh_element(i)%i_therm .eq. 1) then
        nbou =   nbou + 1
        inbou_list(nbou)%i_tri_all = i
        inbou_list(nbou)%vol_norm = mesh_element(i)%V_tri
        vol_sum = vol_sum + inbou_list(nbou)%vol_norm
     endif
  enddo
  print *, nbou, nt_core
  inbou_list(:)%vol_norm = inbou_list(:)%vol_norm/vol_sum

! guess for start profiles:
  do i=1,ntri
     psi_loc = ( mesh_point(mesh_element(i)%i_knot(1))%psi_pol +           &
                 mesh_point(mesh_element(i)%i_knot(2))%psi_pol +           &
                 mesh_point(mesh_element(i)%i_knot(3))%psi_pol )*one3rd
!
     dens_prof=(psi_loc - psimin)/psimax*di0 + d_min
     ddens_dpsi=1.d0/psimax*di0
     temp_prof=(psi_loc - psimin)/psimax*ti0 + t_min
     dtemp_dpsi=1.d0/psimax*ti0
!     ePhi_prof=(psi_loc - psimin)/psimax*ePhi0
     vt2inv=amass/(2.d0*ev2erg*temp_prof)
!
     ePhi_prof=0.d0
     do j=1,3
       k=mesh_element(i)%i_knot(j)
       psi_loc = mesh_point(k)%psi_pol
       ePhi_knot=(psi_loc - psimin)/psimax*ePhi0
       mesh_point(k)%PhiovB2 = ePhi_knot/(mesh_point(k)%b_mod)**2
       ePhi_prof=ePhi_prof+ePhi_knot
     enddo
     ePhi_prof=ePhi_prof*one3rd
     dePhi_dpsi=1.d0/psimax*ePhi0
!
!     mesh_element(i)%D_part(1:nsorts) = (psi_loc - psimin)/psimax*di0 + d_min
!     mesh_element(i)%V_part(1:nsorts) = 0.d0
!     mesh_element(i)%T_part(1:nsorts) = ti0 
     mesh_element(i)%D_part(1:nsorts) = dens_prof
     mesh_element(i)%V_part(1:nsorts) = 0.d0
     mesh_element(i)%T_part(1:nsorts) = temp_prof
     mesh_element(i)%ePhi_tri = ePhi_prof
!
     mesh_element(i)%thermforces(1,:)=ddens_dpsi/dens_prof-1.5d0*dtemp_dpsi/temp_prof !&
!TODO uncomment again                                     +charge*dePhi_dpsi/temp_prof
     mesh_element(i)%thermforces(2,:)=dtemp_dpsi/temp_prof*vt2inv
! Caution: thermodynamic force A_2 has been divided by v_T^2 here
                                      
  enddo

  do j=1, nsorts
     do i=1,ntri
        mesh_element(i)%Dm_part(j) = mesh_element(i)%D_part(j)                       &
                                /sum(mesh_element(:)%D_part(j)*mesh_element(:)%V_tri)
     enddo
  enddo

  n_therm = i_therm_s(ntri)
  allocate(D_therm(n_therm,nsorts), V_therm(n_therm,nsorts), T_therm(n_therm,nsorts))
  allocate(vol_therm(n_therm))

  vol_sum = 0.0d0
  vol_therm = 0.0d0
  do i=1,ntri
     vol_sum = vol_sum + mesh_element(i)%V_tri
     vol_therm(mesh_element(i)%i_therm) = vol_therm(mesh_element(i)%i_therm) + mesh_element(i)%V_tri
  enddo 
  print *,'volumes', sum(vol_therm), vol_sum

! The following works for quadrangles only (tr.s summarized pairwise):
  D_therm = 0.0d0
  V_therm = 0.0d0
  T_therm = 0.0d0
  do j=1, nsorts
     do i=1,ntri
        D_therm(mesh_element(i)%i_therm,j) = D_therm(mesh_element(i)%i_therm,j)                  &
                                               + mesh_element(i)%D_part(j)*mesh_element(i)%V_tri
        V_therm(mesh_element(i)%i_therm,j) = V_therm(mesh_element(i)%i_therm,j)                  &
                                               + mesh_element(i)%V_part(j)*mesh_element(i)%V_tri
        T_therm(mesh_element(i)%i_therm,j) = T_therm(mesh_element(i)%i_therm,j)                  &
                                               + mesh_element(i)%T_part(j)*mesh_element(i)%V_tri
     enddo
     D_therm(:,j) = D_therm(:,j)/vol_therm(:)
     V_therm(:,j) = V_therm(:,j)/vol_therm(:)
     T_therm(:,j) = T_therm(:,j)/vol_therm(:)
  enddo  

! gradients:  
  do i=1,ntri
     r_bc = (mesh_point(mesh_element(i)%i_knot(1))%rcoord + mesh_point(mesh_element(i)%i_knot(2))%rcoord  &
           + mesh_point(mesh_element(i)%i_knot(3))%rcoord)*one3rd
     z_bc = (mesh_point(mesh_element(i)%i_knot(1))%zcoord + mesh_point(mesh_element(i)%i_knot(2))%zcoord  &
           + mesh_point(mesh_element(i)%i_knot(3))%zcoord)*one3rd

     call grads_etc(i)
     psi_loc = 0.d0
     bmod_loc = 0.d0
     do m=1,legs
        psi_loc = psi_loc + mesh_point(mesh_element(i)%i_knot(m))%psi_pol
        bmod_loc = bmod_loc + mesh_point(mesh_element(i)%i_knot(m))%b_mod
     end do
     write(120,204)r_bc, z_bc, psi_loc*one3rd, bmod_loc*one3rd,                       &
          mesh_element(i)%dPsi_dR, mesh_element(i)%dPsi_dZ, mesh_element(i)%dBinv_dR,  &
          mesh_element(i)%dBinv_dZ, mesh_element(i)%ePhi_tri
  enddo  

! interpolate the initial guess profiles to the mesh knots:
  do i=2,npoint ! excluding O-point
     ePhi_knot = 0.d0
     D_knot = 0.d0
     V_knot = 0.d0
     T_knot = 0.d0
     do j=1,mesh_point(i)%n_owners
!        ePhi_knot = ePhi_knot + mesh_point(i)%weight_intp(j)                             &
!             * mesh_element(mesh_point(i)%i_owner_tri(j))%ePhi_tri
        do k=1,nsorts
           D_knot(k) = D_knot(k) + mesh_point(i)%weight_intp(j)                          &
                * D_therm(mesh_element(mesh_point(i)%i_owner_tri(j))%i_therm,k)
           V_knot(k) = V_knot(k) + mesh_point(i)%weight_intp(j)                          &
                * V_therm(mesh_element(mesh_point(i)%i_owner_tri(j))%i_therm,k)
           T_knot(k) = T_knot(k) + mesh_point(i)%weight_intp(j)                          &
                * T_therm(mesh_element(mesh_point(i)%i_owner_tri(j))%i_therm,k)
        enddo
     enddo
!     mesh_point(i)%PhiovB2 = ePhi_knot/(mesh_point(i)%b_mod)**2
     do k=1,nsorts
        mesh_point(i)%pl_parms_knot(1,k) = D_knot(k)
        mesh_point(i)%pl_parms_knot(2,k) = V_knot(k)
        mesh_point(i)%pl_parms_knot(3,k) = T_knot(k)
     enddo
  enddo
! O-point:
  ePhi_knot = 0.d0
  D_knot = 0.d0
  V_knot = 0.d0
  T_knot = 0.d0
!  do j=1, mesh_point(1)%n_owners
!     ePhi_knot = ePhi_knot + mesh_point(1)%weight_intp(j)                             &
!               * mesh_element(mesh_point(1)%i_owner_tri(j))%ePhi_tri
!  enddo
!  mesh_point(1)%PhiovB2 = ePhi_knot/(mesh_point(1)%b_mod)**2
!
  do j=1, nbou
     do k=1,nsorts
        D_knot(k) = D_knot(k) + inbou_list(j)%vol_norm                          &
             * D_therm(mesh_element(inbou_list(j)%i_tri_all)%i_therm,k)
        V_knot(k) = V_knot(k) + inbou_list(j)%vol_norm                          &
             * V_therm(mesh_element(inbou_list(j)%i_tri_all)%i_therm,k)
        T_knot(k) = T_knot(k) + inbou_list(j)%vol_norm                          &
             * T_therm(mesh_element(inbou_list(j)%i_tri_all)%i_therm,k)
     enddo
  enddo
  do k=1,nsorts
     mesh_point(1)%pl_parms_knot(1,k) = D_knot(k)
     mesh_point(1)%pl_parms_knot(2,k) = V_knot(k)
     mesh_point(1)%pl_parms_knot(3,k) = T_knot(k)
  enddo
!
!
!
! Find common edges in neighbouring triangles (entry edges)
!
  do i=1,ntri
    do j=1,3
      k=mesh_element(i)%neighbour(j)
      if(k.gt.0) then
        j1=mesh_element(i)%i_knot(j)
        j2=mesh_element(i)%i_knot(modulo(j,3)+1)
        j3=0
!
        do m=1,3
          if(mesh_element(k)%i_knot(m).eq.j1) then
            j3=m
            exit
          endif
        enddo
!
        if(j3.gt.0) then
          j1=j3
        else
          print *,'common point not found'
        endif
!
        do m=1,3
          if(mesh_element(k)%i_knot(m).eq.j2) then
            j3=m
            exit
          endif
        enddo
!
        if(j3.gt.0) then
          j2=j3
        else
          print *,'common point not found'
        endif
!
        j3=j1+j2
        m=0
!
        select case (j3)
        case(3)  !3=1+2
          m=1
        case(5)  !5=2+3
          m=2
        case(4)  !4=3+1
          m=3
        end select
!
        if(m.gt.0) then
          mesh_element(i)%neighbour_edge(j)=m
        else
          print *,'no edge selected'
        endif
      else 
        mesh_element(i)%neighbour_edge(j)=0
      endif
    enddo
  enddo
!
! Test of edges:
  do i=1,ntri
    do j=1,3
      k=mesh_element(i)%neighbour(j)
      m=mesh_element(i)%neighbour_edge(j)
      if(k.gt.0) then
        if(mesh_element(k)%neighbour(m).ne.i) then
          print *,'wrong edges in triangle ',i
        endif
      endif
    enddo
  enddo
!
!
!
  

  open(1,file='points.dat',form='unformatted')
  write(1) npoint
  write(1) mesh_point
  close(1)

  bphicovar = btf*rtf
! Write value for linear interpolation of bphicovar
  do i=1,ntri
     call gen_dbphicovdpsi(i)
  end do
  
  time_l(:) = 6.0d-1
  fulltime = 0.d0
  open(1,file='triangles.dat',form='unformatted')
  write(1) ntri
  write(1) mesh_element
  write(1) bphicovar, fulltime, time_l
  close(1)
  open(1,file='thermostat.dat',form='unformatted')
  write(1) n_therm
  write(1) D_therm
  write(1) V_therm
  write(1) T_therm
  write(1) vol_therm
  close(1)

  allocate(flux_j(nbou_qq,nsorts), flux_q(nbou_qq,nsorts), tor_mom(nbou_qq,nsorts))
  flux_j = 1.d0
  flux_q = 1.d0
  tor_mom = 1.d0
  open(1,file='wall_loads.dat',form='unformatted')
  write(1) nbou_qq
  write(1) flux_j
  write(1) flux_q
  write(1) tor_mom
  close(1)

  open(1,file='boundary.dat',form='unformatted')
  write(1) nbou 
  write(1) inbou_list
  close(1)

! to start particles according to given profiles (in the core only!):
  allocate(vol_core(nr_core, 2*nt_core), ind_core_tri(nr_core, 2*nt_core)) 
  vol_core(:,:) = 0.d0
  ind_core_tri(:,:) = 0
  do j=1,nr_core-1
     i_tht = 0
     do i=1, ntri
        if(mesh_element(i)%iy_qq .eq. j) then
           i_tht = i_tht + 1
           vol_core(j,i_tht) = mesh_element(i)%V_tri 
           ind_core_tri(j,i_tht) = i
        endif
     enddo
  enddo
  open(1,file='source.dat',form='unformatted')
  write(1) nr_core, nt_core
  write(1) vol_core
  write(1) ind_core_tri
  close(1)

!!$2000  format(a12,f19.5)
2001  format(a17)
2002  format(a28)
204 format(1000(ES15.8E2,1X))

end program readcarre_m
!------------------------------------------------------------------------------------
subroutine grads_etc(icell)
  use from_nrtype
  use mesh_mod,    only : mesh_point, mesh_element
  implicit none
  integer(i4b), intent(in) :: icell
  real(dp) :: r2r1, r3r1, z2z1, z3z1, f2f1, f3f1, Psi_1, Binv_1, PhiovB2_1

  r2r1 = mesh_point(mesh_element(icell)%i_knot(2))%rcoord -                      & 
         mesh_point(mesh_element(icell)%i_knot(1))%rcoord

  r3r1 = mesh_point(mesh_element(icell)%i_knot(3))%rcoord -                      & 
         mesh_point(mesh_element(icell)%i_knot(1))%rcoord

  z2z1 = mesh_point(mesh_element(icell)%i_knot(2))%zcoord -                      & 
         mesh_point(mesh_element(icell)%i_knot(1))%zcoord

  z3z1 = mesh_point(mesh_element(icell)%i_knot(3))%zcoord -                      & 
         mesh_point(mesh_element(icell)%i_knot(1))%zcoord

  Binv_1 = 1.0_dp/mesh_point(mesh_element(icell)%i_knot(1))%b_mod
  f2f1 =  1.0_dp/mesh_point(mesh_element(icell)%i_knot(2))%b_mod - Binv_1
  f3f1 =  1.0_dp/mesh_point(mesh_element(icell)%i_knot(3))%b_mod - Binv_1
  mesh_element(icell)%dBinv_dR = (z3z1*f2f1 - z2z1*f3f1)*mesh_element(icell)%oneoverD3
  mesh_element(icell)%dBinv_dZ = (r2r1*f3f1 - r3r1*f2f1)*mesh_element(icell)%oneoverD3

  Psi_1 = mesh_point(mesh_element(icell)%i_knot(1))%psi_pol
  f2f1 =  mesh_point(mesh_element(icell)%i_knot(2))%psi_pol - Psi_1
  f3f1 =  mesh_point(mesh_element(icell)%i_knot(3))%psi_pol - Psi_1
  mesh_element(icell)%dPsi_dR = (z3z1*f2f1 - z2z1*f3f1)*mesh_element(icell)%oneoverD3
  mesh_element(icell)%dPsi_dZ = (r2r1*f3f1 - r3r1*f2f1)*mesh_element(icell)%oneoverD3

  PhiovB2_1 = mesh_point(mesh_element(icell)%i_knot(1))%PhiovB2
  f2f1 =  mesh_point(mesh_element(icell)%i_knot(2))%PhiovB2 - PhiovB2_1
  f3f1 =  mesh_point(mesh_element(icell)%i_knot(3))%PhiovB2 - PhiovB2_1
  mesh_element(icell)%dPhiovB2_dR = (z3z1*f2f1 - z2z1*f3f1)*mesh_element(icell)%oneoverD3
  mesh_element(icell)%dPhiovB2_dZ = (r2r1*f3f1 - r3r1*f2f1)*mesh_element(icell)%oneoverD3

  return
end subroutine grads_etc

subroutine gen_dbphicovdpsi(ktri)
  ! generate dbphicov/dpsi in a triangle and its neighbor of higher index if it exists
  ! not guaranteed to work outside core plasma
  
  use mesh_mod
  implicit none
  integer, intent(in) :: ktri
  integer :: kno(3), kno2(3)! indices of knots for current and adjacent triangle
  integer :: k, kfar        ! furthest point in psi
  real(dp) :: dtemp, dclose ! closest distance in psi

  integer :: ktri2
  integer :: knear2(2), kfar2, cmatch  ! for matching points
  real(dp) :: absgradpsi1, absgradpsi2 ! absolute values of grad psi
  real(dp) :: deltab, deltapsi1, deltapsi2

  kno(:) = mesh_element(ktri)%i_knot(:)
  
  ! find point for which edges are not on the same flux surface
  do k=1,3
     dtemp = abs(mesh_point(kno(mod(k,3)+1))%psi_pol-mesh_point(kno(k))%psi_pol) 
     if ((dtemp<dclose) .or. (k==1)) then
        dclose = dtemp
        kfar = mod(k+1,3)+1
     end if
  end do

  ! find adjacent triangle of higher index
  cmatch = 0
  do ktri2=ktri+1,ntri
     cmatch = 0
     do k = 1,3
        if (kno(mod(kfar,3)+1) == mesh_element(ktri2)%i_knot(k) .or. &
             kno(mod(kfar+1,3)+1) == mesh_element(ktri2)%i_knot(k)) then
           cmatch=cmatch+1
           knear2(cmatch) = k
        end if
        if (cmatch >= 2) exit
     end do
     if (cmatch >= 2) exit
  end do

  if (cmatch < 2) return

  kno2(:) = mesh_element(ktri2)%i_knot(:)  

  ! find far knot of adjacent triangle
  do kfar2=1,3
     if(.not.(any(knear2==kfar2))) exit
  end do  

  absgradpsi1 = sqrt(mesh_element(ktri)%dPsi_dR**2 + mesh_element(ktri)%dPsi_dZ**2)
  absgradpsi2 = sqrt(mesh_element(ktri2)%dPsi_dR**2 + mesh_element(ktri2)%dPsi_dZ**2)
  
  deltab = bphicovar*(absgradpsi1-absgradpsi2)/(absgradpsi1+absgradpsi2)  
  deltapsi1 = mesh_point(kno(mod(kfar,3)+1))%psi_pol - mesh_point(kno(kfar))%psi_pol
  deltapsi2 = mesh_point(kno2(kfar2))%psi_pol - mesh_point(kno2(mod(kfar2,3)+1))%psi_pol

!print *, ktri, deltab/bphicovar, absgradpsi1, absgradpsi2
  
  mesh_element(ktri)%knot_h  = kfar
  mesh_element(ktri2)%knot_h = kfar2
  mesh_element(ktri)%dbphicovdpsi = deltab/deltapsi1
  mesh_element(ktri2)%dbphicovdpsi = deltab/deltapsi2
!print *, ktri, mesh_element(ktri)%dbphicovdpsi
!print *, ktri2, mesh_element(ktri2)%dbphicovdpsi

end subroutine gen_dbphicovdpsi

