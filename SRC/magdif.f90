module magdif
  use from_nrtype
  use constants, only: clight, pi
  use sparse_mod, only: remap_rc
  use mesh_mod, only: knot, triangle, mesh_point, mesh_element, mesh_element_rmp, ntri,&
       bphicovar
  use gettormode_mod, only: B_Rn,B_Zn,n_tor_out
  implicit none


  integer, parameter :: nt_core=200
  !integer, parameter :: nr_core=75
  !integer, parameter :: nr_max=74
  integer, parameter :: nr_max=74
  

  integer :: n ! harmonic index of perturbation
  real(dp), allocatable :: pres0(:)
  real(dp), allocatable :: psi(:)
  real(dp), allocatable :: q(:)
  real(dp), allocatable :: B2avg(:)
  complex(dp), allocatable :: presn(:)
  complex(dp), allocatable :: curr(:,:)
  complex(dp), allocatable :: Bnflux(:,:), Bnp(:) 

  complex(dp), parameter :: imun = (0d0,1d0)
  
  contains
  
  subroutine assemble_system(nrow, a, b, c, q, d, du, Mq)
    integer, intent(in) :: nrow                ! number of system rows
    complex(dp), intent(in) :: a(:), b(:), c(:) ! system coefficients
    complex(dp), intent(in) :: q(:)             ! right-hand side without mass matrix applied
    complex(dp), intent(out) :: d(nrow)         ! diagonal of stiffness matrix
    complex(dp), intent(out) :: du(nrow)        ! superdiagonal of stiffness matrix + A(n,1)
    complex(dp), intent(out) :: Mq(:)           ! right-hand side with mass matrix applied

    integer :: k

    do k=1,nrow-1
       d(k) = -a(k) + b(k)/2d0
       du(k) = a(k) + b(k)/2d0
       Mq(k) = (c(k+1)*q(k+1)+c(k)*q(k))/2d0
    end do
    
    d(nrow) = -a(nrow) + b(nrow)/2d0
    du(nrow) = a(1) + b(1)/2d0
    Mq(nrow) = (c(1)*q(1)+c(nrow)*q(nrow))/2d0
  end subroutine assemble_system
  
subroutine assemble_system2(nrow, a1, a2, b1, b2, d, du)
    integer, intent(in) :: nrow                ! number of system rows
    complex(dp), intent(in) :: a1(:), a2(:), b1(:), b2(:) ! system coefficients
    complex(dp), intent(out) :: d(nrow)         ! diagonal of stiffness matrix
    complex(dp), intent(out) :: du(nrow)        ! superdiagonal of stiffness matrix + A(n,1)
    
    d = a1 + b1/2d0
    du(1:(nrow-1)) = a2(1:(nrow-1)) + b2(1:(nrow-1))/2d0
    du(nrow) = a2(1) + b2(1)/2d0
  end subroutine assemble_system2
  
  subroutine assemble_sparse(nrow, d, du, nz, irow, icol, aval)
    integer, intent(in)  :: nrow                          ! number of system rows
    complex(dp), intent(in)  :: d(nrow)                   ! diagnonal
    complex(dp), intent(in)  :: du(nrow)                  ! upper diagonal
    integer, intent(out) :: nz                            ! number of system rows
    integer, intent(out) :: irow(2*nrow), icol(2*nrow)    ! matrix index representation
    complex(dp), intent(out) :: aval(2*nrow)              ! matrix values
    
    integer :: k
    
    irow(1) = 1
    icol(1) = 1
    aval(1) = d(1)
 
    irow(2) = nrow
    icol(2) = 1
    aval(2) = du(nrow)

    do k = 2,nrow
       ! off-diagonal 
       irow(2*k-1) = k-1
       icol(2*k-1) = k
       aval(2*k-1) = du(k-1)

       ! diagonal 
       irow(2*k) = k
       icol(2*k) = k
       aval(2*k) = d(k)
    end do

    nz = 2*nrow
  end subroutine assemble_sparse

  subroutine solve_full(n, d, du, alpha, Mq)
    integer, intent(in) :: n                   ! number of system rows
    complex(8), intent(in) :: d(:)             ! diagonal of stiffness matrix
    complex(8), intent(in) :: du(:)            ! superdiagonal of stiffness matrix
    complex(8), intent(in) :: alpha            ! single entry in stiffness matrix to be periodic
    complex(8), intent(in) :: Mq(:)            ! right-hand side with mass matrix applied
    ! TODO: implement
  end subroutine solve_full
  
  subroutine solve_cycl_tridiag(n, d, du, alpha, Mq)
    integer, intent(in) :: n                   ! number of system rows
    complex(8), intent(in) :: d(:)             ! diagonal of stiffness matrix
    complex(8), intent(in) :: du(:)            ! superdiagonal of stiffness matrix
    complex(8), intent(in) :: alpha            ! single entry in stiffness matrix to be periodic
    complex(8), intent(in) :: Mq(:)            ! right-hand side with mass matrix applied
    ! TODO: implement
  end subroutine solve_cycl_tridiag

  subroutine solve_sparse(n, d, du, alpha, Mq)
    integer, intent(in) :: n                   ! number of system rows
    complex(8), intent(in) :: d(:)             ! diagonal of stiffness matrix
    complex(8), intent(in) :: du(:)            ! superdiagonal of stiffness matrix
    complex(8), intent(in) :: alpha            ! single entry in stiffness matrix to be periodic
    complex(8), intent(in) :: Mq(:)            ! right-hand side with mass matrix applied
    ! TODO: implement
  end subroutine solve_sparse

  subroutine common_triangles(knot1, knot2, common_tri)
    type(knot), intent(in) :: knot1, knot2
    integer, intent(out) :: common_tri(2)

    integer :: k, l, kcom
    kcom = 0
    do k=1,knot1%n_owners
       do l=1,knot2%n_owners
          if (knot1%i_owner_tri(k)==knot2%i_owner_tri(l)) then
             kcom = kcom+1
             if (kcom>2) stop "Error: more than two common triangles for knots"
             common_tri(kcom) = knot1%i_owner_tri(k)
          end if
       end do
    end do
  end subroutine common_triangles

  subroutine get_edge_coord(knots, r, z, lr, lz)

    type(knot), intent(in) :: knots(3)
    real(dp), dimension(3) :: r, z, lr, lz

    r(1) = (knots(1)%rcoord + knots(2)%rcoord)/2d0
    r(2) = (knots(2)%rcoord + knots(3)%rcoord)/2d0
    r(3) = (knots(3)%rcoord + knots(1)%rcoord)/2d0
    
    z(1) = (knots(1)%zcoord + knots(2)%zcoord)/2d0
    z(2) = (knots(2)%zcoord + knots(3)%zcoord)/2d0
    z(3) = (knots(3)%zcoord + knots(1)%zcoord)/2d0
    
    lr(1) = knots(2)%rcoord - knots(1)%rcoord
    lr(2) = knots(3)%rcoord - knots(2)%rcoord
    lr(3) = knots(1)%rcoord - knots(3)%rcoord
    
    lz(1) = knots(2)%zcoord - knots(1)%zcoord
    lz(2) = knots(3)%zcoord - knots(2)%zcoord
    lz(3) = knots(1)%zcoord - knots(3)%zcoord
  end subroutine get_edge_coord
  
  subroutine assign_curr3(ktri, p0pr, B2a)
    
    integer, intent(in) :: ktri
    real(dp), intent(in) :: p0pr, B2a
    
    type(triangle) :: elem
    type(knot) :: knots(3)
    
    real(dp), dimension(3) :: r,z,lr,lz
    real(dp) :: Br,Bp,Bz,dBrdR,dBrdp,dBrdZ         &
         ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ
    real(dp) :: B2, j0parB
    
    integer :: k, indl3

    elem = mesh_element(ktri)
    knots = mesh_point(elem%i_knot)
    
    call get_edge_coord(knots, r, z, lr, lz)
    indl3 = mod(elem%knot_h,3)+1

    call field(r(indl3),0d0,z(indl3),Br,Bp,Bz,dBrdR,&
         dBrdp,dBrdZ,dBpdR,dBpdp,dBpdZ,&
         dBzdR,dBzdp,dBzdZ)
    
    curr(ktri, indl3) = clight/(Br**2+Bp**2+Bz**2)*(&
         Bp*r(indl3)*(presn(elem%i_knot(mod(indl3,3)+1))-presn(elem%i_knot(indl3)))&
         -imun*n*(presn(elem%i_knot(mod(indl3,3)+1))+presn(elem%i_knot(indl3)))/2d0*&
         (Br*lr(indl3)+Bz*lz(indl3)))
    
    B2 = Br**2+Bp**2+Bz**2
    j0parB = -clight*p0pr*r(indl3)*(1d0/B2 - 1d0/B2a)
    curr(ktri, indl3) = curr(ktri, indl3) - j0parB*Bnflux(ktri,indl3)/B2

  end subroutine assign_curr3
    
    
  subroutine assign_currents(ktri, I1, I2)
    ! assigns currents in the correct direction
    
    integer, intent(in) :: ktri
    complex(dp), intent(in) :: I1, I2

    type(triangle) :: elem
    type(knot) :: knots(3)

    integer :: indl3

    elem = mesh_element(ktri)
    knots = mesh_point(elem%i_knot)    
    indl3 = mod(elem%knot_h,3)+1

    if (minval(elem%i_knot) == elem%i_knot(elem%knot_h)) then
       ! type 1 triangle
       curr(ktri, elem%knot_h) = I2
       curr(ktri, mod(indl3,3)+1) = -I1
    else       
       ! type 2 triangle
       curr(ktri, elem%knot_h) = -I1
       curr(ktri, mod(indl3,3)+1) = I2
    end if
  end subroutine assign_currents

  function dpsidr(ktri)
    integer, intent(in) :: ktri
    real(8) dpsidr

    type(triangle) :: elem, elem2
    type(knot) :: knots(3), knots2(3)
    real(8) :: rdist
    integer :: neigh
    
    elem = mesh_element(ktri)
    knots = mesh_point(elem%i_knot)
    
    neigh = elem%neighbour(mod(elem%knot_h,3)+1)
    
    if (neigh<0 .OR. elem%knot_h==0) then
       dpsidr = 0
       return
    end if
    
    elem2 = mesh_element(neigh)
    knots2 = mesh_point(elem2%i_knot)
    rdist = sqrt((knots2(elem2%knot_h)%rcoord-knots(elem%knot_h)%rcoord)**2+&
         (knots2(elem2%knot_h)%zcoord-knots(elem%knot_h)%zcoord)**2)
    
    if (minval(elem%i_knot) == elem%i_knot(elem%knot_h)) then
       ! type 1 triangle
       dpsidr = (knots2(elem2%knot_h)%psi_pol - knots(elem%knot_h)%psi_pol)/rdist
    else       
       ! type 2 triangle
       dpsidr = (knots(elem%knot_h)%psi_pol - knots2(elem2%knot_h)%psi_pol)/rdist
    end if
  end function dpsidr

  function get_source(ktri,kl)
    integer, intent(in) :: ktri,kl
    real(dp) :: p0pr
    complex(dp) get_source

    type(triangle) :: elem
    type(knot), dimension(3) :: knots, knots_s

    integer, dimension(3) :: iknot_s
    complex(dp) pres_source, cur_source
    real(dp), dimension(3) :: r, z, lr, lz
    real(dp), dimension(3) :: Br,Bp,Bz,dBrdR,dBrdp,dBrdZ         &
         ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ
    integer :: k

    p0pr = (pres0(kl+1)-pres0(kl))/(psi(kl+1)-psi(kl))

    elem = mesh_element(ktri)
    knots = mesh_point(elem%i_knot)
    
    iknot_s(1) = elem%i_knot(mod(elem%knot_h+1,3)+1)
    iknot_s(2) = elem%i_knot(elem%knot_h)
    iknot_s(3) = elem%i_knot(mod(elem%knot_h,3)+1)
    knots_s = mesh_point(iknot_s)

    call get_edge_coord(knots_s, r, z, lr, lz)

    do k=1,3
       call field(r(k),0d0,z(k),Br(k),Bp(k),Bz(k),dBrdR(k),&
            dBrdp(k),dBrdZ(k),dBpdR(k),dBpdp(k),dBpdZ(k),&
            dBzdR(k),dBzdp(k),dBzdZ(k))
    end do
    
    pres_source = r(1)/2d0*(presn(iknot_s(1))-presn(iknot_s(2)))/&
         (knots_s(1)%psi_pol - knots_s(2)%psi_pol)&
         + r(2)/2d0*(presn(iknot_s(3))-presn(iknot_s(2)))/&
         (knots_s(3)%psi_pol - knots_s(2)%psi_pol)

    pres_source = pres_source - Bnp(ktri)*2d0/(Bp(1)+Bp(2))*p0pr
    pres_source = clight*pres_source

    cur_source = j0phi(r(1),0d0,z(1),p0pr,B2avg(kl))*(Bnp(ktri)/Bp(1)&
         + Bnflux(ktri,mod(elem%knot_h+1,3)+1)/(knots_s(2)%psi_pol-knots_s(1)%psi_pol))&
         + j0phi(r(2),0d0,z(2),p0pr,B2avg(kl))*(Bnp(ktri)/Bp(2)&
         + Bnflux(ktri,elem%knot_h)/(knots_s(3)%psi_pol-knots_s(2)%psi_pol))
!    print *, Bnp(ktri)/Bp(1)
!pres_source = 0d0    
!cur_source = 0d0
    
    get_source = pres_source + cur_source
    get_source = -(0d0,1d0)*n*get_source*elem%det_3/2d0 ! weigh with surface area
!    get_source = 0d0
  end function get_source
  
  function j0phi(r,p,z,p0pr,B2a)
    real(dp), intent(in) :: r,p,z,p0pr,B2a
    real(dp) :: j0phi, j0phipar, j0phiperp
    real(dp) :: Br,Bp,Bz,dBrdR,dBrdp,dBrdZ         &
         ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ
    real(dp) :: B2, Bpol2, Btor2

    call field(r,p,z,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ          &
         ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)

    B2 = (Br**2+Bp**2+Bz**2)
    Bpol2 = Br**2+Bz**2
    Btor2 = Bp**2
    
    j0phipar = clight*p0pr*Btor2*(1d0/B2 - 1d0/B2a)
    j0phiperp = clight*p0pr*Bpol2/B2
!print *, j0phipar, j0phiperp, Btor2, Bpol2
    j0phi = r*(j0phipar + j0phiperp)
  end function j0phi

  subroutine gen_B2avg(kl)
    integer, intent(in) :: kl

    real(dp) :: r,p,z,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ         &
         ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ
    integer :: kt
    type(triangle) :: elem
    type(knot), dimension(3) :: knots

    B2avg(kl) = 0d0
    
    if(kl==1) then
       do kt=1,nt_core
          elem = mesh_element(kt)
          knots = mesh_point(elem%i_knot)
          
          r = (sum(knots%rcoord)+knots(elem%knot_h)%rcoord)/4d0
          z = (sum(knots%zcoord)+knots(elem%knot_h)%zcoord)/4d0
          call field(r,p,z,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ          &
               ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
          B2avg = B2avg + Br**2 + Bp**2 + Bz**2
       end do !kt
       B2avg = B2avg/nt_core
    else
       do kt=1,2*nt_core
          elem = mesh_element(nt_core+(kl-2)*2*nt_core+kt)
          knots = mesh_point(elem%i_knot)
          
          r = (sum(knots%rcoord)+knots(elem%knot_h)%rcoord)/4d0
          z = (sum(knots%zcoord)+knots(elem%knot_h)%zcoord)/4d0
          call field(r,p,z,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ          &
               ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
          B2avg = B2avg + Br**2 + Bp**2 + Bz**2
       end do !kt
       B2avg = B2avg/(2*nt_core)
    end if
  end subroutine gen_B2avg

  subroutine set_Bnflux_notheta(ktri)
    ! generates Bnflux and Bnp without theta component
    ! based on mesh_element_rmp(k)%bnorm_vac
    integer :: ktri
    type(triangle) :: elem
    type(knot), dimension(3) :: knots, knots_s
    integer, dimension(3) :: iknot_s
    complex(dp) :: Bnflux3
    real(dp), dimension(3) :: r, z, lr, lz
    integer :: common_tri(2)
    
    elem = mesh_element(ktri)
    knots = mesh_point(elem%i_knot)

    if(elem%knot_h==0) then
       Bnflux(ktri,:) = 0d0
       Bnp(ktri) = 0d0
       return
    end if
    
    iknot_s(1) = elem%i_knot(mod(elem%knot_h+1,3)+1)
    iknot_s(2) = elem%i_knot(elem%knot_h)
    iknot_s(3) = elem%i_knot(mod(elem%knot_h,3)+1)

    knots_s = mesh_point(iknot_s)
    
    call get_edge_coord(knots_s, r, z, lr, lz)

    call common_triangles(knots_s(1), knots_s(3), common_tri)
    
    Bnflux3 = r(3)*(knots_s(1)%b_mod+knots_s(3)%b_mod)/2d0&
         *1d0/(knots_s(1)%psi_pol-knots_s(2)%psi_pol)*&
         (mesh_element_rmp(common_tri(1))%bnorm_vac+&
         mesh_element_rmp(common_tri(2))%bnorm_vac)/2d0

!    Bnflux(ktri,mod(elem%knot_h+1,3)+1) = r(1)/r(3)*Bnflux3*(lr(1)*lr(3)+lz(1)*lz(3))/sqrt(lr(3)**2+lz(3)**2)
!    Bnflux(ktri,elem%knot_h) = r(2)/r(3)*Bnflux3*(lr(2)*lr(3)+lz(2)*lz(3))/sqrt(lr(3)**2+lz(3)**2)
    Bnflux(ktri,mod(elem%knot_h,3)+1) = Bnflux3
!    Bnflux(ktri,mod(elem%knot_h,3)+1) = 1d0
    Bnflux(ktri,mod(elem%knot_h+1,3)+1) = 0d0
    Bnflux(ktri,elem%knot_h) = 0d0
!if (ktri<200) then
!   print *, knots_s(1)%psi_pol-knots_s(2)%psi_pol
!   print *, knots_s(3)%psi_pol-knots_s(2)%psi_pol
!   print *, knots_s(1)%psi_pol-knots_s(3)%psi_pol
!   print *, '...'
!end if
    
    Bnp(ktri) = -1d0/((0d0,1d0)*n)*sum(Bnflux(ktri,:))/(elem%det_3/2d0)
    
  end subroutine set_Bnflux_notheta
  
  subroutine init_safety_factor
    integer :: kl, kt, k, nt_loop
    type(triangle) :: elem
    type(knot), dimension(3) :: knots
    real(dp) :: r, z, Br,Bp,Bz,dBrdR,dBrdp,dBrdZ         &
         ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ, rinv
    
    q = 0d0
    do kl = 1,nr_max
       do kt=1,nt_loop
          if(kl==1) then
             elem = mesh_element(kt)
          else
             elem = mesh_element(nt_core+(kl-2)*2*nt_core+kt)
          end if
          knots = mesh_point(elem%i_knot)
          r = sum(knots(:)%rcoord)/3d0
          z = sum(knots(:)%zcoord)/3d0
          rinv = sum(1d0/knots(:)%rcoord**2)/3d0
          call field(r,0d0,z,Br,Bp,Bz,dBrdR,&
               dBrdp,dBrdZ,dBpdR,dBpdp,dBpdZ,&
               dBzdR,dBzdp,dBzdZ)
          q(kl) = q(kl) + Bp/(2d0*pi)*elem%det_3/2d0 !bphicovar*rinv*elem%V_tri
       end do
       q(kl) = q(kl)/(psi(kl+1)-psi(kl))
       if(kl==1) nt_loop = 2*nt_core
    end do
  end subroutine init_safety_factor
end module magdif

