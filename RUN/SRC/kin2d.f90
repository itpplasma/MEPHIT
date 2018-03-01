program k2d
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
                            n_tormode, dt0_factor
  use ranstate,      only : idum, idum2, iv, iy
  use hujnja, only  : numstepmax, rel_therm_max !, v3_sum_in, v3_sum_out, v3_in, v3_out,       &
                      !part_num
!, ubila_rr, dobavila_rr, istochnik, stenka, proebal
  use for_mpi, only : mype, npes, mpi_p_root, time_local, Dm_local, Dp_local, Vp_local,        &
                      Tp_local, ephi_local, mype_sym, evolname,              &
                      flux_j, flux_q, tor_mom, Vt2_local, ephi_sheath, dens_phi
  
  use arnoldi_mod, only: ieigen, ngrow, eigvecs, tol, arnoldi_ierr => ierr
#ifdef PARALLEL
  use mpi
!  include 'mpif.h'
#endif
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
  integer(i4b), parameter :: max_openmp_threads=256
  double precision, dimension(0:5), save :: seed, cg_rng, saverandom
  double precision, dimension(0:5,0:max_openmp_threads-1) :: cg_rng_sav
! dispersions
  double precision, dimension(:,:), allocatable :: D_Dm, D_Dp, D_Vp, D_Tp, D_Dt, D_Vt, D_Tt 
  double precision, dimension(:), allocatable :: D_ephi
  double precision, dimension(:,:), allocatable :: D_j, D_q, D_tm
  double precision, dimension(nsorts) :: D_time_l, part_sink, temp_out
  double precision :: dummy_re,dummy_im
  double complex   :: bnorm_vac,cur1,cur2
  double precision :: eps_urelax,delta,dens_reduce_factor
  double precision, dimension(3) :: R_vert,Z_vert
  double complex, dimension(:), allocatable :: cur_R,cur_Z,bnorm_plas,bnorm_plas_prev
  
  double precision :: R_ax=163.d0, Z_ax=15.d0

  integer(i4b), parameter :: nritz = 30  ! for Arnoldi iterations
  double complex, dimension(nritz) :: ritznum
  integer :: info
  double complex, dimension(:), allocatable :: coefren
  integer,        dimension(:), allocatable :: ipiv
  double complex, dimension(:,:), allocatable :: amat, bvec
  
  logical :: do_arnoldi, precond
  
  double complex, dimension(:), allocatable :: cdummy8_1, cdummy8_2
  character(len=8) :: temp
  
  external  ran_u_01
#ifdef PARALLEL
  call MPI_INIT(ierr)
  if(ierr .ne. MPI_SUCCESS) then 
     write(*,*)'MPI init. err.'
     stop 1
  endif
  call MPI_COMM_SIZE(MPI_COMM_WORLD,npes,ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,mype,ierr)
#else
  npes = 1
  mype = 0
#endif

  oneovernp = 1.d0/dble(npes)
  if(max_openmp_threads .lt. npes) stop 'too few threads'

  allocate(cur_R(nsorts),cur_Z(nsorts))
  
  n_tormode=2
  dt0=2d-4
  dt0_factor(1)=1.d0
  dt0_factor(2)=sqrt(amass(2)/amass(1))
  if (mype .eq. mpi_p_root) print *,'dt0 factors: ',dt0_factor
  eps_urelax=1.d0
  dens_reduce_factor=1.d0
     
  if (mype .eq. mpi_p_root) then
! read everything
     open(1,file='kin2d.inp')
     read (1,*) bmod_ref                   ! reference field, G, for
                                           ! Boozer ! $B_{00}$
     read (1,*) nptrace, source0, T_source ! number of particles to trace,
                                           ! source intencity  
     read (1,*) relax_tau                  ! underrelaxation factor 
     read (1,*) wcv(1), wct(1)             ! weights of cells, electrons
     read (1,*) wcv(2), wct(2)             ! weights of cells, ions
     read (1,*) D_perp                     ! anomalous diff.coeff., cm^2/s
     read (1,*) n_evol, n_prfs             ! number of steps between outputs
                                           ! (average vals. and profiles)
     read (1,*) npgroup                    ! number of particles per macrostep
     read (1,*) restart                    ! restart from last output
     close(1)

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
     allocate(mesh_element_rmp(ntri))
!
!
!
     open(1,file='RMP_EXCHANGE/hpsi_vac.dat')
     open(2,file='RMP_EXCHANGE/wsrr_rec_fix.dat')              !<= SRR
     do i=1,ntri
        read (1,*) dummy_re,dummy_im
        mesh_element_rmp(i)%bnorm_vac=cmplx(dummy_re,dummy_im)
mesh_element_rmp(i)%bnorm_vac = 1d0 ! check with constant source
        read (2,*) dummy_re                                    !<= SRR
        mesh_element_rmp(i)%wsrr_rec_fix=dummy_re              !<= SRR
     end do
     close(1)
     close(2)                                                  !<= SRR

!print *,sum(flux_j),sum(flux_q),sum(tor_mom)
!stop
     open(1,file='START_PRFS/thermostat.dat',form='unformatted')
     read(1) n_therm
     allocate(D_therm(n_therm,nsorts),V_therm(n_therm,nsorts),T_therm(n_therm,nsorts))
     allocate(vol_therm(n_therm))
     read(1) D_therm
     read(1) V_therm
     read(1) T_therm
     read(1) vol_therm
     close(1)

!!$allocate(v3_sum_in(n_therm), v3_sum_out(n_therm), v3_in(n_therm), v3_out(n_therm),&
!!$         part_num(n_therm))
!!$v3_sum_in(:) = 0.d0
!!$v3_sum_out(:) = 0.d0
!!$part_num = 0.d0

     open(1,file='START_PRFS/boundary.dat',form='unformatted')
     read(1) ntri_inbou
     allocate(inbou_list(ntri_inbou))
     read(1) inbou_list
!!$     read(1) i_pfz, i_dpr, i_sol, i_dpl, i_inb
!!$     nboucell = i_inb(2) - i_pfz(1) + 1
!!$     allocate(R_bou(nboucell), Z_bou(nboucell), S_bou(nboucell))
!!$     read(1) R_bou, Z_bou, S_bou
     close(1)

     open(1,file='START_PRFS/wall_loads.dat',form='unformatted')
     read(1) nboucell
     allocate(flux_j(nboucell,nsorts), flux_q(nboucell,nsorts), tor_mom(nboucell,nsorts))
     read(1) flux_j
     read(1) flux_q
     read(1) tor_mom
     close(1)

     open(1010,file='RESULTS/errors.dat',status='replace')
     close(1010)
!     open(1001,file='RESULTS/evol.dat',status='replace')
!     close(1001)
     open(1004,file='RESULTS/time_evol.dat',status='unknown')
     close(1004)
!     open(1005,file='RESULTS/bookkeeping.dat',status='unknown')
     !     close(1005)
     open(3210,file='convergence.dat')
     write(3210,*)
     close(3210)
  endif ! (mype .eq. mpi_p_root)

#ifdef PARALLEL
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  call broadcast_all

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif

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

  allocate(Dm_local(ntri,nsorts), Dp_local(ntri,nsorts), Vp_local(ntri,nsorts),&
           Tp_local(ntri,nsorts), ephi_local(ntri), Vt2_local(ntri,nsorts) )
  allocate(D_Dt(n_therm,nsorts), D_Vt(n_therm,nsorts), D_Tt(n_therm,nsorts))
  allocate(D_Dm(ntri,nsorts), D_Dp(ntri,nsorts), D_Vp(ntri,nsorts), &
           D_Tp(ntri,nsorts), D_ephi(ntri))
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
! random numbers:
  if(.not. restart) then
     seed =  17491.d0 ! 12345.d0 !
     do i=1, max_openmp_threads
        call RngStream(seed, cg_rng)
        cg_rng_sav(:,i-1) = cg_rng
     enddo
     cg_rng = cg_rng_sav(:,mype+1)
  else
     if(npes .ne. 1) stop 'wrong restart of RNG'
     open(1,file='START_PRFS/restart_rng.int',form='unformatted')
     read(1) cg_rng
     close(1)
  endif

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

  allocate(sum_tr(ntri,nsorts), time_coll(ntri,nsorts), v_tr(ntri,nsorts),&
           en_mean_tr(ntri,nsorts), t4Dm(ntri,nsorts),   &
           j_tr(nboucell,nsorts), q_tr(nboucell,nsorts), &
           tm_tr(nboucell,nsorts), vt2_tr(ntri,nsorts))

  !

  

do_arnoldi = .False.

!if (mype .ne. mpi_p_root) then   
!   print *, "Slave process running inner loop: ", mype
!end if

!if (mype .eq. mpi_p_root) then
!  print *, "Master process running outer loop: ", mype
!end if
  
  allocate(bnorm_plas_prev(ntri))
  allocate(bnorm_plas(ntri))
  bnorm_plas_prev=(0.d0,0.d0)
  bnorm_plas=(0.d0,0.d0)
  
  if(do_arnoldi) then
     !  Find eigenvalues (Arnoldi)
     if (mype .eq. mpi_p_root) print *, "Doing ", nritz, "Arnoldi iterations..."
     iptrace=1
     ieigen=1
     tol = 0.7d0
     precond = .True.
     call arnoldi(ntri, nritz, ritznum, next_iteration)    
     precond = .False.
     if(ngrow<1) then
        do_arnoldi = .False.
     else if (mype .eq. mpi_p_root) then
        write(1234,*) ritznum

        ! write out eigenvectors
        do i=1,ngrow
          do j=1,ntri
            write (5000+i,*) real(eigvecs(j,i)), aimag(eigvecs(j,i))
          end do
        end do

        allocate(amat(ngrow,ngrow),bvec(ngrow,ngrow),ipiv(ngrow),coefren(ngrow))
        bvec=(0.d0,0.d0)
        do i=1,ngrow
           bvec(i,i)=(1.d0,0.d0)
           do j=1,ngrow
              amat(i,j)=sum(conjg(eigvecs(:,i))*eigvecs(:,j))*(ritznum(j)-(1.d0,0.d0))
           enddo
        enddo
        !
        call zgesv(ngrow,ngrow,amat,ngrow,ipiv,bvec,ngrow,info)
        !
        if(info.ne.0) then
           if(info.gt.0) then
              print *,'iterator: singular matrix in zgesv'
           else
              print *,'iterator: argument ',-info,' has illigal value in zgesv'
           endif
           deallocate(coefren,amat,bvec,ipiv)
           deallocate(eigvecs)
           return
        endif
     endif
  endif

  call MPI_BCAST(ngrow, 1, MPI_INTEGER4, mpi_p_root, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(do_arnoldi, 1, MPI_INTEGER4, mpi_p_root, MPI_COMM_WORLD, ierr)
  
  
     
  !  Main loop:
  !
  !
  
  bnorm_plas_prev=(0.d0,0.d0)
  bnorm_plas=(0.d0,0.d0)

  ! for reading input at first iteration     
  !open(1,file='RMP_EXCHANGE/hpsi_relaxed0.dat')
  !do i=1,ntri
  !   read (1,*) dummy_re,dummy_im
  !   bnorm_plas_prev(i)=cmplx(dummy_re,dummy_im)
  !enddo
  !close(1)


  precond = .True. ! TEST with .True.: always use the same random numbers
  do iptrace=1,macrosteps
     ! BEGIN test vacuum field only
     !     bnorm_plas_prev=(0.d0,0.d0)
     !     bnorm_plas=(0.d0,0.d0)
     ! END test vacuum field only

     call next_iteration(ntri, bnorm_plas_prev, bnorm_plas)
!print *, mype, 'next iteration finished', do_arnoldi
     if (do_arnoldi) then
        if (mype .eq. mpi_p_root) then
!print *, mype, 'root', ngrow
           do j=1,ngrow
              coefren(j)=ritznum(j)*sum(bvec(j,:)                           &
                   *matmul(transpose(conjg(eigvecs(:,1:ngrow))),bnorm_plas-bnorm_plas_prev))
           enddo
           bnorm_plas = bnorm_plas - matmul(eigvecs(:,1:ngrow),coefren(1:ngrow))
!print *, 'bnorm_plas updated'
!else
!print *, mype, 'not root', ngrow
        endif
        call MPI_BCAST(bnorm_plas, ntri, MPI_DOUBLE_COMPLEX, mpi_p_root, MPI_COMM_WORLD, ierr)
     endif
     bnorm_plas_prev = bnorm_plas
  end do
  
#ifdef PARALLEL
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  call MPI_FINALIZE(ierr)
#endif
  stop
! ----------------------------------------------------------------------------------------------
contains
  
  subroutine next_iteration(n, hold, hnew)
    integer, intent(in)      :: n         ! unused, just for arnoldi, use ntri instead
    complex(dp), intent(in)  :: hold(ntri)
    complex(dp), intent(out) :: hnew(ntri)

    if (precond) then
       ! reset RNG for preconditioner
       seed =  17491.d0 ! 12345.d0 !
       do i=1, max_openmp_threads
          call RngStream(seed, cg_rng)
          cg_rng_sav(:,i-1) = cg_rng
       enddo
       cg_rng = cg_rng_sav(:,mype+1)
! TEST: same random numbers for all processes
!cg_rng = cg_rng_sav(:,1)
!print *, mype, cg_rng
    end if
     
    do i=1,ntri
       mesh_element_rmp(i)%currents(:,:)=(0.d0,0.d0)
       mesh_element_rmp(i)%denspert(:)=(0.d0,0.d0)
       mesh_element_rmp(i)%pprespert(:)=(0.d0,0.d0)
       mesh_element_rmp(i)%parcurrpert(:)=(0.d0,0.d0)
       mesh_element_rmp(i)%bnorm_times_thermforces=mesh_element(i)%thermforces &
            *(mesh_element_rmp(i)%bnorm_vac+hold(i))
    end do
    
    if(.false.) then !=====>turn off update<=====
       do i=2,npoint ! excluding O-point
          ePhi_knot = 0.d0
          D_knot = 0.d0
          V_knot = 0.d0
          T_knot = 0.d0
          do j=1,mesh_point(i)%n_owners
             !           ePhi_knot = ePhi_knot + mesh_point(i)%weight_intp(j)                  &
             !                * mesh_element(mesh_point(i)%i_owner_tri(j))%ePhi_tri
             do k=1,nsorts
                D_knot(k) = D_knot(k) + mesh_point(i)%weight_intp(j)                           &
                     * D_therm(mesh_element(mesh_point(i)%i_owner_tri(j))%i_therm,k)
                V_knot(k) = V_knot(k) + mesh_point(i)%weight_intp(j)                           &
                     * V_therm(mesh_element(mesh_point(i)%i_owner_tri(j))%i_therm,k)
                T_knot(k) = T_knot(k) + mesh_point(i)%weight_intp(j)                           &
                     * T_therm(mesh_element(mesh_point(i)%i_owner_tri(j))%i_therm,k)
             enddo
          enddo
          !        mesh_point(i)%PhiovB2 = ePhi_knot/(mesh_point(i)%b_mod)**2
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
       do j=1, ntri_inbou
          !           ePhi_knot = ePhi_knot + mesh_point(i)%weight_intp(j)           &
          !                * mesh_element(mesh_point(i)%i_owner_tri(j))%ePhi_tri
          do k=1,nsorts
             D_knot(k) = D_knot(k) + inbou_list(j)%vol_norm                          &
                  * D_therm(mesh_element(inbou_list(j)%i_tri_all)%i_therm,k)
             V_knot(k) = V_knot(k) + inbou_list(j)%vol_norm                          &
                  * V_therm(mesh_element(inbou_list(j)%i_tri_all)%i_therm,k)
             T_knot(k) = T_knot(k) + inbou_list(j)%vol_norm                          &
                  * T_therm(mesh_element(inbou_list(j)%i_tri_all)%i_therm,k)
          enddo
       enddo
       !        mesh_point(i)%PhiovB2 = ePhi_knot/(mesh_point(i)%b_mod)**2
       do k=1,nsorts
          mesh_point(1)%pl_parms_knot(1,k) = D_knot(k)
          mesh_point(1)%pl_parms_knot(2,k) = V_knot(k)
          mesh_point(1)%pl_parms_knot(3,k) = T_knot(k)
       enddo

!!$     do i=1,ntri
!!$        call grad_PhiovB2(i)
!!$     enddo
    endif !=====>turn off update<=====

    ! copy arrays for MPI:
    do is=1, nsorts
       time_local(is) = time_l(is)
       Dm_local(:,is) = mesh_element(:)%Dm_part(is)
       Dp_local(:,is) = mesh_element(:)%D_part(is)
       Vp_local(:,is) = mesh_element(:)%V_part(is)
       Tp_local(:,is) = mesh_element(:)%T_part(is)
       Vt2_local(:,is) = mesh_element(:)%Vt2_part(is)
    enddo
    ephi_local(:) = mesh_element(:)%ePhi_tri

    !print *,'T min = ', minval(T_therm), mype, '-2'
    !print *, cg_rng, npgroup, iptrace, mype, '-2'
!    call macrostep(cg_rng, npgroup, iptrace)a          !<= usual (collisional) macrostep
    call macrostep_ba(cg_rng, npgroup, iptrace)        !<= new (collisionless) macrostep 
    !print *, cg_rng, npgroup, iptrace, mype, '-1'
    !print *,'T min = ', minval(T_therm), mype, '-1'
#ifdef PARALLEL
!!$     D_time_l = 0.d0
!!$     D_Dm = 0.d0
!!$     D_Dp = 0.d0
!!$     D_Vp = 0.d0
!!$     D_Tp = 0.d0
!!$     D_Dt = 0.d0
!!$     D_Vt = 0.d0
!!$     D_Tt = 0.d0
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    
    write (temp, '(I1.1)') mype

! Plot currents for all MPI jobs (DEBUG)
!    open(1,file='PLOTTING/curr_toplot'//trim(temp)//'.dat')
!       do i=1,ntri
!          R_vert(:) = mesh_point( mesh_element(i)%i_knot(:) )%rcoord
!          Z_vert(:) = mesh_point( mesh_element(i)%i_knot(:) )%zcoord
!          delta=(Z_vert(2)-Z_vert(1))*(R_vert(2)-R_vert(3)) &
!               -(Z_vert(2)-Z_vert(3))*(R_vert(2)-R_vert(1))
!          do isort=1,nsorts
!             cur1=mesh_element_rmp(i)%currents(1,isort)
!             cur2=mesh_element_rmp(i)%currents(2,isort)
!             cur_R(isort)=(cur1*(R_vert(2)-R_vert(3))+cur2*(R_vert(2)-R_vert(1)))/delta
!             cur_Z(isort)=(cur1*(Z_vert(2)-Z_vert(3))+cur2*(Z_vert(2)-Z_vert(1)))/delta
!          enddo
!          write (1,*) real(cur_R),real(cur_Z),dimag(cur_R),dimag(cur_Z)
!       enddo
!       close(1)
!       
!       open(1,file='RMP_EXCHANGE/currents'//trim(temp)//'.dat')
!       do i=1,ntri
!          write(1,*) real(mesh_element_rmp(i)%currents(:,1)              &
!               +mesh_element_rmp(i)%currents(:,2)),            &
!               dimag(mesh_element_rmp(i)%currents(:,1)             &
!               +mesh_element_rmp(i)%currents(:,2))
!       enddo
!       close(1)
!
    call average_mpi()
#else
    do is=1, nsorts
       time_l(is) = time_local(is)
       mesh_element(:)%Dm_part(is) = Dm_local(:,is)
       mesh_element(:)%D_part(is) = Dp_local(:,is)
       mesh_element(:)%V_part(is) = Vp_local(:,is)
       mesh_element(:)%T_part(is) = Tp_local(:,is)
       mesh_element(:)%Vt2_part(is) = Vt2_local(:,is)
    enddo
    mesh_element(:)%ePhi_tri = ephi_local(:)
#endif
!!$     V_therm = 0.0d0
!!$     do is=1, nsorts
!!$        do i=1,ntri
!!$           ind_therm = mesh_element(i)%i_therm
!!$           V_therm(ind_therm,is) = V_therm(ind_therm,is) +       &
!!$                mesh_element(i)%V_part(is)*mesh_element(i)%V_tri/vol_therm(ind_therm)
!!$        enddo
!!$     enddo

!!$!attention, electrons:
!!$D_therm(:,isort1) = D_therm(:,isort2) 
!!$V_therm(:,isort1) = 0.0d0
!!$do i=1,ntri
!!$   ind_therm = mesh_element(i)%i_therm
!!$   V_therm(ind_therm,isort1) = V_therm(ind_therm,isort1) +       &
!!$        mesh_element(i)%V_part(isort2)*mesh_element(i)%V_tri/vol_therm(ind_therm)
!!$enddo
!!$T_therm(:,isort1) = 0.0d0
!!$do i=1,ntri
!!$   ind_therm = mesh_element(i)%i_therm
!!$   T_therm(ind_therm,isort1) = T_therm(ind_therm,isort1) +       &
!!$        mesh_element(i)%T_part(isort2)*mesh_element(i)%V_tri/vol_therm(ind_therm)
!!$enddo


!!$     do is=1, nsorts
!!$        call plates(is)
!!$     enddo
    if (mype .eq. mpi_p_root) then
       if((iptrace/n_prfs)*n_prfs .eq. iptrace)print *,'macrostep # ',iptrace
!        if((iptrace/n_evol)*n_evol .eq. iptrace) then 
       if(.false.) then  !======> this branch is turned off<========
          open(1004,file='RESULTS/time_evol.dat', position='append',status='old')
          dens_tp_avr = 0.d0
          temp_tp_avr = 0.d0
          v2_tp_avr = 0.d0
          v2_avr = 0.d0
          temp_avr = 0.d0
          do i=1,ntri
             dens_tp_avr = dens_tp_avr + mesh_element(i)%D_part(isort2)*mesh_element(i)%V_tri
             temp_tp_avr = temp_tp_avr + mesh_element(i)%T_part(isort2)*mesh_element(i)%V_tri
             !    v2_tp_avr = v2_tp_avr + mesh_element(i)%Vt2_part(isort2)*mesh_element(i)%V_tri
             v2_tp_avr = v2_tp_avr + ( mesh_element(i)%V_part(isort2) )**2                &
                  *amass(isort2)*erg2ev*5.d-1*mesh_element(i)%V_tri
             v2_avr = v2_avr + ( V_therm(mesh_element(i)%i_therm,isort2) )**2    &
                  *amass(isort2)*erg2ev*5.d-1*mesh_element(i)%V_tri
             temp_avr = temp_avr + T_therm(mesh_element(i)%i_therm,isort2)*mesh_element(i)%V_tri
          enddo
          dens_tp_avr = dens_tp_avr/sum(mesh_element(:)%V_tri)
          temp_tp_avr = temp_tp_avr/sum(mesh_element(:)%V_tri)
          v2_tp_avr = v2_tp_avr/sum(mesh_element(:)%V_tri)
          v2_avr = v2_avr/sum(mesh_element(:)%V_tri)
          temp_avr = temp_avr/sum(mesh_element(:)%V_tri)
!!$           part_sink(isort2) = (dp_j_l(isort2) + dp_j_r(isort2) +  j_outer_wall(isort2) +   &
!!$                                j_pfz_wall(isort2))/sum(mesh_element(:)%V_tri)
!!$           temp_out(isort2) = (dp_q_l(isort2) + dp_q_r(isort2) + q_outer_wall(isort2) + q_pfz_wall(isort2))/    &
!!$                              (dp_j_l(isort2) + dp_j_r(isort2) + j_outer_wall(isort2) + j_pfz_wall(isort2) +    &
!!$                               epsilon(1.d0))/1.5d0*erg2ev 
          write(1004,204)fulltime, dens_tp_avr, temp_tp_avr, v2_tp_avr, temp_avr, v2_avr, time_l(isort2),      &
               maxval(V_therm(:,isort2)), minval(V_therm(:,isort2)),                                           &
               maxval(T_therm(:,isort2)), minval(T_therm(:,isort2))

!!$                dp_j_l(isort2),  dp_j_r(isort2),  j_outer_wall(isort2),  j_pfz_wall(isort2),&
!!$                dp_q_l(isort2),  dp_q_r(isort2),  q_outer_wall(isort2),  q_pfz_wall(isort2),&
!!$                dp_tm_l(isort2), dp_tm_r(isort2), tm_outer_wall(isort2), tm_pfz_wall(isort2),&
!!$                part_sink(isort2), temp_out(isort2)

          close(1004)
          !           open(1005,file='RESULTS/bookkeeping.dat', position='append',status='old')
          !           write(1005,205)istochnik, stenka, proebal, ubila_rr, dobavila_rr
          !           close(1005)
       endif
!        if((iptrace/n_prfs)*n_prfs .eq. iptrace) then ! write everything
       if(.false.) then  !======> this branch is turned off<========
          open(1,file='RESULTS/points.dat',form='unformatted')
          write(1) npoint
          write(1) mesh_point
          close(1)
          open(1,file='RESULTS/triangles.dat',form='unformatted')
          write(1) ntri
          write(1) mesh_element
          write(1) bphicovar, fulltime, time_l
          close(1)
          open(1,file='RESULTS/wall_loads.dat',form='unformatted')
          write(1) flux_j
          write(1) flux_q
          write(1) tor_mom
          close(1)
!!$v3_in(:) = 0.d0
!!$v3_out(:) = 0.d0
!!$where (part_num .gt. 0.d0) 
!!$   v3_in(:) = v3_sum_in(:)/part_num(:)
!!$   v3_out(:) = v3_sum_out(:)/part_num(:)
!!$end where


          open(1,file='RESULTS/thermostat.dat',form='unformatted')
          write(1) n_therm
          write(1) D_therm
          write(1) V_therm
          write(1) T_therm
          write(1) vol_therm
          !           write(1) v3_in
          !           write(1) v3_out
          close(1) 

          open(1,file='RESULTS/wpartoverwcell.dat')
          do i=1,ntri
             r_bc = (mesh_point(mesh_element(i)%i_knot(1))%rcoord + mesh_point(mesh_element(i)%i_knot(2))%rcoord  &
                  + mesh_point(mesh_element(i)%i_knot(3))%rcoord)*one3rd
             z_bc = (mesh_point(mesh_element(i)%i_knot(1))%zcoord + mesh_point(mesh_element(i)%i_knot(2))%zcoord  &
                  + mesh_point(mesh_element(i)%i_knot(3))%zcoord)*one3rd
             write(1,204)r_bc, z_bc, rel_therm_max(i)
          enddo
          close(1)
          !

!!$           open(1,file='RESULTS/v3inout.dat')
!!$           do i=1,ntri
!!$              ind_therm = mesh_element(i)%i_therm
!!$              r_bc = (mesh_point(mesh_element(i)%i_knot(1))%rcoord + mesh_point(mesh_element(i)%i_knot(2))%rcoord  &
!!$                   + mesh_point(mesh_element(i)%i_knot(3))%rcoord)*one3rd
!!$              z_bc = (mesh_point(mesh_element(i)%i_knot(1))%zcoord + mesh_point(mesh_element(i)%i_knot(2))%zcoord  &
!!$                   + mesh_point(mesh_element(i)%i_knot(3))%zcoord)*one3rd
!!$              if(part_num(ind_therm) .gt. 0.d0) then
!!$                 write(1,204)r_bc, z_bc, v3_in(ind_therm), mesh_element(i)%V_part(2), mesh_element(i)%D_part(2) &
!!$                      , atan2((z_bc-Z_ax),(r_bc-R_ax))
!!$              endif
!!$           enddo
!!$           close(1)
!!$#ifdef PARALLEL
!!$           open(1,file='RESULTS/dispersion.dat')
!!$           do i=1,ntri
!!$              r_bc = (mesh_point(mesh_element(i)%i_knot(1))%rcoord + mesh_point(mesh_element(i)%i_knot(2))%rcoord  &
!!$                   + mesh_point(mesh_element(i)%i_knot(3))%rcoord)*one3rd
!!$              z_bc = (mesh_point(mesh_element(i)%i_knot(1))%zcoord + mesh_point(mesh_element(i)%i_knot(2))%zcoord  &
!!$                   + mesh_point(mesh_element(i)%i_knot(3))%zcoord)*one3rd
!!$              j = mesh_element(i)%i_therm
!!$              write(1,204)r_bc, z_bc, D_Dm(i,isort), D_Dp(i,isort), D_Vp(i,isort), D_Tp(i,isort), &
!!$                                      D_Dt(j,isort), D_Vt(j,isort), D_Tt(j,isort), D_ephi(i)
!!$           enddo
!!$           close(1)
!!$#endif
       endif
!
       open(1,file='RESULTS/restart_rng.int',form='unformatted')
       write(1) cg_rng
       close(1)
!
       open(1,file='PLOTTING/curr_toplot.dat')
       open(2,file='PLOTTING/dens_toplot.dat')
       open(3,file='PLOTTING/ppres_toplot.dat')
       open(11,file='PLOTTING/parcurr_toplot.dat')
       do i=1,ntri
          R_vert(:) = mesh_point( mesh_element(i)%i_knot(:) )%rcoord
          Z_vert(:) = mesh_point( mesh_element(i)%i_knot(:) )%zcoord
          delta=(Z_vert(2)-Z_vert(1))*(R_vert(2)-R_vert(3)) &
               -(Z_vert(2)-Z_vert(3))*(R_vert(2)-R_vert(1))
          do isort=1,nsorts
             cur1=mesh_element_rmp(i)%currents(1,isort)
             cur2=mesh_element_rmp(i)%currents(2,isort)
             cur_R(isort)=(cur1*(R_vert(2)-R_vert(3))+cur2*(R_vert(2)-R_vert(1)))/delta
             cur_Z(isort)=(cur1*(Z_vert(2)-Z_vert(3))+cur2*(Z_vert(2)-Z_vert(1)))/delta
          enddo
          write (1,*) real(cur_R),real(cur_Z),dimag(cur_R),dimag(cur_Z)
          write (2,*) real(mesh_element_rmp(i)%denspert(:)), &
                      dimag(mesh_element_rmp(i)%denspert(:))
          write (3,*) real(mesh_element_rmp(i)%pprespert(:)), &
                      dimag(mesh_element_rmp(i)%pprespert(:))
          write (11,*) real(mesh_element_rmp(i)%parcurrpert(:)), &
                       dimag(mesh_element_rmp(i)%parcurrpert(:))
       enddo
       close(1)
       close(2)
       close(3)
       close(11)
 stop      
!
       open(1,file='RMP_EXCHANGE/currents.dat')
       do i=1,ntri
          write(1,*) real(mesh_element_rmp(i)%currents(:,1)              &
               +mesh_element_rmp(i)%currents(:,2)),            &
               dimag(mesh_element_rmp(i)%currents(:,1)             &
               +mesh_element_rmp(i)%currents(:,2))
       enddo
       close(1)
!
       open(1,file='RMP_EXCHANGE/currents1.dat')
       do i=1,ntri
          write(1,*) real(mesh_element_rmp(i)%currents(:,1)),            &
               dimag(mesh_element_rmp(i)%currents(:,1))
       enddo
       close(1)
!
       open(1,file='RMP_EXCHANGE/currents2.dat')
       do i=1,ntri
          write(1,*) real(mesh_element_rmp(i)%currents(:,2)),            &
               dimag(mesh_element_rmp(i)%currents(:,2))
       enddo
       close(1)
!
       call update_field(ierr)
       print *,'Maxwell solver exit status = ',ierr
!       
       open(1,file='RMP_EXCHANGE/hpsi.dat')
       open(2,file='PLOTTING/hpsi_relaxed.dat')
       do i=1,ntri
          read (1,*) dummy_re,dummy_im
          hnew(i)=cmplx(dummy_re,dummy_im)*dens_reduce_factor
          hnew(i)=(1.d0-eps_urelax)*hold(i)+eps_urelax*hnew(i)
          write (2,*) real(hnew(i)),dimag(hnew(i))
       enddo
       close(1)
       close(2)
!
       open(3210,file='convergence.dat',position='append')
       write(3210,*) iptrace,sum(abs(hnew))
       close(3210)
!
!
    end if
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call MPI_BCAST(hnew, ntri, MPI_DOUBLE_COMPLEX, mpi_p_root, MPI_COMM_WORLD, ierr)
    !print *,'T min = ', minval(T_therm), mype, '4'
    
204 format(1000(ES12.5E2,1X))
205 format(1000(ES23.16E2,1X))
   end subroutine next_iteration

#ifdef PARALLEL
  subroutine broadcast_all()
    
    implicit none

    integer(i4b) :: nword, len_int, len_dbl
    integer(i4b) :: ierr_c, ierr, k, mpi_triangle, mpi_bound, mpi_mesh_point, mpi_tri_rmp
    integer(kind=mpi_address_kind) :: ext_dbl, ext_lb_dbl, ext_int, ext_lb_int
    integer(i4b), parameter :: nelems_tri = 23 ! 4, 2*legs, 10, 5*nsorts, 1, 1 
    integer(kind=mpi_address_kind), dimension(nelems_tri) :: struct_disp_tri
    integer(i4b), dimension(nelems_tri) :: struct_blen_tri, struct_type_tri
    integer(kind=mpi_address_kind), dimension(2) :: struct_disp_bou
    integer(i4b), dimension(2) :: struct_blen_bou, struct_type_bou
    integer(i4b), parameter :: nelems_poi =  10 !6,legs*nsorts,1,n_owners_max,n_owners_max
    integer(kind=mpi_address_kind), dimension(nelems_poi) :: struct_disp_poi
    integer(i4b), dimension(nelems_poi) :: struct_blen_poi, struct_type_poi
    integer(i4b), parameter :: nelems_trirmp = 3 !  
    integer(kind=mpi_address_kind), dimension(nelems_trirmp) :: struct_disp_trirmp
    integer(i4b), dimension(nelems_trirmp) :: struct_blen_trirmp, struct_type_trirmp

! kin2d.inp:
    call MPI_BCAST(bmod_ref,  1, MPI_DOUBLE_PRECISION, mpi_p_root, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(nptrace,   1, MPI_INTEGER4, mpi_p_root, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(source0,   1, MPI_DOUBLE_PRECISION, mpi_p_root, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(T_source,  1, MPI_DOUBLE_PRECISION, mpi_p_root, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(relax_tau, 1, MPI_DOUBLE_PRECISION, mpi_p_root, MPI_COMM_WORLD, ierr)
    nword = nsorts
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call MPI_BCAST(wcv(1), nword, MPI_DOUBLE_PRECISION, mpi_p_root, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(wct(1), nword, MPI_DOUBLE_PRECISION, mpi_p_root, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(D_perp,     1, MPI_DOUBLE_PRECISION, mpi_p_root, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(n_evol,     1, MPI_INTEGER4, mpi_p_root, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(n_prfs,     1, MPI_INTEGER4, mpi_p_root, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(npgroup,    1, MPI_INTEGER4, mpi_p_root, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(restart,    1, MPI_DOUBLE_PRECISION, mpi_p_root, MPI_COMM_WORLD, ierr)
! triangles.dat:
    call MPI_BCAST(ntri, 1, MPI_INTEGER4, mpi_p_root, MPI_COMM_WORLD, ierr)
    if(mype .ne. mpi_p_root) then
       allocate(mesh_element(ntri))
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
! MPI-type creation
    struct_blen_tri(:) = 1
    struct_blen_tri(3:5) = 3 ! legs
    struct_blen_tri(16:20) = 2 ! nsorts
    struct_blen_tri(23) = 4 ! 2*nsorts
    len_int = 2 + 3 !(*legs)
    len_dbl = 10 + 5 + 1 ! 5(*nsorts) + 1(2*nsorts)

    CALL MPI_GET_ADDRESS(mesh_element(1)%ix_qq,        struct_disp_tri(1), ierr)
    if(ierr .ne. MPI_SUCCESS) print *,'1'
    CALL MPI_GET_ADDRESS(mesh_element(1)%iy_qq,        struct_disp_tri(2), ierr)
    if(ierr .ne. MPI_SUCCESS) print *,'2'
    CALL MPI_GET_ADDRESS(mesh_element(1)%i_knot(1),    struct_disp_tri(3), ierr)
    if(ierr .ne. MPI_SUCCESS) print *,'4'
    CALL MPI_GET_ADDRESS(mesh_element(1)%neighbour(1), struct_disp_tri(4), ierr)
    if(ierr .ne. MPI_SUCCESS) print *,'5'
    CALL MPI_GET_ADDRESS(mesh_element(1)%neighbour_edge(1), struct_disp_tri(5), ierr)
    if(ierr .ne. MPI_SUCCESS) print *,'6'
    CALL MPI_GET_ADDRESS(mesh_element(1)%sizmaxtri,    struct_disp_tri(6), ierr)
    if(ierr .ne. MPI_SUCCESS) print *,'7'
    CALL MPI_GET_ADDRESS(mesh_element(1)%det_3,        struct_disp_tri(7), ierr)
    if(ierr .ne. MPI_SUCCESS) print *,'8'
    CALL MPI_GET_ADDRESS(mesh_element(1)%oneoverD3,    struct_disp_tri(8), ierr)
    if(ierr .ne. MPI_SUCCESS) print *,'9'
    CALL MPI_GET_ADDRESS(mesh_element(1)%V_tri,        struct_disp_tri(9), ierr)
    if(ierr .ne. MPI_SUCCESS) print *,'10'
    CALL MPI_GET_ADDRESS(mesh_element(1)%dBinv_dR,     struct_disp_tri(10), ierr)
    if(ierr .ne. MPI_SUCCESS) print *,'11'
    CALL MPI_GET_ADDRESS(mesh_element(1)%dBinv_dZ,     struct_disp_tri(11), ierr)
    if(ierr .ne. MPI_SUCCESS) print *,'12'
    CALL MPI_GET_ADDRESS(mesh_element(1)%dPsi_dR,      struct_disp_tri(12), ierr)
    if(ierr .ne. MPI_SUCCESS) print *,'13'
    CALL MPI_GET_ADDRESS(mesh_element(1)%dPsi_dZ,      struct_disp_tri(13), ierr)
    if(ierr .ne. MPI_SUCCESS) print *,'14'
    CALL MPI_GET_ADDRESS(mesh_element(1)%dPhiovB2_dR,  struct_disp_tri(14), ierr)
    if(ierr .ne. MPI_SUCCESS) print *,'15'
    CALL MPI_GET_ADDRESS(mesh_element(1)%dPhiovB2_dZ,  struct_disp_tri(15), ierr)
    if(ierr .ne. MPI_SUCCESS) print *,'16'
    CALL MPI_GET_ADDRESS(mesh_element(1)%Dm_part(1),   struct_disp_tri(16), ierr)
    if(ierr .ne. MPI_SUCCESS) print *,'17'
    CALL MPI_GET_ADDRESS(mesh_element(1)%D_part(1),    struct_disp_tri(17), ierr)
    if(ierr .ne. MPI_SUCCESS) print *,'18'
    CALL MPI_GET_ADDRESS(mesh_element(1)%V_part(1),    struct_disp_tri(18), ierr)
    if(ierr .ne. MPI_SUCCESS) print *,'19'
    CALL MPI_GET_ADDRESS(mesh_element(1)%T_part(1),    struct_disp_tri(19), ierr)
    if(ierr .ne. MPI_SUCCESS) print *,'20'
    CALL MPI_GET_ADDRESS(mesh_element(1)%Vt2_part(1),  struct_disp_tri(20), ierr)
    if(ierr .ne. MPI_SUCCESS) print *,'21'
    CALL MPI_GET_ADDRESS(mesh_element(1)%ePhi_tri,     struct_disp_tri(21), ierr)
    if(ierr .ne. MPI_SUCCESS) print *,'22'
    CALL MPI_GET_ADDRESS(mesh_element(1)%i_therm,      struct_disp_tri(22), ierr)
    if(ierr .ne. MPI_SUCCESS) print *,'22'
    CALL MPI_GET_ADDRESS(mesh_element(1)%thermforces(1,1),  struct_disp_tri(23), ierr)
    if(ierr .ne. MPI_SUCCESS) print *,'23'
    do i = 2, nelems_tri
       struct_disp_tri(i) = struct_disp_tri(i) - struct_disp_tri(1)
    enddo
    struct_disp_tri(1) = 0_mpi_address_kind

    struct_type_tri(:) = (/ ( MPI_INTEGER4, k=1,len_int ),                                        &
         ( MPI_DOUBLE_PRECISION, k=1,len_dbl), &
    MPI_INTEGER4, MPI_DOUBLE_PRECISION /)

!!$print *,size(struct_blen_tri), struct_blen_tri, mype
!!$print *,struct_disp_tri, mype
!!$print *,struct_type_tri, ntri, mype
!!$print *,len_int,len_dbl, mype

    call mpi_type_create_struct(nelems_tri, struct_blen_tri(1:nelems_tri),&
         struct_disp_tri(1:nelems_tri),   &
         struct_type_tri(1:nelems_tri), mpi_triangle, ierr_c )
    if(ierr_c .ne. MPI_SUCCESS) print *,'Error in MPI type creation mpi_triangle'
    call mpi_type_commit(mpi_triangle, ierr_c)
    if(ierr_c .ne. MPI_SUCCESS) print *,'Error in MPI type commit mpi_triangle'
    nword = ntri
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    do i = 1, ntri
       call MPI_BCAST(mesh_element(i)%ix_qq, 1, mpi_triangle, mpi_p_root, MPI_COMM_WORLD, ierr)
       if(ierr .ne. MPI_SUCCESS) print *,'Error in BCAST of derived type mpi_triangle'
    end do
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)

!    if(mype.eq.0) print *,'PE #', mype, mesh_element(ntri)
!    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!    if(mype.eq.1) print *,'PE #', mype, mesh_element(ntri)
!    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call MPI_BCAST(bphicovar, 1, MPI_DOUBLE_PRECISION, mpi_p_root, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(fulltime,  1, MPI_DOUBLE_PRECISION, mpi_p_root, MPI_COMM_WORLD, ierr) 
!   but used the version of fulltime from root PE
    nword = nsorts
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call MPI_BCAST(time_l(1), nword, MPI_DOUBLE_PRECISION, mpi_p_root, MPI_COMM_WORLD, ierr)
! thermostat:
    call MPI_BCAST(n_therm, 1, MPI_INTEGER4, mpi_p_root, MPI_COMM_WORLD, ierr)
    if(mype .ne. mpi_p_root) then
       allocate(D_therm(n_therm,nsorts),V_therm(n_therm,nsorts),T_therm(n_therm,nsorts))
    endif
    nword = nsorts*n_therm
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call MPI_BCAST(D_therm(1,1), nword, MPI_DOUBLE_PRECISION, mpi_p_root, MPI_COMM_WORLD, ierr)
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call MPI_BCAST(V_therm(1,1), nword, MPI_DOUBLE_PRECISION, mpi_p_root, MPI_COMM_WORLD, ierr)
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call MPI_BCAST(T_therm(1,1), nword, MPI_DOUBLE_PRECISION, mpi_p_root, MPI_COMM_WORLD, ierr)
    if(mype .ne. mpi_p_root) then
       allocate(vol_therm(n_therm))
    endif
    nword = n_therm


!print *,'T min = ', minval(T_therm), mype, '0'


    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call MPI_BCAST(vol_therm(1), nword, MPI_DOUBLE_PRECISION, mpi_p_root, MPI_COMM_WORLD, ierr)
! inbou.dat:
    call MPI_BCAST(ntri_inbou, 1, MPI_INTEGER4, mpi_p_root, MPI_COMM_WORLD, ierr)
    nword = ntri_inbou
    if(mype .ne. mpi_p_root) then
       allocate(inbou_list(ntri_inbou))
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)

    struct_blen_bou(:) = (/1, 1/)
    CALL MPI_GET_ADDRESS(inbou_list(1)%i_tri_all,        struct_disp_bou(1), ierr)
    CALL MPI_GET_ADDRESS(inbou_list(1)%vol_norm,         struct_disp_bou(2), ierr)
    struct_disp_bou(2) = struct_disp_bou(2) - struct_disp_bou(1)
    struct_disp_bou(1) = 0_mpi_address_kind

    struct_type_bou(:) = (/ MPI_INTEGER4, MPI_DOUBLE_PRECISION /)
    call mpi_type_create_struct( size(struct_blen_bou), struct_blen_bou(:), struct_disp_bou(:), &
                                 struct_type_bou(:), mpi_bound, ierr_c )
    call mpi_type_commit(mpi_bound, ierr_c)
    if(ierr_c .ne. MPI_SUCCESS) print *,'Error in MPI type definition mpi_bound'
    do i = 1, ntri_inbou
       call MPI_BCAST(inbou_list(i)%i_tri_all, 1, mpi_bound, mpi_p_root, MPI_COMM_WORLD, ierr)
       if(ierr .ne. MPI_SUCCESS) print *,'Error in BCAST of derived type mpi_bound'
    enddo
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
! print *,mype,inbou_list(ntri_inbou)%i_tri_all

! mesh_points:
    call MPI_BCAST(npoint, 1, MPI_INTEGER4, mpi_p_root, MPI_COMM_WORLD, ierr)
    if(mype .ne. mpi_p_root) then
       allocate(mesh_point(npoint))
    endif
! MPI-type creation
    len_int = 1 + 1 ! n_owners_max
    len_dbl = 6 + 1 ! 3*nsorts 
    struct_blen_poi(:) = 1
    struct_blen_poi(7) = legs*nsorts  ! 3*2
    struct_blen_poi(9:10) = n_owners_max

    CALL MPI_GET_ADDRESS(mesh_point(1)%rcoord,             struct_disp_poi(1), ierr)
    CALL MPI_GET_ADDRESS(mesh_point(1)%zcoord,             struct_disp_poi(2), ierr)
    CALL MPI_GET_ADDRESS(mesh_point(1)%psi_pol,            struct_disp_poi(3), ierr)
    CALL MPI_GET_ADDRESS(mesh_point(1)%b_mod,              struct_disp_poi(4), ierr)
    CALL MPI_GET_ADDRESS(mesh_point(1)%b_phi,              struct_disp_poi(5), ierr)
    CALL MPI_GET_ADDRESS(mesh_point(1)%PhiovB2,            struct_disp_poi(6), ierr)
    CALL MPI_GET_ADDRESS(mesh_point(1)%pl_parms_knot(1,1), struct_disp_poi(7), ierr)
    CALL MPI_GET_ADDRESS(mesh_point(1)%n_owners,           struct_disp_poi(8), ierr)
    CALL MPI_GET_ADDRESS(mesh_point(1)%i_owner_tri(1),     struct_disp_poi(9), ierr)
    CALL MPI_GET_ADDRESS(mesh_point(1)%weight_intp(1),     struct_disp_poi(10), ierr)

    do i=nelems_poi,1,-1
       struct_disp_poi(i) = struct_disp_poi(i) - struct_disp_poi(1)
    enddo

    struct_type_poi(:) = (/ ( MPI_DOUBLE_PRECISION, k=1,len_dbl),                             &
                            ( MPI_INTEGER4, k=1,len_int ),                                    &
                              MPI_DOUBLE_PRECISION /)
    call mpi_type_create_struct(size(struct_blen_poi), struct_blen_poi(:), struct_disp_poi(:),&
                                struct_type_poi(:), mpi_mesh_point, ierr_c )
    call mpi_type_commit(mpi_mesh_point, ierr_c)
    if(ierr_c .ne. MPI_SUCCESS) print *,'Error in MPI type definition, mpi_mesh_point '

    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    do i = 1, npoint
       call MPI_BCAST(mesh_point(i)%rcoord, 1, mpi_mesh_point, mpi_p_root, MPI_COMM_WORLD,ierr)
       if(ierr .ne. MPI_SUCCESS) print *,'Error in BCAST of derived type mpi_mesh_point'
    end do
!    if(mype.eq.0) print *,'PE #', mype, mesh_point(npoint)%i_owner_tri(:)
!    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!    if(mype.eq.1) print *,'PE #', mype, mesh_point(npoint)%i_owner_tri(:)
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!!$print *,size(struct_blen_poi), struct_blen_poi, mype
!!$print *,struct_disp_poi, mype
!!$print *,struct_type_poi, mype
!!$print *,len_int, len_dbl, mype
!!$stop

! the whole boundary: 
    call MPI_BCAST(i_pfz(1), 2, MPI_INTEGER4, mpi_p_root, MPI_COMM_WORLD, ierr)
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call MPI_BCAST(i_dpr(1), 2, MPI_INTEGER4, mpi_p_root, MPI_COMM_WORLD, ierr)
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!    div_r_count = i_dpr(2) - i_dpr(1) + 1
    call MPI_BCAST(i_sol(1), 2, MPI_INTEGER4, mpi_p_root, MPI_COMM_WORLD, ierr)
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call MPI_BCAST(i_dpl(1), 2, MPI_INTEGER4, mpi_p_root, MPI_COMM_WORLD, ierr)
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!    div_l_count = i_dpl(2) - i_dpl(1) + 1
    call MPI_BCAST(i_inb(1), 2, MPI_INTEGER4, mpi_p_root, MPI_COMM_WORLD, ierr)
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call MPI_BCAST(nboucell, 1, MPI_INTEGER4, mpi_p_root, MPI_COMM_WORLD, ierr)
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    if(mype .ne. mpi_p_root) then
!       allocate(R_bou(nboucell), Z_bou(nboucell), S_bou(nboucell))
       allocate(flux_j(nboucell,nsorts), flux_q(nboucell,nsorts), tor_mom(nboucell,nsorts))
    endif
!    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!    call MPI_BCAST(R_bou(1), nboucell, MPI_DOUBLE_PRECISION, mpi_p_root, MPI_COMM_WORLD, ierr)
!    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!    call MPI_BCAST(Z_bou(1), nboucell, MPI_DOUBLE_PRECISION, mpi_p_root, MPI_COMM_WORLD, ierr)
!    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!    call MPI_BCAST(S_bou(1), nboucell, MPI_DOUBLE_PRECISION, mpi_p_root, MPI_COMM_WORLD, ierr)
!    call MPI_BARRIER(MPI_COMM_WORLD,ierr)

    nword = nboucell*nsorts
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call MPI_BCAST(flux_j(1,1), nword, MPI_DOUBLE_PRECISION, mpi_p_root, MPI_COMM_WORLD, ierr)
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call MPI_BCAST(flux_q(1,1), nword, MPI_DOUBLE_PRECISION, mpi_p_root, MPI_COMM_WORLD, ierr)
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call MPI_BCAST(tor_mom(1,1), nword, MPI_DOUBLE_PRECISION, mpi_p_root, MPI_COMM_WORLD, ierr)
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)

    ! triangle RMP data:
    if(mype .ne. mpi_p_root) then
       allocate(mesh_element_rmp(ntri))
    endif
    
    ! MPI-type creation
    len_dbl = 3 ! 1 + 2*nsorts + legs*nsorts
    struct_blen_trirmp(:) = 1
    struct_blen_trirmp(2) = 2*nsorts  ! 2*2
    struct_blen_trirmp(3) = legs*nsorts  ! 3*2

    CALL MPI_GET_ADDRESS(mesh_element_rmp(1)%bnorm_vac,  struct_disp_trirmp(1), ierr)
    CALL MPI_GET_ADDRESS(mesh_element_rmp(1)%bnorm_times_thermforces,struct_disp_trirmp(2),ierr)
    CALL MPI_GET_ADDRESS(mesh_element_rmp(1)%currents, struct_disp_trirmp(3), ierr)

    do i=nelems_trirmp,1,-1
       struct_disp_trirmp(i) = struct_disp_trirmp(i) - struct_disp_trirmp(1)
    enddo

    struct_type_trirmp(:) = (/ ( MPI_DOUBLE_COMPLEX, k=1,len_dbl) /)
    
    call mpi_type_create_struct(size(struct_blen_trirmp), struct_blen_trirmp(:),&
         struct_disp_trirmp(:), struct_type_trirmp(:), mpi_tri_rmp, ierr_c )
    call mpi_type_commit(mpi_tri_rmp, ierr_c)
    if(ierr_c .ne. MPI_SUCCESS) print *,'Error in MPI type definition, mpi_tri_rmp '

    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    do i = 1, ntri
       call MPI_BCAST(mesh_element_rmp(i)%bnorm_vac, 1, mpi_tri_rmp, mpi_p_root, &
            MPI_COMM_WORLD, ierr)
       if(ierr .ne. MPI_SUCCESS) print *,'Error in BCAST of derived type mpi_tri_rmp'
    end do
    
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    return
  end subroutine broadcast_all
!----------------------------------------------------------------------------------------------
  subroutine average_mpi()

    use for_mpi, only : dummy8, dummy8_2, dummy8_1
    implicit none
    integer(i4b) :: nword, is
    double precision, dimension(nsorts) :: dummy0

! TEST: DO NOT UPDATE PROFILES    
if (.false.) then 
    allocate(dummy8(ntri,nsorts), dummy8_2(ntri,nsorts), dummy8_1(ntri))

    dummy0 = 0.d0
    call MPI_ALLREDUCE(time_local, dummy0, nsorts, MPI_DOUBLE_PRECISION, MPI_SUM, &
         MPI_COMM_WORLD, ierr)
    if(ierr .ne. MPI_SUCCESS) then 
       write(*,*)'time_local MPI error #',ierr
    endif
    time_l(:) = dummy0(:)*oneovernp
    time_local = time_local**2
    dummy0 = 0.d0
    call MPI_ALLREDUCE(time_local, dummy0, nsorts, MPI_DOUBLE_PRECISION, MPI_SUM, &
         MPI_COMM_WORLD, ierr)
    if(ierr .ne. MPI_SUCCESS) then 
       write(*,*)'time_local MPI error #',ierr
    endif
    time_local(:) = dummy0(:)*oneovernp
    D_time_l(:) = sqrt(time_local(:) - time_l(:)**2)

    nword = ntri*nsorts
    dummy8 = 0.d0
    call MPI_ALLREDUCE(Dm_local(1,1), dummy8(1,1), nword, MPI_DOUBLE_PRECISION, MPI_SUM, &
         MPI_COMM_WORLD, ierr)
    if(ierr .ne. MPI_SUCCESS) then 
       write(*,*)'Dm_local MPI error #',ierr
    endif
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    do is=1,nsorts
       mesh_element(:)%Dm_part(is) = dummy8(:,is)*oneovernp
    enddo
    Dm_local = Dm_local**2
    dummy8 = 0.d0
    call MPI_ALLREDUCE(Dm_local(1,1), dummy8(1,1), nword, MPI_DOUBLE_PRECISION, MPI_SUM, &
         MPI_COMM_WORLD, ierr)
    if(ierr .ne. MPI_SUCCESS) then 
       write(*,*)'Dm_local MPI error #',ierr
    endif
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    do is=1,nsorts
       Dm_local(:,is) = dummy8(:,is)*oneovernp
       D_Dm(:,is) = sqrt(max((Dm_local(:,is) - mesh_element(:)%Dm_part(is)**2),0.0d0))
    enddo
    dummy8 = 0.d0
    call MPI_ALLREDUCE(Dp_local(1,1), dummy8(1,1), nword, MPI_DOUBLE_PRECISION, MPI_SUM, &
         MPI_COMM_WORLD, ierr)
    if(ierr .ne. MPI_SUCCESS) then 
       write(*,*)'Dp_local MPI error #',ierr
    endif
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    do is=1,nsorts
       mesh_element(:)%D_part(is) = dummy8(:,is)*oneovernp
    enddo
    Dp_local = Dp_local**2
    dummy8 = 0.d0
    call MPI_ALLREDUCE(Dp_local(1,1), dummy8(1,1), nword, MPI_DOUBLE_PRECISION, MPI_SUM, &
         MPI_COMM_WORLD, ierr)
    if(ierr .ne. MPI_SUCCESS) then 
       write(*,*)'Dp_local MPI error #',ierr
    endif
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    do is=1,nsorts
       Dp_local(:,is) = dummy8(:,is)*oneovernp
       D_Dp(:,is) = sqrt(max((Dp_local(:,is) - mesh_element(:)%D_part(is)**2),0.0d0))
    enddo
    dummy8 = 0.d0
    call MPI_ALLREDUCE(Vp_local(1,1), dummy8(1,1), nword, MPI_DOUBLE_PRECISION, MPI_SUM, &
         MPI_COMM_WORLD, ierr)
    if(ierr .ne. MPI_SUCCESS) then 
       write(*,*)'Vp_local MPI error #',ierr
    endif
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    do is=1,nsorts
       mesh_element(:)%V_part(is) = dummy8(:,is)*oneovernp
    enddo
    Vp_local = Vp_local**2
    dummy8 = 0.d0
    call MPI_ALLREDUCE(Vp_local(1,1), dummy8(1,1), nword, MPI_DOUBLE_PRECISION, MPI_SUM, &
         MPI_COMM_WORLD, ierr)
    if(ierr .ne. MPI_SUCCESS) then 
       write(*,*)'D_Vp MPI error #',ierr
    endif
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    do is=1,nsorts
    Vp_local(:,is) = dummy8(:,is)*oneovernp
    D_Vp(:,is) = sqrt(max((Vp_local(:,is) - mesh_element(:)%V_part(is)**2),0.0d0))
    enddo
    dummy8 = 0.d0
    call MPI_ALLREDUCE(Tp_local(1,1), dummy8(1,1), nword, MPI_DOUBLE_PRECISION, MPI_SUM, &
         MPI_COMM_WORLD, ierr)
    if(ierr .ne. MPI_SUCCESS) then 
       write(*,*)'Tp_local MPI error #',ierr
    endif
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    do is=1,nsorts
       mesh_element(:)%T_part(is) = dummy8(:,is)*oneovernp
    enddo
    Tp_local = Tp_local**2
    dummy8 = 0.d0
    call MPI_ALLREDUCE(Tp_local(1,1), dummy8(1,1), nword, MPI_DOUBLE_PRECISION, MPI_SUM, &
         MPI_COMM_WORLD, ierr)
    if(ierr .ne. MPI_SUCCESS) then 
       write(*,*)'D_Tp MPI error #',ierr
    endif
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    do is=1,nsorts
       Tp_local(:,is) = dummy8(:,is)*oneovernp
       D_Tp(:,is) = sqrt(max((Tp_local(:,is) - mesh_element(:)%T_part(is)**2),0.0d0))
    enddo

    dummy8 = 0.d0
    call MPI_ALLREDUCE(ephi_local(1), dummy8_1(1), ntri, MPI_DOUBLE_PRECISION, MPI_SUM, &
         MPI_COMM_WORLD, ierr)
    if(ierr .ne. MPI_SUCCESS) then 
       write(*,*)'ephi_local MPI error #',ierr
    endif
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    mesh_element(:)%ePhi_tri = dummy8_1(:)*oneovernp
    ephi_local = ephi_local**2
    dummy8 = 0.d0
    call MPI_ALLREDUCE(ephi_local(1), dummy8_1(1), ntri, MPI_DOUBLE_PRECISION, MPI_SUM, &
         MPI_COMM_WORLD, ierr)
    if(ierr .ne. MPI_SUCCESS) then 
       write(*,*)'D_ephi MPI error #',ierr
    endif
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    ephi_local(:) = dummy8_1(:)*oneovernp
    D_ephi(:) = sqrt(max((ephi_local(:) - mesh_element(:)%ePhi_tri**2),0.0d0))

    deallocate(dummy8, dummy8_2)
    allocate(dummy8(n_therm,nsorts), dummy8_2(n_therm,nsorts))

    dummy8 = 0.d0
    dummy8_2 = D_therm
    nword = n_therm*nsorts
    call MPI_ALLREDUCE( dummy8_2(1,1), dummy8(1,1), nword, MPI_DOUBLE_PRECISION, MPI_SUM, &
         MPI_COMM_WORLD, ierr)
    if(ierr .ne. MPI_SUCCESS) then 
       write(*,*)'D_therm MPI error #',ierr
    endif
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    do is=1,nsorts
       D_therm(1:n_therm,is) = dummy8(1:n_therm,is)*oneovernp
    enddo
    dummy8_2 = dummy8_2**2
    dummy8 = 0.d0
    call MPI_ALLREDUCE( dummy8_2(1,1), dummy8(1,1), nword, MPI_DOUBLE_PRECISION, MPI_SUM, &
         MPI_COMM_WORLD, ierr)
    if(ierr .ne. MPI_SUCCESS) then 
       write(*,*)'D_Dt MPI error #',ierr
    endif
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    do is=1,nsorts
       dummy8_2(1:n_therm,is) = dummy8(1:n_therm,is)*oneovernp
       D_Dt(:,is) = sqrt(max((dummy8_2(:,is) - D_therm(1:n_therm,is)**2),0.0d0))
    enddo
    dummy8 = 0.d0
    dummy8_2 = V_therm
    call MPI_ALLREDUCE( dummy8_2(1,1), dummy8(1,1), nword, MPI_DOUBLE_PRECISION, MPI_SUM, &
         MPI_COMM_WORLD, ierr)
    if(ierr .ne. MPI_SUCCESS) then 
       write(*,*)'V_therm MPI error #',ierr
    endif


!print *,'T min = ', minval(T_therm), mype, '1'


    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    do is=1,nsorts
       V_therm(1:n_therm,is) = dummy8(1:n_therm,is)*oneovernp
    enddo
    dummy8_2 = dummy8_2**2
    do is=1,nsorts
       T_therm(1:n_therm,is) = T_therm(1:n_therm,is) +          &
            (dummy8_2(1:n_therm,is)    - V_therm(1:n_therm,is)**2)    & ! **2 above
            *erg2eV*amass(is)*one3rd*w_cell_V(1:n_therm,is)/w_cell_T(1:n_therm,is)
    enddo
    dummy8 = 0.d0
    call MPI_ALLREDUCE( dummy8_2(1,1), dummy8(1,1), nword, MPI_DOUBLE_PRECISION, MPI_SUM, &
         MPI_COMM_WORLD, ierr)
    if(ierr .ne. MPI_SUCCESS) then 
       write(*,*)'D_Vt MPI error #',ierr
    endif


!print *,'T min = ', minval(T_therm), mype, '2'

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    do is=1,nsorts
       dummy8_2(1:n_therm,is) = dummy8(1:n_therm,is)*oneovernp
       D_Vt(:,is) = sqrt(max((dummy8_2(:,is) - V_therm(1:n_therm,is)**2),0.0d0))
    enddo
    dummy8 = 0.d0
    dummy8_2 = T_therm
    call MPI_ALLREDUCE( dummy8_2(1,1), dummy8(1,1), nword, MPI_DOUBLE_PRECISION, MPI_SUM, &
         MPI_COMM_WORLD, ierr)
    if(ierr .ne. MPI_SUCCESS) then 
       write(*,*)'T_therm MPI error #',ierr
    endif
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    do is=1,nsorts
       T_therm(1:n_therm,is) = dummy8(1:n_therm,is)*oneovernp
    enddo
    dummy8_2 = dummy8_2**2
    dummy8 = 0.d0
    call MPI_ALLREDUCE( dummy8_2(1,1), dummy8(1,1), nword, MPI_DOUBLE_PRECISION, MPI_SUM, &
         MPI_COMM_WORLD, ierr)
    if(ierr .ne. MPI_SUCCESS) then 
       write(*,*)'D_Vt MPI error #',ierr
    endif
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    do is=1,nsorts
       dummy8_2(1:n_therm,is) = dummy8(1:n_therm,is)*oneovernp
       D_Tt(1:n_therm,is) = sqrt(max((dummy8_2(1:n_therm,is) - T_therm(1:n_therm,is)**2),0.0d0))
    enddo
    deallocate(dummy8, dummy8_2, dummy8_1)
endif ! TEST: NO PROFILE UPDATES
    
    allocate(cdummy8_1(ntri))
    allocate(cdummy8_2(ntri))

    do i=1,legs
       do is=1,nsorts
          cdummy8_1 = mesh_element_rmp(:)%currents(i,is)
          cdummy8_2 = (0.d0,0.d0)
          call MPI_ALLREDUCE(cdummy8_1(1), cdummy8_2(1), ntri, MPI_DOUBLE_COMPLEX, MPI_SUM, &
               MPI_COMM_WORLD, ierr)
          mesh_element_rmp(:)%currents(i,is) = cdummy8_2*oneovernp
       end do
    end do
    
    deallocate(cdummy8_1, cdummy8_2)
    
    
!print *,'T min = ', minval(T_therm), mype, '3'

!!$    dummy8 = 0.d0
!!$    dummy8_2(1:ntri) = flux_j(1:ntri,isort)
!!$    call MPI_ALLREDUCE( dummy8_2(1), dummy8(1), ntri, MPI_DOUBLE_PRECISION, MPI_SUM, &
!!$         MPI_COMM_WORLD, ierr)
!!$    if(ierr .ne. MPI_SUCCESS) then 
!!$       write(*,*)'flux_j MPI error #',ierr
!!$    endif
!!$    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
!!$    flux_j(1:ntri,isort) = dummy8(1:ntri)*oneovernp
!!$    dummy8_2(1:ntri) = dummy8_2(1:ntri)**2
!!$    dummy8 = 0.d0
!!$    call MPI_ALLREDUCE( dummy8_2(1), dummy8(1), ntri, MPI_DOUBLE_PRECISION, MPI_SUM, &
!!$         MPI_COMM_WORLD, ierr)
!!$    if(ierr .ne. MPI_SUCCESS) then 
!!$       write(*,*)'D_j MPI error #',ierr
!!$    endif
!!$    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
!!$    dummy8_2(1:ntri) = dummy8(1:ntri)*oneovernp
!!$    D_j(1:ntri,isort) = sqrt(max((dummy8_2(1:ntri) - flux_j(1:ntri,isort)**2),0.0d0))
!!$
!!$    dummy8 = 0.d0
!!$    dummy8_2(1:ntri) = flux_q(1:ntri,isort)
!!$    call MPI_ALLREDUCE( dummy8_2(1), dummy8(1), ntri, MPI_DOUBLE_PRECISION, MPI_SUM, &
!!$         MPI_COMM_WORLD, ierr)
!!$    if(ierr .ne. MPI_SUCCESS) then 
!!$       write(*,*)'flux_q MPI error #',ierr
!!$    endif
!!$    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
!!$    flux_q(1:ntri,isort) = dummy8(1:ntri)*oneovernp
!!$    dummy8_2(1:ntri) = dummy8_2(1:ntri)**2
!!$    dummy8 = 0.d0
!!$    call MPI_ALLREDUCE( dummy8_2(1), dummy8(1), ntri, MPI_DOUBLE_PRECISION, MPI_SUM, &
!!$         MPI_COMM_WORLD, ierr)
!!$    if(ierr .ne. MPI_SUCCESS) then 
!!$       write(*,*)'D_q MPI error #',ierr
!!$    endif
!!$    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
!!$    dummy8_2(1:ntri) = dummy8(1:ntri)*oneovernp
!!$    D_q(1:ntri,isort) = sqrt(max((dummy8_2(1:ntri) - flux_q(1:ntri,isort)**2),0.0d0))
!!$
!!$    dummy8 = 0.d0
!!$    dummy8_2(1:ntri) = tor_mom(1:ntri,isort)
!!$    call MPI_ALLREDUCE( dummy8_2(1), dummy8(1), ntri, MPI_DOUBLE_PRECISION, MPI_SUM, &
!!$         MPI_COMM_WORLD, ierr)
!!$    if(ierr .ne. MPI_SUCCESS) then 
!!$       write(*,*)'tor_mom MPI error #',ierr
!!$    endif
!!$    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
!!$    tor_mom(1:ntri,isort) = dummy8(1:ntri)*oneovernp
!!$    dummy8_2(1:ntri) = dummy8_2(1:ntri)**2
!!$    dummy8 = 0.d0
!!$    call MPI_ALLREDUCE( dummy8_2(1), dummy8(1), ntri, MPI_DOUBLE_PRECISION, MPI_SUM, &
!!$         MPI_COMM_WORLD, ierr)
!!$    if(ierr .ne. MPI_SUCCESS) then 
!!$       write(*,*)'D_Vt MPI error #',ierr
!!$    endif
!!$    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
!!$    dummy8_2(1:ntri) = dummy8(1:ntri)*oneovernp
!!$    D_tm(1:ntri,isort) = sqrt(max((dummy8_2(1:ntri) - tor_mom(1:ntri,isort)**2),0.0d0))
  end subroutine average_mpi
#endif
!!$!--------------------------------------------------------------------------------------------
!!$  subroutine plates(isort)
!!$   integer :: l, m, n_bou, isort
!!$   double precision, dimension(:), allocatable :: j_div_l, q_div_l, j_div_r, q_div_r
!!$   double precision, dimension(:), allocatable :: tm_div_l, tm_div_r
!!$   character*2 sort
!!$   character*30 filename
!!$
!!$   allocate(j_div_l(nboucell), q_div_l(nboucell), tm_div_l(nboucell))
!!$   allocate(j_div_r(nboucell), q_div_r(nboucell), tm_div_r(nboucell))
!!$
!!$   j_div_l(:) = 0.d0
!!$   q_div_l(:) = 0.d0
!!$   tm_div_l(:) = 0.d0
!!$   dp_j_l(isort) = 0.d0
!!$   dp_q_l(isort) = 0.d0
!!$   dp_tm_l(isort) = 0.d0
!!$   j_div_r(:) = 0.d0
!!$   q_div_r(:) = 0.d0
!!$   tm_div_r(:) = 0.d0
!!$   dp_j_r(isort) = 0.d0
!!$   dp_q_r(isort) = 0.d0
!!$   dp_tm_r(isort) = 0.d0
!!$   j_outer_wall(isort) = 0.d0
!!$   j_pfz_wall(isort) = 0.d0
!!$   q_outer_wall(isort) = 0.d0
!!$   q_pfz_wall(isort) = 0.d0
!!$   tm_outer_wall(isort) = 0.d0
!!$   tm_pfz_wall(isort) = 0.d0
!!$!print *,'flux_j*tau = ', sum(flux_j)*time_l(isort)
!!$   do l=1, ntri
!!$      n_bou = mesh_element(l)%iqq_gc
!!$      if(n_bou.ge.i_sol(1) .and. n_bou.le.i_sol(2)) then
!!$         j_outer_wall(isort) = j_outer_wall(isort) + flux_j(n_bou,isort)
!!$         q_outer_wall(isort) = q_outer_wall(isort) + flux_q(n_bou,isort)*5.d-1*amass(isort)
!!$         tm_outer_wall(isort) = tm_outer_wall(isort) + tor_mom(n_bou,isort)*amass(isort)
!!$      else if(n_bou.ge.i_pfz(1) .and. n_bou.le.i_pfz(2)) then
!!$         j_pfz_wall(isort) = j_pfz_wall(isort) + flux_j(n_bou,isort)
!!$         q_pfz_wall(isort) = q_pfz_wall(isort) + flux_q(n_bou,isort)*5.d-1*amass(isort)
!!$         tm_pfz_wall(isort) = tm_pfz_wall(isort) + tor_mom(n_bou,isort)*amass(isort)
!!$      else if(n_bou.ge.i_dpl(1) .and. n_bou.le.i_dpl(2)) then
!!$         dp_j_l(isort) = dp_j_l(isort) + flux_j(n_bou,isort)
!!$         dp_q_l(isort) = dp_q_l(isort) + flux_q(n_bou,isort)*5.d-1*amass(isort)
!!$         dp_tm_l(isort) =    dp_tm_l(isort) + tor_mom(n_bou,isort)*amass(isort)
!!$         j_div_l(n_bou) = j_div_l(n_bou) + flux_j(n_bou,isort)/S_bou(n_bou)
!!$         q_div_l(n_bou) = q_div_l(n_bou) + flux_q(n_bou,isort)/S_bou(n_bou)*5.d-1*amass(isort)
!!$         tm_div_l(n_bou) = tm_div_l(n_bou) + tor_mom(n_bou,isort)/S_bou(n_bou)*amass(isort)
!!$      else if(n_bou.ge.i_dpr(1) .and. n_bou.le.i_dpr(2)) then
!!$         dp_j_r(isort) = dp_j_r(isort) + flux_j(n_bou,isort)
!!$         dp_q_r(isort) = dp_q_r(isort) + flux_q(n_bou,isort)*5.d-1*amass(isort)
!!$         dp_tm_r(isort) =    dp_tm_r(isort) + tor_mom(n_bou,isort)*amass(isort)
!!$         j_div_r(n_bou) = j_div_r(n_bou) + flux_j(n_bou,isort)/S_bou(n_bou)
!!$         q_div_r(n_bou) = q_div_r(n_bou) + flux_q(n_bou,isort)/S_bou(n_bou)*5.d-1*amass(isort)
!!$         tm_div_r(n_bou) = tm_div_r(n_bou) + tor_mom(n_bou,isort)/S_bou(n_bou)*amass(isort)
!!$      end if
!!$   enddo
!!$
!!$   write(sort,'(i2)') isort
!!$   if(isort .le. 9) then
!!$      filename = 'RESULTS/leftplate.0' // trim(adjustl(sort))
!!$   else
!!$      filename = 'RESULTS/leftplate.'  // trim(adjustl(sort))
!!$   endif
!!$   open(1,file=filename)
!!$
!!$ !  open(1,file='RESULTS/leftplate.dat')
!!$   do l = i_dpl(1), i_dpl(2)
!!$      write(1,204) R_bou(l), Z_bou(l), j_div_l(l), q_div_l(l), tm_div_l(l)
!!$   enddo
!!$   close(1)
!!$
!!$   if(isort .le. 9) then
!!$      filename = 'RESULTS/rightplate.0' // trim(adjustl(sort))
!!$   else
!!$      filename = 'RESULTS/rightplate.'  // trim(adjustl(sort))
!!$   endif
!!$   open(1,file=filename)
!!$!   open(1,file='RESULTS/rightplate.dat')
!!$   do l = i_dpr(1), i_dpr(2)
!!$      write(1,204) R_bou(l), Z_bou(l), j_div_r(l), q_div_r(l), tm_div_r(l)
!!$   enddo
!!$   close(1)
!!$204 format(1000(ES12.5E2,1X))
!!$   return
!!$  end subroutine plates

!
  subroutine update_field(stat)
!
    !use from_nrtype
    !use mesh_mod, only: ntri, mesh_element, mesh_element_rmp
!
    integer, intent(out) :: stat
    call execute_command_line (&
         "cd RMP_EXCHANGE; PATH=/temp/ert/local/bin:$PATH /temp/ert/local/bin/FreeFem++ "//&
         "maxwell.edp maxwell.msh currents.dat 2 > freefem.out 2>&1; cd ..",&
    exitstat=stat)
  end subroutine update_field


end program k2d
