!
  implicit none
!
  integer :: npoi,i
!
  integer :: mpol,ntor,isw_Phi_m
  integer :: nskip
  double precision :: am_i,z_i,Rtor,Btor,rplas,sepval,x,skip_radius,ev
  double precision, dimension(:), allocatable :: psi,qsaf,bcovar_phi,Phi_0,avR2nabpsi2, &
                                                 dens_e,dens_i,temp_e,temp_i,anu_e,anu_i,r
  double complex,   dimension(:), allocatable :: bpsi_over_bphi,parcur_over_b0,Phi_m
!
  double precision, dimension(:), allocatable :: dPhi_0_dpsi
  double complex,   dimension(:), allocatable :: Phi_m_aligned
!
! Switch for potential perturbation in j_parallel (0 - skip, 1 -use):
  isw_Phi_m=1
! poloidal mode number:
  mpol=-6
! toroidal mode number:
  ntor=2
! mass number of ions:
  am_i=2.d0
! charge number of ions:
  z_i=1.d0
!
  open(1,file='equil_r_q_psi.dat')
  read(1,*)
  read(1,*)
  read(1,*)
  npoi=0
  do
    read(1,*,end=1) x
    npoi=npoi+1
  enddo
1 close(1)
!
  allocate(psi(npoi),qsaf(npoi),bcovar_phi(npoi),Phi_0(npoi),avR2nabpsi2(npoi),         &
           dens_e(npoi),dens_i(npoi),temp_e(npoi),temp_i(npoi),anu_e(npoi),anu_i(npoi), &
           bpsi_over_bphi(npoi),parcur_over_b0(npoi),Phi_m(npoi),r(npoi))
!
  open(1,file='equil_r_q_psi.dat')
  read(1,*) 
  read(1,*)
  read(1,*)
  do i=1,npoi
    read(1,*) r(i),qsaf(i),psi(i)
  enddo
  close(1)
!
  open(1,file='btor_rbig.dat')
  read(1,*) Btor,Rtor
  close(1)
!
  bcovar_phi=Btor*Rtor
!
  call first_deriv(npoi,r,psi,avR2nabpsi2)
!
  avR2nabpsi2=(Rtor*avR2nabpsi2)**2
!
  skip_radius=10.d0 !make profiles parabolic within this radius
  nskip=0
!
  open(1,file='profs.dat')
  do i=1,npoi
    read(1,*) x,dens_e(i),temp_e(i),temp_i(i),Phi_0(i)
    if(nskip.eq.0.and.x.gt.skip_radius) nskip=i
  enddo
  close(1)
!
  dens_e(1:nskip)=dens_e(nskip)                   &
                 +(dens_e(nskip+1)-dens_e(nskip)) &
                 /(r(nskip+1)**2-r(nskip)**2)     &
                 *(r(1:nskip)**2-r(nskip)**2)
  temp_e(1:nskip)=temp_e(nskip)                   &
                 +(temp_e(nskip+1)-temp_e(nskip)) &
                 /(r(nskip+1)**2-r(nskip)**2)     &
                 *(r(1:nskip)**2-r(nskip)**2)
  temp_i(1:nskip)=temp_i(nskip)                   &
                 +(temp_i(nskip+1)-temp_i(nskip)) &
                 /(r(nskip+1)**2-r(nskip)**2)     &
                 *(r(1:nskip)**2-r(nskip)**2)
  Phi_0(1:nskip) =Phi_0(nskip)                    &
                 +(Phi_0(nskip+1)-Phi_0(nskip))   &
                 /(r(nskip+1)**2-r(nskip)**2)     &
                 *(r(1:nskip)**2-r(nskip)**2)
  
! collision frequencies:
  ev=1.6022d-12
  anu_e = 7.7d-6*(1.d0+Z_i)*dens_e/(temp_e/ev)**1.5d0                  &
        * (24.d0-log(sqrt(dens_e)*ev/temp_e))
  anu_i = 1.8d-7*Z_i**3/sqrt(am_i)*dens_e/(temp_i/ev)**1.5d0           &
        * (23.d0-log(Z_i**2.5d0*sqrt(2.d0*dens_e)*(ev/temp_i)**1.5d0))
!
! radial magnetic field perturbation:
  bpsi_over_bphi=(1.d0,0.d0)
!
!
  call response_current(isw_Phi_m,mpol,ntor,npoi,am_i,z_i,Rtor,     &
                        psi,qsaf,bcovar_phi,Phi_0,avR2nabpsi2,      &
                        dens_e,temp_e,temp_i,anu_e,anu_i,           &
                        bpsi_over_bphi,parcur_over_b0,Phi_m)
!
!
! compute potential aligned with perturbed flux surfaces:
  allocate(dPhi_0_dpsi(npoi),Phi_m_aligned(npoi))
!
  call first_deriv(npoi,psi,Phi_0,dPhi_0_dpsi)
!
  Phi_m_aligned=(0.d0,1.d0)*bpsi_over_bphi*qsaf*dPhi_0_dpsi/(mpol+qsaf*ntor)
!
  open(1,file='pert_potential.dat')
  open(2,file='pert_current.dat')
  open(3,file='collision_frecs.dat')
  do i=1,npoi
!
! pertubation of the potentiali (solution to Eq.(95)), aligned potential (Eq.(27)):
    write(1,*) r(i),real(Phi_m(i)),aimag(Phi_m(i)),real(Phi_m_aligned(i)),aimag(Phi_m_aligned(i)),psi(i)
!
! paralel current density divided by mod-B (Eq.(97)):
    write(2,*) r(i),real(parcur_over_b0(i)),aimag(parcur_over_b0(i)),psi(i)
!
! collision frequencies:
    write(3,*) r(i),anu_e(i),anu_i(i),psi(i)
!
  enddo
  close(1)
  close(2)
  close(3)
!
  end
