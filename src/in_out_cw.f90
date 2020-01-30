module check_cw
  interface
     subroutine in_out_cw(R0, Z0, R, Z, out_of_cw, r_intsec, z_intsec, theta)
       use const, only : twopi
       implicit none
       real(kind=8), intent(in) :: R0, Z0 ! polar center
       real(kind=8), intent(inout) :: R, Z   ! point to be checked
       real(kind=8), intent(in), optional :: theta ! for call with prescribed angle
       real(kind=8), intent(out) :: r_intsec, z_intsec ! point of intersection (R_o,Z_o,R,Z) with CW
       logical, intent(out) :: out_of_cw
     end subroutine in_out_cw
  end interface
end module check_cw
! ---------------------------------------------------------------------------
subroutine in_out_cw(R0, Z0, R, Z, out_of_cw, r_intsec, z_intsec, theta)
! checks if a point inside of CW, finds coords of point of intersection
! of the line theta=const with CW
  use input_files, only : iunit, convexfile
  use const, only : twopi
  implicit none
  real(kind=8), intent(in) :: R0, Z0 ! polar center
  real(kind=8), intent(inout) :: R, Z   ! point to be checked, irrelevant if theta prescribed
  real(kind=8), intent(in), optional :: theta ! for call with prescribed angle
  real(kind=8), intent(out) :: r_intsec, z_intsec ! point of intersection (R_o,Z_o,R,Z) with CW
  logical, intent(out) :: out_of_cw
  integer i, j, numint, nrz ! number of points "convex wall" in input file
  integer, parameter :: nrzmx=1000 ! possible max. of nrz
  real(kind=8), dimension(:), allocatable :: r_cw, z_cw, rad_cw, zet_cw, tht_cw 
  logical :: firstcall=.true.
  real(kind=8) :: rho, tht, tht_1, rho_intsec, drox, dzox, rzox, zrox, drii, dzii, alph
  save firstcall, nrz, rad_cw, zet_cw, tht_cw, tht_1, drox, dzox, rzox, zrox

  if(firstcall) then
     firstcall = .false.
     allocate(r_cw(nrzmx), z_cw(nrzmx))
     nrz = 0
     r_cw = 0.
     z_cw = 0.
     open(iunit,file=trim(convexfile))
     do i=1,nrzmx
        read(iunit,*,END=10)r_cw(i),z_cw(i)
        nrz = nrz + 1
     enddo
10   continue
     close(iunit)

     allocate(rad_cw(1:nrz), zet_cw(1:nrz), tht_cw(1:nrz))

     rad_cw(1:nrz) = r_cw(1:nrz)
     zet_cw(1:nrz) = z_cw(1:nrz)
     deallocate(r_cw, z_cw)
  endif
  tht_1 =  atan2((zet_cw(1)-Z0),(rad_cw(1)-R0))
  do i=1,nrz-1
     tht_cw(i) = atan2((zet_cw(i)-Z0),(rad_cw(i)-R0)) - tht_1
     if (tht_cw(i) .ge. twopi) tht_cw(i) = tht_cw(i) - (int(tht_cw(i)/twopi))*twopi
     if (tht_cw(i) .lt. 0.) tht_cw(i) = tht_cw(i) + (int(abs(tht_cw(i))/twopi) +1)*twopi
  enddo
  tht_cw(1) = 0.d0
  tht_cw(nrz) = twopi   

  if(present(theta)) then
     tht = theta - tht_1
  else
     tht = atan2((Z-Z0),(R-R0)) - tht_1
  endif
     if (tht .ge. twopi) tht = tht - (int(tht/twopi))*twopi
     if (tht .lt. 0.) tht = tht + (int(abs(tht)/twopi) +1)*twopi

  do i=1,nrz-1
     if( tht.ge.tht_cw(i) .and. tht.le.tht_cw(i+1) ) then
        numint = i
        go to 1
     endif
  enddo
  stop 'in_out_cw: no theta interval found'
1 continue

  if(present(theta)) then
     rho = sqrt( (rad_cw(numint)-R0)**2 + (zet_cw(numint)-Z0)**2 )
     R = R0 + rho*cos(theta)
     Z = Z0 + rho*sin(theta)
!print *, R, Z, r0, z0, rho, sin(theta), cos(theta)
!pause
  endif

  drox = R0 - R
  dzox = Z0 - Z
  rzox = drox/dzox
  zrox = dzox/drox
  rho = sqrt( (R-R0)**2 + (Z-Z0)**2 )

  drii = rad_cw(numint+1) - rad_cw(numint)
  dzii = zet_cw(numint+1) - zet_cw(numint)
  if (abs(drii) .ge. abs(dzii)) then
!print *,'1'
     alph = dzii/drii
     r_intsec = (zet_cw(numint) - Z + R*zrox - rad_cw(numint)*alph)/(zrox - alph)
     z_intsec = Z + (r_intsec - R)*zrox
  else
!print *,'2'
     alph = drii/dzii
     z_intsec = (rad_cw(numint) - R + Z*rzox - zet_cw(numint)*alph)/(rzox - alph)
     r_intsec = R + (z_intsec - Z)*rzox           
  endif
  rho_intsec =  sqrt( (r_intsec-R0)**2 + (z_intsec-Z0)**2 )
!print *, 'in_out_cw', rho, rho_intsec, alph,i

  if(rho.ge.rho_intsec .and. (.not.present(theta))) then
     out_of_cw = .true.
  else
     out_of_cw = .false.
  endif
  return

end subroutine in_out_cw
