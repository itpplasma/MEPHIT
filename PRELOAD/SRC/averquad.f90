!
  implicit none
!
  integer :: i
  double precision, dimension(8) :: parcurr_e,parcurr_o
!
  open(1,file='Bn_flux.dat')
  open(2,file='Bn_flux_avquad.dat')
!
  do
    read (1,*,end=1) parcurr_e
    read (1,*,end=1) parcurr_o
    write (2,*) 0.5d0*(parcurr_e+parcurr_o) 
    write (2,*) 0.5d0*(parcurr_e+parcurr_o) 
  enddo
!
1 close(1)
  close(2)
  end
