module field_interface
  use from_nrtype
  use mesh_mod, only: ntri, mesh_element, mesh_element_rmp

  contains

  subroutine read_field
    double precision :: dummy_re,dummy_im
    double complex   :: bnorm_vac,bnorm_plas
    
    open(1,file='RMP_EXCHANGE/hpsi_vac.dat')
    open(2,file='RMP_EXCHANGE/hpsi.dat')
    do i=1,ntri
      mesh_element_rmp(i)%currents(:,:) = (0.d0,0.d0)
      read (1,*) dummy_re,dummy_im
      bnorm_vac=cmplx(dummy_re,dummy_im)
      read (2,*) dummy_re,dummy_im
      bnorm_plas=cmplx(dummy_re,dummy_im)
      !bnorm_plas=(0.d0,0.d0)
      mesh_element_rmp(i)%bnorm_times_thermforces=mesh_element(i)%thermforces &
&       *(bnorm_vac+bnorm_plas)
    enddo
    close(1)
    close(2)
  end subroutine read_field

  subroutine update_field(stat)
    integer, intent(out) :: stat
    call execute_command_line (&
         "cd RMP_EXCHANGE; /temp/ert/local/bin/FreeFem++ "//&
         "maxwell.edp maxwell.msh currents.dat 2 > freefem.out 2>&1; cd ..",&
    exitstat=stat)
  end subroutine update_field
end module field_interface
