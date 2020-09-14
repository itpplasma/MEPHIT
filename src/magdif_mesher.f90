program magdif_preprocess

  use hdf5_tools, only: h5_init, h5_deinit
  use magdif_conf, only: conf, magdif_config_read, conf_arr, log, magdif_log
  use magdif_mesh, only: generate_mesh
  use magdif_pert, only: compute_kilca_vac_coeff, compute_kilca_vacuum, &
       check_kilca_vacuum, check_RT0
  use mesh_mod, only: mesh_point, mesh_element, mesh_element_rmp

  implicit none

  character(len = 1024) :: config_file, unprocessed_geqdsk

  if (command_argument_count() >= 1) then
     call get_command_argument(1, config_file)
     if (command_argument_count() >= 2) then
        call get_command_argument(2, unprocessed_geqdsk)
     else
        error stop 'expected path to unprocessed G EQDSK file as second parameter'
     end if
  else
     error stop 'expected path to magdif config file as first parameter'
  endif

  call h5_init
  call magdif_config_read(conf, config_file)
  call conf_arr%read(config_file, conf%m_min, conf%m_max)
  log = magdif_log('-', conf%log_level, conf%quiet)
  call generate_mesh(unprocessed_geqdsk)
  call h5_deinit
  ! TODO: proper cleanup
  if (allocated(mesh_element_rmp)) deallocate(mesh_element)
  if (allocated(mesh_element)) deallocate(mesh_element)
  if (allocated(mesh_point)) deallocate(mesh_point)
end program magdif_preprocess
