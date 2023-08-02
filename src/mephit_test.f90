program mephit_test

  use magdata_in_symfluxcoor_mod, only: load_magdata_in_symfluxcoord
  use hdf5_tools, only: h5_init, h5overwrite
  use mephit_conf, only: conf, config_read, conf_arr, logger, datafile
  use mephit_util, only: init_field, geqdsk_import_hdf5
  use mephit_mesh, only: equil, mesh, mesh_read, read_cache, check_mesh, &
       write_illustration_data, flux_func_cache_check, check_safety_factor, check_curr0
  use mephit_iter, only: mephit_deinit

  implicit none
  character(len = 1024) :: config_filename

  if (command_argument_count() >= 1) then
     call get_command_argument(1, config_filename)
  else
     error stop 'expected path to config file as first parameter'
  end if
  call config_read(conf, config_filename)
  call logger%init('-', conf%log_level, conf%quiet)
  call h5_init
  h5overwrite = .true.
  call geqdsk_import_hdf5(equil, datafile, 'equil')
  call init_field(equil)
  call mesh_read(mesh, datafile, 'mesh')
  call read_cache
  call load_magdata_in_symfluxcoord
  call conf_arr%read(conf%config_file, mesh%m_res_min, mesh%m_res_max)

  if (conf%check_mesh) call check_mesh
  call write_illustration_data(5, 8, 256, 256)
  call flux_func_cache_check
  call check_safety_factor
  call check_curr0

  call mephit_deinit
       
end program mephit_test
