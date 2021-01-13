program vacfield

  use iso_fortran_env, only: dp => real64
  use hdf5_tools, only: h5_init, h5_deinit, h5overwrite
  use magdif_conf, only: conf, conf_arr, magdif_config_read, log, magdif_log, datafile
  use magdif_mesh, only: mesh, read_mesh_cache, equil
  use magdif_util, only: get_field_filenames, init_field
  use magdif_pert, only: compute_Bn_nonres, compute_kilca_vac_coeff, compute_kilca_vacuum, &
       check_kilca_vacuum, RT0_check_redundant_edges, RT0_check_div_free, RT0_t, RT0_init, &
       RT0_deinit, RT0_write, check_RT0, debug_fouriermodes, debug_Bnvac_rectplot, &
       debug_Bmnvac, compute_Bnvac

  implicit none

  character(len = 1024) :: config_file, gfile, pfile, convexfile
  type(RT0_t) :: Bn

  if (command_argument_count() >= 1) then
     call get_command_argument(1, config_file)
  else
     error stop 'expected path to config file as first parameter'
  endif
  call magdif_config_read(conf, config_file)
  call conf_arr%read(config_file, conf%m_min, conf%m_max)
  log = magdif_log('-', conf%log_level, conf%quiet)

  call h5_init
  h5overwrite = .true.
  call read_mesh_cache
  call RT0_init(Bn, mesh%ntri)
  if (conf%kilca_scale_factor /= 0) then
     call compute_kilca_vac_coeff
     call compute_kilca_vacuum(Bn)
     call check_kilca_vacuum
     call check_RT0(Bn)
  else
     if (conf%nonres) then
        call compute_Bn_nonres(Bn)
     else
        call get_field_filenames(gfile, pfile, convexfile)
        call init_field(gfile, pfile, convexfile)
        call equil%read(trim(gfile))
        call equil%classify
        call compute_Bnvac(Bn)
        call debug_Bnvac_rectplot
        call debug_Bmnvac
        call debug_fouriermodes
     end if
  end if
  call RT0_check_redundant_edges(Bn, 'Bnvac')
  call RT0_check_div_free(Bn, mesh%n, 1d-9, 'Bnvac')
  call RT0_write(Bn, datafile, 'Bnvac', 'magnetic field (vacuum)', 'G', 1)
  call RT0_deinit(Bn)
  call h5_deinit

end program vacfield
