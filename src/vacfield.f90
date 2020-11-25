program vacfield

  use iso_fortran_env, only: dp => real64
  use hdf5_tools, only: h5_init, h5_deinit
  use magdif_conf, only: conf, conf_arr, magdif_config_read, log, magdif_log, datafile
  use magdif_mesh, only: mesh, read_mesh_cache
  use magdif_util, only: get_field_filenames, init_field
  use magdif_pert, only: compute_Bn_nonres, compute_kilca_vac_coeff, compute_kilca_vacuum, &
       check_kilca_vacuum, RT0_check_redundant_edges, RT0_check_div_free, RT0_t, RT0_init, &
       RT0_deinit, RT0_write, check_RT0

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
        call compute_Bnvac(Bn)
     end if
  end if
  call RT0_check_redundant_edges(Bn, 'Bnvac')
  call RT0_check_div_free(Bn, mesh%n, 1d-9, 'Bnvac')
  call RT0_write(Bn, datafile, 'Bnvac', 'magnetic field (vacuum)', 'G', 1)
  call RT0_deinit(Bn)
  call h5_deinit

contains

  subroutine compute_Bnvac(Bn)
    use iso_c_binding, only: c_long
    use field_mod, only: ipert, iequil
    use field_c_mod, only: icall_c
    use magdif_conf, only: conf
    use magdif_util, only: gauss_legendre_unit_interval, imun
    use magdif_mesh, only: mesh
    type(RT0_t), intent(inout) :: Bn
    integer, parameter :: order = 2
    integer :: ktri, ke, k
    real(dp) :: dum, R, Z, edge_R, edge_Z, node_R(4), node_Z(4)
    real(dp), dimension(order) :: points, weights
    complex(dp) :: B_R, B_Z

    call gauss_legendre_unit_interval(order, points, weights)
    ! initialize vacuum field
    ipert = 1
    iequil = 0
    icall_c = -1
    call field_c(0d0, 0d0, 0d0, dum, dum, dum, dum, dum, dum, dum, dum, dum, dum, dum, dum)
    do ktri = 1, mesh%ntri
       node_R = mesh%node_R([mesh%tri_node(:, ktri), mesh%tri_node(1, ktri)])
       node_Z = mesh%node_Z([mesh%tri_node(:, ktri), mesh%tri_node(1, ktri)])
       do ke = 1, 3
          edge_R = node_R(ke + 1) - node_R(ke)
          edge_Z = node_Z(ke + 1) - node_Z(ke)
          do k = 1, order
             R = node_R(ke) * points(k) + node_R(ke + 1) * points(order - k + 1)
             Z = node_Z(ke) * points(k) + node_Z(ke + 1) * points(order - k + 1)
             call spline_bpol_n(conf%n, R, Z, B_R, B_Z)
             Bn%DOF(ke, ktri) = Bn%DOF(ke, ktri) + &
                  weights(k) * (B_R * edge_Z - B_Z * edge_R) * R
          end do
       end do
       ! toroidal flux via zero divergence
       Bn%comp_phi(ktri) = imun / mesh%n * sum(Bn%DOF(:, ktri)) / mesh%area(ktri)
    end do
    ! reset to equilibrium field
    ipert = 0
    iequil = 1
  end subroutine compute_Bnvac
end program vacfield
