program vacfield

  use iso_fortran_env, only: dp => real64
  use hdf5_tools, only: h5_init, h5_deinit
  use magdif_conf, only: conf, conf_arr, magdif_config_read, log, magdif_log, &
       decorate_filename
  use magdif_mesh, only: mesh, read_mesh
  use magdif_util, only: initialize_globals
  use magdif_pert, only: compute_Bn_nonres, compute_kilca_vac_coeff, compute_kilca_vacuum, &
       check_kilca_vacuum, check_redundant_edges, check_div_free, &
       write_vector_dof, write_vector_plot

  implicit none

  character(len = 1024) :: config_file
  complex(dp), allocatable :: Bnflux(:, :), Bnphi(:)

  if (command_argument_count() >= 1) then
     call get_command_argument(1, config_file)
  else
     error stop 'expected path to config file as first parameter'
  endif
  call magdif_config_read(conf, config_file)
  call conf_arr%read(config_file, conf%m_min, conf%m_max)
  log = magdif_log('-', conf%log_level, conf%quiet)

  call h5_init
  call read_mesh
  allocate(Bnflux(mesh%ntri, 3))
  allocate(Bnphi(mesh%ntri))
  if (conf%kilca_scale_factor /= 0) then
     call compute_kilca_vac_coeff
     call compute_kilca_vacuum(Bnflux, Bnphi)
     call check_kilca_vacuum
  else
     if (conf%nonres) then
        call compute_Bn_nonres(Bnflux, Bnphi)
     else
        call initialize_globals(mesh%R_O, mesh%Z_O)
        call compute_Bnvac(Bnflux, Bnphi)
     end if
  end if
  call check_redundant_edges(Bnflux, .false., 'vacuum B_n')
  call check_div_free(Bnflux, Bnphi, mesh%n, 1d-9, 'vacuum B_n')
  call write_vector_dof(Bnflux, Bnphi, conf%Bn_vac_file)
  call write_vector_plot(Bnflux, Bnphi, decorate_filename(conf%Bn_vac_file, 'plot_', ''))
  call h5_deinit
  if (allocated(Bnflux)) deallocate(Bnflux)
  if (allocated(Bnphi)) deallocate(Bnphi)

contains

  subroutine compute_Bnvac(Bnflux, Bnphi)
    use iso_c_binding, only: c_long
    use magdif_conf, only: conf
    use magdif_util, only: gauss_legendre_unit_interval, imun
    use magdif_mesh, only: mesh
    complex(dp), intent(out) :: Bnflux(:, :), Bnphi(:)
    integer, parameter :: order = 2
    integer :: ktri, ke, k
    real(dp) :: R, Z, edge_R, edge_Z, node_R(4), node_Z(4)
    real(dp), dimension(order) :: points, weights
    complex(dp) :: B_R, B_Z

    Bnflux = (0d0, 0d0)
    Bnphi = (0d0, 0d0)
    call gauss_legendre_unit_interval(order, points, weights)
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
             Bnflux(ktri, ke) = Bnflux(ktri, ke) + &
                  weights(k) * (B_R * edge_Z - B_Z * edge_R) * R
          end do
       end do
       ! toroidal flux via zero divergence
       Bnphi(ktri) = imun / mesh%n * sum(Bnflux(ktri, :)) / mesh%area(ktri)
    end do
  end subroutine compute_Bnvac
end program vacfield
