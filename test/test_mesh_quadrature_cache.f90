program test_mesh_quadrature_cache
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use mephit_mesh, only: mesh_t, mesh_read, mesh_write, mesh_deinit
    use hdf5_tools, only: HID_T, h5_create, h5_close
    implicit none
    type(mesh_t) :: original, restored
    integer(HID_T) :: handle
    real(dp) :: expected_R(3), expected_Z(3)

    ! One physical triangle has three edges and three interior quadrature nodes.
    ! Its vertices are (10,0), (12,0), (10,2), in centimeters.
    original%npoint = 3
    original%ntri = 1
    original%nedge = 3
    original%nflux = 1
    original%n = 3
    original%m_res_min = 4
    original%m_res_max = 4
    original%GL_order = 2
    original%GL2_order = 3
    original%R_O = 10.0_dp
    original%Z_O = 0.0_dp
    original%R_X = 0.0_dp
    original%Z_X = 0.0_dp
    original%R_min = 10.0_dp
    original%R_max = 12.0_dp
    original%Z_min = 0.0_dp
    original%Z_max = 2.0_dp
    allocate (original%node_R(3), original%node_Z(3), &
        original%node_theta_flux(3), original%node_theta_geom(3), &
        original%tri_node(3, 1), original%tri_node_F(1), &
        original%tri_theta2pi(1), original%tri_theta_extent(2, 1), &
        original%tri_RZ_extent(2, 2, 1), original%orient(1), &
        original%edge_node(2, 3), original%edge_tri(2, 3), original%tri_edge(3, 1), &
        original%mid_R(3), original%mid_Z(3), original%edge_R(3), original%edge_Z(3), &
        original%area(1), original%cntr_R(1), original%cntr_Z(1), &
        original%res_modes(1), original%kp_max(1), original%kp_low(1), &
        original%kt_max(1), original%kt_low(1), original%gpec_jacfac(1), &
        original%GL_weights(2), original%GL_R(2, 3), original%GL_Z(2, 3), &
        original%GL2_weights(3), original%GL2_R(3, 1), original%GL2_Z(3, 1))
    original%node_R = [10.0_dp, 12.0_dp, 10.0_dp]
    original%node_Z = [0.0_dp, 0.0_dp, 2.0_dp]
    original%node_theta_flux = [0.0_dp, 0.0_dp, 0.0_dp]
    original%node_theta_geom = [0.0_dp, 0.0_dp, 0.0_dp]
    original%tri_node = reshape([1, 2, 3], [3, 1])
    original%tri_node_F = [1]
    original%tri_theta2pi = [.false.]
    original%tri_theta_extent = reshape([0.0_dp, 0.0_dp], [2, 1])
    original%tri_RZ_extent = reshape([10.0_dp, 12.0_dp, 0.0_dp, 2.0_dp], [2, 2, 1])
    original%orient = [.true.]
    original%edge_node = reshape([1, 2, 2, 3, 3, 1], [2, 3])
    original%edge_tri = reshape([1, 0, 1, 0, 1, 0], [2, 3])
    original%tri_edge = reshape([1, 2, 3], [3, 1])
    original%mid_R = [11.0_dp, 11.0_dp, 10.0_dp]
    original%mid_Z = [0.0_dp, 1.0_dp, 1.0_dp]
    original%edge_R = [2.0_dp, -2.0_dp, 0.0_dp]
    original%edge_Z = [0.0_dp, 2.0_dp, -2.0_dp]
    original%area = [2.0_dp]
    original%cntr_R = [32.0_dp/3.0_dp]
    original%cntr_Z = [2.0_dp/3.0_dp]
    allocate (original%res_ind(4:4), original%psi_res(4:4), &
        original%rad_norm_res(4:4), original%Delta_psi_res_curr(4:4))
    original%res_ind = 1
    original%psi_res = 1.0_dp
    original%rad_norm_res = 1.0_dp
    original%res_modes = [4]
    original%Delta_psi_res_curr = 0.0_dp
    allocate (original%damping(0:1), original%avg_R2gradpsi2(0:1))
    original%damping = 0.0_dp
    original%avg_R2gradpsi2 = 0.0_dp
    original%kp_max = [2]
    original%kp_low = [1]
    original%kt_max = [1]
    original%kt_low = [1]
    original%gpec_jacfac = [1.0_dp]
    original%GL_weights = [0.5_dp, 0.5_dp]
    original%GL_R = reshape([10.5_dp, 11.5_dp, 11.5_dp, 10.5_dp, &
        10.0_dp, 10.0_dp], [2, 3])
    original%GL_Z = reshape([0.0_dp, 0.0_dp, 0.5_dp, 1.5_dp, &
        1.5_dp, 0.5_dp], [2, 3])
    ! Degree-two triangle rule: barycentric nodes (2/3,1/6,1/6) and permutations.
    expected_R = [31.0_dp/3.0_dp, 34.0_dp/3.0_dp, 31.0_dp/3.0_dp]
    expected_Z = [1.0_dp/3.0_dp, 1.0_dp/3.0_dp, 4.0_dp/3.0_dp]
    original%GL2_weights = [1.0_dp/3.0_dp, 1.0_dp/3.0_dp, 1.0_dp/3.0_dp]
    original%GL2_R = reshape(expected_R, [3, 1])
    original%GL2_Z = reshape(expected_Z, [3, 1])
    call h5_create('triangle-cache.h5', handle)
    call h5_close(handle)
    call mesh_write(original, 'triangle-cache.h5', 'mesh')
    call mesh_read(restored, 'triangle-cache.h5', 'mesh')
    if (restored%ntri /= 1) error stop 'Triangle count changed'
    if (restored%nedge /= 3) error stop 'Edge count changed'
    if (any(shape(restored%GL2_R) /= [3, 1])) error stop 'Triangle R shape'
    if (any(shape(restored%GL2_Z) /= [3, 1])) error stop 'Triangle Z shape'
    if (any(shape(restored%GL_R) /= [2, 3])) error stop 'Edge R shape'
    if (any(shape(restored%GL_Z) /= [2, 3])) error stop 'Edge Z shape'
    if (maxval(abs(restored%GL2_R(:, 1) - expected_R)) > 1e-12_dp) &
        error stop 'Interior triangle R nodes changed'
    if (maxval(abs(restored%GL2_Z(:, 1) - expected_Z)) > 1e-12_dp) &
        error stop 'Interior triangle Z nodes changed'
    if (maxval(abs(restored%GL_R - original%GL_R)) > 1e-12_dp) &
        error stop 'Edge R nodes changed'
    if (maxval(abs(restored%GL_Z - original%GL_Z)) > 1e-12_dp) &
        error stop 'Edge Z nodes changed'
    ! The normalized rule integrates the triangle centroid independently.
    if (abs(sum(restored%GL2_weights*restored%GL2_R(:, 1)) - 32.0_dp/3.0_dp) &
        > 1e-12_dp) error stop 'Triangle centroid R integral changed'
    if (abs(sum(restored%GL2_weights*restored%GL2_Z(:, 1)) - 2.0_dp/3.0_dp) &
        > 1e-12_dp) error stop 'Triangle centroid Z integral changed'
    call mesh_deinit(original)
    call mesh_deinit(restored)
    print *, 'Triangle and edge quadrature survive native mesh write/read'
end program test_mesh_quadrature_cache
