program test_maxwell_boundary
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use maxwell_boundary_m, only: validate_outer_polygon, read_mesh_exterior_scale
    use mephit_mesh, only: mesh_t, mesh_read
    use mephit_conf, only: conf
    use hdf5_tools, only: HID_T, h5_create, h5_add, h5_close, h5_open
    implicit none
    real(dp) :: inner_R(128), inner_Z(128), outer_R(32), outer_Z(32)
    real(dp) :: angle, scale, effective
    integer :: i, level
    integer(HID_T) :: h5id
    character(len=32) :: probe
    type(mesh_t) :: cached_mesh

    call get_command_argument(1, probe)
    do i = 1, 128
        angle = 2.0_dp*acos(-1.0_dp)*real(i-1, dp)/128.0_dp
        inner_R(i) = 620.0_dp+62.0_dp*cos(angle)
        inner_Z(i) = 17.0_dp+62.0_dp*sin(angle)
    end do
    do level = 2, 6
        scale = real(level, dp)
        if (probe == 'near') scale = 1.001_dp
        if (probe == 'axis') scale = 11.0_dp
        if (probe == 'nan') scale = ieee_value(1.0_dp, ieee_quiet_nan)
        do i = 1, 32
            angle = 2.0_dp*acos(-1.0_dp)*real(i-1, dp)/32.0_dp
            outer_R(i) = 620.0_dp+62.0_dp*scale*cos(angle)
            outer_Z(i) = 17.0_dp+62.0_dp*scale*sin(angle)
        end do
        call validate_outer_polygon(inner_R, inner_Z, outer_R, outer_Z, scale)
    end do
    call h5_create('maxwell-scale-replay.h5', h5id)
    ! A conflicting requested config must never overwrite the mesh provenance.
    call h5_add(h5id, 'config/maxwell_outer_scale', 3.0_dp)
    if (probe == 'native-present') call h5_add(h5id, 'mesh/maxwell_outer_scale', 2.0_dp)
    if (probe == 'native-new') call h5_add(h5id, 'mesh/maxwell_outer_scale', 3.0_dp)
    call h5_close(h5id)
    if (index(probe, 'native-') == 1) then
        conf%maxwell_outer_scale = 3.0_dp
        if (probe == 'native-new') conf%maxwell_outer_scale = 2.0_dp
        call mesh_read(cached_mesh, 'maxwell-scale-replay.h5', 'mesh')
        error stop 'Mismatched native cached mesh was accepted'
    end if
    call h5_open('maxwell-scale-replay.h5', h5id)
    effective = 9.0_dp
    call read_mesh_exterior_scale(h5id, 'mesh', 2.0_dp, effective)
    if (effective /= 2.0_dp) error stop 'Legacy scale must reset to two'
    call h5_close(h5id)
    call h5_create('maxwell-scale-two.h5', h5id)
    call h5_add(h5id, 'mesh/maxwell_outer_scale', 2.0_dp)
    call h5_close(h5id)
    call h5_open('maxwell-scale-two.h5', h5id)
    call read_mesh_exterior_scale(h5id, 'mesh', 2.0_dp, effective)
    if (effective /= 2.0_dp) error stop 'Stored default scale was not retained'
    call h5_close(h5id)
    call h5_create('maxwell-scale-present.h5', h5id)
    call h5_add(h5id, 'mesh/maxwell_outer_scale', 3.0_dp)
    call h5_close(h5id)
    call h5_open('maxwell-scale-present.h5', h5id)
    call read_mesh_exterior_scale(h5id, 'mesh', 3.0_dp, effective)
    if (effective /= 3.0_dp) error stop 'Matching stored scale was not retained'
    call h5_close(h5id)
    print *, 'Outer polygons and HDF5 mesh-scale replay passed'
end program test_maxwell_boundary
