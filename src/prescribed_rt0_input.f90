module prescribed_rt0_input
    use iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use hdf5_tools, only: HID_T, HSIZE_T, SIZE_T, h5_open, h5_close, h5_get, &
        h5_get_dataset_info, h5_open_rw, h5_add, h5_create_parent_groups
    implicit none
    private
    public :: read_prescribed_rt0, check_prescribed_toroidal, write_prescribed_source
    character(len=*), parameter :: group = 'prescribed_rt0/'
contains
    subroutine read_prescribed_rt0(file, n, r, z, en, tn, te, orient, area, &
            flux, phi, bound)
        character(len=*), intent(in) :: file
        integer, intent(in) :: n, en(:, :), tn(:, :), te(:, :)
        logical, intent(in) :: orient(:)
        real(dp), intent(in) :: r(:), z(:), area(:)
        complex(dp), intent(out) :: flux(:), phi(:)
        real(dp), intent(out) :: bound(:)
        integer(HID_T) :: fid
        integer :: value
        integer, allocatable :: imat(:, :), flags(:)
        real(dp), allocatable :: values(:), re(:), im(:)
        character(len=1024) :: label

        if (n == 0) error stop 'Prescribed RT0 requires nonzero signed n'
        if (size(r) /= size(z)) error stop 'Prescribed RT0 mesh coordinate shape'
        if (size(en, 1) /= 2) error stop 'Prescribed RT0 edge topology shape'
        if (size(tn, 1) /= 3) error stop 'Prescribed RT0 triangle topology shape'
        if (any(shape(tn) /= shape(te))) error stop 'Prescribed RT0 triangle shapes'
        if (size(area) /= size(tn, 2)) error stop 'Prescribed RT0 area shape'
        if (size(orient) /= size(area)) error stop 'Prescribed RT0 orientation shape'
        if (size(flux) /= size(en, 2)) error stop 'Prescribed RT0 flux shape'
        if (size(phi) /= size(area)) error stop 'Prescribed RT0 toroidal shape'
        if (size(bound) /= size(area)) error stop 'Prescribed RT0 bound shape'
        if (.not. all(ieee_is_finite(r))) error stop 'Nonfinite mesh R'
        if (.not. all(ieee_is_finite(z))) error stop 'Nonfinite mesh Z'
        if (.not. all(ieee_is_finite(area))) error stop 'Nonfinite mesh area'
        if (any(r <= 0) .or. any(area <= 0)) error stop 'Invalid prescribed RT0 mesh'
        call h5_open(file, fid)
        call expect_integer('schema_version', 1)
        call expect_integer('signed_n', n)
        call expect_integer('nnode', size(r))
        call expect_integer('nedge', size(en, 2))
        call expect_integer('ntri', size(area))
        call expect_integer('index_base', 1)
        call expect_integer('phase_sign', 1)
        call expect_integer('real_field_multiplier', 2)
        call expect_label('coordinate_unit', 'cm')
        call expect_label('edge_flux_unit', 'G*cm^2')
        call expect_label('field_unit', 'G')
        call require_label('mesh_image_sha256')
        call require_label('source_input_sha256')
        call require_label('generator_version_or_sha256')
        call require_label('evaluation_receipt_sha256')
        call require_label('source_description')
        call check_real('node_R_cm', r)
        call check_real('node_Z_cm', z)
        call check_real('area_cm2', area)
        call check_int('edge_node', en)
        call check_int('tri_node', tn)
        call check_int('tri_edge', te)
        allocate(flags(size(orient)))
        call require_shape(fid, 'orient', [size(flags)])
        call h5_get(fid, group // 'orient', flags)
        if (any(flags /= merge(1, 0, orient))) error stop 'Prescribed RT0 orientation mismatch'
        allocate(re(size(flux)), im(size(flux)))
        call read_real('edge_flux_real', re)
        call read_real('edge_flux_imag', im)
        flux = cmplx(re, im, dp)
        deallocate(re, im)
        allocate(re(size(phi)), im(size(phi)))
        call read_real('phi_area_reference_real', re)
        call read_real('phi_area_reference_imag', im)
        phi = cmplx(re, im, dp)
        call read_real('phi_error_bound_G', bound)
        if (any(bound < 0)) error stop 'Negative prescribed RT0 error bound'
        call h5_close(fid)
    contains
        subroutine expect_integer(name, expected)
            character(len=*), intent(in) :: name
            integer, intent(in) :: expected
            call h5_get(fid, group // name, value)
            if (value /= expected) error stop 'Prescribed RT0 integer metadata mismatch'
        end subroutine expect_integer
        subroutine expect_label(name, expected)
            character(len=*), intent(in) :: name, expected
            call h5_get(fid, group // name, label)
            if (trim(label) /= expected) error stop 'Prescribed RT0 unit mismatch'
        end subroutine expect_label
        subroutine require_label(name)
            character(len=*), intent(in) :: name
            call h5_get(fid, group // name, label)
            if (len_trim(label) == 0) error stop 'Missing prescribed RT0 provenance'
        end subroutine require_label
        subroutine read_real(name, out)
            character(len=*), intent(in) :: name
            real(dp), intent(out) :: out(:)
            call require_shape(fid, name, [size(out)])
            call h5_get(fid, group // name, out)
            if (.not. all(ieee_is_finite(out))) error stop 'Nonfinite prescribed RT0 data'
        end subroutine read_real
        subroutine check_real(name, expected)
            character(len=*), intent(in) :: name
            real(dp), intent(in) :: expected(:)
            allocate(values(size(expected)))
            call read_real(name, values)
            if (any(values /= expected)) error stop 'Prescribed RT0 geometry mismatch'
            deallocate(values)
        end subroutine check_real
        subroutine check_int(name, expected)
            character(len=*), intent(in) :: name
            integer, intent(in) :: expected(:, :)
            call require_shape(fid, name, shape(expected))
            allocate(imat(size(expected, 1), size(expected, 2)))
            call h5_get(fid, group // name, imat)
            if (any(imat /= expected)) error stop 'Prescribed RT0 topology mismatch'
            deallocate(imat)
        end subroutine check_int
    end subroutine read_prescribed_rt0

    subroutine require_shape(fid, name, expected)
        integer(HID_T), intent(in) :: fid
        character(len=*), intent(in) :: name
        integer, intent(in) :: expected(:)
        integer(HSIZE_T) :: dimensions(size(expected))
        integer(SIZE_T) :: type_size
        integer :: kind, ierr
        call h5_get_dataset_info(fid, group // name, dimensions, kind, type_size, ierr)
        if (ierr /= 0) error stop 'Prescribed RT0 missing or malformed dataset'
        if (any(dimensions /= int(expected, HSIZE_T))) &
            error stop 'Prescribed RT0 dataset shape mismatch'
    end subroutine require_shape

    subroutine check_prescribed_toroidal(actual, reference, bound)
        complex(dp), intent(in) :: actual(:), reference(:)
        real(dp), intent(in) :: bound(:)
        if (size(actual) /= size(reference)) error stop 'Toroidal reference shape'
        if (size(bound) /= size(actual)) error stop 'Toroidal bound shape'
        if (.not. all(ieee_is_finite(real(actual)))) error stop 'Nonfinite RT0 toroidal real'
        if (.not. all(ieee_is_finite(aimag(actual)))) error stop 'Nonfinite RT0 toroidal imag'
        if (any(abs(actual - reference) > bound)) &
            error stop 'Prescribed RT0 independent toroidal reference failed'
    end subroutine check_prescribed_toroidal

    subroutine write_prescribed_source(output, input, defect, bound)
        character(len=*), intent(in) :: output, input
        real(dp), intent(in) :: defect(:), bound(:)
        integer(HID_T) :: fid
        call h5_open_rw(output, fid)
        call h5_create_parent_groups(fid, 'vac/source/')
        call h5_add(fid, 'vac/source/schema_version', 1)
        call h5_add(fid, 'vac/source/input_file', input)
        call h5_add(fid, 'vac/source/description', &
            'Prescribed magnetic RT0 field; discrete consistency does not certify vacuum curl')
        call h5_add(fid, 'vac/source/toroidal_defect_G', defect)
        call h5_add(fid, 'vac/source/toroidal_error_bound_G', bound)
        call h5_close(fid)
    end subroutine write_prescribed_source
end module prescribed_rt0_input
