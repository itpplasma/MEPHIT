!> A resonance with damping switched on and zero current width must be reported.
program test_damping_warning
    use mephit_conf, only: logger
    use mephit_mesh, only: warn_undamped_resonances
    implicit none
    character(len=*), parameter :: logfile = 'test_damping_warning.log'
    character(len=1024) :: line
    integer :: fid, ios, nwarn
    logical :: ok

    call logger%init(logfile, 2, .true.)
    ! m = 2 and 4 have zero width
    call warn_undamped_resonances(.true., 2, [0d0, 1d0, 0d0])
    ! no damping requested, no warning
    call warn_undamped_resonances(.false., 2, [0d0, 0d0])
    ! all widths set, no warning
    call warn_undamped_resonances(.true., 2, [1d0, 2d0])
    call logger%deinit

    ok = .true.
    nwarn = 0
    open (newunit=fid, file=logfile, status='old', action='read')
    do
        read (fid, '(a)', iostat=ios) line
        if (ios /= 0) exit
        if (index(line, 'Warning: damp') == 0) cycle
        nwarn = nwarn + 1
        if (nwarn == 1) ok = ok .and. index(line, 'm = 2,') > 0
        if (nwarn == 2) ok = ok .and. index(line, 'm = 4,') > 0
    end do
    close (fid, status='delete')
    ok = ok .and. nwarn == 2
    if (.not. ok) error stop 'damping warning test failed'
end program test_damping_warning
