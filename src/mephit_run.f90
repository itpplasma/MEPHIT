program mephit_run

  use mephit_iter, only: mephit_main
  use hdf5_tools, only: HID_T, h5_init, h5_deinit, h5overwrite, h5_defer_close

  implicit none

  character(len = 1024) :: str_runmode, config, suffix
  integer :: runmode

  call process_command_arguments(1, str_runmode, 'numeric runmode', .false.)
  read (str_runmode, *) runmode
  call process_command_arguments(2, config, 'path to config file', .false.)
  call process_command_arguments(3, suffix, 'file basename suffix', .true.)

  call h5_init
  h5overwrite = .true.
  ! MEPHIT updates one output image through many close/reopen calls.  Keep
  ! that image in Fortio memory and flush it once during mephit_deinit.
  h5_defer_close = .true.
  call mephit_main(runmode, config, suffix)
  call h5_deinit
  h5_defer_close = .false.

contains

  subroutine process_command_arguments(pos, val, name, accept_empty)
    use iso_fortran_env, only: error_unit
    integer, intent(in) :: pos
    character(len=*), intent(inout) :: val
    character(len=*), intent(in) :: name
    logical, intent(in) :: accept_empty

    if (command_argument_count() < pos) then
      write (error_unit, '("Expected ", a, " at argument position ", i0)') name, pos
      error stop
    end if
    call get_command_argument(pos, val)
    if (.not. accept_empty .and. len_trim(val) == 0) then
      write (error_unit, '("Empty string for ", a, " at argument position ", i0)') name, pos
      error stop
    end if
  end subroutine process_command_arguments

end program mephit_run
