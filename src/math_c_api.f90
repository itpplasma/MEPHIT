module math_c_api
  ! C-callable wrappers around the libneo math kit, replacing GSL calls in
  ! the C sources. No module state; reentrant as long as the C callback is.

  use iso_fortran_env, only: dp => real64
  use iso_c_binding, only: c_int, c_double, c_ptr, c_funptr, c_f_procpointer

  implicit none

  private

  abstract interface
    function c_integrand(x, params) bind(C)
      import :: c_double, c_ptr
      real(c_double), intent(in), value :: x
      type(c_ptr), intent(in), value :: params
      real(c_double) :: c_integrand
    end function c_integrand
  end interface

contains

  !> Adaptive Gauss-Kronrod 21-point quadrature of \p f over [\p a, \p b] with
  !> QUADPACK QAG semantics, matching gsl_integration_qag with GSL_INTEG_GAUSS21.
  subroutine mephit_integrate_gk(f, params, a, b, epsabs, epsrel, limit, &
    result, abserr, ierr) bind(C, name = 'mephit_integrate_gk')
    use neo_gauss_kronrod, only: integrate_gk
    type(c_funptr), intent(in), value :: f
    type(c_ptr), intent(in), value :: params
    real(c_double), intent(in), value :: a, b, epsabs, epsrel
    integer(c_int), intent(in), value :: limit
    real(c_double), intent(out) :: result, abserr
    integer(c_int), intent(out) :: ierr
    procedure(c_integrand), pointer :: f_c
    integer :: ierr_f

    call c_f_procpointer(f, f_c)
    call integrate_gk(wrapped_integrand, a, b, epsabs, epsrel, result, abserr, &
      ierr_f, key = 21, limit = int(limit))
    ierr = int(ierr_f, c_int)

  contains

    function wrapped_integrand(x) result(fx)
      real(dp), intent(in) :: x
      real(dp) :: fx

      fx = f_c(real(x, c_double), params)
    end function wrapped_integrand
  end subroutine mephit_integrate_gk

  !> Gauss-Legendre nodes and weights on [\p a, \p b] in ascending node order,
  !> matching gsl_integration_glfixed_point.
  subroutine mephit_gauss_legendre(n, a, b, x, w) &
    bind(C, name = 'mephit_gauss_legendre')
    use neo_gauss_quadrature, only: gauss_legendre_ab
    integer(c_int), intent(in), value :: n
    real(c_double), intent(in), value :: a, b
    real(c_double), intent(out) :: x(n), w(n)

    call gauss_legendre_ab(int(n), a, b, x, w)
  end subroutine mephit_gauss_legendre

end module math_c_api
