module f_regtests
  implicit none
  private

  interface compare
     module procedure compare_l, compare_i, compare_f, compare_d, compare_s
     module procedure compare_ai, compare_af, compare_ad
  end interface compare

  public :: run, compare, verify, ensure_error
contains
  subroutine run(func)
    use yaml_output
    use dictionaries
    implicit none
    interface
       subroutine func(test)
         character(len = *), intent(out) :: test
       end subroutine func
    end interface

    character(len = max_field_length) :: test
    type(dictionary), pointer :: errors

    nullify(errors)
    call f_err_open_try()
    call func(test)
    call f_err_close_try(errors)

    call yaml_sequence(advance = "no")
    if (associated(errors)) then
       call yaml_mapping_open(test)
       call yaml_map("status", "failed")
       call yaml_map("error", errors)
       call yaml_mapping_close()
       call dict_free(errors)
    else
       call yaml_map(test, "succeed")
    end if
  end subroutine run

  subroutine verify(val, label)
    use dictionaries
    implicit none
    logical, intent(in) :: val
    character(len = *), intent(in) :: label

    if (.not. val) call f_err_throw(label // ": failed")
  end subroutine verify

  subroutine ensure_error(id, name)
    use dictionaries
    use yaml_strings
    implicit none
    integer, intent(in), optional :: id
    character(len = *), intent(in), optional :: name

    integer :: errid

    if (.not. f_err_check(id, name)) then
       if (present(name)) then
          call f_err_throw("error " // name // " expected")
       else if (present(id)) then
          call f_err_throw("error (" // trim(yaml_toa(id)) // ") expected")
       else
          call f_err_throw("error expected")
       end if
    else
       errid = f_err_pop(id, name)
    end if
  end subroutine ensure_error
  
  subroutine compare_l(val, expected, label)
    use dictionaries
    use yaml_strings
    implicit none
    logical, intent(in) :: val, expected
    character(len = *), intent(in), optional :: label

    if (val .neqv. expected) then
       if (present(label)) then
          call f_err_throw(label // ": expected " // yaml_toa(expected) // ", value " // yaml_toa(val))
       else
          call f_err_throw("expected " // yaml_toa(expected) // ", value " // yaml_toa(val))
       end if
    end if
  end subroutine compare_l

  subroutine compare_i(val, expected, label)
    use dictionaries
    use f_precisions, only: f_integer
    use yaml_strings
    implicit none
    integer(f_integer), intent(in) :: val, expected
    include 'compare-scalar-inc.f90'
  end subroutine compare_i

  subroutine compare_f(val, expected, label, tol)
    use dictionaries
    use f_precisions, only: f_simple
    use yaml_strings
    implicit none
    real(f_simple), intent(in) :: val, expected
    real, intent(in), optional :: tol
    include 'compare-tol-scalar-inc.f90'
  end subroutine compare_f

  subroutine compare_d(val, expected, label, tol)
    use dictionaries
    use f_precisions, only: f_double
    use yaml_strings
    implicit none
    real(f_double), intent(in) :: val, expected
    real, intent(in), optional :: tol
    include 'compare-tol-scalar-inc.f90'
  end subroutine compare_d

  subroutine compare_s(val, expected, label)
    use dictionaries
    use yaml_strings
    implicit none
    character(len = *), intent(in) :: val, expected
    include 'compare-scalar-inc.f90'
  end subroutine compare_s

  subroutine compare_ai(val, expected, label)
    use dictionaries
    use f_precisions, only: f_integer
    use yaml_strings
    implicit none
    integer(f_integer), dimension(:), intent(in) :: val, expected
    include 'compare-array-inc.f90'
  end subroutine compare_ai

  subroutine compare_af(val, expected, label, tol)
    use dictionaries
    use f_precisions, only: f_simple
    use yaml_strings
    implicit none
    real(f_simple), dimension(:), intent(in) :: val, expected
    real, intent(in), optional :: tol
    include 'compare-tol-array-inc.f90'
  end subroutine compare_af

  subroutine compare_ad(val, expected, label, tol)
    use dictionaries
    use f_precisions, only: f_double
    use yaml_strings
    implicit none
    real(f_double), dimension(:), intent(in) :: val, expected
    real, intent(in), optional :: tol
    include 'compare-tol-array-inc.f90'
  end subroutine compare_ad
end module f_regtests
