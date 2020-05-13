    character(len = *), intent(in), optional :: label
    logical :: failure

    if (f_err_raise(size(val) /= size(expected), "size mismatch")) return
    if (present(tol)) then
       failure = any(abs(val - expected) >= tol)
    else
       failure = any(val /= expected)
    end if
    if (failure) then
       if (present(label)) then
          call f_err_throw(label // ": maximum difference of " // yaml_toa(maxval(abs(val - expected))))
       else
          call f_err_throw("maximum difference of " // yaml_toa(maxval(abs(val - expected))))
       end if
    end if
