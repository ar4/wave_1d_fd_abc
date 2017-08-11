module oneway_plain

  implicit none

contains

  subroutine step(f1, f2, model_padded, dt, dx, sources, sources_x,    &
      num_steps, pad_width)

    real, intent (in out), dimension (:) :: f1
    real, intent (in out), dimension (:) :: f2
    real, intent (in), dimension (:) :: model_padded
    real, intent (in) :: dt
    real, intent (in) :: dx
    real, intent (in), dimension (:, :) :: sources
    integer, intent (in), dimension (:) :: sources_x
    integer, intent (in) :: num_steps
    integer, intent (in) :: pad_width

    integer :: step_idx
    logical :: even

    do step_idx = 1, num_steps
    even = (mod (step_idx, 2) == 0)
    if (even) then
      call one_step(f2, f1, model_padded, dt, dx, sources, sources_x,  &
        step_idx, pad_width)
    else
      call one_step(f1, f2, model_padded, dt, dx, sources, sources_x,  &
        step_idx, pad_width)
    end if
    end do

  end subroutine step


  subroutine one_step(f, fp, model_padded, dt, dx, sources, sources_x, &
      step_idx, pad_width)

    real, intent (in), dimension (:) :: f
    real, intent (in out), dimension (:) :: fp
    real, intent (in), dimension (:) :: model_padded
    real, intent (in) :: dt
    real, intent (in) :: dx
    real, intent (in), dimension (:, :)  :: sources
    integer, intent (in), dimension (:) :: sources_x
    integer, intent (in) :: step_idx
    integer, intent (in) :: pad_width

    integer :: i
    integer :: nx_padded
    integer :: num_sources

    nx_padded = size(f)
    num_sources = size(sources, dim=1)

    do i = pad_width + 1, nx_padded - pad_width
    call fd_interior(f, fp, model_padded, dt, dx, i)
    end do

    ! source term
    do i = 1, num_sources
    call add_source(fp, model_padded, dt, sources(i, step_idx),        &
      sources_x(i), pad_width)
    end do

  end subroutine one_step


  subroutine fd_interior(f, fp, model_padded, dt, dx, i)

    real, intent (in), dimension (:) :: f
    real, intent (in out), dimension (:) :: fp
    real, intent (in), dimension (:) :: model_padded
    real, intent (in) :: dt
    real, intent (in) :: dx
    integer, intent (in) :: i

    real :: f_x

    f_x = first_x_deriv(f, i, dx, 1)

    fp(i) = -model_padded(i) * dt * f_x + f(i)

  end subroutine fd_interior


  subroutine add_source(fp, model_padded, dt, source, source_x,        &
      total_pad)

    real, intent (in out), dimension (:) :: fp
    real, intent (in), dimension (:) :: model_padded
    real, intent (in) :: dt
    real, intent (in)  :: source
    integer, intent (in) :: source_x
    integer, intent (in) :: total_pad

    integer :: sx

    sx = source_x + total_pad + 1;
    fp(sx) = fp(sx) + (model_padded(sx) * dt * source)

  end subroutine add_source


  pure function first_x_deriv(f, i, dx, direction)

    real, intent (in), dimension (:) :: f
    integer, intent (in) :: i
    real, intent (in) :: dx
    integer, intent (in) :: direction

    real :: first_x_deriv

    first_x_deriv = 0.0

    if (direction == 1) then
      first_x_deriv = (f(i) - f(i-1))/dx
    else if (direction == -1) then
      first_x_deriv = (f(i) - f(i+1))/dx
    end if

  end function first_x_deriv

end module oneway_plain