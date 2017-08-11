module oneway

        implicit none

contains

        subroutine step(f1, f2, model_padded, dt, dx, sources,         &
                        sources_x, num_steps, pml_width, pad_width)

                real, intent (in out), dimension (:) :: f1
                real, intent (in out), dimension (:) :: f2
                real, intent (in), dimension (:) :: model_padded
                real, intent (in) :: dt
                real, intent (in) :: dx
                real, intent (in), dimension (:, :) :: sources
                integer, intent (in), dimension (:) :: sources_x
                integer, intent (in) :: num_steps
                integer, intent (in) :: pml_width
                integer, intent (in) :: pad_width

                integer :: step_idx
                logical :: even

                do step_idx = 1, num_steps
                even = (mod (step_idx, 2) == 0)
                if (even) then
                        call one_step(f2, f1, model_padded, dt, dx,    &
                                sources, sources_x, step_idx,          &
                                pml_width, pad_width)
                else
                        call one_step(f1, f2, model_padded, dt, dx,    &
                                sources, sources_x, step_idx,          &
                                pml_width, pad_width)
                end if
                end do

        end subroutine step


        subroutine one_step(f, fp, model_padded, dt, dx, sources,      &
                        sources_x, step_idx, abc_width, pad_width)

                real, intent (in), dimension (:) :: f
                real, intent (in out), dimension (:) :: fp
                real, intent (in), dimension (:) :: model_padded
                real, intent (in) :: dt
                real, intent (in) :: dx
                real, intent (in), dimension (:, :)  :: sources
                integer, intent (in), dimension (:) :: sources_x
                integer, intent (in) :: step_idx
                integer, intent (in) :: abc_width
                integer, intent (in) :: pad_width

                integer :: i
                integer :: nx_padded
                integer :: num_sources
                real :: lambda

                nx_padded = size(f)
                num_sources = size(sources, dim=1)

                ! left absorbing boundary
                do i = pad_width + 1, pad_width + 1 + abc_width
                lambda = max(real(i - pad_width - 1-4) / (abc_width-4),&
                        0.0)
                call fd_abc(f, fp, model_padded, dt, dx, i, -1, lambda)
                end do

                ! interior (no absorption)
                do i = pad_width + abc_width + 2,                      &
                        nx_padded - abc_width - pad_width
                call fd_interior(f, fp, model_padded, dt, dx, i)
                end do

                ! right absorbing boundary
                do i = nx_padded - abc_width - pad_width + 1,          &
                        nx_padded - pad_width
                lambda = max(real(nx_padded - pad_width - i-4) /       &
                        (abc_width-4), 0.0)
                call fd_abc(f, fp, model_padded, dt, dx, i, +1, lambda)
                end do

                ! source term
                do i = 1, num_sources
                call add_source(fp, model_padded, dt,                  &
                        sources(i, step_idx), sources_x(i),            &
                        abc_width + pad_width)
                end do

        end subroutine one_step


        subroutine fd_interior(f, fp, model_padded, dt, dx, i)

                real, intent (in), dimension (:) :: f
                real, intent (in out), dimension (:) :: fp
                real, intent (in), dimension (:) :: model_padded
                real, intent (in) :: dt
                real, intent (in) :: dx
                integer, intent (in) :: i

                real :: f_xx

                f_xx = second_x_deriv(f, i, dx)
                fp(i) = model_padded(i)**2 * dt**2 * f_xx +            &
                        2 * f(i) - fp(i)

        end subroutine fd_interior


        subroutine fd_abc(f, fp, model_padded, dt, dx, i, direction,   &
                        lambda)

                real, intent (in), dimension (:) :: f
                real, intent (in out), dimension (:) :: fp
                real, intent (in), dimension (:) :: model_padded
                real, intent (in) :: dt
                real, intent (in) :: dx
                integer, intent (in) :: i
                integer, intent (in) :: direction
                real, intent (in) :: lambda

                real :: f_xx
                real :: f_x
                real :: two_way
                real :: one_way

                f_xx = second_x_deriv(f, i, dx)
                f_x = first_x_deriv(f, i, dx, direction)

                two_way = model_padded(i)**2 * dt**2 * f_xx +          &
                        2 * f(i) - fp(i)
                one_way = -model_padded(i) * dt * f_x + f(i)

                fp(i) = two_way * lambda + one_way * (1.0 - lambda)

        end subroutine fd_abc


        subroutine add_source(fp, model_padded, dt, source, source_x,  &
                        total_pad)

                real, intent (in out), dimension (:) :: fp
                real, intent (in), dimension (:) :: model_padded
                real, intent (in) :: dt
                real, intent (in)  :: source
                integer, intent (in) :: source_x
                integer, intent (in) :: total_pad

                integer :: sx

                sx = source_x + total_pad + 1;
                fp(sx) = fp(sx) + (model_padded(sx)**2 * dt**2 * source)

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


        pure function second_x_deriv(f, i, dx)

                real, intent (in), dimension (:) :: f
                integer, intent (in) :: i
                real, intent (in) :: dx

                real :: second_x_deriv

                !second_x_deriv = (                                     &
                !  -735*f(i-8)+15360*f(i-7)                             &
                !  -156800*f(i-6)+1053696*f(i-5)                        & 
                !  -5350800*f(i-4)+22830080*f(i-3)                      & 
                !  -94174080*f(i-2)+538137600*f(i-1)                    & 
                !  -924708642*f(i+0)                                    & 
                !  +538137600*f(i+1)-94174080*f(i+2)                    & 
                !  +22830080*f(i+3)-5350800*f(i+4)                      & 
                !  +1053696*f(i+5)-156800*f(i+6)                        & 
                !  +15360*f(i+7)-735*f(i+8))/(302702400*dx**2)
                second_x_deriv = (f(i-1) - 2*f(i) + f(i+1))/dx**2


        end function second_x_deriv

end module oneway
